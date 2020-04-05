library(lubridate)
library(rgdal)
library(dplyr)
library(Hmisc)
library(plotly)
library(interp) #library(akima)
library(sp)
library(pracma)
library(gwtools)
source("CVHM_budget_functions.R")

mile <- 1609.34
sqml <- mile*mile
cvhm_tm <- seq.Date(from = as.Date(paste0(1961, "/", 4, "/1")), to = as.Date(paste0(2003, "/", 9, "/1")),by = "month")

startTime <-  as.Date(paste0(1992, "/", 10, "/1")) # 1999/10
endTime <-  as.Date(paste0(2003, "/", 9, "/1")) # 2003/9

# Create a data frame with the options
prefix <- "cvhm"
timestring <- paste0(year(startTime),'m',  month(startTime), '_', year(endTime),'m',  month(endTime))
simFolder <- paste0('Sim_', timestring, '/')
cbcf_path <- '/media/giorgk/FIONA/UBU_BACKUP/Documents/UCDAVIS/CVHM_DATA/ClaudiaRun/'
gis_path  <- 'gis_data/'
std_Htol <- 0.5
LU <- 'LU_2000'
RCH_data <- "CVHM_rch.RData"
STRM_data <- "CVHM_strm.RData"
WELL_data <- "CVHM_wells.RData"
# This is a coefficient that determines how the budget should be equalized. 
# 0 -> The volume of water will be reduced to make the budget 0
# 1 -> The volume of water will be increased to make the budget 0
# Any value in between will adjust the volumes by reducing and increasing proportionally 
Budget_coef <- 0.0 
opt <- data.frame(prefix, timestring, simFolder, cbcf_path, std_Htol, LU, RCH_data, STRM_data, WELL_data, Budget_coef)
rm(prefix, timestring, simFolder, cbcf_path, gis_path, std_Htol, LU, RCH_data, STRM_data, WELL_data, Budget_coef)



# prepare the zones variable
Basin_dwr_ids <- vector(mode = "list", length = 3)
Basin_dwr_ids[[1]] <- 1:7
Basin_dwr_ids[[2]] <- 8:13
Basin_dwr_ids[[3]] <- 14:21
farms <- readOGR(dsn = "../gis_data/", layer = "BAS_active_3310")

group_zones <- array(data = NA, dim = c(441,98))
for (i in 1:length(Basin_dwr_ids)) {
  for (j in 1:length(Basin_dwr_ids[[i]])) {
    ids <- which(farms$dwr_sbrgns == Basin_dwr_ids[[i]][j])
    group_zones[ cbind(farms$ROW[ids], farms$COLUMN_[ids]) ] <- i
  }
}

# Create the simulation folder if doesnt exist
if (!dir.exists(as.character(opt$simFolder))){
  dir.create(as.character(opt$simFolder))
}

{ # Run budget and generate wells
  # Average and equalize budgets
  sumBUD <- AverageBudgets(startTime, endTime, cvhm_tm, group_zones, opt)
  # Generate wells
  Wells <- GenerateWells(sumBUD$sumWELL/sumBUD$totdays)
  # Save Bugdets and wells to a file
  save(list = c("sumBUD", "Wells"), file = paste0(opt$simFolder,"AdjustedBudget_", opt$Budget_coef,".RData"))
}
{ # Or load them
  load(file  = paste0(opt$simFolder,"AdjustedBudget_", opt$Budget_coef,".RData"))
}


# Write the Recharge function
load("chvm2DinterpNodes.RData") # generated in PreProcessData.R
rch_fl <- paste0(opt$simFolder, opt$prefix, '_', opt$timestring, '_bf', opt$Budget_coef, '_RCH.npsat')
write2DscatterFile(rch_fl, (sumBUD$sumRCH/sumBUD$totdays)/sqml, bas_bary_coords, buffer_nodes, opt)

# Write boundary conditions and initial top elevation
load(file = "CVHM_headTop.RData") # generated in ReadCVHMfiles.R
load(file = "CVHM_LNS_IJ.RData") # generated in PreProcessData.R
headAv <- writeHeadBC(startTime, endTime, cvhm_tm, headTop, PNTS, LNS_IJ, opt )
init_top_fl <- paste0(opt$simFolder, opt$prefix, '_', opt$timestring, '_initTop_h_' ,opt$std_Htol, '.npsat')
write2DscatterFile(init_top_fl, headAv, bas_bary_coords, buffer_nodes, opt)

# prepare Rivers
# First run the section "Preprocess Streams" from the PreProcessData.R script
# This will generate a file this is needed from the function below.
# The following function is going to create a subfolder tempHOUStreams and
# write one file for each river.
Nfiles <- writeStreams4Houdini(sumBUD$sumSTRM/sumBUD$totdays, opt)
# Then play the animation in the Houdini scene to generate the polygonal river files
# in the same folder
# After the #####_OUT.dat files have been generated run the following to create the stream input file
writeStreamsNPSATinputv1(53, sumBUD$sumSTRM/sumBUD$totdays,opt)
writeStreamsNPSATinput(Nfiles, opt)

# Make sure the bottom elevation is below the top elevation
top <- gwtools::npsat.ReadScattered(init_top_fl)
bot <- gwtools::npsat.ReadScattered("CVHM_Bottom_3310.npsat")
# make sure that the two files are written with the same order. The following dhould be 0
sum(sqrt((top$X - bot$X)^2 + (top$Y - bot$Y)^2))
neg <- which(top$V < bot$V) # find those that the bottom is higher
min(top$V[neg]) # Find the minimum top of those
bot$V[bot$V > -30] <- -30 # Set everything above a threshold a constant value
gwtools::npsat.WriteScattered("CVHM_Bottom1_3310.npsat", 2, "HOR", "SIMPLE", bot)

# write wells to file
top <- gwtools::npsat.ReadScattered(init_top_fl)
bot <- gwtools::npsat.ReadScattered("CVHM_Bottom1_3310.npsat")
writeWells(Wells, top, bot, opt, 10, 10)


{
  p <- plot_ly()
  p <- add_trace(p,x=Wells$X,y=Wells$Y, mode = 'markers')
  p <- add_trace(p,x=top$X,y=top$Y, mode = 'markers')
  p <- add_trace(p,x=t$x[ii],y=t$y[ii], mode = 'markers')
  p
}
{
  p <- plot_ly()
  p <- add_trace(p, x = cvhm_mesh_coords[,1], y = cvhm_mesh_coords[,2], mode = 'markers',marker = list(size = 10))
  p <- add_trace(p, x =buffer_nodes[,1], y = buffer_nodes[,2], mode = 'markers')
  p
}

t <- interp::interp(buffer_nodes[,1], buffer_nodes[,2], rand(length(cvhm_mesh_coords[,1]),1), 
               Wells$X, Wells$Y, linear = F, extrap = T, duplicate = "mean", output = "points")
which(is.na(t$z))

{# Test stream budget
  strm <- gwtools::npsat.ReadStreams("Sim_1992m10_2003m9/cvhm_1992m10_2003m9_bf0_STREAMS.npsat")
  sum(strm$Q)
}

{# Test wells
  wells <- gwtools::npsat.ReadWells("Sim_1992m10_2003m9/cvhm_1992m10_2003m9_bf0_WELL.npsat")
  sum(wells$Q)
}

