library(lubridate)
library(rgdal)
library(dplyr)
library(Hmisc)
library(plotly)
library(akima)
library(sp)
library(pracma)
source("~/Documents/CODES/npsat_engine/Rwrkspc/npsat_Rfun.R")
source("CVHM_budget_functions.R")

mile <- 1609.34
sqml <- mile*mile
cvhm_tm <- seq.Date(from = as.Date(paste0(1961, "/", 4, "/1")), to = as.Date(paste0(2003, "/", 9, "/1")),by = "month")

startTime <-  as.Date(paste0(1995, "/", 7, "/1")) # 1999/10
endTime <-  as.Date(paste0(1999, "/", 12, "/1")) # 2003/9

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
Budget_coef <- 0.5 
opt <- data.frame(prefix, timestring, simFolder, cbcf_path, std_Htol, LU, RCH_data, STRM_data, WELL_data, Budget_coef)
rm(prefix, timestring, simFolder, cbcf_path, gis_path, std_Htol, LU, RCH_data, STRM_data, WELL_data, Budget_coef)



# prepare the zones variable
Basin_dwr_ids <- vector(mode = "list", length = 3)
Basin_dwr_ids[[1]] <- 1:7
Basin_dwr_ids[[2]] <- 8:13
Basin_dwr_ids[[3]] <- 14:21
farms <- readOGR(dsn = "../gis_data/", layer = "FMP")

group_zones <- array(data = NA, dim = c(441,98))
for (i in 1:length(Basin_dwr_ids)) {
  for (j in 1:length(Basin_dwr_ids[[i]])) {
    ids <- which(farms$dwr_sbrgns == Basin_dwr_ids[[i]][j])
    group_zones[ cbind(farms$ROW[ids], farms$COLUMN_[ids]) ] <- i
  }
}

# Average and equalize budgets
sumBUD <- AverageBudgets(startTime, endTime, cvhm_tm, group_zones, opt)

# Generate wells
Wells <- GenerateWells(sumBUD[[3]]/sumBUD[[4]])



# Write files ------------
if (!dir.exists(as.character(opt$simFolder))){
  dir.create(as.character(opt$simFolder))
}

# Write Recharge function
load("chvm2DinterpNodes.RData")
rch_fl <- paste0(opt$simFolder, opt$prefix, '_', opt$timestring, '_bf' ,opt$Budget_coef, '_RCH.npsat')
write2DscatterFile(rch_fl, (sumBUD[[1]]/sumBUD[[4]])/sqml, cvhm_mesh_coords, buffer_nodes, opt)

# Write boundary conditions and initial top elevation
load(file = "CVHM_headTop.RData")
load(file = "CVHM_LNS_IJ.RData")
headAv <- writeHeadBC(startTime, endTime, cvhm_tm, headTop, PNTS, LNS_IJ, opt )
init_top_fl <- paste0(opt$simFolder, opt$prefix, '_', opt$timestring, '_initTop_h_' ,opt$std_Htol, '.npsat')
write2DscatterFile(init_top_fl, headAv, cvhm_mesh_coords, buffer_nodes, opt)





