source("../../C2VsimCG/Rwrkspc/c2vsim_io.R")
# This is a complamentary file to the CVHM2NPSAT markdown documentation file

############### LOAD DATA ############################
# From the ReadCVHMfiles.R run the first section that loads the budget file.
# The BUD variable is a list of 510 time steps of the CVHM simulation.
# In each time step there is a list of 11 budget terms.
# The storage is in the 1st and 7 th element of the inner list

# For the C2VSIM Fine grid we load saved data from another project
load(file = "../../C2VSIM_FG_OR/C2Vsim_FG_v2/NPSAT/C2VSIM_STORAGE_TS.RData")
load(file = "../../C2VsimCG/Rwrkspc/cumul_GWBUD_base.RData")
c2vsimGWBUD_FG <- c2vsim.readGWBUD(filename =  "../../C2VSIM_FG_OR/C2Vsim_FG_v2/c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Results/C2VSimFG_GW_Budget.bud",
                                NtimeSteps = 504, Nskip = c(8,9*ones(20,1)),  CG = F)
c2vsimGWBUD_FG <- c2vsim.cumGWBUD(c2vsimGWBUD_FG)


# sumarize the to storage components into one time series
cvhm_tm <- seq.Date(from = as.Date(paste0(1961, "/", 4, "/1")), to = as.Date(paste0(2003, "/", 9, "/1")),by = "month")
ndays <- monthDays(cvhm_tm)
cvhm_storage_monthly <- vector(mode = "numeric", length = length(BUD))
for (i in 1:length(BUD)) {
  cvhm_storage_monthly[i] <- cvhm_storage_monthly[i] + sum(BUD[[i]][[1]][[4]])*ndays[i]
  cvhm_storage_monthly[i] <- cvhm_storage_monthly[i] + sum(BUD[[i]][[7]][[4]])*ndays[i]
}

# Calculate yearly storage

cvhm_tm_year <- seq.Date(from = as.Date(paste0(1962, "/", 9, "/1")), to = as.Date(paste0(2003, "/", 9, "/1")),by = "year")
cvhm_storage_yearly <- matrix(data = cvhm_storage_monthly[7:510], nrow = 12, ncol = length(seq(7, 510))/12)
cvhm_storage_yearly <- apply(cvhm_storage_yearly, 2, sum)

### Calculate the fine grid storage
c2vsim_storage_yearly <- matrix(data = c2vsim_storage_monthly, nrow = 12, ncol = 42)
c2vsim_storage_yearly <- apply(c2vsim_storage_yearly, 2, sum)
c2vsim_tm_year <- seq.Date(from = as.Date(paste0(1973, "/", 9, "/1")), to = as.Date(paste0(2015, "/", 9, "/1")),by = "year")
# verification 
c2vsim_storage_yearly1 <- matrix(data = c2vsimGWBUD_FG$BS-c2vsimGWBUD_FG$ES, nrow = 12, ncol = 42)
c2vsim_storage_yearly1 <- apply(c2vsim_storage_yearly1, 2, sum)


### Calculate the Coarse grid storage
c2vsimCG_storage_yearly <- matrix(data = cumGWBUD_CG$BS - cumGWBUD_CG$ES, nrow = 12, ncol = 88)
c2vsimCG_storage_yearly <- apply(c2vsimCG_storage_yearly, 2, sum)
c2vsimCG_tm_year <- seq.Date(from = as.Date(paste0(1922, "/", 9, "/1")), to = as.Date(paste0(2009, "/", 9, "/1")),by = "year")


{
  p <- plot_ly()
  p <- add_trace(p, x = cvhm_tm_year[12:42], y = -cumsum(c(0, cvhm_storage_yearly[13:42]))/1233.48/1000000, 
                 type = 'scatter', mode = 'lines', name = 'CVHM')
  p <- add_trace(p, x = c2vsim_tm_year, y = -cumsum(c(0,c2vsim_storage_yearly))/43559.9/1000000, 
                 type = 'scatter', mode = 'lines', name = 'C2VSim FG' )
  #p <- add_trace(p, x = c2vsim_tm_year, y = -cumsum(c(0,c2vsim_storage_yearly1))/1000000, type = 'scatter', mode = 'lines' )
  p <- add_trace(p, x = c2vsimCG_tm_year[52:88], y = -cumsum(c(0, c2vsimCG_storage_yearly[53:88]))/1000000,
                 type = 'scatter', mode = 'lines', name = 'C2VSim CG' )
  
  p %>%
    layout(xaxis = list(title = "Water year"), 
           yaxis = list(title = "Cumulative change in groundwater storage [MAF]"),
           legend = list(x = 0.02, y = 0.1))
}


### Create a shapefile containing the BAS active cells and add the adjusted and raw budget terms
bas_active <- readOGR(dsn = "../gis_data/", layer = "BAS_active")
bas_data <-  bas_active@data
bas_data <- bas_data[,1:4]
bas_data$rch_adj <- 0
bas_data$well_adj <- 0
bas_data$strm_adj <- 0
bas_data$rch_raw <- 0
bas_data$well_raw <- 0
bas_data$strm_raw <- 0
for (i in 1:length(bas_active)) {
  bas_data$rch_adj[i] <- sumBUD$sumRCH[bas_data$ROW[i], bas_data$COLUMN_[i]]/sumBUD$totdays
  bas_data$well_adj[i] <- sumBUD$sumWELL[bas_data$ROW[i], bas_data$COLUMN_[i]]/sumBUD$totdays
  bas_data$strm_adj[i] <- sumBUD$sumSTRM[bas_data$ROW[i], bas_data$COLUMN_[i]]/sumBUD$totdays
  bas_data$rch_raw[i] <- sumBUD$sumRCHraw[bas_data$ROW[i], bas_data$COLUMN_[i]]/sumBUD$totdays
  bas_data$well_raw[i] <- sumBUD$sumWELLraw[bas_data$ROW[i], bas_data$COLUMN_[i]]/sumBUD$totdays
  bas_data$strm_raw[i] <- sumBUD$sumSTRMraw[bas_data$ROW[i], bas_data$COLUMN_[i]]/sumBUD$totdays
}
bas_active@data <- bas_data
writeOGR(obj = bas_active, dsn = "../gis_data/", layer = "bud_92_03",driver = ogrDrivers()[18,1])

# Write the average stresses to shapefile
# read the bas active file
bas_active <- readOGR(dsn = "../gis_data/", layer = "BAS_active")
# and the shapefile that contains the farm ids
fmp <- readOGR(dsn = "../gis_data/", layer = "FMP")
# Identify the active cells in the farm shape file
rc_active <- sub2ind(bas_active$ROW, bas_active$COLUMN_, 441)
rc_fmp <- sub2ind(fmp$ROW, fmp$COLUMN_, 441)
# Remove the inactive cells
fmp$active <- 0
tmp <- match(rc_fmp, rc_active)
fmp$active[which(!is.na(tmp))] <- 1
fmp <- subset(fmp, active == 1) # THIS is doing the trick
# keep only the usefull fields
fmp_data <- fmp@data
fmp_data <- fmp_data[c(-1, -5, -(7:13))]
fmp@data <- fmp_data
# Add the recharge wells and stresses
rc_ind <- sub2ind(fmp$ROW, fmp$COLUMN_, 441)
fmp$Rch_m <- (sumBUD$sumRCH[rc_ind]/sumBUD$totdays)/(1609.34*1609.34)
fmp$Well_m3 <- (sumBUD$sumWELL[rc_ind]/sumBUD$totdays)
fmp$Strm_m <- (sumBUD$sumSTRM[rc_ind]/sumBUD$totdays)/(1609.34*1609.34)
writeOGR(obj = fmp, dsn = "Sim_1992m10_2003m9/", layer = "Average_bud_92_03",driver = ogrDrivers()[18,1])
