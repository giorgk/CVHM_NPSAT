library(rgdal)
library(Hmisc) # for the monthDays function

source("modflow_functions.R")
source("myUtils.R")

mile <- 1609.34

# CVHM time steps:
# 4/1961 - 9/2003
cvhm_tm <- seq.Date(from = as.Date(paste0(1961, "/", 4, "/1")), to = as.Date(paste0(2003, "/", 9, "/1")),by = "month")

# READ ENTIRE Budget from CVHM -----
BUD <- modflow.readCBC(filename = "../../CVHM_DATA/1996_cap/1996_cap/cbcf.out")

# Extract the Recharge component from the budget and convert it to 2D matrices ----
RCH <- modflow.gather(BUD,11)
for (i in 1:length(RCH)) {
  RCH[[i]] <- aperm(array(RCH[[i]][,2],c(98,441)),c(2,1))
}
save(RCH, file = "CVHM_rch.RData")

# Gather Stream recharge
STRM <- modflow.gather(BUD, 8)
save(STRM, file = "CVHM_strm.RData")

# Gather Wells recharge
MNW <- modflow.gather(BUD, 9)
FRW <- modflow.gather(BUD, 10) 
save(MNW, FRW, file = "CVHM_wells.RData")

HEAD <- modflow.readArrayASCII(filename = "../../CVHM_DATA/1996_cap/1996_cap/headsout.txt",nlay = 10)
headTop <- modflow.extractTopLayInfo(D = HEAD, inactive = HEAD[[1]][1,1,1])
save(headTop, file = "CVHM_headTop.RData")


# month ids to average recharge for year: e.g 2003-> 499:510
id_m <- 499:510
RCH_av <- matrix(data = 0, nrow = dim(RCH[[1]])[1], ncol = dim(RCH[[1]])[2])
ndays <- 0
for (i in id_m) {
  RCH_av <- RCH_av + RCH[[i]]*monthDays(cvhm_tm[i])
  ndays <- ndays + monthDays(cvhm_tm[i])
}
RCH_av <- RCH_av/ndays
RCH_av <- RCH_av/(mile*mile)

# Read Bas shapefile from the CVHM database
# The bas file was actually converted to a 3310 projection.
# First it was uploaded to EE, exported from EE as bas6EE and converted to 3310 in QGIS 
cvhm_bas <- readOGR(dsn = "../../CVHM_NPSAT/gis_data/bas6EE_3310.shp")
bas_data <- cvhm_bas@data
# add an extra row to remember the order
bas_data$order <- 1:dim(bas_data)[1]
# find the active and inactive cells
View(colnames(bas_data))
# the ids of the the columns with active layer indication
lay_id <- c(1,2,4,5,7,9,14,16,23)
bas_data$ACTIVE <- vector(mode = "numeric", length = dim(bas_data)[1])
bas_data$ACTIVE[which(rowSums(bas_data[,lay_id]) != 0)] <- 1

bas_data <- bas_data[,c(10,19,22,6,24,25)]

bas_data <- bas_data[with(bas_data, order(COLUMN_,ROW)),]
bas_data$RCH <- RCH_av[1:dim(bas_data)[1]]
bas_data <- bas_data[with(bas_data, order(order)),]
cvhm_bas@data <- bas_data[-5]
writeOGR(obj =  cvhm_bas, dsn = "../../SWAT_Model/rch_av_2003.shp",layer = "../../SWAT_Model/rch_av_2003",driver = "ESRI Shapefile")

# Attempt to remove inactive cells in R
cvhm_bas_active <- subset(cvhm_bas, cvhm_bas$ACTIVE == 1)
cvhm_bas_active@data <- cvhm_bas_active@data[-5]
writeOGR(obj =  cvhm_bas_active, dsn = "../../SWAT_Model/rch_av_2003.shp",layer = "../../SWAT_Model/rch_av_2003",driver = "ESRI Shapefile")

cvarea <- 0
for (i in 1:length(cvhm_bas_active)) {
  cvarea <- cvarea +  cvhm_bas_active@polygons[[i]]@area
}


# PREPARE Interpolation vertices in 2D ----------
bas_active <- readOGR(dsn = "../gis_data/", layer = "BAS_active")
cvhm_mesh_buffer <- readOGR(dsn = "../gis_data/", layer = "CVHM_Mesh_outline_buffer")
cvhm_mesh_coords <- matrix(data = NA, nrow = length(bas_active), ncol = 4)
# make a list with the centroids of the cell
for (i in 1:length(bas_active)) {
  if (length(bas_active@polygons[[i]]@Polygons) == 1){
    if (dim(bas_active@polygons[[i]]@Polygons[[1]]@coords)[1] == 5){
      cvhm_mesh_coords[i,] <-  c(colSums(bas_active@polygons[[i]]@Polygons[[1]]@coords[-5,])/4, bas_active$ROW[i], bas_active$COLUMN_[i])
    }
    else{
      print(paste(i, "doesn't have 5 coords"))
    }
  }
  else{
    print(paste(i, "doesn't have only 1 polygon"))
  }
}

# read the buffer nodes
buffer_nodes <- cvhm_mesh_buffer@polygons[[1]]@Polygons[[1]]@coords
save(cvhm_mesh_coords, buffer_nodes, file = "chvm2DinterpNodes.RData")


# for the CVHM the function works with the inputs reversed
# r-> is the column, c-> the row and m-> the number of columns
myutils.sub2ind(95,2,98)


myutils.writeRasterAscii(data = t,filename = "temp.txt",llcoord = c(0,0), nodata = 0,cellsize = 1609.34)
