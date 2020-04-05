library(rgdal)
library(plotly)
library(gwtools)
library(pracma)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# CREATE MESH INPUT FILE -----------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load the mesh shapefile and write it to 
cvhm_mesh_shp <- readOGR(dsn = "../gis_data/", layer = "CVHM_mesh_3310")
# create a unique list of nodes and mesh elents that point of the nodelist
ND <- matrix(data = NA, nrow = 50000, ncol = 2)
MSH <- matrix(data = NA, nrow = length(cvhm_mesh_shp), ncol = 4)
cnt_nd <- 0
for (i in 1:length(cvhm_mesh_shp)) {
  if (length(cvhm_mesh_shp@polygons[[i]]@Polygons) != 1){
    print(paste("There are more polygons in", i))
    break
  }
  if (dim(cvhm_mesh_shp@polygons[[i]]@Polygons[[1]]@coords)[1] != 5){
    print(paste("The nodes are not 5 in", i))
    break
  }
  
  for (j in 1:4) {
    x <- cvhm_mesh_shp@polygons[[i]]@Polygons[[1]]@coords[j,1]
    y <- cvhm_mesh_shp@polygons[[i]]@Polygons[[1]]@coords[j,2]
    if (cnt_nd == 0){
      cnt_nd <- cnt_nd + 1
      ND[cnt_nd,] <- c(x,y)
      MSH[i,j] <- cnt_nd
    }
    else{
      dst <- sqrt((ND[1:cnt_nd,1] - x)^2 + (ND[1:cnt_nd,2] - y)^2)
      id <- which(dst < 0.1)
      if (length(id) == 0){
        cnt_nd <- cnt_nd + 1
        ND[cnt_nd,] <- c(x,y)
        MSH[i,j] <- cnt_nd
      }
      else{
        if (length(id) > 1){
          print(paste("More than one nodes have the same coordinates in element", i))
          break
        }
        MSH[i,j] <- id
      }
    }
  }
}
ND <- ND[1:cnt_nd,]
ND <- cbind(ND, zeros(n = dim(ND)[1],m = 1))
# Check the area of the elements
a <- vector(mode = "numeric", length = dim(MSH)[1])
for (i in 1:dim(MSH)[1]) {
  a[i] <- polyarea(ND[MSH[i,],1],ND[MSH[i,],2])
  if (a[i] < 0){
    MSH[i,] <- MSH[i,c(1,4,3,2)]
  }
}

if (!file.exists("CVHM_msh_3310.npsat") == TRUE){
  gwtools::npsat.writeMesh(filename = "CVHM_msh_3310.npsat", nd = ND, msh = MSH-1)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# CREATE BUFFER NODES -------------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Create a buffer layer of nodes to make sure that all interpolation will return a value
cvhm_buffer <- readOGR(dsn =  "../gis_data/", layer = "CVHM_mesh_buffer_3310")
ND_buffer <- matrix(data = NA, nrow = 5000, ncol = 2)
cnt_nd <- 0
for (i in 1:dim(cvhm_buffer@polygons[[1]]@Polygons[[1]]@coords)[1]) {
  if (cnt_nd == 0){
    cnt_nd <- cnt_nd + 1
    ND_buffer[cnt_nd,] <- cvhm_buffer@polygons[[1]]@Polygons[[1]]@coords[i,]
  }
  else{
    x <- cvhm_buffer@polygons[[1]]@Polygons[[1]]@coords[i,1]
    y <- cvhm_buffer@polygons[[1]]@Polygons[[1]]@coords[i,2]
    dst <- sqrt((ND_buffer[1:cnt_nd,1] - x)^2 + (ND_buffer[1:cnt_nd,2] - y)^2)
    id <- which(dst < 0.1)
    if (length(id) == 0){
      cnt_nd <- cnt_nd + 1
      ND_buffer[cnt_nd,] <- c(x,y)
    }
  }
}
ND_buffer <- ND_buffer[1:cnt_nd,]


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# WRITE THE BOTTOM ELEVATION FILE--------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load the information for the bottom of the CV
bas_active <- readOGR(dsn =  "../gis_data/", layer = "BAS_active_3310")
cvhm_bottom <- readOGR(dsn =  "../gis_data/", layer = "DIS_3310")
active_lin_ind <- gwtools::sub2ind(bas_active$ROW, bas_active$COLUMN_, 441)
dis_lin_ind <- gwtools::sub2ind(cvhm_bottom$ROW, cvhm_bottom$COLUMN_,441)
tmp <- match(dis_lin_ind, active_lin_ind)
dis_act_id <- which(!is.na(tmp))
View(tmp[dis_act_id])
cvhm_bottom <- subset(cvhm_bottom, !is.na(tmp))
bottom_points <- matrix(data = NA, nrow = length(cvhm_bottom), ncol = 3)
for (i in 1:length(cvhm_bottom)) {
  bottom_points[i,] <- c(
  mean(cvhm_bottom@polygons[[i]]@Polygons[[1]]@coords[1:4,1]),
  mean(cvhm_bottom@polygons[[i]]@Polygons[[1]]@coords[1:4,2]),
  cvhm_bottom@data$cvr2lay10b[i]
  )
}

# assign elevations to the buffer nodes
buff_Bot_elev <- matrix(data = NA, nrow = dim(ND_buffer)[1], ncol = 3)
for (i in 1:dim(ND_buffer)[1]) {
  dst <- sqrt((ND_buffer[i,1]-bottom_points[,1])^2 + (ND_buffer[i,2]-bottom_points[,2])^2)
  buff_Bot_elev[i,] <- c(
    ND_buffer[i,1], ND_buffer[i,2],
    bottom_points[which.min(dst),3]
  )
}

gwtools::npsat.WriteScattered(filename = "CVHM_Bottom_3310.npsat",PDIM = 2,
                              TYPE = "HOR", MODE = "SIMPLE",
                              DATA = rbind(bottom_points,buff_Bot_elev))



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# WRITE HYDRAULIC CONDUCTIVITY -------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
HK <- vector(mode = "list", length = 10)
VK <- vector(mode = "list", length = 10)
for (i in 1:10) {
  HK[[i]] <- as.matrix(read.table(file = paste0("../hyd_cond/HK_lyr", i, ".txt"), header = F))
  VK[[i]] <- as.matrix(read.table(file = paste0("../hyd_cond/VK_lyr", i, ".txt"), header = F))
}
HK_mat <- matrix(data = NA, nrow = dim(bottom_points)[1], ncol = 10)
VK_mat <- matrix(data = NA, nrow = dim(bottom_points)[1], ncol = 10)
lin_ind <- gwtools::sub2ind(cvhm_bottom$ROW, cvhm_bottom$COLUMN_, 441)
for (i in 1:10) {
  HK_mat[,i] <- HK[[i]][lin_ind] 
  VK_mat[,i] <- VK[[i]][lin_ind]
}

# interpolate conductivity and elevation to the buffer nodes

interLayElev <- cvhm_bottom@data[,6:14]
interLayElev_buff <- matrix(data = NA, nrow = dim(ND_buffer)[1], ncol = dim(interLayElev)[2])
HK_buff <- matrix(data = NA, nrow = dim(ND_buffer)[1], ncol = 10)
VK_buff <- matrix(data = NA, nrow = dim(ND_buffer)[1], ncol = 10)
for (i in 1:dim(ND_buffer)[1]) {
  dst <- sqrt((ND_buffer[i,1]-bottom_points[,1])^2 + (ND_buffer[i,2]-bottom_points[,2])^2)
  interLayElev_buff[i,] <- as.numeric(interLayElev[which.min(dst),])
  HK_buff[i,] <- as.numeric(HK_mat[which.min(dst),])
  VK_buff[i,] <- as.numeric(VK_mat[which.min(dst),])
}


HKdata <- cbind(rbind(bottom_points[,1:2],buff_Bot_elev[,1:2]), c(HK_mat[,1], HK_buff[,1]))
VKdata <- cbind(rbind(bottom_points[,1:2],buff_Bot_elev[,1:2]), c(VK_mat[,1], VK_buff[,1]))
for (i in 2:10) {
  HKdata <- cbind(HKdata, c(interLayElev[,i-1], interLayElev_buff[,i-1]))
  HKdata <- cbind(HKdata, c(HK_mat[,i], HK_buff[,i]))
  VKdata <- cbind(VKdata, c(interLayElev[,i-1], interLayElev_buff[,i-1]))
  VKdata <- cbind(VKdata, c(VK_mat[,i], VK_buff[,i]))
}
# Isolate only the Conductivity values
HKtemp <- HKdata[, seq(3,21,2) ]
VKtemp <- VKdata[, seq(3,21,2) ]
# There are in here and there zero conductivity values. 
# We will replace those values with the nonzero value of the above layer
for (i in 1:dim(HKtemp)[2]) {
  zr_id <- which(HKtemp[,i] == 0)
  if (i == 0 & length(zr_id) > 0){
    print(paste("The first layer has zero HK values"))
  }
  else{
    if (length(zr_id) == 0){
      next
    }
    HKtemp[zr_id,i] <- HKtemp[zr_id,i-1]
  }
  zr_id <- which(VKtemp[,i] == 0)
  if (i == 0 & length(zr_id) > 0){
    print(paste("The first layer has zero VK values"))
  }
  else{
    if (length(zr_id) == 0){
      next
    }
    VKtemp[zr_id,i] <- VKtemp[zr_id,i-1]
  }
}
# Substitute the modified values
HKdata[, seq(3,21,2) ] <- HKtemp
VKdata[, seq(3,21,2) ] <- VKtemp


gwtools::npsat.WriteScattered(filename = "CVHM_HK_3310.npsat",
                              PDIM = 2,TYPE = "FULL",MODE = "STRATIFIED", DATA = HKdata)
gwtools::npsat.WriteScattered(filename = "CVHM_VK_3310.npsat",
                              PDIM = 2,TYPE = "FULL",MODE = "STRATIFIED", DATA = VKdata)

{# Test HK
  HK <- gwtools::npsat.ReadScattered("CVHM_HK_3310.npsat")
  
}
