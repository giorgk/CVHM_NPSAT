library(rgdal)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ------ Identify the CVHM boundaries -------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
mesh_outline <- readOGR(dsn = "../gis_data", layer = "CVHM_Mesh_outline_modif")
mesh <- readOGR(dsn = "../gis_data", layer = "CVHM_mesh")
PNTS <- mesh_outline@polygons[[1]]@Polygons[[1]]@coords
#PNTS <- PNTS[-dim(PNTS)[1],]
LNS_IJ <- matrix(data = NA, nrow = dim(PNTS)[1], ncol = 2)

for (i in 1:length(mesh)) {
  if (i %% 1000 == 0){
    print(i)
  } 
  if (length(mesh@polygons[[i]]@Polygons) > 1){
    warning(paste("Polygon ", i, "seems more complex"))
  }
  poly_coord <- mesh@polygons[[i]]@Polygons[[1]]@coords
  for (j in 1:(dim(poly_coord)[1]-1)) {
    pa <- poly_coord[j,]
    dst <- sqrt((PNTS[,1] - pa[1])^2 + (PNTS[,2] - pa[2])^2)
    ida <- which(dst < 0.01)
    
    if (length(ida) > 0){
      pb <- poly_coord[j+1,]
      dst <- sqrt((PNTS[,1] - pb[1])^2 + (PNTS[,2] - pb[2])^2)
      idb <- which(dst < 0.01)
      
      if (length(idb) > 0){
        for (kk in 1:length(ida)) {
          for (kkk in 1:length(idb)) {
            idmn <- min(ida[kk], idb[kkk])
            idmx <- max(ida[kk], idb[kkk])
            if (idmx - idmn == 1){
              LNS_IJ[idmn,] <- c(mesh$R[i], mesh$C[i])
            }
          }
        }
      }
    }
  }
}
LNS_IJ <- LNS_IJ[-dim(LNS_IJ)[1],]
PNTS <- PNTS[-dim(PNTS)[1],]
save(PNTS, LNS_IJ, file = "CVHM_LNS_IJ.RData")





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ------ Preprocess Streams -------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cvhm_streams <- readOGR(dsn = "../gis_data/", layer = "CVHM_streams")

# group streams segments per stream 
streamNames <- unique(cvhm_streams@data$NAME)
for (i in 1:length(streamNames)) {
  segids <-  which(cvhm_streams$NAME == streamNames[i])
  # this is a list of nodes that are contained in this river
  uniqNodes <- array(data = NA, dim = c(0,2))
  # This is the ids of the segment nodes
  segNodeIds <- matrix(data = NA, nrow = 0, ncol = 2)
  for (j in 1:length(segids)) {
    for (k in 1:length(cvhm_streams@lines[[segids[j]]]@Lines)) {
      line_seg <- cvhm_streams@lines[[segids[j]]]@Lines[[k]]@coords
      for (ii in 1:(dim(line_seg)[1]-1)) {
        p1 <- line_seg[ii,]
        p2 <- line_seg[ii+1,]
        if (dim(uniqNodes)[1] == 0){
          uniqNodes <- rbind(uniqNodes,c(p1))
          uniqNodes <- rbind(uniqNodes,c(p2))
          segNodeIds <- rbind(segNodeIds, c(1,2))
        }
        else{
          id1 <- which(sqrt((uniqNodes[,1] - p1[1])^2 + (uniqNodes[,2] - p1[2])^2) < 0.1)
          if (length(id1) == 0){
            uniqNodes <- rbind(uniqNodes,c(p1))
            id1 <- dim(uniqNodes)[1]
          }
          id2 <- which(sqrt((uniqNodes[,1] - p2[1])^2 + (uniqNodes[,2] - p2[2])^2) < 0.1)
          if (length(id2) == 0){
            uniqNodes <- rbind(uniqNodes,c(p2))
            id2 <- dim(uniqNodes)[1]
          }
          segNodeIds <- rbind(segNodeIds, c(id1,id2))
        }
      }
      
    }
    
    
  }
  
}
