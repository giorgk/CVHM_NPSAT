library(rgdal)
library(plotly)
library(igraph)

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
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ------ Preprocess Streams -------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cvhm_streams <- readOGR(dsn = "../gis_data/", layer = "CVHM_streams")

# Split streams per name
streamNames <- unique(cvhm_streams@data$NAME)
STREAMS <- vector(mode = "list", length = length(streamNames))
for (i in 1:length(streamNames)) {
  segids <-  which(cvhm_streams$NAME == streamNames[i])
  # this is a list of nodes that are contained in this river
  uniqNodes <- array(data = NA, dim = c(0,2))
  # This is the ids of the segment nodes
  segNodeIds <- matrix(data = NA, nrow = 0, ncol = 2)
  segRC <- matrix(data = NA, nrow = 0, ncol = 2)
  segCell <- matrix(data = NA, nrow = 0, ncol = 1)
  for (j in 1:length(segids)) {
    for (k in 1:length(cvhm_streams@lines[[segids[j]]]@Lines)) {
      line_seg <- cvhm_streams@lines[[segids[j]]]@Lines[[k]]@coords
      for (ii in 1:(dim(line_seg)[1]-1)) {
        p1 <- line_seg[ii,]
        p2 <- line_seg[ii+1,]
        if (dim(uniqNodes)[1] == 0){
          uniqNodes <- rbind(uniqNodes,c(p1))
          uniqNodes <- rbind(uniqNodes,c(p2))
          id1 <- 1
          id2 <- 2
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
        }
        segNodeIds <- rbind(segNodeIds, c(id1,id2))
        segRC <- rbind(segRC, c(cvhm_streams$ROW[segids[j]], cvhm_streams$COLUMN_[segids[j]]))
        segCell <- rbind(segCell, c(cvhm_streams$CELLNUM[segids[j]]))
      }
    }
  }
  STREAMS[[i]] <- list(NAME = streamNames[i], segids = segids, XY = uniqNodes, 
                       ID = segNodeIds, RC = segRC, CNUM = segCell)
}

# Plot specific stream
ax <- list(zeroline = F,
           showline = F,
           showticklabels = F,
           showgrid = F,
           scaleanchor = "x")
istream <- 5
pl <- plot_ly()
for (i in 1:dim(STREAMS[[istream]]$ID)[1]) {
  pl <- add_trace(pl,x = STREAMS[[istream]]$XY[STREAMS[[istream]]$ID[i,],1], 
                  y = STREAMS[[istream]]$XY[STREAMS[[istream]]$ID[i,],2], 
                  type = "scatter", mode='lines',line = list(color = 'rgb(50, 50, 250)'))
}
pl %>%
  layout(xaxis = ax, yaxis = ax,showlegend=F)

# order streams
g <- graph_from_edgelist(STREAMS[[1]]$ID)

g <- make_empty_graph(n=4) %>%
  add_edges(c(1,2,2,4))
# Find the number of connections per node

N_ends <- vector(mode = "numeric", length = length(STREAMS))
for (i in 1:length(STREAMS)) {
  print(i)
  flush.console()
  g <- graph_from_edgelist(STREAMS[[i]]$ID,directed = F)
  n_conn <- vector(mode = "numeric", length = dim(STREAMS[[i]]$XY)[1])
  for (j in 1:dim(STREAMS[[i]]$XY)[1]) {
    #n_conn[j] <- length(as_adj_edge_list(g)[[j]])
    n_conn[j] <- length(neighbors(g,j) )
  }
  STREAMS[[i]]$N_conn <- n_conn
  STREAMS[[i]]$start_end_id <-  which(n_conn == 1)
  N_ends[i] <- length(STREAMS[[i]]$start_end_id)
  
  rem_nodes <- STREAMS[[i]]$start_end_id
  if (N_ends[i] > 2){
    # find longest path
    dst <- 0
    for (ii in 1:(N_ends[i]-1)) {
      for (jj in (ii+1):N_ends[i]) {
        dd <- distances(g, STREAMS[[i]]$start_end_id[ii], STREAMS[[i]]$start_end_id[jj])
        if (dd > dst ){
          istart <- STREAMS[[i]]$start_end_id[ii]
          iend <- STREAMS[[i]]$start_end_id[jj]
          dst <- dd
        }
      }
    }
  }
  else{
    istart <- STREAMS[[i]]$start_end_id[1]
    iend <- STREAMS[[i]]$start_end_id[2]
  }
  rem_nodes <- rem_nodes[-which(rem_nodes == istart)]
  rem_nodes <- rem_nodes[-which(rem_nodes == iend)]
  STREAMS[[i]]$IDordered <- as.numeric( shortest_paths(g,istart)$vpath[[iend]])
  if (length(rem_nodes) > 0){
    # find the tributaries
  }
}


  
