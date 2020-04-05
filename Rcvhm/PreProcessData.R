library(rgdal)
library(plotly)
library(igraph)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ------ Prepare Mesh and buffer nodes ------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load the basic file
bas6_3310 <- readOGR(dsn = "../gis_data/", layer = 'BAS6_3310')
# load the farm layer
fmp_3310 <- readOGR(dsn = "../gis_data/", layer = 'FMP_3310')
# join the databases 
id_active <- apply(bas6_3310@data[, c(5:13)], 1, sum) != 0

bas_active <- bas6_3310
bas_active@data <- cbind(bas6_3310@data, fmp_3310@data)
# isolate the active cells
bas_active <- subset(bas_active, id_active)
# keep only the usefull links
bas_active@data <- bas_active@data[c(2:4,29:35)]
writeOGR(obj = bas_active, dsn = "../gis_data/", layer = 'BAS_active_3310', driver = ogrDrivers()[18,1])


# PREPARE and SAVE 2D interpolation data  ----------
bas6_3310 <- readOGR(dsn = "../gis_data/", layer = 'BAS_active_3310')
cvhm_mesh_buffer <- readOGR(dsn = "../gis_data/", layer = "CVHM_mesh_buffer_3310")
bas_bary_coords <- matrix(data = NA, nrow = length(bas6_3310), ncol = 4)
# make a list with the centroids of the cell
for (i in 1:length(bas6_3310)) {
  if (length(bas6_3310@polygons[[i]]@Polygons) == 1){
    if (dim(bas6_3310@polygons[[i]]@Polygons[[1]]@coords)[1] == 5){
      bas_bary_coords[i,] <-  c(colSums(bas6_3310@polygons[[i]]@Polygons[[1]]@coords[-5,])/4, bas6_3310$ROW[i], bas6_3310$COLUMN_[i])
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
buffer_nodes <- buffer_nodes[-dim(buffer_nodes)[1],]
save(bas_bary_coords, buffer_nodes, file = "chvm2DinterpNodes.RData")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ------ Identify the CVHM boundaries -------------
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
mesh_outline <- readOGR(dsn = "../gis_data", layer = "CVHM_mesh_outline_3310")
mesh <- readOGR(dsn = "../gis_data", layer = "CVHM_mesh_3310")
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
cvhm_streams <- readOGR(dsn = "../gis_data/", layer = "CVHM_streams_3310")

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

# Fields of STREAMS list
# segids : are the ids in the cvhm_streams shapefile that are part of the same river
# XY : are the coordinates of the river nodes
# ID : are the indices of the segment in the XY talbe
# RC : is the row and column of the segment
# CNUM is the column number.
# Therefore the rate of the segment between the river nodes ID[i,1] - ID[i,2] will be assigned by the RC[i,] row column cell.

{
  { # Tests_______________________________
    istream <- 1
    # Plot specific stream
    ax <- list(zeroline = F,
               showline = F,
               showticklabels = F,
               showgrid = F,
               scaleanchor = "x")
    pl <- plot_ly()
    for (i in 1:dim(STREAMS[[istream]]$ID)[1]) {
      pl <- add_trace(pl,x = STREAMS[[istream]]$XY[STREAMS[[istream]]$ID[i,],1], 
                      y = STREAMS[[istream]]$XY[STREAMS[[istream]]$ID[i,],2], 
                      type = "scatter", mode='lines',line = list(color = 'rgb(50, 50, 250)'))
    }
    pl %>%
      layout(xaxis = ax, yaxis = ax,showlegend=F)
  }
  
  # order streams
  g <- graph_from_edgelist(STREAMS[[1]]$ID)
  
  g <- make_empty_graph(n=4) %>%
    add_edges(c(1,2,2,4))
  # Tests______________________________
}

# Find the number of connections per node
# For each stream order the vertices of the main stream and identify and order the tributaries
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
    tributs <- vector(mode = "list", length = length(rem_nodes))
    for (j in 1:length(rem_nodes)) {
      tributary <- rem_nodes[j]
      while (TRUE){
        ilast <- tributary[length(tributary)]
        nb <- as.numeric(neighbors(g,ilast))
        if (length(nb) == 1){
          tributary <- c(tributary, nb)
        }
        else{
          # remove the one that is connected with the previous
          dlt <- c()
          for (ii in 1:length(nb)) {
            if (length(which(tributary == nb[ii])) == 1){
              dlt <- c(dlt, ii)
            }
          }
          nb <- nb[-dlt]
          if (length(nb) == 1){
            tributary <- c(tributary, nb)
          }
          else{
            print(paste("WHAT???  " , i))
          }
        }
        
        if (length(which(STREAMS[[i]]$IDordered == tributary[length(tributary)])) == 1){
          break
        }
      }
      tributs[[j]] <- list(IDordered = tributary)
      
    }
    STREAMS[[i]]$Tributs <- tributs 
  }
  else{
    STREAMS[[i]]$Tributs <- NA
  }
}
# IDordered are the ids of XY ordered along the stream

for (i in 1:length(STREAMS)) {
  IDorderRC <- vector(mode = "numeric", length = length(STREAMS[[i]]$IDordered)-1)
  for (j in seq(1, length(STREAMS[[i]]$IDordered)-1, 1) ) {
    id1 <- STREAMS[[i]]$IDordered[j]
    id2 <- STREAMS[[i]]$IDordered[j+1]
    idrc <- which(STREAMS[[i]]$ID[,1] == id1 & STREAMS[[i]]$ID[,2] == id2)
    if (length(idrc) == 0){
      idrc <- which(STREAMS[[i]]$ID[,1] == id2 & STREAMS[[i]]$ID[,2] == id1)
    }
    if (length(idrc) == 0){
      print(paste("Cant find the segment", id1, "-", id2, "of stream",i))
    }
    else{
      IDorderRC[j] <- idrc
    }
  }
  STREAMS[[i]]$IDorderRC <- IDorderRC
  
  if (length(STREAMS[[i]]$Tributs)){
    if (is.na(STREAMS[[i]]$Tributs[[1]])){
      next
    }
  }
  
  
  for (k in 1:length(STREAMS[[i]]$Tributs)) {
    IDorderRC <- vector(mode = "numeric", length = length(STREAMS[[i]]$Tributs[[k]]$IDordered)-1)
    for (j in seq(1, length(STREAMS[[i]]$Tributs[[k]]$IDordered)-1, 1)) {
      id1 <- STREAMS[[i]]$Tributs[[k]]$IDordered[j]
      id2 <- STREAMS[[i]]$Tributs[[k]]$IDordered[j+1]
      idrc <- which(STREAMS[[i]]$ID[,1] == id1 & STREAMS[[i]]$ID[,2] == id2)
      if (length(idrc) == 0){
        idrc <- which(STREAMS[[i]]$ID[,1] == id2 & STREAMS[[i]]$ID[,2] == id1)
      }
      if (length(idrc) == 0){
        print(paste("Cant find the segment", id1, "-", id2, "of stream",i))
      }
      else{
        IDorderRC[j] <- idrc
      }
    }
    STREAMS[[i]]$Tributs[[k]]$IDorderRC <- IDorderRC
  }
}



# Make a unique list of cells and calculate the total river length per cell
stream_cell_unique <- matrix(data = NA, nrow = 0, ncol = 2)
cell_riv_len <- c()
for (i in 1:length(STREAMS)) {
  for (j in 1:dim(STREAMS[[i]]$RC)[1]) {
    rr <- STREAMS[[i]]$RC[j,1]
    cc <- STREAMS[[i]]$RC[j,2]
    if (dim(stream_cell_unique)[1] == 0){
      stream_cell_unique <- rbind(stream_cell_unique,c(rr,cc))
      cell_riv_len <- c(cell_riv_len,0)
      id <- 1
    }
    else{
      id <- which(stream_cell_unique[,1] == rr & stream_cell_unique[,2] == cc)
      if (length(id) == 0){
        stream_cell_unique <- rbind(stream_cell_unique,c(rr,cc))
        cell_riv_len <- c(cell_riv_len,0)
        id <- dim(stream_cell_unique)[1]
      }
    }
    p1 <- STREAMS[[i]]$ID[j,1]
    p2 <- STREAMS[[i]]$ID[j,2]
    cell_riv_len[id] <- cell_riv_len[id] + sqrt(sum((STREAMS[[i]]$XY[p1,] - STREAMS[[i]]$XY[p2,])^2))
  }
}
save(STREAMS, stream_cell_unique,cell_riv_len, file = "CVHM_STREAM_Preproces.RData")

  
