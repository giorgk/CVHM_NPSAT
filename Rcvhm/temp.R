temp <- function(){
  cvhm_streams <- readOGR(dsn = "../gis_data/", layer = "CVHM_streams")
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
  
  return(STREAMS)
}

temp1 <- function(STREAMS){
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
              print(paste("WHAT???"))
            }
          }
          
          if (length(which(STREAMS[[i]]$IDordered == tributary[length(tributary)])) == 1){
            break
          }
        }
        tributs[[j]] <- tributary
        
      }
      STREAMS[[i]]$Tributs <- tributs
    }
  }
}

temp2 <- function(STREAMS){
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
}

temp3 <- function(filename, nd, msh){
  write(c(dim(nd)[1], dim(msh)[1]), file = filename, append = FALSE)
  utils::write.table(nd, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  utils::write.table(msh, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}


temp4 <- function(filename, maxChar = 1000){
  allLines <- readLines(filename)
  # read the number of polygons
  ind <- 1
  tmp <- strsplit(substr(allLines[ind], 1, maxChar)[[1]], split = " ")[[1]]
  tmp <- tmp[which(tmp != "")]
  Npoly <- as.numeric(tmp)
  streamCoords <- vector(mode = "list", length = Npoly)
  Area <- vector(mode = "numeric", length = Npoly)
  Rate <- vector(mode = "numeric", length = Npoly)
  Q <- vector(mode = "numeric", length = Npoly)
  for (i in 1:Npoly) {
    ind <- ind + 1
    tmp <- strsplit(substr(allLines[ind], 1, maxChar)[[1]], split = " ")[[1]]
    tmp <- tmp[which(tmp != "")]
    tmp <- as.numeric(tmp)
    Npnts <- tmp[1]
    Rate[i] <- tmp[2]
    coords <- matrix(data = NA, nrow = Npnts, ncol = 2)
    for (j in 1:Npnts) {
      ind <- ind + 1
      tmp <- strsplit(substr(allLines[ind], 1, maxChar)[[1]], split = " ")[[1]]
      tmp <- tmp[which(tmp != "")]
      tmp <- as.numeric(tmp)
      coords[j,] <- tmp
    }
    Area[i] <- abs(pracma::polyarea(coords[,1], coords[,2]))
    Q[i] <- Area[i]*Rate[i]
    streamCoords[[i]] <- coords 
  }
  return(list(Coords = streamCoords, Rate = Rate, Area = Area, Q = Q))
}

qqq <- vector(mode = "numeric", length = length(STREAMS))
n_unique <- vector(mode = "numeric", length = length(STREAMS))
ccc <- c()
for (i in 1:length(STREAMS)) {
  ucels <- unique(gwtools::sub2ind(STREAMS[[i]]$RC[,1],STREAMS[[i]]$RC[,2],441))
  ccc <- c(ccc,ucels)
  n_unique[i] <- length(ucels)
  qqq[i] <- sum(sumBUD$sumSTRM[ucels])/sumBUD$totdays
}

sum(sumBUD$sumSTRM[gwtools::sub2ind(stream_cell_unique[,1],stream_cell_unique[,2],441)])/sumBUD$totdays

temp5 <- function(filename, PDIM = 2){
  # Read header
  tmp <- readLines(filename, n = 4)
  N <- as.numeric(strsplit(tmp[4], " ")[[1]])
  if (PDIM == 1){
    cnames <- c("X")
  }
  else if (PDIM == 2){
    cnames <- c("X", "Y", "V")
  }
  if(N[2] == 3){
    cnames <- c(cnames, "V")
  }
  else if (N[2] > 3){
    cnames <- c(cnames, "V1")
    for (i in seq(1, (N[2]-1)/2, 1)) {
      cnames <- c(cnames, paste0("Z", i), paste0("V", i+1))
    }
    
  }
  DATA <- utils::read.table(file = filename,
                            header = FALSE, sep = "", skip = 4, nrows = N[1],
                            quote = "",fill = TRUE,
                            col.names = cnames)
  
  return(DATA)
}
