temp <- function(){
  streamNames <- unique(cvhm_streams@data$NAME)
  STREAMS <- vector(mode = "list", length = length(streamNames))
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
    STREAMS[[i]] <- list(NAME = streamNames[i], segids = segids, XY = uniqNodes, ID = segNodeIds)
  }
  
  return(STREAMS)
}