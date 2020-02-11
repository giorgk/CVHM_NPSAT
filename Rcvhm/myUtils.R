myutils.writeRasterAscii <- function(data, filename, llcoord, cellsize, nodata){
  write(paste("NCOLS", dim(data)[2]), file = filename, append = F)
  write(paste("NROWS", dim(data)[1]), file = filename, append = T)
  write(paste("XLLCORNER", llcoord[1]), file = filename, append = T)
  write(paste("YLLCORNER", llcoord[2]), file = filename, append = T)
  write(paste("CELLSIZE", cellsize), file = filename, append = T)
  write(paste("NODATA_VALUE", nodata), file = filename, append = T)
  write.table(data, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

#' myutils.sub2ind is the equivelant of matlab sub2ind. 
#' Given the row r and column m of a matrix with m rows calculate the linear index
#' 
#' https://stackoverflow.com/questions/4452039/rs-equivalent-to-ind2sub-sub2ind-in-matlab
#' https://cran.r-project.org/doc/contrib/Hiebeler-matlabR.pdf
#'
#' @param r row
#' @param c column
#' @param m number of rows
#'
#' @return index
#' @export
#'
#' @examples
myutils.sub2ind <- function(r,c,m){
  return((c-1)*m + r)
}

myutils.cvhm_sub2ind <- function(r,c,m){
  return((c-1)*m + r)
}