AverageBudgets <- function(startTime, endTime, cvhm_tm, zones, opt){
  out <- vector(mode = "list", length = 4)
  # load the data
  load(file = as.character(opt$RCH_data))
  load(file = as.character(opt$STRM_data))
  load(file = as.character(opt$WELL_data))

  ndays <- monthDays(cvhm_tm)
  sy <- year(startTime)
  sm <- month(startTime)
  ey <- year(endTime)
  em <- month(endTime)
  inperiod <-  F
  sumRCH <- array(data = 0, dim = dim(RCH[[1]]))
  sumWELL <- array(data = 0, dim = dim(RCH[[1]]))
  sumSTRM <- array(data = 0, dim = dim(RCH[[1]]))
  totdays <- 0
  for (i in 1:length(cvhm_tm)) {

    y <- year(cvhm_tm[i])
    m <- month(cvhm_tm[i])

    if(y == sy & m == sm){
      inperiod <- T
    }

    if(inperiod){
      print(cvhm_tm[i])
      flush.console()
      # Recharge
      sumRCH <- sumRCH + RCH[[i]]*ndays[i]

      # Streams
      ij <- STRM[[i]][,2]*1000 + STRM[[i]][,3]
      tmp <- split(seq_along(ij),ij)
      unique_cells <- as.numeric(names(tmp))
      cell_values <- vector(mode = "numeric", length = length(unique_cells))
      for (j in 1:length(unique_cells)) {
        cell_values[j] <- sum(STRM[[i]][tmp[[j]],4])
      }
      ij <- cbind(floor(unique_cells/1000), unique_cells - floor(unique_cells/1000)*1000 )
      sumSTRM[ij] <- sumSTRM[ij] + cell_values*ndays[i]

      #Multi Node Wells
      ij <- MNW[[i]][,2]*1000 + MNW[[i]][,3]
      tmp <- split(seq_along(ij),ij)
      unique_cells <- as.numeric(names(tmp))
      cell_values <- vector(mode = "numeric", length = length(unique_cells))
      for (j in 1:length(unique_cells)) {
        cell_values[j] <- sum(MNW[[i]][tmp[[j]],4])
      }
      ij <- cbind(floor(unique_cells/1000), unique_cells - floor(unique_cells/1000)*1000 )
      sumWELL[ij] <- sumWELL[ij] + cell_values*ndays[i]

      # Farm wells
      ij <- FRW[[i]][,2]*1000 + FRW[[i]][,3]
      tmp <- split(seq_along(ij),ij)
      unique_cells <- as.numeric(names(tmp))
      cell_values <- vector(mode = "numeric", length = length(unique_cells))
      for (j in 1:length(unique_cells)) {
        cell_values[j] <- sum(FRW[[i]][tmp[[j]],4])
      }
      ij <- cbind(floor(unique_cells/1000), unique_cells - floor(unique_cells/1000)*1000 )
      sumWELL[ij] <- sumWELL[ij] + cell_values*ndays[i]
      totdays <- totdays + ndays[i]
    }

    if (y == ey & m == em){
      inperiod <- F
    }
  }

  # Calculate the budgets and the correct ratio terms per zone
  zone_ids <- unique(array(data = group_zones, dim = c(441*98,1)))
  zone_ids <- zone_ids[-which(is.na(zone_ids))]
  posFlows <- array(data = 0, dim = dim(RCH[[1]]))
  negFlows <- array(data = 0, dim = dim(RCH[[1]]))
  posFlows[sumRCH > 0] <- posFlows[sumRCH > 0] + sumRCH[sumRCH > 0]
  negFlows[sumRCH < 0] <- negFlows[sumRCH < 0] + sumRCH[sumRCH < 0]

  posFlows[sumSTRM > 0] <- posFlows[sumSTRM > 0] + sumSTRM[sumSTRM > 0]
  negFlows[sumSTRM < 0] <- negFlows[sumSTRM < 0] + sumSTRM[sumSTRM < 0]

  posFlows[sumWELL > 0] <- posFlows[sumWELL > 0] + sumWELL[sumWELL > 0]
  negFlows[sumWELL < 0] <- negFlows[sumWELL < 0] + sumWELL[sumWELL < 0]

  #err <- sum(sumRCH) + sum(sumSTRM) + sum(sumWELL)
  err_ratio <- vector(mode = "numeric", length = length(zone_ids))
  for (i in 1:length(zone_ids)) {
    zn_ind <- which(group_zones == zone_ids[i])
    #err <- sum(posFlows[zn_ind]) + sum(negFlows[zn_ind])
    a <- sum(posFlows[zn_ind])
    b <- abs(sum(negFlows[zn_ind]))
    ua <- 1 - ((a-b)*opt$Budget_coef)/a
    ub <- 1 + ((a-b)*(1-opt$Budget_coef))/b
    err_ratio[i] <- a/b

    znRCH <- sumRCH[zn_ind]
    znSTRM <- sumSTRM[zn_ind]
    znWELL <- sumWELL[zn_ind]

    znRCH[znRCH > 0] <- znRCH[znRCH > 0]*ua
    znRCH[znRCH < 0] <- znRCH[znRCH < 0]*ub

    znSTRM[znSTRM > 0] <- znSTRM[znSTRM > 0]*ua
    znSTRM[znSTRM < 0] <- znSTRM[znSTRM < 0]*ub

    znWELL[znWELL > 0] <- znWELL[znWELL > 0]*ua
    znWELL[znWELL < 0] <- znWELL[znWELL < 0]*ub

    sumRCH[zn_ind] <- znRCH
    sumWELL[zn_ind] <- znWELL
    sumSTRM[zn_ind] <- znSTRM
  }
  out[[1]] <- sumRCH
  out[[2]] <- sumSTRM
  out[[3]] <- sumWELL
  out[[4]] <- totdays
  out[[5]] <- err_ratio
  return(out)
}

write2DscatterFile <- function(filename, Dmat, mesh_coords, buff_pnt, opt){
  D <-  cbind(mesh_coords[,1], mesh_coords[,2],Dmat[mesh_coords[,3:4]])

  DBuff <- interpp(x = D[,1], y = D[,2], z = D[,3], linear = F, extrap = T, duplicate = "mean",
                     xo = buff_pnt[,1], yo = buff_pnt[,2])
  Dbufdata <- cbind(DBuff[[1]], DBuff[[2]], DBuff[[3]])

  npsat.input.WriteScattered(filename = filename, PDIM = 2, TYPE = "HOR", MODE = "SIMPLE",
                             DATA = rbind(D, Dbufdata))
}

writeHeadBC <- function(startTime, endTime, cvhm_tm, headTop, PNTS, LNS_IJ, opt){
  ndays <- monthDays(cvhm_tm)
  sy <- year(startTime)
  sm <- month(startTime)
  ey <- year(endTime)
  em <- month(endTime)
  sy_id <- NA
  ey_id <- NA
  for (i in 1:length(cvhm_tm)){
    if(year(cvhm_tm[i]) == sy & month(cvhm_tm[i]) == sm){
      sy_id <- i
    }
    if(year(cvhm_tm[i]) == ey & month(cvhm_tm[i]) == em){
      ey_id <- i
    }
  }

  # Calculate mean and standard deviation of the heads
  head_tmp <- array(data = NA, dim = c(dim(headTop[[1]])[1], dim(headTop[[1]])[2], length(sy_id:ey_id)))
  ind <- 1
  for (i in sy_id:ey_id) {
    head_tmp[,,ind] <- headTop[[i]]
    ind <- ind + 1
  }

  headAV <- apply(head_tmp, c(1,2), mean)
  headStd <- apply(head_tmp, c(1,2), sd)

  # Assign head values to segments
  LNS_HeadInfo <- matrix(data = NA, nrow = dim(LNS_IJ)[1], ncol = 2)
  for (i in 1:dim(LNS_IJ)[1]) {
    LNS_HeadInfo[i,] <- c(headAV[LNS_IJ[i,1],LNS_IJ[i,2]], headStd[LNS_IJ[i,1],LNS_IJ[i,2]])
  }

  # Interpolate the head and standard deviations values on the nodes
  PNTSheadInfo <- matrix(data = NA, nrow = dim(PNTS)[1], ncol = 2)
  for (i in 1:dim(LNS_IJ)[1]) {
    if (i == 1){
      PNTSheadInfo[i,] <- c((LNS_HeadInfo[i, 1] + LNS_HeadInfo[dim(LNS_HeadInfo)[1], 1])/2,
                            (LNS_HeadInfo[i, 2] + LNS_HeadInfo[dim(LNS_HeadInfo)[1], 2])/2)
    }
    else{
      PNTSheadInfo[i,] <- c((LNS_HeadInfo[i, 1] + LNS_HeadInfo[i-1, 1])/2,
                            (LNS_HeadInfo[i, 2] + LNS_HeadInfo[i-1, 2])/2)
    }
  }

  # find the number of boundary functions baes on the head standard deviation
  cnt_bnd <- 0
  next_line <- 1
  BND_LINES <- vector(mode = "list", length = 0)
  bnd_id <- vector(mode = "numeric", length = 0)
  for (i in 1:dim(PNTS)[1]) {
    if (PNTSheadInfo[i,2] < opt$std_Htol){
      if (next_line == 1){
        cnt_bnd <- cnt_bnd + 1
        next_line <- 0
      }
      bnd_id <- c(bnd_id, i)
    }
    else{
      next_line <- 1
      if (length(bnd_id) > 1){
        BND_LINES <- c(BND_LINES, list(bnd_id))
      }
      bnd_id <- vector(mode = "numeric", length = 0)
    }
  }

  if (!dir.exists(paste0(as.character(opt$simFolder),'BC_files'))){
    dir.create(paste0(as.character(opt$simFolder),'BC_files'))
  }

  # write the file using boundary function
  filename <- paste0(opt$simFolder, opt$prefix, '_', opt$timestring, '_BC_h_' ,opt$std_Htol, '.npsat')
  write(length(BND_LINES), file = filename, append = FALSE)
  for (i in 1:length(BND_LINES)) {
    file_bnd <- paste0("bnd_", opt$timestring, '_h', opt$std_Htol, '_',i, '.npsat')
    write(paste0("EDGETOP 0 BC_files/", file_bnd), file = filename, append = T)
    path_file_bnd <- paste0(opt$simFolder, 'BC_files/', file_bnd)
    write("BOUNDARY_LINE", file = path_file_bnd, append = F)
    txt <- paste(length(BND_LINES[[i]]), 1, 1) # Npnts Ndata tolerance
    write(txt, file = path_file_bnd, append = T)
    tmp_data <- cbind(PNTS[BND_LINES[[i]],], PNTSheadInfo[BND_LINES[[i]],1])
    write.table(tmp_data, file = path_file_bnd, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  }

  return(headAV)

}

GenerateWells <- function(wellAv){
  welldata <- readOGR(dsn = "../../CVHM_DATA/WellData",layer = "CVwells_30y_proj")
  farms <- readOGR(dsn = "../gis_data", layer = "FMP")
  farms_poly <- readOGR(dsn = "../gis_data",layer = "FARMS_poly")
  # Find the volume of water that is used for Urban and Ag for each farm
  dwr_id <- unique(farms$dwr_sbrgns)
  dwr_id <- dwr_id[-which(dwr_id == 0)]
  AgUrbPump <- matrix(data = NA, nrow = length(dwr_id), ncol = 2)
  for (i in 1:length(dwr_id)) {
    dwr_cells <- which(farms$dwr_sbrgns == dwr_id[i])
    r <- farms$ROW[dwr_cells]
    c <- farms$COLUMN_[dwr_cells]
    lu <- farms$LU_2000[dwr_cells]
    farmQ <- wellAv[cbind(r,c)]
    AgUrbPump[i,] <- c(sum(farmQ[which(lu != 2)]), sum(farmQ[which(lu == 2)]))
  }
  farms_without_streams <- readOGR(dsn = "../gis_data/", layer = "FARMS_poly_exclSTREAMS")

  # Generate random pumpings
  Ql <- 0.15 # lower limit when sampling from the distribution
  Qu <- 0.925 # upper limit when sampling from the distribution
  WELLS_df <- data.frame(matrix(data = NA, nrow = 50000, ncol = 6))
  colnames(WELLS_df) <- c("X", "Y", "Q", "D", "SL", "Type")
  cnt_wells <- 1

  for (i in 1:length(dwr_id)) {
    print(paste("farm", dwr_id[i]))
    id_farm_pl <- which(farms_poly$dwr_sbrgns==dwr_id[i])

    pl <- plotPolygon(farms_poly@polygons[[id_farm_pl]]@Polygons)
    #makePlot(pl)

    in_sample_wells <- inMultiPolygon(
      mpoly = farms_poly@polygons[[id_farm_pl]]@Polygons,
      x = welldata@coords[,1], y = welldata@coords[,2])
    tmp_well_data <- welldata[in_sample_wells,]
    ipb <- which(tmp_well_data$type == 1)
    iag <- which(tmp_well_data$type == 0)

    tmp_well_data_pb <- tmp_well_data[ipb,]
    tmp_well_data_ag <- tmp_well_data[iag,]

    AG_pdf <- Calc_2D_PDF(tmp_well_data_ag@coords[,1], tmp_well_data_ag@coords[,2], 15)
    PB_pdf <- Calc_2D_PDF(tmp_well_data_pb@coords[,1], tmp_well_data_pb@coords[,2], 15)

    # Pumping distribution (This is different for Ag and Urban wells unless there are
    # very few urban wells in the farm.
    Qag <- tmp_well_data_ag$Q
    Qpb <- tmp_well_data_pb$Q
    Qag <- Qag[-which(is.na(Qag) | Qag <= 0)]
    Qpb <- Qpb[-which(is.na(Qpb) | Qpb <= 0)]
    # In some farms there are very few public supply wells
    if (length(Qpb) < 5){
      Qpb <- c(Qpb, Qag)
    }
    # we have to sample from a log10 distribution
    Qag <- log10(Qag)
    Qpb <- log10(Qpb)

    #qtmp <- vector(mode = "numeric", length = 1000)
    #for(i in 1:1000){
    #  qtmp[i] <- sample(Qag[-which(Qag < 2.6437)],1)
    #}


    # Pumping - depth distribution
    Qtmp <- tmp_well_data$Q
    Dtmp <- tmp_well_data$depth
    id <- which(!is.na(Qtmp) & !is.na(Dtmp))
    Qtmp <- Qtmp[id]
    Dtmp <- Dtmp[id]
    id <- which(Qtmp > 0 & Dtmp > 0)
    Qtmp <- Qtmp[id]
    Dtmp <- Dtmp[id]
    QD_pdf <- Calc_2D_PDF(log10(Qtmp), log10(Dtmp), 15)

    # Depth - screen length
    Dtmp <- tmp_well_data$depth
    SLtmp <- tmp_well_data$bot - tmp_well_data$top
    id <- which(!is.na(Dtmp) & !is.na(SLtmp))
    Dtmp <- Dtmp[id]
    SLtmp <- SLtmp[id]
    id <- which(Dtmp > 0 & SLtmp > 0)
    Dtmp <- Dtmp[id]
    SLtmp <- SLtmp[id]
    DS_pdf <- Calc_2D_PDF(log10(Dtmp), log10(SLtmp), 15)

    Qag_farm <- 0
    Qpb_farm <- 0

    #% The data in Rich well data set are in gpm while the units of CVHM are
    # in m^3/day. To convert the pumping to m^3/day Q = Q*5.450992992.
    # However this correpond to the design pumping. The argument is that
    # each year only six months operate with this rate therefore we divide
    # the pumping by 6 months Q = Q*5.451/6.
    # As we dont want rates lower than 400 m^3/day this would be
    # log10(400*(6/5.451)) = 2.6437 gpm.
    # 10^2.6437*5.451/6

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # AG WELLS

    while (Qag_farm < abs(AgUrbPump[i,1])){
      # generate a valid well location

      BB <- BBoxMultiPoly(farms_without_streams@polygons[[id_farm_pl]]@Polygons)
      xr <- BB$X[1] + (BB$X[2] - BB$X[1])*rand(1)[1,1]
      yr <- BB$Y[1] + (BB$Y[2] - BB$Y[1])*rand(1)[1,1]
      is_in <- inMultiPolygon(farms_without_streams@polygons[[id_farm_pl]]@Polygons, xr, yr)
      if (length(is_in) == 0){next}
      if (cnt_wells > 1){
        # check if the generated point is closer to a previously generated well
        if (min(sqrt((WELLS_df$X - xr)^2 + (WELLS_df$Y - yr)^2), na.rm = T) < 400){next}
      }
      # accept the point with certain probability
      r_accept <- rand(1)[1,1]
      if (r_accept > sample2Dpdf(AG_pdf, xr, yr)){next}


      # assign a random Pumping with minimum 400 m^3/day
      Qr <- sample(Qag[-which(Qag < 2.6437)],1)
      if (abs(AgUrbPump[i,1]) - Qag_farm + (10^Qr)*5.450992992/6 < 400){
        Qr <- log10((AgUrbPump[i,1] - Qag_farm)*6/5.450992992)
      }

      # assign Depth and screen length
      DSr <- AssignQDS_v2(QD_pdf,DS_pdf,Qr)
      if (DSr$V == 0){next}
      WELLS_df[cnt_wells,] <- c(xr, yr, 10^Qr*5.450992992/6, 10^DSr$D, 10^DSr$SL, 0)
      Qag_farm <- Qag_farm + 10^Qr*5.450992992/6
      cnt_wells <- cnt_wells + 1
      pl <- plotwell(pl, xr, yr, 'rgba(255, 182, 193, .9)')
      #if (cnt_wells %% 1000 == 0){
      #  makePlot(pl)
      #}
      print(paste("Nwells: ", cnt_wells-1, " | Farm:",  dwr_id[i], " AG ", Qag_farm, " out of ", abs(AgUrbPump[i,1]) ))
    }

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # URBAN WELLS
    
    while (Qpb_farm < abs(AgUrbPump[i,2])){

      BB <- BBoxMultiPoly(farms_without_streams@polygons[[id_farm_pl]]@Polygons)
      xr <- BB$X[1] + (BB$X[2] - BB$X[1])*rand(1)[1,1]
      yr <- BB$Y[1] + (BB$Y[2] - BB$Y[1])*rand(1)[1,1]
      is_in <- inMultiPolygon(farms_without_streams@polygons[[id_farm_pl]]@Polygons, xr, yr)
      if (length(is_in) == 0){next}
      if (cnt_wells > 0){
        # check if the generated point is closer to a previously generated well
        if (min(sqrt((WELLS_df$X - xr)^2 + (WELLS_df$Y - yr)^2),na.rm = T) < 400){next}
      }
      # accept the point with certain probability
      r_accept <- rand(1)[1,1]
      if (r_accept > sample2Dpdf(PB_pdf, xr, yr)){next}


      # assign a random Pumping with minimum 400 m^3/day
      Qr <- sample(Qpb[-which(Qpb < 2.6437)],1)
      if (abs(AgUrbPump[i,2]) - Qpb_farm + (10^Qr)*5.450992992/6 < 400){
        Qr <- log10((AgUrbPump[i,2] - Qpb_farm)*6/5.450992992)
      }

      # assign Depth and screen length
      DSr <- AssignQDS_v2(QD_pdf,DS_pdf,Qr)
      if (DSr$V == 0){next}
      WELLS_df[cnt_wells,] <- c(xr, yr, 10^Qr*5.450992992/6, 10^DSr$D, 10^DSr$SL, 1)
      Qpb_farm <- Qpb_farm + 10^Qr*5.450992992/6
      cnt_wells <- cnt_wells + 1
      pl <- plotwell(pl, xr, yr, 'rgba(9, 148, 65, .5)')
      #if (cnt_wells %% 1000 == 0){
      #  makePlot(pl)
      #}
      print(paste("Nwells: ", cnt_wells-1, " | Farm:",  dwr_id[i], " PB ", Qpb_farm, " out of ", abs(AgUrbPump[i,2]) ))
    }
    pppp <- T
  }
  return(WELLS_df[-(cnt_wells:50000),])
}

# Return the indices of the points that are in the polygon
inMultiPolygon <- function(mpoly,x,y){
  out <- vector(mode = "numeric", length = length(x))
  for (i in 1:length(mpoly)) {
    out <- out + point.in.polygon(point.x = x,point.y = y,
                     pol.x = mpoly[[i]]@coords[,1], pol.y = mpoly[[i]]@coords[,2])
  }
  return(which(out !=0))
}

BBoxMultiPoly <- function(mpoly){
  x <- c(99999999999999, -99999999999999)
  y <- c(99999999999999, -99999999999999)
  for (i in 1:length(mpoly)) {
    xmn <- min(mpoly[[i]]@coords[,1])
    xmx <- max(mpoly[[i]]@coords[,1])
    ymn <- min(mpoly[[i]]@coords[,2])
    ymx <- max(mpoly[[i]]@coords[,2])
    if (xmn < x[1]){x[1] <- xmn}
    if (xmx > x[2]){x[2] <- xmx}
    if (ymn < y[1]){y[1] <- ymn}
    if (ymx > y[2]){y[2] <- ymx}
  }
  return(list(X = x, Y = y))
}

Calc_2D_PDF <- function(x, y, s){
  xlm <- c(min(x), max(x))
  ylm <- c(min(y), max(y))

  dx <- linspace(xlm[1], xlm[2], 101)
  dy <- linspace(ylm[1], ylm[2], 101)
  dx_n <- linspace(0, 100, 101)
  dy_n <- linspace(0, 100, 101)

  x_n <- interp1(x = dx, y = dx_n, xi = x)
  y_n <- interp1(x = dy, y = dy_n, xi = y)
  XYg <- meshgrid(0:100)

  V <- array(data = 0, dim = dim(XYg$X))
  for (i in 1:length(x_n)) {
    xm <- x_n[i]
    ym <- y_n[i]
    V <- V + phiRBF(nrmRBF(XYg,xm, ym), s)
  }

  V <- V/max(V)
  #DXYg <- meshgrid(dx, dy)
  return(list(X = dx, Y = dy, V = V))
}

nrmRBF <- function(XYgrid, xm, ym){
  return(sqrt((XYgrid$X - xm)^2 + (XYgrid$Y - ym)^2))
}

phiRBF <- function(r,s){
  return(exp(-(r/s)^2))
}

sample2Dpdf <- function(P, x, y){
  nx <- length(P$X)
  if (x < P$X[1]){
    ix <- 1
    iix <- 1
    ux <- 1
  }
  else if(x > P$X[nx]){
    ix <- nx
    iix <- nx
    ux <- 1
  }
  else{
    tx <- interp1(P$X, seq_along(P$X), x)
    ix <- floor(tx)
    ux <- tx - ix
    iix <- min(ix + 1, dim(P$V)[2])
  }

  ny <- length(P$Y)
  if (y < P$Y[1]){
    iy <- 1
    iiy <- 1
    uy <- 1
  }
  else if(y > P$Y[ny]){
    iy <- ny
    iiy <- ny
    uy <- 1
  }
  else{
    ty <- interp1(P$Y, seq_along(P$Y), y)
    iy <- floor(ty)
    uy <- ty - iy
    iiy <- min(iy + 1, dim(P$V)[1])
  }
  if (ix < 1 | ix > dim(P$V)[2]){
    wtf <- T
  }
  if (iix < 1 | iix > dim(P$V)[2]){
    wtf <- T
  }
  if (iy < 1 | iy > dim(P$V)[1]){
    wtf <- T
  }
  if (iiy < 1 | iiy > dim(P$V)[1]){
    wtf <- T
  }
  
  v1 <- P$V[iy,ix]*(1-ux) + P$V[iy,iix]* ux
  v2 <- P$V[iiy,ix]*(1-ux) + P$V[iiy,iix]* ux
  

  v <- v1*(1-uy) + v2*uy

  return(v)
}

AssignQDS_v2 <- function(QD, DS, Qr, minD = 50, minSL = 50, cnt_lim = 500){
  dmin <- max( c(min(QD$Y), min(DS$X),1.699) )
  dmax <- min(max(QD$Y), max(DS$X))
  slmin <- max(min(DS$Y), 1.699)
  slmax <- min(max(DS$Y), 3.1)

  cnt <- 0;
  while(T){
    cnt <- cnt + 1
    Dr <- dmin + (dmax - dmin)*rand(1)[1,1]
    Sr <- slmin + (slmax - slmin)*rand(1)[1,1]
    r1 <- rand(1)[1,1]
    r2 <- rand(1)[1,1]
    R1 <- sample2Dpdf(QD, Qr, Dr)
    R2 <- sample2Dpdf(DS, Dr, Sr)
    if (r1 < R1 & r2 < R2){
      return(list(V=1,D = Dr, SL = Sr))
      break
    }
    if (cnt > cnt_lim){
      return(list(V=0))
      break
    }
  }
}

plotPolygon <- function(poly){
  p <- plot_ly()
  for (i in 1:length(poly)) {
    p <- add_trace(p, x = poly[[i]]@coords[,1], y = poly[[i]]@coords[,2], type = 'scatter', mode = "lines",
                   line = list(color = 'rgb(50, 50, 250)'))
  }
  #p %>%
  #  layout(xaxis = ax, yaxis = ax,showlegend=F)

  return(p)
}

plotwell <- function(p, xw, yw , clr){
  #color <- 'rgba(255, 182, 193, .9)'
  p <- add_trace(p, x = xw, y = yw, type = 'scatter', mode = "markers",
                 marker = list(color = clr) )
  return(p)
}

makePlot <- function(p){
  ax <- list(zeroline = F,
             showline = F,
             showticklabels = F,
             showgrid = F,
             scaleanchor = "x")
  p %>%
    layout(xaxis = ax, yaxis = ax,showlegend=F)

}


















