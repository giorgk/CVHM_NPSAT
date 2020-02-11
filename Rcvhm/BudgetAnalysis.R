# Plot the main budget terms 
# Recharge 
RCH <- modflow.gather(BUD,11)
for (i in 1:length(RCH)) {
  RCH[[i]] <- aperm(array(RCH[[i]][,2],c(98,441)),c(2,1))
}


# Wells & Streams
MNW <- modflow.gather(BUD, 9)
FRW <- modflow.gather(BUD, 10)
STRM <- modflow.gather(BUD, 8)
well_array <- array(data = 0, dim = c(dim(RCH[[1]])[1], dim(RCH[[1]])[2], length(RCH)))
strm_array <- array(data = 0, dim = c(dim(RCH[[1]])[1], dim(RCH[[1]])[2], length(RCH)))
for (i in 1:length(FRW)){
  for (j in 1:dim(FRW[[i]])[1]) {
    r <- FRW[[i]][1,2]
    c <- FRW[[i]][1,3]
    well_array[r,c,i] <- well_array[r,c,i] + FRW[[i]][j,4]
  }
  for (j in 1:dim(MNW[[i]])[1]) {
    r <- MNW[[i]][1,2]
    c <- MNW[[i]][1,3]
    well_array[r,c,i] <- well_array[r,c,i] + MNW[[i]][j,4]
  }
  for (j in 1:dim(STRM[[i]])[1]) {
    r <- STRM[[i]][j,2]
    c <- STRM[[i]][j,3]
    strm_array[r,c,i] <- strm_array[r,c,i] + STRM[[i]][j,4]
  }
}

well_monthly <- vector(mode = "numeric", length = length(cvhm_tm))
rch_monthly <- vector(mode = "numeric", length = length(cvhm_tm))
strm_monthly <- vector(mode = "numeric", length = length(cvhm_tm))
ndays <- monthDays(cvhm_tm)

for (i in 1:length(cvhm_tm)) {
  well_monthly[i] <- well_monthly[i] + sum(FRW[[i]][,4])*ndays[i]
  well_monthly[i] <- well_monthly[i] + sum(MNW[[i]][,4])*ndays[i]
  rch_monthly[i] <- rch_monthly[i] + sum(RCH[[i]])*ndays[i]
  strm_monthly[i] <- strm_monthly[i] + sum(STRM[[i]])*ndays[i]
}

well_yearly <- array(data =  well_monthly[7:510], dim = c(12,42))
rch_yearly <- array(data =  rch_monthly[7:510], dim = c(12,42))
strm_yearly <- array(data =  strm_monthly[7:510], dim = c(12,42))


Time <- rep(1962:2003,3)
Type <- c(rep("Well",42), rep("Rch",42), rep("Strm",42))
value <- c(-(colSums(well_yearly)/1.23348e+9), (colSums(rch_yearly)/1.23348e+9), (colSums(strm_yearly)/1.23348e+9))
plot_df <- data.frame(Time,Type,value)

# plot budgets are bar plot
ggplot(plot_df, aes(fill=Type, y=value, x=Time)) + 
  theme_minimal()+
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(name = "MAF", breaks=seq(0,16,2), labels=as.character(seq(0,16,2)))+
  scale_x_continuous(name = "Year", breaks=seq(1960,2000,5), labels=as.character(seq(1960,2000,5)))+
  scale_fill_discrete(labels = c("Recharge","Stream leackage","Pumping"))+
  theme(axis.title.x = element_text(size = 12, hjust = 0.5, vjust = 1.12),
        axis.text = element_text(size = 12),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
  
ggsave(filename = "CVHM_Budget_annual.png", plot = g, device="png")

# plot budgets are step line
tm <- 1962:2003
Well <- -(colSums(well_yearly)/1.23348e+9)
Rch <- (colSums(rch_yearly)/1.23348e+9)
Strm <- (colSums(strm_yearly)/1.23348e+9)
CVHMBUD_MAF <- data.frame(tm,Well,Rch,Strm)

ggplot(CVHMBUD_MAF, aes(x=tm)) + 
  geom_step(aes(y = Well, color = "Pumping"), size = 0.8)+
  geom_step(aes(y = Rch, color = "Recharge"), size = 0.8)+
  geom_step(aes(y = Strm, color = "Stream leackage"), size = 0.8)+
  scale_color_manual("", breaks = c("Pumping", "Recharge","Stream leackage"),  
                        values = c('#619cff','#f8766d','#00ba38'))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 16, hjust = 0.5, vjust = 2.12 , color = "black"),
        legend.position = c(0.20, 0.77),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.5,"line"),
        legend.key=element_blank(),
        legend.justification = "center")+
  labs(x="Year", y="MAF")

save(CVHMBUD_MAF, file = "CVHMBUD_MAF.RData")

# ------- Caclulate recharge rate per year  -------
# first convert it to m^3/month
for (i in 1:length(RCH)) {
  RCH[[i]] <- RCH[[i]]*ndays[i]
}

# Calculate it to m^3/year
RCH_YR <- vector(mode = "list", length = length(1962:2003))
sind <- 7
for (i in 1:length(RCH_YR)) {
  tmp <- array(data = 0, dim = dim(RCH[[1]]))
  eind <- sind + 11
  for (j in sind:eind) {
    tmp <- tmp +RCH[[j]]
  }
  RCH_YR[[i]] <- tmp
  sind <- sind + 12
}

# Now convert it to mm/year
for (i in 1:length(RCH_YR)) {
  RCH_YR[[i]] <- 1000*(RCH_YR[[i]]/(mile*mile))
}

tmp <- hist(RCH_YR[[1]][which(RCH_YR[[1]] != 0)],nclass = 200)
df <- data.frame(tmp$mids,tmp$counts)
colnames(df) <- c("mids","counts")

p <- plot_ly(data = df)
for (i in 1:length(RCH_YR)) {
  non_zero_rch <- RCH_YR[[i]][which(RCH_YR[[i]] != 0)]
  tmp <- hist(non_zero_rch,nclass = 75)
  df <- data.frame(tmp$mids,tmp$counts)
  colnames(df) <- c("mids","counts")
  df$counts <- 100*(tmp$counts /sum(tmp$counts))
  p <- add_trace(p, data = df, x=~mids,y=~counts,
                 type = 'scatter', mode = 'none', fill = 'tonexty')
  #p <- add_trace(p, data = df, x=rep(1961+i,length(tmp$mids)),y=~mids,z=~counts,
  #               type = 'scatter3d', mode = 'lines', fill = 'tonexty')
}

p1 <- p %>%
  layout(title = 'CVHM recharge distributions',
         xaxis = list(title = "Recharge [mm/year]", titlefont = list(size = 18), zeroline = F),
         yaxis = list(title = "Percentage of grid cells [%]", titlefont = list(size = 18)), 
         showlegend = FALSE)
p1
  

plot(tmp$mids,tmp$counts,"s")  