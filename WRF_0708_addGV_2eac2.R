######### Intial stuff
rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing')
library(abind)
library("R.matlab")
library(fields)
library(maps)
figdir="~/Documents/ECLs/WRFruns/0708/EACpaper"

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")
c=3 ## Which version do I want? 3 is the default



############# Doing Fei's GV
###
###

## First, make 9-point smoothed running average GV

wrfdirs=c("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN_2eac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_2eac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN_2eac/out/")
dom="d02"
cat="rad2_p100"

wnames=c("R1","R2","R3",
         "R1_notopo","R2_notopo","R3_notopo",
         "R1_BRAN","R2_BRAN","R3_BRAN",
         "R1_BRAN_noeac","R2_BRAN_noeac","R3_BRAN_noeac",
         "R1_BRAN_2eac","R2_BRAN_2eac","R3_BRAN_2eac")

fixes_GV<-events_GV<-list()
dates=seq.POSIXt(from=as.POSIXct("2007010100",format="%Y%m%d%H",tz="GMT"),
                 to=as.POSIXct("2008123118",format="%Y%m%d%H",tz="GMT"),
                 by="6 hours")

dom="d02"
#for(cat in c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5"))
for(w in 13:length(wnames))
{
  tmp=read.csv(paste(wrfdirs[w],"Analysis/GV_6hrly_timeseries.txt",sep=""),header=F)
  GV=rep(NaN,length(tmp[,1]))
  for(i in 5:(length(GV)-4)) GV[i]=mean(tmp[(i-4):(i+4),1])
  
  ## Need a GV thresh
  
  dayGV=aggregate(GV,by=list(day = cut(dates, "days", right = TRUE)),max)
  GVthresh=quantile(dayGV[,2],0.9,na.rm=T)

  wrf=read.csv(paste("ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_v2_typing_impactsC.csv",sep=""),
               stringsAsFactors = F)
  wrf=wrf[,-1]
  wrf$Date2=as.POSIXct(paste(as.character(wrf$Date),substr(wrf$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  I=match(wrf$Date2,dates)
  wrf$GV=GV[I]
  
  wrf_E=read.csv(paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_v2_typing_impactsC.csv",sep=""),
                 stringsAsFactors = F)
  wrf_E=wrf_E[,-1]
  wrf_E$GV=NaN
  
  for(i in 1:length(wrf_E$ID))
  {
    I=which(wrf$ID==wrf_E$ID[i] & wrf$Location==1 & !is.na(wrf$GV))
    if(length(I)>0) wrf_E$GV[i]=max(wrf$GV[I])
  }
  
  wrf_E$GVthresh=as.numeric(wrf_E$GV>=GVthresh)
  
  write.csv(wrf,paste("ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_v2_typing_impactsC_GV.csv",sep=""))
  write.csv(wrf_E,paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_v2_typing_impactsC_GV.csv",sep=""))
}