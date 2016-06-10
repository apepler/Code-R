rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')
figdir="~/Documents/ECLs/WRFruns/0708/EACpaper"

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->Useful
readMat('~/Documents/GDI/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")
c=3
dir='/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing/'
events<-fixes<-comp<-list()
comp_rain<-comp_wind<-array(0,c(101,101,3,5))
tlist=c("","_notopo","_BRAN","_BRAN_noeac","_BRAN_2eac")
dir1=c(36,36,37,37,45)

n=1
dom="d02"
for(type in 1:5)
    for(r in 1:3)
    {
      if(type<5)
      {
      events[[n]]=read.csv(paste(dir,"ECLevents_",dom,"_0708_R",r,tlist[type],"_",cat[c],"_typing_impacts.csv",sep=""))
      fixes[[n]]=read.csv(paste(dir,"ECLfixes_",dom,"_0708_R",r,tlist[type],"_",cat[c],"_typing_impacts.csv",sep=""))
      } else {
        events[[n]]=read.csv(paste(dir,"ECLevents_",dom,"_0708_R",r,tlist[type],"_",cat[c],"_v2_typing_impacts.csv",sep=""))
        fixes[[n]]=read.csv(paste(dir,"ECLfixes_",dom,"_0708_R",r,tlist[type],"_",cat[c],"_v2_typing_impacts.csv",sep=""))
      }
      
      events[[n]]$Year=floor(events[[n]]$Date1/10000)
      events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
      fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
      fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
      
      fixes[[n]]$Location2<-0
      I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
      fixes[[n]]$Location2[I]<-1
      I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
      fixes[[n]]$Location2[I]<-1
      I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
      fixes[[n]]$Location2[I]<-1
      
      ### Rain stuff
      
      if(type<5) a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLrain_d02_0708_",cat[c],"_centred.nc",sep="")) else
        a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLrain_d02_0708_",cat[c],"_v2_centred.nc",sep=""))
      tmp=var.get.nc(a,"ECLrain")
      I=which(fixes[[n]]$Location2==1)
      fixes[[n]]$MeanRainS=apply(tmp[,1:50,],3,mean,na.rm=T)
      fixes[[n]]$MeanRainN=apply(tmp[,52:101,],3,mean,na.rm=T)
      comp_rain[,,r,type]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
      
      if(type<5) a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLwind_d02_0708_",cat[c],".nc",sep="")) else
        a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLwind_d02_0708_",cat[c],"_v2.nc",sep=""))
      tmp=var.get.nc(a,"ECL_WS10")
      I=which(fixes[[n]]$Location2==1)
      comp_wind[,,r,type]=apply(tmp[,,I],c(1,2),mean,na.rm=T)      
      n=n+1     
    }

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
#bb2=c(-100,-4,-3,-2,-1,0,1,2,3,4,100)
bb2=c(-100,-5,-2,-1,0,1,2,5,100)

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = brks[seq(1, length(brks), subsampleg)])
}

pal <- color.palette(c("white","cyan","blue","black"), c(20,20,20))
pdf(file=paste(figdir,"/ECL_rainC_0708_",dom,"_",cat[c],"_SST_location2.pdf",sep=""),width=11.5,height=4,pointsize=12)
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))
par(mar=c(3,3,4,1))
image(apply(comp_rain[,,,4],c(1,2),mean),breaks=c(seq(0,14,1),100),col=pal(15),main="NoEAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
image(apply(comp_rain[,,,3],c(1,2),mean),breaks=c(seq(0,14,1),100),col=pal(15),main="BRAN",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
image(apply(comp_rain[,,,5],c(1,2),mean),breaks=c(seq(0,14,1),100),col=pal(15),main="2EAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
ColorBar(c(seq(0,14,1),100),pal(15),subsampleg=2)
dev.off()

pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
#bb2=c(-100,-4,-3,-2,-1,0,1,2,3,4,100)
bb2=c(-100,-5,-2,-1,0,1,2,5,100)
pdf(file=paste(figdir,"/ECL_rainC_0708_",dom,"_",cat[c],"_SST_location2_change.pdf",sep=""),width=8,height=4,pointsize=12)
layout(cbind(1,2,3),width=c(1,1,0.3))
par(mar=c(3,3,4,1))
image(apply(comp_rain[,,,4]-comp_rain[,,,3],c(1,2),mean),breaks=seq(-5,5,1),col=pal(10),main="NoEAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
image(apply(comp_rain[,,,5]-comp_rain[,,,3],c(1,2),mean),breaks=seq(-5,5,1),col=pal(10),main="2EAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
ColorBar(seq(-5,5),pal(10),subsampleg=1)
dev.off()


pal <- color.palette(c("white","cyan","blue","black"), c(20,20,20))
pdf(file=paste(figdir,"/ECL_wind_0708_",dom,"_",cat[c],"_SST_location.pdf",sep=""),width=11.5,height=4,pointsize=12)
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))
par(mar=c(3,3,4,1))
image(apply(comp_wind[,,,4],c(1,2),mean),breaks=c(seq(0,14,1),100),col=pal(15),main="NoEAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
image(apply(comp_wind[,,,3],c(1,2),mean),breaks=c(seq(0,14,1),100),col=pal(15),main="BRAN",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
image(apply(comp_wind[,,,5],c(1,2),mean),breaks=c(seq(0,14,1),100),col=pal(15),main="2EAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
ColorBar(c(seq(0,14,1),100),pal(15),subsampleg=2)
dev.off()

pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
#bb2=c(-100,-4,-3,-2,-1,0,1,2,3,4,100)
bb2=c(-100,-5,-2,-1,0,1,2,5,100)
pdf(file=paste(figdir,"/ECL_wind_0708_",dom,"_",cat[c],"_SST_location_change.pdf",sep=""),width=8,height=4,pointsize=12)
layout(cbind(1,2,3),width=c(1,1,0.3))
par(mar=c(3,3,4,1))
image(apply(comp_wind[,,,4]-comp_wind[,,,3],c(1,2),mean),breaks=seq(-2.5,2.5,0.5),col=pal(10),main="NoEAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
image(apply(comp_wind[,,,5]-comp_wind[,,,3],c(1,2),mean),breaks=seq(-2.5,2.5,0.5),col=pal(10),main="2EAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
ColorBar(seq(-2.5,2.5,0.5),pal(10),subsampleg=1)
dev.off()


pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
#bb2=c(-100,-4,-3,-2,-1,0,1,2,3,4,100)
bb2=c(-100,-5,-2,-1,0,1,2,5,100)
pdf(file=paste(figdir,"/ECL_0708_",dom,"_",cat[c],"_notopo_location2_change.pdf",sep=""),width=9,height=4,pointsize=12)
layout(cbind(1,2,3,4),width=c(1,0.3,1,0.3))
par(mar=c(3,3,4,1))
image(apply(comp_rain[,,,2]-comp_rain[,,,1],c(1,2),mean),breaks=seq(-5,5,1),col=pal(10),main="Rain",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
ColorBar(seq(-5,-5,1),pal(10),subsampleg=1)
image(apply(comp_wind[,,,2]-comp_wind[,,,1],c(1,2),mean),breaks=seq(-2.5,2.5,0.5),col=pal(10),main="Wind",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
ColorBar(seq(-2.5,2.5,0.5),pal(10),subsampleg=1)
dev.off()


################ SLP/GV composites :)

rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')
figdir="~/Documents/ECLs/WRFruns/0708/EACpaper"

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->Useful
readMat('~/Documents/GDI/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")
c=3
dir='/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing/'
events<-fixes<-comp<-list()
comp_slp<-comp_gv<-array(0,c(21,21,3,5))
tlist=c("","_notopo","_BRAN","_BRAN_noeac","_BRAN_2eac")
dir1=c(36,36,37,37,45)

n=1
dom="d01"
for(type in 1:5)
  for(r in 1:3)
  {
    if(type<5)
    {
      events[[n]]=read.csv(paste(dir,"ECLevents_",dom,"_0708_R",r,tlist[type],"_",cat[c],"_typing_impacts.csv",sep=""))
      fixes[[n]]=read.csv(paste(dir,"ECLfixes_",dom,"_0708_R",r,tlist[type],"_",cat[c],"_typing_impacts.csv",sep=""))
    } else {
      events[[n]]=read.csv(paste(dir,"ECLevents_",dom,"_0708_R",r,tlist[type],"_",cat[c],"_v2_typing_impacts.csv",sep=""))
      fixes[[n]]=read.csv(paste(dir,"ECLfixes_",dom,"_0708_R",r,tlist[type],"_",cat[c],"_v2_typing_impacts.csv",sep=""))
    }
    
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    
    fixes[[n]]$Location2<-0
    I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    
    ### Rain stuff
    
    if(type<5) a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLslp_0708_",cat[c],".nc",sep="")) else
      a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLslp_0708_",cat[c],"_v2.nc",sep=""))
    tmp=var.get.nc(a,"ECL_slp")
    I=which(fixes[[n]]$Location==1)
    comp_slp[,,r,type]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
    tmp=var.get.nc(a,"ECL_gv")
    comp_gv[,,r,type]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
 
    n=n+1     
  }

pdf(file=paste(figdir,"/ECL_slp_0708_",dom,"_",cat[c],"_SST_location.pdf",sep=""),width=10.5,height=4,pointsize=12)
layout(cbind(1,2,3))
par(mar=c(3,3,4,1))
contour(apply(comp_slp[,,,4],c(1,2),mean),main="NoEAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(comp_slp[,,,3],c(1,2),mean),main="BRAN",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(comp_slp[,,,5],c(1,2),mean),main="2EAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()

pdf(file=paste(figdir,"/ECL_gv_0708_",dom,"_",cat[c],"_SST_location.pdf",sep=""),width=10.5,height=4,pointsize=12)
layout(cbind(1,2,3))
par(mar=c(3,3,4,1))
contour(apply(comp_gv[,,,4],c(1,2),mean),main="NoEAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(comp_gv[,,,3],c(1,2),mean),main="BRAN",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(comp_gv[,,,5],c(1,2),mean),main="2EAC",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()


