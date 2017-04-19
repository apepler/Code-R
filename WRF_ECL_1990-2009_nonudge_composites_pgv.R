rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
library(maps)

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->Useful
readMat('~/Documents/GDI/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

cat="p100_rad2cv1"
cat2="rad2_p100_cv1.0"
slp<-gv<-fixes<-list()
comp_slp_d01<-comp_gv_d01<-array(0,c(21,21,2))
comp_gv_d02<-comp_slp_d02<-array(0,c(101,101,2))

tlist=c("ERA-nonudge","ERA-nonudge_notopo")
dirs=c("/srv/ccrc/data34/z3478332/WRF/ERA-nonudge/","/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nonudging_notopo/out/impact/")

end2=c("_2","")

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}
ColorBar2 <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(1, length(brks)-1, subsampleg)])
}

n=1
for(dom in c("d01","d02"))
  for(r in 1:2)
  {
    fixes[[n]]<-read.csv(paste("outputUM_",tlist[r],"/",cat,end2[r],"/",dom,"/ECLfixes_",tlist[r],"_",dom,"_",cat2,"_typing_impactsC2.csv",sep=""))
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    fixes[[n]]$Location2<-0
    #       I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    #       fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    
    I=which(fixes[[n]]$Location==1)
    fixes[[n]]=fixes[[n]][I,]
    
    a=open.nc(paste(dirs[r],"ECLcomposites_",tlist[r],"_",dom,"_",cat2,"_slpgv.nc",sep=""))
    slp[[n]]=var.get.nc(a,"ECL_slp")
    gv[[n]]=var.get.nc(a,"ECL_gv")

    if(dom=="d01")
    {
      comp_slp_d01[,,r]=apply(slp[[n]],c(1,2),mean,na.rm=T)
      comp_gv_d01[,,r]=apply(gv[[n]],c(1,2),mean,na.rm=T)
    } else {
      comp_slp_d02[,,r]=apply(slp[[n]],c(1,2),mean,na.rm=T)
      comp_gv_d02[,,r]=apply(gv[[n]],c(1,2),mean,na.rm=T)
    }
    n=n+1     
  }

pdf(file=paste(figdir,"ECLcomposite_slppv_d02_",cat,"_NoNudgeNoTopo_location.pdf",sep=""),width=7,height=4,pointsize=12)
layout(cbind(1,2))
par(mar=c(3,3,4,1))
contour(comp_slp_d02[,,1],main="a) Control",axes=F,cex.main=2,lwd=2)
contour(comp_gv_d02[,,1],levels=seq(-1,1,0.1),add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(comp_slp_d02[,,2],main="b) NoTopography",axes=F,cex.main=2,lwd=2)
contour( comp_gv_d02[,,2],levels=seq(-1,1,0.1),add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()