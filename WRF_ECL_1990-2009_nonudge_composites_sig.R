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
events<-fixes<-list()
comp_rain_d01<-comp_wind_d01<-array(0,c(21,21,4,2))
comp_rain_d02<-comp_wind_d02<-array(0,c(101,101,4,2))

ECLrain<-ECLwind<-list()

tlist=c("ERA-nonudge","ERA-nonudge_notopo","ERA-nudge","ERA-nudge_notopo")
dirs=c("/srv/ccrc/data34/z3478332/WRF/ERA-nonudge/","/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nonudging_notopo/out/impact/",
       "/srv/ccrc/data34/z3478332/WRF/ERA-nudge/","/srv/ccrc/data34/z3478332/WRF/output/ERAI_R2_nudging_notopo_19902009/out/impact/")

end2=c("_2","","_2","")

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
      events[[n]]<-read.csv(paste("outputUM_",tlist[r],"/",cat,end2[r],"/",dom,"/ECLevents_",tlist[r],"_",dom,"_",cat2,"_typing_impactsC2.csv",sep=""))
      events[[n]]$Year=floor(events[[n]]$Date1/10000)
      events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
      events[[n]][events[[n]]==-Inf]=NaN
      
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
      
      ### Rain stuff
      
      a=open.nc(paste(dirs[r],"ECLrain_",tlist[r],"_",dom,"_",cat2,"_centred",end2[r],".nc",sep=""))
      ECLrain[[n]]<-var.get.nc(a,"ECLrain")
      a=open.nc(paste(dirs[r],"ECLwind_",tlist[r],"_",dom,"_",cat2,end2[r],".nc",sep=""))
      ECLwind[[n]]<-var.get.nc(a,"ECL_WS10")
      
      if(dom=="d01")
      {
        I=which(fixes[[n]]$Location==1)
        comp_rain_d01[,,r,1]=apply(ECLrain[[n]][,,I],c(1,2),mean,na.rm=T)
        comp_wind_d01[,,r,1]=apply(ECLwind[[n]][,,I],c(1,2),mean,na.rm=T)
        
        I=which(fixes[[n]]$Location2==1)
        comp_rain_d01[,,r,2]=apply(ECLrain[[n]][,,I],c(1,2),mean,na.rm=T)
        comp_wind_d01[,,r,2]=apply(ECLwind[[n]][,,I],c(1,2),mean,na.rm=T)
      } else if(r<=2) {
        I=which(fixes[[n]]$Location==1)
        comp_rain_d02[,,r,1]=apply(ECLrain[[n]][,,I],c(1,2),mean,na.rm=T)
        comp_wind_d02[,,r,1]=apply(ECLwind[[n]][,,I],c(1,2),mean,na.rm=T)
        
        I=which(fixes[[n]]$Location2==1)
        comp_rain_d02[,,r,2]=apply(ECLrain[[n]][,,I],c(1,2),mean,na.rm=T)
        comp_wind_d02[,,r,2]=apply(ECLwind[[n]][,,I],c(1,2),mean,na.rm=T)      
        }
      n=n+1     
    }

pal <- color.palette(c("white","cyan","blue","black"), c(20,20,20))
pdf(file=paste(figdir,"ECLwind_composite_d02_NoNudgevNoTopo_Location2.pdf",sep=""),width=7.5,height=4,pointsize=12)
layout(cbind(1,2,3),width=c(1,1,0.3))
par(mar=c(3,3,4,1))
for(r in 1:2)
{
image(comp_wind_d02[,,r,2],breaks=c(seq(0,14,1),100),col=pal(15),main=tlist[r],axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
}
ColorBar(c(seq(0,14,1),100),pal(15),subsampleg=2)
dev.off()

pdf(file=paste(figdir,"ECLwind_composite_NoNudgevNoTopo_Location2a_change_vres_sig2.pdf",sep=""),width=12,height=4,pointsize=12)
layout(cbind(1,2,3,4,5),width=c(1,0.3,1,1,0.3))
pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))
pal2 <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))

par(mar=c(3,3,4,1), mgp = c(3, 1, 0),cex=1)
image(comp_wind_d02[,,1,2],breaks=c(seq(0,14,1),100),col=pal1(15),main="Average windfall",axes=F)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))

ColorBar(c(seq(0,14,1),100),pal1(15),subsampleg=2)

sig_d01=matrix(0,21,21)
for(i in 1:21)
  for(j in 1:21)
  {
    a=t.test(ECLwind[[1]][i,j,fixes[[1]]$Location2==1],ECLwind[[2]][i,j,fixes[[2]]$Location2==1])
    sig_d01[i,j]=a$p.value
  }

par(mar=c(3,3,4,1), mgp = c(3, 1, 0),cex=1)
image(comp_wind_d01[,,2,2]-comp_wind_d01[,,1,2],breaks=c(-5,seq(-3,3.5,0.5)),col=pal2(14),main="50km resolution",axes=F)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
# contour(sig_d01<=0.05,add=T,drawlabels=F)

lon<-lat<-seq(0,1,length.out=21)
sigmask=which(sig_d01<=0.05,arr.ind=T)
sigmask2=cbind(lon[sigmask[,1]],lat[sigmask[,2]])
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.2)


sig_d02=matrix(0,101,101)
for(i in 1:101)
  for(j in 1:101)
  {
    a=t.test(ECLwind[[3]][i,j,fixes[[3]]$Location2==1],ECLwind[[4]][i,j,fixes[[4]]$Location2==1])
    sig_d02[i,j]=a$p.value
  }

image(comp_wind_d02[,,2,2]-comp_wind_d02[,,1,2],breaks=c(-5,seq(-3,3.5,0.5)),col=pal2(14),main="10km resolution",axes=F)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
# contour(sig_d02<=0.05,add=T,drawlabels=F)
lon<-lat<-seq(0,1,length.out=101)
sigmask=which(sig_d02<=0.05,arr.ind=T)
sigmask2=cbind(lon[sigmask[,1]],lat[sigmask[,2]])

I=which(sigmask2[,1]%%0.02==0 & sigmask2[,2]%%0.02==0)
points(lon[sigmask[I,1]],lat[sigmask[I,2]],col="black",pch=19,cex=0.2)

ColorBar(seq(-3.5,3.5,0.5),pal2(14),subsampleg=2)
dev.off()


pdf(file=paste(figdir,"ECLrain_composite_NoNudgeNoTopo_Location2a_vres.pdf",sep=""),width=8,height=4,pointsize=12)
layout(cbind(1,2,3),width=c(1,1,0.3))
pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))
par(mar=c(3,3,4,1), mgp = c(3, 1, 0),cex=1)
image(comp_rain_d01[,,2,2],breaks=c(seq(0,14,1),100),col=pal1(15),main="50 km resolution",axes=F)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
par(mar=c(3,3,4,1), mgp = c(3, 1, 0),cex=1)
image(comp_rain_d02[,,2,2],breaks=c(seq(0,14,1),100),col=pal1(15),main="10 km resolution",axes=F)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
ColorBar(c(seq(0,14,1),100),pal1(15),subsampleg=2)
dev.off()

pdf(file=paste(figdir,"ECLwind_composite_NoNudgevNoTopo_Location2a_change_vres.pdf",sep=""),width=12,height=4,pointsize=12)
layout(cbind(1,2,3,4,5),width=c(1,0.3,1,1,0.3))
pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))
pal2 <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))

par(mar=c(3,3,4,1), mgp = c(3, 1, 0),cex=1)
image(comp_wind_d02[,,1,2],breaks=c(seq(0,14,1),100),col=pal1(15),main="Average wind speed (m/s)",axes=F)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))

ColorBar(c(seq(0,14,1),100),pal1(15),subsampleg=2)

par(mar=c(3,3,4,1), mgp = c(3, 1, 0),cex=1)
image(comp_wind_d01[,,2,2]-comp_wind_d01[,,1,2],breaks=c(-5,seq(-3,3.5,0.5)),col=pal2(14),main="50km resolution",axes=F)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
image(comp_wind_d02[,,2,2]-comp_wind_d02[,,1,2],breaks=c(-5,seq(-3,3.5,0.5)),col=pal2(14),main="10km resolution",axes=F)
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))

ColorBar(seq(-3.5,3.5,0.5),pal2(14),subsampleg=2)
dev.off()




