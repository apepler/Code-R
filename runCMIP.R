rm(list=ls())
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
source('~/Documents/R/CMIP5_play.R')
for(i in 1:37) CMIP_wind(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N") 

f1=open.nc('/srv/ccrc/data23/z3478332/NCEP/Monthly/uwnd.mon.mean.nc')
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")
timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
##If using southern region
#NCEP=apply(var.get.nc(f1,'uwnd',c(61,49,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)
##If using northern region
NCEP=apply(var.get.nc(f1,'uwnd',c(61,47,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)

Ucomp=matrix(0,12,39)
for(i in 1:12)
{
  Ucomp[i,1]=i
  I=which(timeN[,1]<=2005 & timeN[,2]==i)
  Ucomp[i,2]=mean(NCEP[I])
}

coolC<-warmC<-array(0,c(144,73,38))
for(i in 1:37)
{
a=CMIP_data(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N") 
#Ucomp[,i+2]=unlist(a[1])
coolC[,,i]=unlist(a[2])
warmC[,,i]=unlist(a[3])
}



##Mean, SD, quartiles
coolC2<-warmC2<-array(0,c(144,73,5))
coolC2[,,1]=apply(coolC,c(1,2),mean,na.rm=T)
coolC2[,,2]=apply(coolC,c(1,2),sd,na.rm=T)
coolC2[,,3]=apply(coolC,c(1,2),quantile,0.25,na.rm=T)
coolC2[,,4]=apply(coolC,c(1,2),quantile,0.5,na.rm=T)
coolC2[,,5]=apply(coolC,c(1,2),quantile,0.7,na.rm=T)
warmC2[,,1]=apply(warmC,c(1,2),mean,na.rm=T)
warmC2[,,2]=apply(warmC,c(1,2),sd,na.rm=T)
warmC2[,,3]=apply(warmC,c(1,2),quantile,0.25,na.rm=T)
warmC2[,,4]=apply(warmC,c(1,2),quantile,0.5,na.rm=T)
warmC2[,,5]=apply(warmC,c(1,2),quantile,0.7,na.rm=T)
Ucomp2=cbind(Ucomp[,1:2],matrix(0,12,5))
Ucomp2[,3]=apply(Ucomp[,3:39],1,mean,na.rm=T)
Ucomp2[,4]=apply(Ucomp[,3:39],1,sd,na.rm=T)
Ucomp2[,5]=apply(Ucomp[,3:39],1,quantile,0.25,na.rm=T)
Ucomp2[,6]=apply(Ucomp[,3:39],1,quantile,0.5,na.rm=T)
Ucomp2[,7]=apply(Ucomp[,3:39],1,quantile,0.75,na.rm=T)

##Plotting for Uwind - two versions - quartiles & std

library(hydroGOF)
tiff(file="~/Documents/Data/CMIP5/Uwind_meansd.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp2[,3],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(0,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(0,10),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","CMIP5 Average"),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
plotbandsonly(Ucomp2[,3]-Ucomp2[,4],Ucomp2[,3]+Ucomp2[,4], bands.col=rgb(0.9,0.9,0.9),ylim=c(0,10),alpha=0.5)
lines(seq(1,12),Ucomp2[,2],col="blue",lwd=2)
lines(seq(1,12),Ucomp2[,3],col="red",lwd=2)
dev.off()
tiff(file="~/Documents/Data/CMIP5/Uwind_quart.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp2[,6],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(0,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(0,10),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","CMIP5 Median"),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
plotbandsonly(Ucomp2[,5],Ucomp2[,7], bands.col=rgb(0.9,0.9,0.9),ylim=c(0,10),alpha=0.5)
lines(seq(1,12),Ucomp2[,2],col="blue",lwd=2)
lines(seq(1,12),Ucomp2[,6],col="red",lwd=2)
dev.off()

##Plotting for Corrs - Lots of plots - mean, q25,median,q75
coolC2=coolC2[,,-2]
warmC2=warmC2[,,-2]
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(14)
cm[7:8]="white"
coolC2[coolC2>=0.69]=0.69
coolC2[coolC2<=(-0.69)]=-0.69
coolC2[abs(coolC2)<=0.27]=NaN
warmC2[warmC2>=0.69]=0.69
warmC2[warmC2<=(-0.69)]=-0.69
warmC2[abs(warmC2)<=0.27]=NaN
names=c("mean","Q25","Q50","Q75")
Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",list[1,1],"_historical_r1i1p1_",list[1,2],"-",list[1,3],"_2.5deg.nc",sep="")
f1=open.nc(Rfile)
lat=var.get.nc(f1,"lat")
lon=var.get.nc(f1,"lon")

for(i in 1:4)
{
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_CMIP5_",names[i],"_coola.tiff",sep=""), height=500, width=800)
  image.plot(lon,lat,coolC2[,,i],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_CMIP5_",names[i],"_warma.tiff",sep=""), height=500, width=800)
  image.plot(lon,lat,warmC2[,,i],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}


###############AWAP COMPARISON############

rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")

f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")
timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(144,73,1,length(timeN[,1])),unpack=T)
lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
Uind1=apply(uwnd[61,47:51,],2,mean) ##25-35
Uind2=apply(uwnd[61,49:53,],2,mean) ##30-40
years=cbind(seq(1950,2005),matrix(0,56,2))
for(i in 1:length(years[,1]))
{
  I=which(timeN[,1]==years[i,1] & timeN[,2]>=5 & timeN[,2]<=10)
  years[i,2]=mean(Uind1[I])
  I=which((timeN[,1]==years[i,1] & timeN[,2]<=4) | (timeN[,1]==years[i,1]-1 & timeN[,2]>=11))
  years[i,3]=mean(Uind1[I])
}
load('~/Documents/Data/AWAPrain.RData')
yearsA=seq(1900,2012)
I=which(yearsA>=1950 & yearsA<=2005)
cool=cool[,,I]
yearsB=seq(1901,2012)
I=which(yearsB>=1950 & yearsB<=2005)
warm=warm[,,I] #So has the year ending in 1950 - year ending in 2005

library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

library(fields)
cool1<-warm1<-cool2<-warm2<-array(0,c(144,73,56))
for(i in 1:56)
{
   a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(cool[,,i])),list(x=lon,y=lat))
   cool1[,,i]=a$z
   a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(warm[,,i])),list(x=lon,y=lat))
   warm1[,,i]=a$z
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(cool[,,i]*Useful$mask)),list(x=lon,y=lat))
  cool2[,,i]=a$z
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(warm[,,i]*Useful$mask)),list(x=lon,y=lat))
  warm2[,,i]=a$z
}

coolC<-warmC<-matrix(0,144,73)
for(i in 1:144)
  for(j in 1:73)
  {
    coolC[i,j]=cor(cool1[i,j,],years[,2])
    warmC[i,j]=cor(warm1[i,j,],years[,3])
  }

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(14)
cm[7:8]="white"
coolC[coolC>=0.69]=0.69
coolC[coolC<=(-0.69)]=-0.69
warmC[warmC>=0.69]=0.69
warmC[warmC<=(-0.69)]=-0.69

tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_AWAPregrid_cool.tiff", height=500, width=600,pointsize=20)
image.plot(lon,lat[73:1],coolC[,73:1],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_AWAPregrid_warm.tiff", height=500, width=600,pointsize=20)
image.plot(lon,lat[73:1],warmC[,73:1],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

coolC[abs(coolC)<=0.27]=NaN
warmC[abs(warmC)<=0.27]=NaN
tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_sig_AWAPregrid_cool.tiff", height=500, width=600,pointsize=20)
image.plot(lon,lat[73:1],coolC[,73:1],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_sig_AWAPregrid_warm.tiff", height=500, width=600,pointsize=20)
image.plot(lon,lat[73:1],warmC[,73:1],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

coolC=apply(cool1,c(1,2),mean,na.rm=T)
warmC=apply(warm1,c(1,2),mean,na.rm=T)
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
tiff(file="~/Documents/Data/CMIP5/Rain_AWAPregrid_warma.tiff", height=500, width=800)
par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
image(lon,lat[73:1],warmC[,73:1],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
dev.off()
tiff(file="~/Documents/Data/CMIP5/Rain_AWAPregrid_coola.tiff", height=500, width=800)
par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
image(lon,lat[73:1],coolC[,73:1],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
dev.off()


###Finally, need to do the same regridding for WRF
rm(list=ls())
load("~/Documents/Data/JE_WRF/WRF_25deg.RData")
# library(akima)
# cool2<-array(0,c(144,73,31))
# warm2<-array(0,c(144,73,30))
# for(i in 1:31)
# {
#   a=as.vector(coolR[,,i])
#   a[which(is.na(a))]=0
#   b=interp(lont,latt,a,lon,lat)
#   cool2[,,i]=b$z
#   a=as.vector(warmR[,,i])
#   a[which(is.na(a))]=0
#   b=interp(lont,latt,a,lon,lat)
#   warm2[,,i]=b$z
# }
# 
# f1=open.nc('~/Documents/Data/WRF_meanwind.nc')
# U=var.get.nc(f1,'U')
# U2=array(0,c(144,73,372))
# 
# for(i in 1:372)
# {
#   a=as.vector(U[,,i])
#   a[which(is.na(a))]=0
#   b=interp(lont,latt,a,lon,lat)
#   U2[,,i]=b$z
# }
Uind1=apply(U2[61,47:51,],2,mean) ##25-35
Uind2=apply(U2[61,49:53,],2,mean) ##30-40

years=cbind(seq(1979,2009),matrix(NaN,31,2))
for(i in 1:length(years[,1]))
{
  I=which(time2[,1]==years[i,1] & time2[,2]>=5 & time2[,2]<=10)
  years[i,2]=mean(Uind1[I])
  I=which((time2[,1]==years[i,1]+1 & time2[,2]<=4) | (time2[,1]==years[i,1] & time2[,2]>=11))
  years[i,3]=mean(Uind1[I])
}
years[31,3]=NaN

coolC<-warmC<-matrix(0,144,73)
for(i in 1:144) 
  for(j in 1:73) 
  {
    coolC[i,j]=cor(cool2[i,j,],years[,2])
    warmC[i,j]=cor(warm2[i,j,],years[1:30,3])
  }

library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
library(fields)
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(14)
cm[7:8]="white"
coolC[coolC>=0.69]=0.69
coolC[coolC<=(-0.69)]=-0.69
warmC[warmC>=0.69]=0.69
warmC[warmC<=(-0.69)]=-0.69
coolC[abs(coolC)<=0.27]=NaN
warmC[abs(warmC)<=0.27]=NaN
tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_sig_WRFregrid_cool.tiff", height=500, width=800)
image.plot(lon,lat[73:1],coolC[,73:1],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_sig_WRFregrid_warm.tiff", height=500, width=800)
image.plot(lon,lat[73:1],warmC[,73:1],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

###Choosing a location for testing

list$inland=0
list$coast=0
i=which(lat==-35)
j=which(lon==145)
for(k in 1:37) list[k,4]=coolC[j,i,k]
i=which(lat==-30)
j=which(lon==152.5)
for(k in 1:37) list[k,5]=coolC[j,i,k]


list[38,]=c("AWAP","190001","201312",0,0)
i=which(lat==-35)
j=which(lon==145)
list[38,4]=coolC[j,i]
i=which(lat==-30)
j=which(lon==152.5)
list[38,5]=coolC[j,i]

##Rain averages
rm(list=ls())
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
source('~/Documents/R/CMIP5_play.R')
coolM<-warmM<-array(0,c(144,73,37))
for(i in 1:37)
{
  a=CMIP_rain(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3])) 
  coolM[,,i]=unlist(a[1])
  warmM[,,i]=unlist(a[2])
}

##Mean, SD, quartiles
coolC2<-warmC2<-array(0,c(144,73,4))
coolC2[,,1]=apply(coolM,c(1,2),mean,na.rm=T)
coolC2[,,2]=apply(coolM,c(1,2),quantile,0.25,na.rm=T)
coolC2[,,3]=apply(coolM,c(1,2),quantile,0.5,na.rm=T)
coolC2[,,4]=apply(coolM,c(1,2),quantile,0.7,na.rm=T)
warmC2[,,1]=apply(warmM,c(1,2),mean,na.rm=T)
warmC2[,,2]=apply(warmM,c(1,2),quantile,0.25,na.rm=T)
warmC2[,,3]=apply(warmM,c(1,2),quantile,0.5,na.rm=T)
warmC2[,,4]=apply(warmM,c(1,2),quantile,0.7,na.rm=T)

library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
names=c("mean","Q25","Q50","Q75")
Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",list[1,1],"_historical_r1i1p1_",list[1,2],"-",list[1,3],"_2.5deg.nc",sep="")
f1=open.nc(Rfile)
lat=var.get.nc(f1,"lat")
lon=var.get.nc(f1,"lon")

for(i in 1:4)
{
  tiff(file=paste("~/Documents/Data/CMIP5/AveRain/Rain_CMIP5_",names[i],"_cool.tiff",sep=""), height=500, width=800)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,coolC2[,,i],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
  dev.off()
  tiff(file=paste("~/Documents/Data/CMIP5/AveRain/Rain_CMIP5_",names[i],"_warm.tiff",sep=""), height=500, width=800)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,warmC2[,,i],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
  dev.off()  
}


###Playing with the list of correlations

data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/25 to 35/comp.csv")
tiff(file="~/Documents/Data/CMIP5/Corrs_scatter.tiff",height=500,width=400)
plot(data[2:39,6],data[2:39,5],type="p",col="red",xlab="Coastal correlation (30S,152.5E)",ylab="Inland correlation (35S,145E)",xlim=c(-0.6,0.2),ylim=c(-0.6,0.6))
abline(h=0,col="grey")
abline(v=0,col="grey")
points(data[2:39,6],data[2:39,5],col="red")
##Actually, want: Black for AWAP, grey for both right, red for only coast right, unfilled red where neither
I=which(data[,6]<=(-0.25) & data[,5]<0.25)
points(data[I,6],data[I,5],col="red",pch=16)
I=which(data[,6]<=(-0.25) & data[,5]>=0.25)
points(data[I,6],data[I,5],col="blue",pch=19)
points(data[1,6],data[1,5],col="black",pch=17,cex=2)
points(median(data[2:39,6]),median(data[2:39,5]),col="red",pch=17,cex=2)
dev.off()


##Testing dates
rm(list=ls())
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list

dates=data.frame(names=list[,1],d1=list[,2],d1a=rep(0,39),d2=list[,3],d2a=rep(0,39),diff=rep(0,39))
library(RNetCDF)
for(i in 1:39)
  {
  name=as.character(list[i,1])
  date1=as.character(list[i,2])
  date2=as.character(list[i,3])
  Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  f1=open.nc(Rfile)
  time=var.get.nc(f1,"time")
  ori=att.get.nc(f1,"time","units")
   ori=unlist(strsplit(ori, split=" "))[3]
   time2=as.Date(time,origin=ori)
   dates[i,3]=as.numeric(format(time2[1],"%Y"))*100+as.numeric(format(time2[1],"%m"))
   dates[i,5]=as.numeric(format(time2[length(time)],"%Y"))*100+as.numeric(format(time2[length(time)],"%m"))
  dd=as.numeric(format(time2,"%Y"))*100+as.numeric(format(time2,"%m"))
  dates[i,6]=length(dd)-length(unique(dd))
}



