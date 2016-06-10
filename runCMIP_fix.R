
##Version 2.0

rm(list=ls())
library(RNetCDF)
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
source('~/Documents/R/CMIP5_play.R')
f1=open.nc('/srv/ccrc/data34/z3478332/NCEP/Monthly/uwnd.mon.mean.nc')
lon=var.get.nc(f1,"lon")
lat=var.get.nc(f1,"lat")
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")
timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
##If using northern region
NCEP=apply(var.get.nc(f1,'uwnd',c(61,47,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)

f1=open.nc('/srv/ccrc/data34/z3478332/ERAI/erai_uwnd_150.nc')
timeE=var.get.nc(f1,'time')
timeE=as.Date((timeE/24-2),origin="1900-01-15")
timeE=cbind(as.numeric(format(timeE,"%Y")),as.numeric(format(timeE,"%m")))
latE=var.get.nc(f1,'latitude')
I=which(latE>=-36 & latE<=-24)
lonE=var.get.nc(f1,'longitude')
J=which(lonE>=148.5 & lonE<=151.5)
##If using northern region
ERAI=apply(var.get.nc(f1,'u',c(J[1],I[1],1),c(length(J),length(I),length(timeE[,1])),unpack=T),3,mean,na.rm=T)
Ucomp=matrix(0,12,39)
for(i in 1:12)
{
  Ucomp[i,1]=i
  I=which(timeN[,1]<=2005 & timeN[,1]>=1950 &timeN[,2]==i)
  Ucomp[i,2]=mean(NCEP[I])
}

gdi1=read.csv("~/Documents/Timeseries/gdi.csv")
GDI=rep(0,788)
for(i in 1:length(GDI))
{
  I=which(gdi1[,1]==timeN[i,1])
  GDI[i]=gdi1[I,timeN[i,2]+1]
}

I=which(timeN[,1]>=1979 & timeN[,1]<=2005)
J=which(timeE[,1]>=1979 & timeE[,1]<=2005)
cor(GDI[I],NCEP[I])


Ucomp2=matrix(0,13,4)
Ucomp2[,1]=seq(0,12)
I=which(timeN[,1]<=2005 & timeN[,1]>=1979)
J=which(timeE[,1]<=2012 & timeE[,1]>=1979)
Ucomp2[1,2]=mean(NCEP[I])
Ucomp2[1,3]=mean(ERAI[J])
Ucomp2[1,4]=cor(NCEP[I],ERAI[J])
for(i in 1:12)
{
  I=which(timeN[,1]<=2012 & timeN[,1]>=1979 & timeN[,2]==i)
  J=which(timeE[,1]<=2012 & timeE[,1]>=1979 & timeE[,2]==i)
  Ucomp2[i+1,2]=mean(NCEP[I])
  Ucomp2[i+1,3]=mean(ERAI[J])
  Ucomp2[i+1,4]=cor(NCEP[I],ERAI[J])
}
names(Ucomp2)=c("Month","NCEP","ERAI","Cor")
write.csv(Ucomp2,file="Ucomp.csv")

for(i in 1:37) Ucomp[,i+2]=CMIP_wind_v2(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N") 
Ucomp=cbind(Ucomp,seq(1:12))

list$Mean=apply(Ucomp[,4:40],2,mean)
list$CoolMean=apply(Ucomp[5:10,4:40],2,mean)

coolC<-warmC<-array(0,c(144,73,37))
for(i in 1:37)
{
  a=CMIP_v2(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N") 
  Ucomp[,i+2]=unlist(a[1])
  coolC[,,i]=unlist(a[2])
  warmC[,,i]=unlist(a[3])
}

save(Ucomp,coolC,warmC,file="CMIP5corrs.RData")

Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",list[1,1],"_historical_r1i1p1_",list[1,2],"-",list[1,3],"_2.5deg.nc",sep="")
f1=open.nc(Rfile)
lat=var.get.nc(f1,"lat")
lon=var.get.nc(f1,"lon")
list$inland=0
list$coast=0
i=which(lat==-35)
j=which(lon==145)
for(k in 1:37) list[k,4]=coolC[j,i,k]
i=which(lat==-30)
j=which(lon==152.5)
for(k in 1:37) list[k,5]=coolC[j,i,k]
write.csv(list,file="comp.csv")

##Mean, SD, quartiles
coolC2<-warmC2<-array(0,c(144,73,5))
coolC2[,,1]=apply(coolC,c(1,2),mean,na.rm=T)
coolC2[,,2]=apply(coolC,c(1,2),quantile,0.25,na.rm=T)
coolC2[,,3]=apply(coolC,c(1,2),quantile,0.5,na.rm=T)
coolC2[,,4]=apply(coolC,c(1,2),quantile,0.75,na.rm=T)
warmC2[,,1]=apply(warmC,c(1,2),mean,na.rm=T)
warmC2[,,2]=apply(warmC,c(1,2),quantile,0.25,na.rm=T)
warmC2[,,3]=apply(warmC,c(1,2),quantile,0.5,na.rm=T)
warmC2[,,4]=apply(warmC,c(1,2),quantile,0.75,na.rm=T)
Ucomp2=cbind(Ucomp[,1:2],matrix(0,12,5))
Ucomp2[,3]=apply(Ucomp[,3:39],1,mean,na.rm=T)
Ucomp2[,4]=apply(Ucomp[,3:39],1,sd,na.rm=T)
Ucomp2[,5]=apply(Ucomp[,3:39],1,quantile,0.25,na.rm=T)
Ucomp2[,6]=apply(Ucomp[,3:39],1,quantile,0.5,na.rm=T)
Ucomp2[,7]=apply(Ucomp[,3:39],1,quantile,0.75,na.rm=T)

##Plotting for Uwind - two versions - quartiles & std
tiff(file="~/Documents/Data/CMIP5/Uwind_all.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp2[,3],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-5,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(-5,10,5),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","CMIP5 Average"),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
for(i in 1:37) lines(seq(1:12),Ucomp[,i+2],col="gray")
lines(seq(1,12),Ucomp2[,2],col="blue",lwd=2)
lines(seq(1,12),Ucomp2[,3],col="red",lwd=2)
abline(h=0,col="black",lwd=2,lty=2)
dev.off()

tiff(file="~/Documents/Data/CMIP5/Uwind_all_erai.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp2[,3],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-5,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(-5,10,5),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","ERAI","CMIP5 Average"),col=c("blue","purple","red"),lwd=2,bty="n",cex=1.5)
for(i in 1:37) lines(seq(1:12),Ucomp[,i+2],col="gray")
lines(seq(1,12),Ucomp2[,2],col="blue",lwd=2)
lines(seq(1,12),Uerai[,2],col="purple",lwd=2)
lines(seq(1,12),Ucomp2[,3],col="red",lwd=2)
dev.off()




library(hydroGOF)
tiff(file="~/Documents/Data/CMIP5/Uwind_meansd.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp2[,3],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-5,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(-5,10,5),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","CMIP5 Average"),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
plotbandsonly(Ucomp2[,3]-Ucomp2[,4],Ucomp2[,3]+Ucomp2[,4], bands.col=rgb(0.9,0.9,0.9),ylim=c(-5,10),alpha=0.5)
lines(seq(1,12),Ucomp2[,2],col="blue",lwd=2)
lines(seq(1,12),Ucomp2[,3],col="red",lwd=2)
dev.off()

tiff(file="~/Documents/Data/CMIP5/Uwind_quart.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp2[,6],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-5,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(-5,10,5),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","CMIP5 Median"),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
plotbandsonly(Ucomp2[,5],Ucomp2[,7], bands.col=rgb(0.9,0.9,0.9),ylim=c(-5,10),alpha=0.5)
lines(seq(1,12),Ucomp2[,2],col="blue",lwd=2)
lines(seq(1,12),Ucomp2[,6],col="red",lwd=2)
dev.off()



##Plotting for Corrs - Lots of plots - mean, q25,median,q75
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
#coolC2[abs(coolC2)<=0.27]=NaN
warmC2[warmC2>=0.69]=0.69
warmC2[warmC2<=(-0.69)]=-0.69
#warmC2[abs(warmC2)<=0.27]=NaN
names=c("mean","Q25","Q50","Q75")

for(i in 1:4)
{
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_CMIP5_",names[i],"_coola_small.tiff",sep=""), height=500, width=600,pointsize=20)
  image.plot(lon,lat,coolC2[,,i],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_CMIP5_",names[i],"_warma_small.tiff",sep=""), height=500, width=600,pointsize=20)
  image.plot(lon,lat,warmC2[,,i],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

coolC2[abs(coolC2)<=0.27]=NaN
warmC2[abs(warmC2)<=0.27]=NaN

for(i in 1:4)
{
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_CMIP5_",names[i],"_coola_small.tiff",sep=""), height=500, width=600,pointsize=20)
  image.plot(lon,lat,coolC2[,,i],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_CMIP5_",names[i],"_warma_small.tiff",sep=""), height=500, width=600,pointsize=20)
  image.plot(lon,lat,warmC2[,,i],xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/25-35/comp.csv")
tiff(file="~/Documents/Data/CMIP5/Corrs_scatter_v3.tiff",height=600,width=600,pointsize=20)
plot(data[2:39,6],data[2:39,5],type="p",col="grey",xlab="Coastal correlation (30S,152.5E)",ylab="Inland correlation (35S,145E)",xlim=c(-0.6,0.2),ylim=c(-0.6,0.6))
abline(h=0,col="grey")
abline(v=0,col="grey")
##Actually, want: Black for AWAP, grey for both right, red for only coast right, unfilled red where neither
I=which(data[,6]<=(-0.265) & data[,5]<0.265)
points(data[I,6],data[I,5],col="darkgrey",pch=16)
I=which(data[,6]<=(-0.265) & data[,5]>=0.265)
points(data[I,6],data[I,5],col="darkgrey",pch=18,cex=2)
points(data[1,6],data[1,5],col="black",pch=17,cex=2)
points(median(data[2:39,6]),median(data[2:39,5]),col="darkgrey",pch=17,cex=2)
dev.off()

data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/20c3m/comp.csv")
tiff(file="~/Documents/Data/CMIP5/CMIP3/20c3m/Corrs_scatter_CMIP3_10km.tiff",height=600,width=600,pointsize=20)
plot(data[2:5,5],data[2:5,4],type="p",col="red",xlab="Coastal correlation (30S,152.5E)",ylab="Inland correlation (35S,145E)",xlim=c(-0.8,0.2),ylim=c(-0.6,0.6))
points(median(data[2:5,5]),median(data[2:5,4]),col="red",pch=17,cex=2)
abline(h=0,col="grey")
abline(v=0,col="grey")
points(data[6:9,5],data[6:9,4],col="blue")
points(median(data[6:9,5]),median(data[6:9,4]),col="blue",pch=17,cex=2)
points(data[10:13,5],data[10:13,4],col="purple")
points(median(data[10:13,5]),median(data[10:13,4]),col="purple",pch=17,cex=2)
points(data[1,5],data[1,4],col="black",pch=17,cex=2)
legend("bottomright",pch=c(17,17,17),col=c("black","red","blue","purple"),legend=c("AWAP","CMIP3","WRF50","WRF10"))
dev.off()


###Actually, also need to do mean annual/cool/warm rainfall

##Version 2.0

rm(list=ls())
library(RNetCDF)
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
source('~/Documents/R/CMIP5_play.R')

coolC<-warmC<-annC<-array(0,c(144,73,37))
for(i in 1:37)
{
  a=CMIP_rain_v2(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N")
  annC[,,i]=unlist(a[1])
  coolC[,,i]=unlist(a[2])
  warmC[,,i]=unlist(a[3])
}

##Mean, SD, quartiles
coolC2<-warmC2<-annC2<-array(0,c(144,73,5))
# annC2[,,1]=apply(annC,c(1,2),mean,na.rm=T)
# annC2[,,2]=apply(annC,c(1,2),quantile,0.25,na.rm=T)
# annC2[,,3]=apply(annC,c(1,2),quantile,0.5,na.rm=T)
# annC2[,,4]=apply(annC,c(1,2),quantile,0.75,na.rm=T)
coolC2[,,1]=apply(coolC,c(1,2),mean,na.rm=T)
coolC2[,,2]=apply(coolC,c(1,2),quantile,0.25,na.rm=T)
coolC2[,,3]=apply(coolC,c(1,2),quantile,0.5,na.rm=T)
coolC2[,,4]=apply(coolC,c(1,2),quantile,0.75,na.rm=T)
warmC2[,,1]=apply(warmC,c(1,2),mean,na.rm=T)
warmC2[,,2]=apply(warmC,c(1,2),quantile,0.25,na.rm=T)
warmC2[,,3]=apply(warmC,c(1,2),quantile,0.5,na.rm=T)
warmC2[,,4]=apply(warmC,c(1,2),quantile,0.75,na.rm=T)

##Plotting for Corrs - Lots of plots - mean, q25,median,q75
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
#   tiff(file=paste("~/Documents/Data/CMIP5/Rain_CMIP5_",names[i],"_ann.tiff",sep=""), height=500, width=800)
#   par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
#   image(lon,lat,annC2[,,i],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
#   image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
#   dev.off()
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_CMIP5_",names[i],"_cool.tiff",sep=""), height=500, width=800)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,coolC2[,,i],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
  dev.off()
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_CMIP5_",names[i],"_warm.tiff",sep=""), height=500, width=800)
  par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,warmC2[,,i],xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
  dev.off()
}

###Need AWAP at native resolution

load('~/Documents/Data/AWAPrain.RData')
yearsA=seq(1900,2012)
I=which(yearsA>=1950 & yearsA<=2005)
cool=apply(cool[,,I],c(1,2),mean)
yearsB=seq(1901,2012)
I=which(yearsB>=1950 & yearsB<=2005) #So has the year ending in 1950 - year ending in 2005
warm=apply(warm[,,I],c(1,2),mean) #So has the year ending in 1950 - year ending in 2005

library(fields)
 a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(cool)),list(x=lon,y=lat))
 cool1=a$z
 a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(warm)),list(x=lon,y=lat))
 warm1=a$z
 a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(cool*Useful$mask)),list(x=lon,y=lat))
 cool2=a$z
 a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(warm*Useful$mask)),list(x=lon,y=lat))
 warm2=a$z

library(RNetCDF)



par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
image(lon,lat,warm2,xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))


source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
for(i in 1:4)
{
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_CMIP5anomPC_",names[i],"_cool.tiff",sep=""), height=500, width=600,pointsize=20)
  par(plt = c(0.1,0.8,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,(coolC2[,,i]-cool2)/cool2,xlab="",ylab="",breaks=c(-5,seq(-1,1,0.2),5),col=pal(12),zlim=c(-5,5),,xlim=c(110,160),ylim=c(-45,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("","-100%","-80%","-60%","-40%","-20%","","+20%","+40%","+60%","+80%","+100%",""),col=pal(12),zlim=c(0,12))
  dev.off()
  tiff(file=paste("~/Documents/Data/CMIP5/Rain_CMIP5anomPC_",names[i],"_warm.tiff",sep=""), height=500, width=600,pointsize=20)
  par(plt = c(0.1,0.8,0.15,0.85),las = 1,cex.axis = 1)
  image(lon,lat,(warmC2[,,i]-warm2)/warm2,xlab="",ylab="",breaks=c(-5,seq(-1,1,0.2),5),col=pal(12),zlim=c(-5,5),xlim=c(110,160),ylim=c(-45,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
  image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("","-100%","-80%","-60%","-40%","-20%","","+20%","+40%","+60%","+80%","+100%",""),col=pal(12),zlim=c(0,12))
  dev.off()
}

####Proportion of models that overestimate?
coolOver<-warmOver<-array(0,c(144,73))
for(i in 1:37)
{
  I=which(coolC[,,i]>cool2)
  coolOver[I]=coolOver[I]+1
  I=which(warmC[,,i]>warm2)
  warmOver[I]=warmOver[I]+1
}
coolOver=coolOver/37
coolOver[is.na(cool2)==1]=NaN
warmOver=warmOver/37
warmOver[is.na(cool2)==1]=NaN

tiff(file=paste("~/Documents/Data/CMIP5/Rain_CMIP5_propOver_cool.tiff",sep=""), height=500, width=600,pointsize=20)
par(plt = c(0.1,0.8,0.15,0.85),las = 1,cex.axis = 1)
image(lon,lat,coolOver,xlab="",ylab="",breaks=seq(0,1,0.1),col=pal(10),zlim=c(0,1),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
image.plot(legend.only=TRUE,breaks=seq(0,10),lab.breaks=c("0%","","20%","","40%","","60%","","80%","","100%"),col=pal(10),zlim=c(0,10))
dev.off()
tiff(file=paste("~/Documents/Data/CMIP5/Rain_CMIP5_propOver_warm.tiff",sep=""), height=500, width=600,pointsize=20)
par(plt = c(0.1,0.8,0.15,0.85),las = 1,cex.axis = 1)
image(lon,lat,warmOver,xlab="",ylab="",breaks=seq(0,1,0.1),col=pal(10),zlim=c(0,1),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
image.plot(legend.only=TRUE,breaks=seq(0,10),lab.breaks=c("0%","","20%","","40%","","60%","","80%","","100%"),col=pal(10),zlim=c(0,10))
dev.off()


###Using a different region for AWAP?
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(14)
cm[7:8]="white"

load('~/Documents/Data/AWAPrain.RData')
yearsA=seq(1900,2012)
I=which(yearsA>=1950 & yearsA<=2005)
cool=cool[,,I]
yearsB=seq(1901,2012)
I=which(yearsB>=1950 & yearsB<=2005) #So has the year ending in 1950 - year ending in 2005
warm=warm[,,I] 

lat2=seq(-50,0,0.5)
lon2=seq(105,180,0.5)
years=cbind(seq(1950,2005),matrix(0,56,2))
cool2=array(NaN,dim=c(length(lon),length(lat),length(years[,1])))
cool3=array(NaN,dim=c(length(lon2),length(lat2),length(years[,1])))

for(i in 1:length(years[,1]))
{
  I=which(timeN[,1]==years[i,1] & timeN[,2]>=5 & timeN[,2]<=10)
  years[i,2]=mean(NCEP[I])
  I=which((timeN[,1]==years[i,1] & timeN[,2]<=4) | (timeN[,1]==years[i,1]-1 & timeN[,2]>=11))
  years[i,3]=mean(NCEP[I])

  
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(cool[,,i])),list(x=lon,y=lat))
  cool2[,,i]=a$z
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(cool[,,i])),list(x=lon2,y=lat2))
  cool3[,,i]=a$z
}

coolC<-warmC<-matrix(0,691,886)
for(i in 1:691)
  for(j in 1:886)
  {
     coolC[i,j]=cor(cool[i,j,],years[,2],use="pairwise.complete.obs")
     warmC[i,j]=cor(warm[i,j,],years[,3],use="pairwise.complete.obs")
  }

coolC2<-matrix(0,144,73)
for(i in 1:144)
  for(j in 1:73)
     coolC2[i,j]=cor(cool2[i,j,],years[,2],use="pairwise.complete.obs")

coolC3<-matrix(0,151,101)
for(i in 1:151)
  for(j in 1:101)
     coolC3[i,j]=cor(cool3[i,j,],years[,2],use="pairwise.complete.obs")

###Full resolution
maskS<-maskN<-matrix(NaN,886,691)
J=which(Useful$y>(-33.75) & Useful$y<(-23.75))
I=which(Useful$x>151.25 & Useful$x<153.75)
maskN[I,J]=1

J=which(Useful$y>(-38.75) & Useful$y<(-31.25))
for(i in 1:length(J))
{
  l1=131.25-(10/7.5)*(31.25+Useful$y[J[i]])
  if(Useful$y[J[i]]<(-36.25)) l2=146.75 else l2=138.75-(8/5)*(31.25+Useful$y[J[i]])
  I=which(Useful$x<l2 & Useful$x>l1)
  maskS[I,J[i]]=1
}


rainN=matrix(0,length(years[,1]),8)
for(i in 1:length(years[,1]))
{
  rainN[i,1]=mean(cool[,,i]*Useful$mask*t(maskS),na.rm=T)
  rainN[i,2]=mean(cool[,,i]*Useful$mask*t(maskN),na.rm=T)
  rainN[i,3]=mean(cool[,,i]*t(maskS),na.rm=T)
  rainN[i,4]=mean(cool[,,i]*t(maskN),na.rm=T)
}

mean(coolC*t(maskN),na.rm=T)
mean(coolC*t(maskS),na.rm=T)
mean(coolC*t(maskN)*Useful$mask,na.rm=T)
mean(coolC*t(maskS)*Useful$mask,na.rm=T)

I=which(Useful$x>152.4 & Useful$x<152.6)
J=which(Useful$y>(-30.1) & Useful$y<(-29.9))
mean(coolC[J,I],na.rm=T)
I=which(Useful$x>144.9 & Useful$x<145.1)
J=which(Useful$y>(-35.1) & Useful$y<(-34.9))
mean(coolC[J,I],na.rm=T)

coolC[coolC>=0.69]=0.69
coolC[coolC<=(-0.69)]=-0.69
coolC[abs(coolC)<=0.27]=NaN

reglats=c(rep(-32.5,4),rep(-35,4),rep(-37.5,3))
reglons=c(seq(135,142.5,2.5),seq(137.5,145,2.5),seq(140,145,2.5))
reglats2=seq(-25,-32.5,-2.5)

tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_sig_AWAP_cool_corrmask_nobar_v2.tiff", height=1000, width=1100,pointsize=20)
#par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
image(as.vector(Useful$x),as.vector(Useful$y),t(coolC*Useful$mask),xlab="",
      breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),cex.axis=2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(reglons,reglats,pch=4,col="gray",lwd=4,cex=2)
polygon(x=c(131.25,141.25,146.75,146.75,138.75,131.25),y=c(-31.25,-31.25,-36.25,-38.75,-38.75,-31.25),lwd=8,border="gray")
points(rep(152.5,4),reglats2,pch=4,col="black",lwd=4,cex=2)
rect(151.25,-33.75,153.75,-23.75,lwd=8)
# par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
# filled.legend(Useful$x,Useful$y,t(coolC*Useful$mask),lev=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)) 
dev.off()

warmC[warmC>=0.69]=0.69
warmC[warmC<=(-0.69)]=-0.69
warmC[abs(warmC)<=0.27]=NaN
tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_sig_AWAP_warm_nobar_v3.tiff", height=1000, width=1100,pointsize=20)
#par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
image(as.vector(Useful$x),as.vector(Useful$y),t(warmC*Useful$mask),xlab="",
      breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),cex.axis=2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
# filled.legend(Useful$x,Useful$y,t(coolC*Useful$mask),lev=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)) 
dev.off()


sigmaskC=which(abs(coolC)*Useful$mask>=0.45,arr.ind=T)
I=which(Useful$x[sigmaskC[,2]]>155 | Useful$x[sigmaskC[,2]]<135)
if(length(I)>1) sigmaskC=sigmaskC[-I,]
I=which(Useful$y[sigmaskC[,1]]>(-23) | Useful$y[sigmaskC[,1]]<(-39))
if(length(I)>1) sigmaskC=sigmaskC[-I,]
sigmaskW=which(abs(warmC)*Useful$mask>=0.45,arr.ind=T)
I=which(Useful$x[sigmaskW[,2]]>155 | Useful$x[sigmaskW[,2]]<135)
if(length(I)>1) sigmaskW=sigmaskW[-I,]
I=which(Useful$y[sigmaskW[,1]]>(-23) | Useful$y[sigmaskW[,1]]<(-39))
if(length(I)>1) sigmaskW=sigmaskW[-I,]
tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_sigdot_AWAP_cool.tiff", height=1000, width=1100,pointsize=20)
filled.contour3(Useful$x,Useful$y,t(coolC*Useful$mask),xlab="",levels=seq(-0.7,0.7,0.1),ylab="",col=cm,
                plot.axes={
                  axis(1,cex.axis=2)
                  axis(2,cex.axis=2)
                },zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
par(xpd = NA)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
points(Useful$x[sigmaskC[,2]],Useful$y[sigmaskC[,1]],col="black",pch=".",cex=2,xlim=c(135,155),ylim=c(-39,-23))
dev.off()
tiff(file="~/Documents/Data/CMIP5/Rain_corrUwind_sigdot_AWAP_warm.tiff", height=1000, width=1100,pointsize=20)
filled.contour3(Useful$x,Useful$y,t(warmC*Useful$mask),xlab="",levels=seq(-0.7,0.7,0.1),ylab="",
                plot.axes={
                  axis(1,cex.axis=2)
                  axis(2,cex.axis=2)
                },col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
par(xpd = NA)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
points(Useful$x[sigmaskW[,2]],Useful$y[sigmaskW[,1]],col="black",pch=".",cex=2,xlim=c(135,155),ylim=c(-39,-23))
dev.off()



####50 km resolution
maskS3<-maskN3<-matrix(NaN,151,101)
J=which(lat2>(-33.75) & lat2<(-23.75))
I=which(lon2>151.25 & lon2<153.75)
maskN3[I,J]=1

J=which(lat2>(-38.75) & lat2<(-31.25))
for(i in 1:length(J))
{
  l1=131.25-(10/7.5)*(31.25+lat2[J[i]])
  if(lat2[J[i]]<(-36.25)) l2=146.75 else l2=138.75-(8/5)*(31.25+lat2[J[i]])
  I=which(lon2<l2 & lon2>l1)
  maskS3[I,J[i]]=1
}

for(i in 1:length(years[,1]))
{
  rainN[i,5]=mean(cool3[,,i]*maskS3,na.rm=T)
  rainN[i,6]=mean(cool3[,,i]*maskN3,na.rm=T)
}

mean(coolC3*maskN3,na.rm=T)
mean(coolC3*maskS3,na.rm=T)
coolC3[lon2==152.5,lat2==-30]
coolC3[lon2==145,lat2==-35]

coolC3[coolC3>=0.69]=0.69
coolC3[coolC3<=(-0.69)]=-0.69
tiff(file="Rain_corrUwind_AWAP_cool_postregrid50.tiff", height=500, width=800)
image.plot(lon2,lat2,coolC3,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

coolC3[abs(coolC3)<=0.27]=NaN
tiff(file="Rain_corrUwind_sig_AWAP_cool_postregrid50.tiff", height=500, width=800)
image.plot(lon2,lat2,coolC3,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

##250 km resolution
maskS2<-maskN2<-matrix(NaN,144,73)
reglats2=seq(-25,-32.5,-2.5)
for(i in 1:4)
{
  I=which(lon==152.5)
  J=which(lat==reglats2[i])
  maskN2[I,J]=1
}

reglats=c(rep(-32.5,4),rep(-35,4),rep(-37.5,3))
reglons=c(seq(135,142.5,2.5),seq(137.5,145,2.5),seq(140,145,2.5))
for(i in 1:length(reglats))
{
  I=which(lat==reglats[i])
  J=which(lon==reglons[i])
  maskS2[J,I]=1
}

for(i in 1:length(years[,1]))
{
  rainN[i,7]=mean(cool2[,,i]*maskS2,na.rm=T)
  rainN[i,8]=mean(cool2[,,i]*maskN2,na.rm=T)
}

mean(coolC2*maskN2,na.rm=T)
mean(coolC2*maskS2,na.rm=T)
coolC2[lon==152.5,lat==-30]
coolC2[lon==145,lat==-35]

##Corr ave rain

coolM=matrix(0,56,8)
for(i in 1:56)
{
coolM[i,1]=mean(cool[,,i]*t(maskS),na.rm=T)
coolM[i,2]=mean(cool[,,i]*t(maskN),na.rm=T)
coolM[i,3]=mean(cool[,,i]*Useful$mask*t(maskS),na.rm=T)
coolM[i,4]=mean(cool[,,i]*Useful$mask*t(maskN),na.rm=T)
coolM[i,5]=mean(cool3[,,i]*maskS3,na.rm=T)
coolM[i,6]=mean(cool3[,,i]*maskN3,na.rm=T)
coolM[i,7]=mean(cool2[,,i]*maskS2,na.rm=T)
coolM[i,8]=mean(cool2[,,i]*maskN2,na.rm=T)
}
for(i in 1:8) print(cor(coolM[,i],years[,2],use="pairwise.complete.obs"))

#############Okay, time to try

rm(list=ls())
library(RNetCDF)
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
source('~/Documents/R/CMIP5_play.R')
list=list[c(-24,-14),]
data=matrix(0,37,4)
for(i in 1:37) data[i,]=CMIP_corrs_v2(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N") 

data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/25-35/newcomp_v2.csv")
tiff(file="~/Documents/Data/CMIP5/Corrs_scatter_avecorr_reg_v2.tiff",height=600,width=600,pointsize=20)
a=6
b=5
plot(data[2:38,a],data[2:38,b],type="p",col="grey",xlab="Coastal correlation",ylab="Inland correlation",xlim=c(-0.7,0.2),ylim=c(-0.7,0.7))
abline(h=0,col="grey")
abline(v=0,col="grey")
##Actually, want: Black for AWAP, grey for both right, red for only coast right, unfilled red where neither
I=which(data[,a]<=(-0.265) & data[,b]<0.265)
points(data[I,a],data[I,b],col="darkgrey",pch=16)
I=which(data[,a]<=(-0.265) & data[,b]>=0.265)
points(data[I,a],data[I,b],col="darkgrey",pch=18,cex=2)
points(data[1,a],data[1,b],col="black",pch=4,cex=2,lwd=4)
points(median(data[2:38,a]),median(data[2:38,b]),col="black",pch=1,cex=2,lwd=4)
dev.off()

data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/25-35/newcomp_v3.csv")
##Only want ONE AWAP version - currently have 2.5 deg (1), 0.5 deg (2), 0.05 deg (3) and 0.05 deg landmask (4)
##Corrs are basically the same, except for 2.5 SW slightly weaker (by ~0.05)
##Let's stick with 2.5 for now
data=data[-c(1,2,3),]
tiff(file="~/Documents/Data/CMIP5/Corrs_scatter_avecorr_reg_v5.tiff",height=600,width=600,pointsize=20)
a=6
b=5
plot(data[2:38,a],data[2:38,b],type="p",col="grey",xlab="East Coast",ylab="South Coast",xlim=c(-0.7,0.2),ylim=c(-0.7,0.7))
abline(h=0,col="grey")
abline(v=0,col="grey")
##Actually, want: Black for AWAP, grey for both right, red for only coast right, unfilled red where neither
I=which(data[,a]<=(-0.265) & data[,b]<0.265)
points(data[I,a],data[I,b],col="darkgrey",pch=16)
I=which(data[,a]<=(-0.265) & data[,b]>=0.265)
points(data[I[2:length(I)],a],data[I[2:length(I)],b],col="darkgrey",pch=18,cex=2)
points(data[1,a],data[1,b],col="black",pch=4,cex=2,lwd=4)
points(median(data[2:38,a]),median(data[2:38,b]),col="black",pch=1,cex=2,lwd=4)
dev.off()

###CMIP5 subset for CMIP3 comparison

list=list[c(7,13,32,33,34,35),]
Ucomp=matrix(0,12,9)
for(i in 1:12)
{
  Ucomp[i,1]=i
  I=which(timeN[,1]<=2005 & timeN[,2]==i)
  Ucomp[i,2]=mean(NCEP[I])
}

coolC<-warmC<-array(0,c(144,73,6))
for(i in 1:6)
{
  a=CMIP_v2(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N") 
  Ucomp[,i+2]=unlist(a[1])
  coolC[,,i]=unlist(a[2])
  warmC[,,i]=unlist(a[3])
}

coolC2<-apply(coolC,c(1,2),median,na.rm=T)
warmC2<-apply(warmC,c(1,2),median,na.rm=T)
Ucomp[,9]<-apply(Ucomp[,3:8],1,median,na.rm=T)

tiff(file="~/Documents/Data/CMIP5/Uwind_CMIP3subset.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp[,2],type="l",col="blue",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-5,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(-5,10,5),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","CMIP5 Average"),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
for(i in 3:8) lines(seq(1:12),Ucomp[,i],col="gray")
lines(seq(1,12),Ucomp[,2],col="blue",lwd=2)
lines(seq(1,12),Ucomp[,9],col="red",lwd=2)
dev.off()

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

fig2="~/Documents/Data/CMIP5/Rain_corrUwind_sig_CMIP5subset_median_cool.tiff"
fig3="~/Documents/Data/CMIP5/Rain_corrUwind_sig_CMIP5subset_median_warm.tiff"
fig2a="~/Documents/Data/CMIP5/Rain_corrUwind_CMIP5subset_median_cool.tiff"
fig3a="~/Documents/Data/CMIP5/Rain_corrUwind_CMIP5subset_median_warm.tiff"

coolC2[coolC2>=0.69]=0.69
coolC2[coolC2<=(-0.69)]=-0.69
warmC2[warmC2>=0.69]=0.69
warmC2[warmC2<=(-0.69)]=-0.69
tiff(file=fig2a, height=500, width=800)
image.plot(lon,lat[73:1],coolC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3a, height=500, width=800)
image.plot(lon,lat[73:1],warmC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

coolC2[abs(coolC2)<=0.27]=NaN
warmC2[abs(warmC2)<=0.27]=NaN
tiff(file=fig2, height=500, width=800)
image.plot(lon,lat[73:1],coolC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=800)
image.plot(lon,lat[73:1],warmC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()


###"Good" models
read.csv("Modellist.csv")->list
list=list[c(23,26,28,29,30,36,37),]
Ucomp=matrix(0,12,9)
for(i in 1:12)
{
  Ucomp[i,1]=i
  I=which(timeN[,1]<=2005 & timeN[,2]==i)
  Ucomp[i,2]=mean(NCEP[I])
}

coolC<-warmC<-array(0,c(144,73,6))
for(i in 1:6)
{
  a=CMIP_v2(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N") 
  Ucomp[,i+2]=unlist(a[1])
  coolC[,,i]=unlist(a[2])
  warmC[,,i]=unlist(a[3])
}

coolC2<-apply(coolC,c(1,2),median,na.rm=T)
warmC2<-apply(warmC,c(1,2),median,na.rm=T)
Ucomp[,9]<-apply(Ucomp[,3:8],1,median,na.rm=T)

tiff(file="~/Documents/Data/CMIP5/Uwind_Good.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp[,2],type="l",col="blue",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-5,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(-5,10,5),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","CMIP5 Average"),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
for(i in 3:8) lines(seq(1:12),Ucomp[,i],col="gray")
lines(seq(1,12),Ucomp[,2],col="blue",lwd=2)
lines(seq(1,12),Ucomp[,9],col="red",lwd=2)
dev.off()

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

fig2="~/Documents/Data/CMIP5/Rain_corrUwind_sig_Good_median_cool.tiff"
fig3="~/Documents/Data/CMIP5/Rain_corrUwind_sig_Good_median_warm.tiff"
fig2a="~/Documents/Data/CMIP5/Rain_corrUwind_Good_median_cool.tiff"
fig3a="~/Documents/Data/CMIP5/Rain_corrUwind_Good_median_warm.tiff"

coolC2[coolC2>=0.69]=0.69
coolC2[coolC2<=(-0.69)]=-0.69
warmC2[warmC2>=0.69]=0.69
warmC2[warmC2<=(-0.69)]=-0.69
tiff(file=fig2a, height=500, width=800)
image.plot(lon,lat,coolC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3a, height=500, width=800)
image.plot(lon,lat,warmC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

coolC2[abs(coolC2)<=0.27]=NaN
warmC2[abs(warmC2)<=0.27]=NaN
tiff(file=fig2, height=500, width=800)
image.plot(lon,lat,coolC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=800)
image.plot(lon,lat,warmC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

fig2="~/Documents/Data/CMIP5/Rain_corrUwind_sig_Good_median_cool_small.tiff"
fig3="~/Documents/Data/CMIP5/Rain_corrUwind_sig_Good_median_warm_small.tiff"
tiff(file=fig2, height=500, width=600,pointsize=20)
image.plot(lon,lat,coolC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=600,pointsize=20)
image.plot(lon,lat,warmC2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()


###Mean rain?


rm(list=ls())
library(RNetCDF)
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
source('~/Documents/R/CMIP5_play.R')
Rmean=matrix(0,37,2)
for(i in 1:37) Rmean[i,]=CMIP_mean_rain(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N") 


####
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

good=c(23,26,28,29,30,36,37)-2
good=c(7,13,32,33,34,35)

medC=apply(coolC[,,good],c(1,2),median,na.rm=T)
medW=apply(warmC[,,good],c(1,2),median,na.rm=T)
medC[medC>=0.69]=0.69
medC[medC<=(-0.69)]=-0.69
medW[medW>=0.69]=0.69
medW[medW<=(-0.69)]=-0.69

sigmaskC=which(abs(medC)>=0.265,arr.ind=T)
I=which(lon[sigmaskC[,1]]>160 | lon[sigmaskC[,1]]<110)
if(length(I)>1) sigmaskC=sigmaskC[-I,]
I=which(lat[sigmaskC[,2]]>(-10) | lat[sigmaskC[,2]]<(-45))
if(length(I)>1) sigmaskC=sigmaskC[-I,]
sigmaskW=which(abs(medW)>=0.265,arr.ind=T)
I=which(lon[sigmaskW[,1]]>160 | lon[sigmaskW[,1]]<110)
if(length(I)>1) sigmaskW=sigmaskW[-I,]
I=which(lat[sigmaskW[,2]]>(-10) | lat[sigmaskW[,2]]<(-45))
if(length(I)>1) sigmaskW=sigmaskW[-I,]

lat=seq(-90,90,2.5)
fig2="~/Documents/Data/CMIP5/Rain_corrUwind_sigdot_CMIP3_median_cool_small.tiff"
fig3="~/Documents/Data/CMIP5/Rain_corrUwind_sigdot_CMIP3_median_warm_small.tiff"
tiff(file=fig2, height=1000, width=1100,pointsize=20)
image(lon,lat,medC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),cex.axis=2,xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon[sigmaskC[,1]],lat[sigmaskC[,2]],col="black",pch=19,cex=0.75)
dev.off()
tiff(file=fig3, height=1000, width=1100,pointsize=20)
image(lon,lat,medW,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),cex.axis=2,xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon[sigmaskW[,1]],lat[sigmaskW[,2]],col="black",pch=19,cex=0.75)
dev.off()


#########Average rain comparison
rm(list=ls())
library(RNetCDF)
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
source('~/Documents/R/CMIP5_play.R')



library("ncdf")

# Define netcdf dimensions: name, units and value
dim1 <- dim.def.ncdf("lon", "degreesE", lon)
dim2 <- dim.def.ncdf("lat", "degreesN", lat)
dim3 <- dim.def.ncdf("quart","",c(0,0.25,0.5,0.75))
dim4 <- dim.def.ncdf("masked","",c(0,1))
# Define objets of class var.def.ncdf for the netcdf file with its name, units, dimension(s) and missing value
var1 <- var.def.ncdf("cool_CMIP", "mm", list(dim1,dim2,dim3), -999)
var2 <- var.def.ncdf("warm_CMIP", "mm", list(dim1,dim2,dim3), -999)
var3 <- var.def.ncdf("cool_AWAP", "mm", list(dim1,dim2,dim4), -999)
var4 <- var.def.ncdf("warm_AWAP", "mm", list(dim1,dim2,dim4), -999)

# Create a netcdf file with its variables var1, var2 and var3
filename <- "meanraintest2.nc"
fid <- create.ncdf(filename, list(var1, var2,var3,var4))
put.var.ncdf(fid, var1, coolC2[,,1:4])
put.var.ncdf(fid, var2, warmC2[,,1:4])

a=abind(cool1,cool3,along=3)
a[is.na(a)]=-999
put.var.ncdf(fid, var3, a)
a=abind(warm1,warm3,along=3)
a[is.na(a)]=-999
put.var.ncdf(fid, var4, a)
# Close netcdf file access
close.ncdf(fid)


######Different mask?
a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=mask),list(x=lon,y=lat))
mask1=a$z
mask1[mask1==0]=NA

mask2=array(NaN,dim=dim(mask1))
for(i in 1:length(lat))
  for(j in 1:length(lon))
  {
    I=which(Useful$x>(lon[j]-1.25) & Useful$x<(lon[j]+1.25))
    J=which(Useful$y>(lat[i]-1.25) & Useful$y<(lat[i]+1.25))
    
    if(length(I)>0 & length(J)>0)     mask2[j,i]=mean(mask[I,J])
  }
mask2[mask2<0]=NA

mask3=mask2
mask3[mask2<0.2]=NA
mask3[mask2>=0.2]=1
image(lon,lat,mask3,xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

cool3=cool1*mask3
warm3=warm1*mask3
