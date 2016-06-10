
##Currently edited to produce figs for WRF10 not WRF50, but quick to change refs.
rm(list=ls())
setwd("~/Documents/Data/CMIP5/CMIP3")
library(RNetCDF)
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

a=23
b=47
f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")
timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
NCEP=apply(var.get.nc(f1,'uwnd',c(61,b,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)

Ucomp=matrix(0,12,5)
for(i in 1:12)
{
  Ucomp[i,1]=i
  I=which(timeN[,1]>=1950 & timeN[,1]<=2005 & timeN[,2]==i)
  Ucomp[i,2]=mean(NCEP[I])
}

###50 km mask??

lat=seq(-50,0,0.5)
lon=seq(105,180,0.5)

maskS<-maskN<-matrix(NaN,151,101)
J=which(lat>(-33.75) & lat<(-23.75))
I=which(lon>151.25 & lon<153.75)
maskN[I,J]=1

J=which(lat>(-38.75) & lat<(-31.25))
for(i in 1:length(J))
{
  l1=131.25-(10/7.5)*(31.25+lat[J[i]])
  if(lat[J[i]]<(-36.25)) l2=146.75 else l2=138.75-(8/5)*(31.25+lat[J[i]])
  I=which(lon<l2 & lon>l1)
  maskS[I,J[i]]=1
}

corrlist=data.frame(Names=rep("aaa",3),Avecorr_SW=rep(NaN,3),Avecorr_NE=rep(NaN,3),
                    Corrave_SW=rep(NaN,3),Corrave_NE=rep(NaN,3),
                    Point_SW=rep(NaN,3),Point_NE=rep(NaN,3), stringsAsFactors=FALSE)
coolA<-warmA<-array(NaN,dim=c(151,101,3))

UWRF=matrix(0,672,4)
I=which(timeN[,1]>=1950 & timeN[,1]<=2005)
UWRF[,1]=NCEP[I]
for(r in 1:3)
{
  W_ufile=paste("/srv/ccrc/data34/z3478332/WRF/WRF_R",r,"_ncep1_uwnd.nc",sep="")  
  f1=open.nc(W_ufile)
  time=var.get.nc(f1,"time")
  year=floor(time/100)
  month=time%%100
  time=cbind(year,month)
  Uave=apply(var.get.nc(f1,'uwnd250',c(61,a,1),c(1,5,length(time[,1])),unpack=T),2,mean,na.rm=T)

  I=which(time[,1]>=1950 & time[,1]<=2005)
  UWRF[,r+1]=Uave[I]
}



count=1
for(r in 1:3)
{
W_ufile=paste("/srv/ccrc/data34/z3478332/WRF/WRF_R",r,"_ncep1_uwnd.nc",sep="")
W_pfile=paste("/srv/ccrc/data34/z3478332/WRF/WRF_R",r,"_ncep1_precip.nc",sep="")
corrlist[count,1]=paste("R",r,sep="")

fig1=paste("Uwnd_NCEP_WRFR",r,"_cool.tiff",sep="")
fig2a=paste("Rain_corrUwind_sig_NCEP_WRFR",r,"_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_NCEP_WRFR",r,"_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_NCEP_WRFR",r,"_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_NCEP_WRFR",r,"_warm.tiff",sep="")

f1=open.nc(W_ufile)
time=var.get.nc(f1,"time")
year=floor(time/100)
month=time%%100
time=cbind(year,month)
Uave=apply(var.get.nc(f1,'uwnd250',c(61,a,1),c(1,5,length(time[,1])),unpack=T),2,mean,na.rm=T)

f1=open.nc(W_pfile)
lat=var.get.nc(f1,"lat50")
lon=var.get.nc(f1,"lon50")
rain=var.get.nc(f1,"prcp50")

I=which(time[,1]>=1950 & time[,1]<=2005)
time=time[I,]
Uave=Uave[I]
rain=rain[,,I]

for(i in 1:12)
{
  I=which(time[,2]==i)
  Ucomp[i,2+r]=mean(Uave[I])
}

#Finally, let's do the two correlation plots
years=seq(min(time[,1]),max(time[,1]))
years=cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
coolM=matrix(0,length(years[,1]),2)
coolR<-warmR<-array(0,c(151,101,length(years[,1])))
for(i in 1:length(years[,1]))
{
  I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
  years[i,2]=mean(Uave[I])
  coolR[,,i]=apply(rain[,,I],c(1,2),sum)
  coolM[i,1]=mean(coolR[,,i]*maskS,na.rm=T)
  coolM[i,2]=mean(coolR[,,i]*maskN,na.rm=T)
  I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
  years[i,3]=mean(Uave[I])
  warmR[,,i]=apply(rain[,,I],c(1,2),sum)
}

coolC<-warmC<-matrix(0,151,101)
for(i in 1:151)
  for(j in 1:101)
  {
    coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
    warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
  }
coolA[,,count]=coolC
warmA[,,count]=warmC
count=count+1

# corrlist[count,2]=mean(coolC*maskS,na.rm=T)
# corrlist[count,3]=mean(coolC*maskN,na.rm=T)
# corrlist[count,4]=cor(coolM[,1],years[,2],use="pairwise.complete.obs")
# corrlist[count,5]=cor(coolM[,2],years[,2],use="pairwise.complete.obs")
# 
# i=which(lat==-35)
# j=which(lon==145)
# corrlist[count,6]=coolC[j,i]
# i=which(lat==-30)
# j=which(lon==152.5)
# corrlist[count,7]=coolC[j,i]
# count=count+1
# 
# coolC[coolC>=0.69]=0.69
# coolC[coolC<=(-0.69)]=-0.69
# warmC[warmC>=0.69]=0.69
# warmC[warmC<=(-0.69)]=-0.69
# 
# tiff(file=fig2, height=500, width=800)
# image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# tiff(file=fig3, height=500, width=800)
# image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# 
# coolC[abs(coolC)<=0.27]=NaN
# warmC[abs(warmC)<=0.27]=NaN
# 
# tiff(file=fig2a, height=500, width=800)
# image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# tiff(file=fig3a, height=500, width=800)
# image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()

}

cool2=apply(coolA,c(1,2),median)
warm2=apply(warmA,c(1,2),median)
fig2a=paste("Rain_corrUwind_sig_NCEP_WRF_median_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_NCEP_WRF_median_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_NCEP_WRF_median_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_NCEP_WRF_median_warm.tiff",sep="")

cool2[cool2>=0.69]=0.69
cool2[cool2<=(-0.69)]=-0.69
warm2[warm2>=0.69]=0.69
warm2[warm2<=(-0.69)]=-0.69

tiff(file=fig2, height=500, width=600,pointsize=20)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=600,pointsize=20)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

cool2[abs(cool2)<=0.27]=NaN
warm2[abs(warm2)<=0.27]=NaN

tiff(file=fig2a, height=500, width=600,pointsize=20)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3a, height=500, width=600,pointsize=20)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

tiff(file="~/Documents/Data/CMIP5/CMIP3/Uwind_ncep.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp[,2],type="l",col="blue",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-5,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(-5,10,5),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","WRF Average"),col=c("blue","purple"),lwd=2,bty="n",cex=1.5)
for(i in 1:3) lines(seq(1:12),Ucomp[,i+2],col="gray")
lines(seq(1,12),Ucomp[,2],col="blue",lwd=2)
lines(seq(1,12),apply(Ucomp[,3:5],1,mean),col="purple",lwd=2)
dev.off()
