
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

name=c("MPI_ECHAM5","CSIRO_Mk3.0","MIROC3.2","CGCM3.1")
name2=c("echam5","csiromk3","miroc3.2","cccma")
name3=c("mpi-echam5","csiromk3","miroc3.2","cgcm3.1")
date1a=c(191001,193101,185001,185001)
date2a=c(200912,200012,200012,200012)

a=23
b=47
f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")
timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
NCEP=apply(var.get.nc(f1,'uwnd',c(61,b,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)

yy=seq(1850,2012)
timeREF=matrix(0,length(yy)*12,3)
n=1
for(i in 1:163)
  for(j in 1:12)
  {
    timeREF[n,1]=yy[i]
    timeREF[n,2]=j
    timeREF[n,3]=yy[i]*100+j
    n=n+1
  }

r=1
n=1

##Uwind plot
for(n in 1:4)
{  
  C_ufile=paste("/srv/ccrc/data34/z3478332/CMIP3/uw850/hus_850hPa_20c3m_",name3[n],"_2.5deg.nc",sep="")
  fig1=paste("Uwind_",name[n],"_WRF10_Rcomp.tiff",sep="")
  
  f1=open.nc(C_ufile)
  u=var.get.nc(f1,"ua")
  CMIP=apply(u[61,a:(a+4),],2,mean)
  timeC=timeREF[(timeREF[,3]>=date1a[n] & timeREF[,3]<=date2a[n]),1:2]
  
  Ucomp=matrix(0,12,9)
  for(i in 1:12)
  {
    Ucomp[i,1]=i
    I=which(timeN[,1]<=2009 & timeN[,1]>=1950 & timeN[,2]==i)
    Ucomp[i,2]=mean(NCEP[I])
    I=which(timeC[,1]<=2009 & timeC[,1]>=1950 & timeC[,2]==i)
    Ucomp[i,3]=mean(CMIP[I])
  }  
  
  for(r in 1:3)
  {
    W_ufile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF_R",r,"_",name2[n],"_uwnd.nc",sep="")
    f1=open.nc(W_ufile)
    time=var.get.nc(f1,"time")
    year=floor(time/100)
    month=time%%100
    time=cbind(year,month)
    Uave=apply(var.get.nc(f1,'uwnd250',c(61,a,1),c(1,5,length(time[,1])),unpack=T),2,mean,na.rm=T)
    
    for(i in 1:12)
    {
      I=which(time[,2]==i)
      Ucomp[i,3+r]=mean(Uave[I])
    }     
  }
  
  for(r in 1:3)
  {
    W_ufile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF10_R",r,"_",name2[n],"_uwnd.nc",sep="")
    f1=open.nc(W_ufile)
    time=var.get.nc(f1,"time")
    year=floor(time/100)
    month=time%%100
    time=cbind(year,month)
    Uave=apply(var.get.nc(f1,'uwnd250',c(9,5,1),c(1,5,length(time[,1])),unpack=T),2,mean,na.rm=T)
    
    for(i in 1:12)
    {
      I=which(time[,2]==i)
      Ucomp[i,6+r]=mean(Uave[I])
    }     
  }
  
  mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  yl=c(-10,10)
  tiff(file=fig1, height=500, width=800)
  plot(seq(1:12),Ucomp[,2],ty="l",lwd=4,col="blue",xlab="",ylab="",axes=F,ylim=yl)
  axis(1, at = seq(1:12), labels = c(mnames),cex.axis=1.5) 
  axis(2, at = seq(yl[1],yl[2]),cex.axis=1.5)
  title(main="Average Uwind (m/s)",cex.main=2)
  lines(seq(1:12),Ucomp[,3],lwd=4,col="red")
  for(i in 4:6) lines(seq(1:12),Ucomp[,i],col="magenta3")
  lines(seq(1:12),apply(Ucomp[,4:6],1,mean),ty="l",lwd=4,col="magenta3")
  for(i in 7:9) lines(seq(1:12),Ucomp[,i],col="purple4")
  lines(seq(1:12),apply(Ucomp[,7:9],1,mean),ty="l",lwd=4,col="purple4")
  abline(h=0,col="gray",lty=2)
  legend("topleft",legend=c("NCEP",name[n],"WRF 50km","WRF 10km"),
         col=c("blue","red","magenta3","purple4"),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()
}

##Uwind plot - ALL
Ucomp1=matrix(0,12,7)
for(i in 1:12)
{
  Ucomp1[i,1]=i
  I=which(timeN[,1]<=2009 & timeN[,1]>=1950 & timeN[,2]==i)
  Ucomp1[i,2]=mean(NCEP[I])
}

Ucomp2<-Ucomp3<-matrix(0,12,13)
count=1

for(n in 1:4)
{  
  C_ufile=paste("/srv/ccrc/data34/z3478332/CMIP3/uw850/hus_850hPa_20c3m_",name3[n],"_2.5deg.nc",sep="")
  f1=open.nc(C_ufile)
  u=var.get.nc(f1,"ua")
  CMIP=apply(u[61,a:(a+4),],2,mean)
  timeC=timeREF[(timeREF[,3]>=date1a[n] & timeREF[,3]<=date2a[n]),1:2]
  
  for(i in 1:12)
  {
    I=which(timeC[,1]<=2009 & timeC[,1]>=1950 & timeC[,2]==i)
    Ucomp1[i,2+n]=mean(CMIP[I])
  }
  
  for(r in 1:3)
  {
    W_ufile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF_R",r,"_",name2[n],"_uwnd.nc",sep="")
    f1=open.nc(W_ufile)
    time=var.get.nc(f1,"time")
    year=floor(time/100)
    month=time%%100
    time=cbind(year,month)
    Uave=apply(var.get.nc(f1,'uwnd250',c(61,a,1),c(1,5,length(time[,1])),unpack=T),2,mean,na.rm=T)
    
    for(i in 1:12)
    {
      I=which(time[,2]==i)
      Ucomp2[i,count]=mean(Uave[I])
    }
    
    W_ufile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF10_R",r,"_",name2[n],"_uwnd.nc",sep="")
    f1=open.nc(W_ufile)
    time=var.get.nc(f1,"time")
    year=floor(time/100)
    month=time%%100
    time=cbind(year,month)
    Uave=apply(var.get.nc(f1,'uwnd250',c(9,5,1),c(1,5,length(time[,1])),unpack=T),2,mean,na.rm=T)
    
    for(i in 1:12)
    {
      I=which(time[,2]==i)
      Ucomp3[i,count]=mean(Uave[I])
    }     
    count=count+1
  }
}

Ucomp1[,7]=apply(Ucomp1[,3:6],1,median)
Ucomp2[,13]=apply(Ucomp2[,1:12],1,median)
Ucomp3[,13]=apply(Ucomp3[,1:12],1,median)
  
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
yl=c(-10,10)
tiff(file="Uwind_CMIP3_WRF10_mean.tiff", height=500, width=800)
  plot(seq(1:12),Ucomp1[,2],ty="l",lwd=4,col="blue",xlab="",ylab="",axes=F,ylim=yl)
  axis(1, at = seq(1:12), labels = c(mnames),cex.axis=1.5) 
  axis(2, at = seq(yl[1],yl[2]),cex.axis=1.5)
  title(main="Average Uwind (m/s)",cex.main=2)
abline(h=0,col="gray",lty=2)
#for(i in 3:6) lines(seq(1:12),Ucomp1[,i],ty="l",col="red",lty=2)
#for(i in 1:12) lines(seq(1:12),Ucomp2[,i],ty="l",col="magenta3",lty=2)
#for(i in 1:12) lines(seq(1:12),Ucomp3[,i],ty="l",col="purple4",lty=2)
lines(seq(1:12),Ucomp1[,2],ty="l",lwd=4,col="blue")
lines(seq(1:12),Ucomp1[,7],ty="l",lwd=4,col="red")
lines(seq(1:12),Ucomp2[,13],ty="l",lwd=4,col="magenta3")
lines(seq(1:12),Ucomp3[,13],ty="l",lwd=4,col="purple4")
legend("topleft",legend=c("NCEP","CMIP3 median","WRF50 median","WRF10 median"),
         col=c("blue","red","magenta3","purple4"),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()


###10 km mask??

lat=seq(-45,-20,0.1)
lon=seq(130,160,0.1)

maskS<-maskN<-matrix(NaN,301,251)
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

corrlist=data.frame(Names=rep("aaa",12),Avecorr_SW=rep(NaN,12),Avecorr_NE=rep(NaN,12),
                    Corrave_SW=rep(NaN,12),Corrave_NE=rep(NaN,12),
                    Point_SW=rep(NaN,12),Point_NE=rep(NaN,12), stringsAsFactors=FALSE)
coolA<-warmA<-array(NaN,dim=c(301,251,12))

count=1
for(r in 1:3)
for(n in 1:4)
{
W_ufile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF10_R",r,"_",name2[n],"_uwnd.nc",sep="")
W_pfile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF10_R",r,"_",name2[n],"_precip.nc",sep="")
corrlist[count,1]=paste(name2[n]," R",r,sep="")

fig2a=paste("Rain_corrUwind_sig_",name[n],"_WRF10R",r,"_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_",name[n],"_WRF10R",r,"_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_",name[n],"_WRF10R",r,"_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_",name[n],"_WRF10R",r,"_warm.tiff",sep="")
print(name[n])

f1=open.nc(W_ufile)
time=var.get.nc(f1,"time")
year=floor(time/100)
month=time%%100
time=cbind(year,month)
Uave=apply(var.get.nc(f1,'uwnd250',c(9,5,1),c(1,5,length(time[,1])),unpack=T),2,mean,na.rm=T)

f1=open.nc(W_pfile)
lat=var.get.nc(f1,"lat10")
lon=var.get.nc(f1,"lon10")
rain=var.get.nc(f1,"prcp10")

#Finally, let's do the two correlation plots
years=seq(min(time[,1]),max(time[,1]))
years=cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
coolM=matrix(0,length(years[,1]),2)
coolR<-warmR<-array(0,c(301,251,length(years[,1])))
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

coolC<-warmC<-matrix(0,301,251)
for(i in 1:301)
  for(j in 1:251)
  {
    coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
    warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
  }
coolA[,,count]=coolC
warmA[,,count]=warmC

corrlist[count,2]=mean(coolC*maskS,na.rm=T)
corrlist[count,3]=mean(coolC*maskN,na.rm=T)
corrlist[count,4]=cor(coolM[,1],years[,2],use="pairwise.complete.obs")
corrlist[count,5]=cor(coolM[,2],years[,2],use="pairwise.complete.obs")

i=which(lat==-35)
j=which(lon==145)
corrlist[count,6]=coolC[j,i]
i=which(lat==-30)
j=which(lon==152.5)
corrlist[count,7]=coolC[j,i]
count=count+1

coolC[coolC>=0.69]=0.69
coolC[coolC<=(-0.69)]=-0.69
warmC[warmC>=0.69]=0.69
warmC[warmC<=(-0.69)]=-0.69

tiff(file=fig2, height=500, width=800)
image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=800)
image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

coolC[abs(coolC)<=0.27]=NaN
warmC[abs(warmC)<=0.27]=NaN

tiff(file=fig2a, height=500, width=800)
image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3a, height=500, width=800)
image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

}

### Plot various medians

Cnames=data.frame(CMIP=rep("aaa",12),R=rep("aaa",12), stringsAsFactors=FALSE)
for(i in 1:12){
  x=strsplit(corrlist[i,1]," ")
  Cnames[i,]=x[[1]]
}

cool2=apply(coolA,c(1,2),median)
warm2=apply(warmA,c(1,2),median)
fig2a=paste("Rain_corrUwind_sig_CMIP3_WRF10_median_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_CMIP3_WRF10_median_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_CMIP3_WRF10_median_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_CMIP3_WRF10_median_warm.tiff",sep="")

cool2[cool2>=0.69]=0.69
cool2[cool2<=(-0.69)]=-0.69
warm2[warm2>=0.69]=0.69
warm2[warm2<=(-0.69)]=-0.69

tiff(file=fig2, height=500, width=600, pointsize=20)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=600, pointsize=20)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

cool2[abs(cool2)<=0.27]=NaN
warm2[abs(warm2)<=0.27]=NaN

tiff(file=fig2a, height=500, width=600, pointsize=20)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3a, height=500, width=600, pointsize=20)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

for(n in 1:4)
{
  I=which(Cnames[,1]==name2[n])

cool2=apply(coolA[,,I],c(1,2),median)
warm2=apply(warmA[,,I],c(1,2),median)
fig2a=paste("Rain_corrUwind_sig_",name[n],"_WRF10_median_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_",name[n],"_WRF10_median_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_",name[n],"_WRF10_median_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_",name[n],"_WRF10_median_warm.tiff",sep="")

cool2[cool2>=0.69]=0.69
cool2[cool2<=(-0.69)]=-0.69
warm2[warm2>=0.69]=0.69
warm2[warm2<=(-0.69)]=-0.69

tiff(file=fig2, height=500, width=600)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=600)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

coolC[abs(coolC)<=0.27]=NaN
warmC[abs(warmC)<=0.27]=NaN

tiff(file=fig2a, height=500, width=600)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3a, height=500, width=600)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
}

rlist=c("R1","R2","R3")
for(r in 1:3)
{
  I=which(Cnames[,2]==rlist[r])
  
  cool2=apply(coolA[,,I],c(1,2),median)
  warm2=apply(warmA[,,I],c(1,2),median)
  fig2a=paste("Rain_corrUwind_sig_CMIP3_WRF10",rlist[r],"_median_cool.tiff",sep="")
  fig3a=paste("Rain_corrUwind_sig_CMIP3_WRF10",rlist[r],"_median_warm.tiff",sep="")
  fig2=paste("Rain_corrUwind_CMIP3_WRF10",rlist[r],"_median_cool.tiff",sep="")
  fig3=paste("Rain_corrUwind_CMIP3_WRF10",rlist[r],"_median_warm.tiff",sep="")
  
  cool2[cool2>=0.69]=0.69
  cool2[cool2<=(-0.69)]=-0.69
  warm2[warm2>=0.69]=0.69
  warm2[warm2<=(-0.69)]=-0.69
  
  tiff(file=fig2, height=500, width=800)
  image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=fig3, height=500, width=800)
  image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  
  coolC[abs(coolC)<=0.27]=NaN
  warmC[abs(warmC)<=0.27]=NaN
  
  tiff(file=fig2a, height=500, width=800)
  image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=fig3a, height=500, width=800)
  image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}


##Our favourite corr plots - but more complex now
data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/25-35/newcomp_v3.csv")
##Only want ONE AWAP version - currently have 2.5 deg (1), 0.5 deg (2), 0.05 deg (3) and 0.05 deg landmask (4)
##Corrs are basically the same, except for 2.5 SW slightly weaker (by ~0.05)
##Let's use full 5 km landmask
data=data[-c(1,2,3),]
cmip=apply(data[2:38,5:10],2,median)

data1=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/20c3m/newcomp_wrf50_v3.csv")
data1=data1[-c(1,2,3),]
data2=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/20c3m/newcomp_wrf10_v3.csv")
data2=data2[-c(1,2,3),]

tiff(file="~/Documents/Data/CMIP5/CMIP3/Corrs_scatter_corr_point_v4_both.tiff",height=600,width=600,pointsize=20)
a=9
b=8
plot(data1[2:5,a],data1[2:5,b],type="p",col="gray50",pch=19,xlab="Coastal correlation",ylab="Inland correlation",xlim=c(-0.8,0.2),ylim=c(-0.8,0.8))
abline(h=0,col="grey")
abline(v=0,col="grey")
points(data1[6:17,a],data1[6:17,b],col="gray65",pch=17)
points(data2[6:17,a],data2[6:17,b],col="gray80",pch=15)
points(data1[1,a],data1[1,b],col="black",pch=4,cex=2,lwd=4)
points(median(data1[2:5,a]),median(data1[2:5,b]),col="gray50",pch=1,cex=2,lwd=4)
points(median(data1[6:17,a]),median(data1[6:17,b]),col="gray65",pch=2,cex=2,lwd=4)
points(median(data2[6:17,a]),median(data2[6:17,b]),col="gray80",pch=0,cex=2,lwd=4)
points(cmip[2],cmip[1],col="black",pch=1,cex=2,lwd=4)
legend("bottomright",c("AWAP","CMIP5","CMIP3","WRF50","WRF10"),
       pch=c(4,1,1,2,0),col=c("black","black","gray50","gray65","gray80"),pt.lwd=2)
dev.off()


