##Notes - 
##ECHAM5 precip has time in MODEL years, need to subtract 400 years
##MIROC3.2 has years 1-100, rather than an actual date. Year 1 should be 1900.

loc="N"

if(loc=="S")
{
  a=21
  b=49
} else {
  a=23
  b=47
}

# #MPI - have to fix for wrong dates
# name="MPI_ECHAM5"
# name1="mpi-echam5"

##CSIRO - Have to fix for wrong dates
#name="CSIRO Mk3.0"
#name1="csiromk3"

#MIROC3.2 - Have to fix for wrong dates
name="MIROC3.2"
name1="miroc3.2"

# #CCCMA CGCM3.1
# name="CGCM3.1"
# name1="cgcm3.1"

setwd("~/Documents/Data/CMIP5/CMIP3")
Rfile=paste("/srv/ccrc/data23/z3478332/CMIP3/precip/prcp_picntrl_",name1,"_2.5deg.nc",sep="")
Ufile=paste("/srv/ccrc/data23/z3478332/CMIP3/uw850/hus_850hPa_picntrl_",name1,"_2.5deg.nc",sep="")
fig1=paste("Uwind_",name,".tiff",sep="")
fig2a=paste("Rain_corrUwind_sig_",name,"_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_",name,"_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_",name,"_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_",name,"_warm.tiff",sep="")

library(RNetCDF)
f1=open.nc(Rfile)
time=var.get.nc(f1,"time")
ori=att.get.nc(f1,"time","units")
ori=unlist(strsplit(ori, split=" "))[3]
time2=as.Date(time,origin=ori)
time=data.frame(Year=as.numeric(format(time2,"%Y")),Month=as.numeric(format(time2,"%m")))
lat=var.get.nc(f1,"lat")
lon=var.get.nc(f1,"lon")
rain=var.get.nc(f1,"pr")
##Fix for CSIRO
# yy=seq(1871,2250)
# mm=seq(1,12)
# n=1
# for(i in 1:380)
#   for(j in 1:12)
#   {
#     time[n,1]=yy[i]
#     time[n,2]=mm[j]
#     n=n+1
#   }

# #Fix for MPI
# time[,1]=time[,1]-400

# #Fix for MIROC
# time[,1]=time[,1]+1900

#Restrict to 1950-2005
I=which(time[,1]>=1950 & time[,1]<=2005)
rain=rain[,,I]
time=time[I,]

I=which(time[,2] %in% c(1,3,5,7,8,10,12))
rain[,,I]=rain[,,I]*31*60*60*24
I=which(time[,2] %in% c(4,6,9,11))
rain[,,I]=rain[,,I]*30*60*60*24
I=which(time[,2]==2)
rain[,,I]=rain[,,I]*28*60*60*24

##Uwind - is orig on different time period
f1=open.nc(Ufile)
u=var.get.nc(f1,"ua")
Uave=apply(u[61,a:(a+4),],2,mean)
time=var.get.nc(f1,"time")
ori=att.get.nc(f1,"time","units")
ori=unlist(strsplit(ori, split=" "))[3]
time2=as.Date(time,origin=ori)
time=data.frame(Year=as.numeric(format(time2,"%Y")),Month=as.numeric(format(time2,"%m")))
# # Fix for CSIRO
# yy=seq(1931,2030)
# mm=seq(1,12)
# n=1
# for(i in 1:100)
#   for(j in 1:12)
#   {
#     time[n,1]=yy[i]
#     time[n,2]=mm[j]
#     n=n+1
#   }

# #Fix for MIROC3.2
# time[,1]=time[,1]+1900

I=which(time[,1]>=1950 & time[,1]<=2005)
Uave=Uave[I]
time=time[I,]

#Finally, let's do the two correlation plots
years=seq(min(time[,1]),max(time[,1]))
years=cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
coolR<-warmR<-array(0,c(144,73,length(years[,1])))
for(i in 1:length(years[,1]))
{
  I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
  years[i,2]=mean(Uave[I])
  coolR[,,i]=apply(rain[,,I],c(1,2),sum)
  I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
  years[i,3]=mean(Uave[I])
  warmR[,,i]=apply(rain[,,I],c(1,2),sum)
}

coolC<-warmC<-matrix(0,144,73)
for(i in 1:144)
  for(j in 1:73)
  {
    coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
    warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
  }
i=which(lat==-35)
j=which(lon==145)
coolC[j,i]
i=which(lat==-30)
j=which(lon==152.5)
coolC[j,i]

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


##Need to have NCEP for comparison! 
f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")

timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(144,73,1,length(timeN[,1])),unpack=T)
NCEP=apply(var.get.nc(f1,'uwnd',c(61,b,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)

##So, for colocated years 1948-2005, want the mean wind for each month & plot
Ucomp=matrix(0,12,3)
for(i in 1:12)
{
  Ucomp[i,1]=i
  I=which(timeN[,1]<=2005 & timeN[,2]==i)
  Ucomp[i,2]=mean(NCEP[I])
  I=which(time[,1]>=1948 & time[,1]<=2005 & time[,2]==i)
  Ucomp[i,3]=mean(Uave[I])
  }

mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
yl=c(floor(min(Ucomp[,2:3])),ceiling(max(Ucomp[,2:3])))
tiff(file=fig1, height=500, width=800)
plot(seq(1:12),Ucomp[,2],ty="l",lwd=2,col="blue",xlab="",ylab="",axes=F,ylim=yl)
axis(1, at = seq(1:12), labels = c(mnames),cex.axis=1.5) 
axis(2, at = seq(yl[1],yl[2]),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
lines(seq(1:12),Ucomp[,3],ty="l",lwd=2,col="red")
abline(h=0,col="gray",lty=2)
legend("topleft",legend=c("NCEP",name),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
dev.off()


