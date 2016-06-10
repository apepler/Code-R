##Notes - 
##ECHAM5 precip has time in MODEL years, need to subtract 400 years
##MIROC3.2 has years 1-100, rather than an actual date. Year 1 should be 1900.
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

loc="N"

if(loc=="S")
{
  a=21
  b=49
} else {
  a=23
  b=47
}

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

#MPI - have to fix for wrong dates
name=c("MPI_ECHAM5","CSIRO Mk3.0","MIROC3.2","CGCM3.1")
name1=c("mpi-echam5","csiromk3","miroc3.2","cgcm3.1")
Wfile=c("/srv/ccrc/data23/z3478332/CMIP3/WRF_echam5_monthly.nc",
        "/srv/ccrc/data23/z3478332/CMIP3/WRF_csiromk3_monthly.nc",
        "/srv/ccrc/data23/z3478332/CMIP3/WRF_miroc3.2_monthly.nc",
        "/srv/ccrc/data23/z3478332/CMIP3/WRF_cccma_monthly.nc")
date1=c(186001,187101,185001,185001)
date2=c(210012,200012,200012,200012)
date1a=c(191001,193101,185001,185001)
date2a=c(200912,200012,200012,200012)

#coolC<-warmC<-array(NaN,c(144,73,4))


lat=seq(-90,90,2.5)
lon=seq(0,357.5,2.5)
#reglats=c(rep(-32.5,7),rep(-35,5),rep(-37.5,3))
#reglons=c(seq(130,145,2.5),seq(135,145,2.5),seq(140,145,2.5))
reglats=c(rep(-32.5,4),rep(-35,4),rep(-37.5,3))
reglons=c(seq(135,142.5,2.5),seq(137.5,145,2.5),seq(140,145,2.5))

maskS<-maskN<-matrix(NaN,144,73)
for(i in 1:15)
{
  I=which(lat==reglats[i])
  J=which(lon==reglons[i])
  maskS[J,I]=1
}
reglats2=seq(-25,-32.5,-2.5)
for(i in 1:4)
{
  I=which(lon==152.5)
  J=which(lat==reglats2[i])
  maskN[I,J]=1
}

data=matrix(0,4,4)

coolC<-warmC<-array(0,c(144,73,4))
for(n in 1:4)
{
  #Rfile=paste("/srv/ccrc/data34/z3478332/CMIP3/precip/prcp_picntrl_",name1[n],"_2.5deg.nc",sep="")
  #Ufile=paste("/srv/ccrc/data34/z3478332/CMIP3/uw850/hus_850hPa_picntrl_",name1[n],"_2.5deg.nc",sep="")
  Rfile=paste("/srv/ccrc/data34/z3478332/CMIP3/precip/prcp_20c3m_",name1[n],"_2.5deg.nc",sep="")
  Ufile=paste("/srv/ccrc/data34/z3478332/CMIP3/uw850/hus_850hPa_20c3m_",name1[n],"_2.5deg.nc",sep="")
  # fig1=paste("Uwind_",name[n],".tiff",sep="")
  # fig2a=paste("Rain_corrUwind_sig_",name[n],"_cool.tiff",sep="")
  # fig3a=paste("Rain_corrUwind_sig_",name[n],"_warm.tiff",sep="")
  # fig2=paste("Rain_corrUwind_",name[n],"_cool.tiff",sep="")
  # fig3=paste("Rain_corrUwind_",name[n],"_warm.tiff",sep="")
  
  f1=open.nc(Rfile)
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  rain=var.get.nc(f1,"pr")
  time=timeREF[(timeREF[,3]>=date1[n] & timeREF[,3]<=date2[n]),1:2]
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
  time=timeREF[(timeREF[,3]>=date1a[n] & timeREF[,3]<=date2a[n]),1:2]
  I=which(time[,1]>=1950 & time[,1]<=2005)
  Uave=Uave[I]
  time=time[I,]
  
  #Finally, let's do the two correlation plots
  years=seq(min(time[,1]),max(time[,1]))
  coolM<-years<-cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
  coolR<-warmR<-array(0,c(144,73,length(years[,1])))
  for(i in 1:length(years[,1]))
  {
    I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
    years[i,2]=mean(Uave[I])
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    #  coolM[i,1]=mean(coolR[,,i]*maskS,na.rm=T)
    #  coolM[i,2]=mean(coolR[,,i]*maskN,na.rm=T)
    I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
    years[i,3]=mean(Uave[I])
    warmR[,,i]=apply(rain[,,I],c(1,2),sum)
  }
  
  #coolC=matrix(0,144,73)
  for(i in 1:144)
    for(j in 1:73)
    {
      coolC[i,j,n]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
      warmC[i,j,n]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
    }
  
  # data[n,1]=mean(coolC*maskS,na.rm=T)
  # data[n,2]=mean(coolC*maskN,na.rm=T)
  # data[n,3]=cor(coolM[,1],years[,2],use="pairwise.complete.obs")
  # data[n,4]=cor(coolM[,2],years[,2],use="pairwise.complete.obs")
  
}

# coolC<-warmC<-matrix(0,144,73)
# for(i in 1:144)
#   for(j in 1:73)
#   {
#     coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
#     warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
#   }
# i=which(lat==-35)
# j=which(lon==145)
# coolC[j,i]
# i=which(lat==-30)
# j=which(lon==152.5)
# coolC[j,i]
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
# 
# 
# ##Need to have NCEP for comparison! 
# f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
# timeN=var.get.nc(f1,'time')
# timeN=as.Date((timeN/24-2),origin="0001-01-01")
# 
# timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
# uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(144,73,1,length(timeN[,1])),unpack=T)
# NCEP=apply(var.get.nc(f1,'uwnd',c(61,b,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)
# 
# ##So, for colocated years 1948-2005, want the mean wind for each month & plot
# f1=open.nc(Wfile)
# uwind=var.get.nc(f1,"uw850")
# timeW=var.get.nc(f1,"time")
# year=floor(timeW/100)
# month=timeW%%100
# timeW=cbind(year,month)
# UaveW=apply(uwind[61,a:(a+4),],2,mean)
# 
# Ucomp=matrix(0,12,4)
# for(i in 1:12)
# {
#   Ucomp[i,1]=i
#   I=which(timeN[,1]<=2005 & timeN[,2]==i)
#   Ucomp[i,2]=mean(NCEP[I])
#   I=which(time[,1]>=1948 & time[,1]<=2005 & time[,2]==i)
#   Ucomp[i,3]=mean(Uave[I])
#   I=which(timeW[,2]==i)
#   Ucomp[i,4]=mean(UaveW[I])
#   }
# 
# mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# yl=c(floor(min(Ucomp[,2:3])),ceiling(max(Ucomp[,2:3])))
# tiff(file=fig1, height=500, width=800)
# plot(seq(1:12),Ucomp[,2],ty="l",lwd=2,col="blue",xlab="",ylab="",axes=F,ylim=yl)
# axis(1, at = seq(1:12), labels = c(mnames),cex.axis=1.5) 
# axis(2, at = seq(yl[1],yl[2]),cex.axis=1.5)
# title(main="Average Uwind (m/s)",cex.main=2)
# lines(seq(1:12),Ucomp[,3],ty="l",lwd=2,col="red")
# lines(seq(1:12),Ucomp[,4],ty="l",lwd=2,col="darkgreen")
# abline(h=0,col="gray",lty=2)
# legend("topleft",legend=c("NCEP",name,"WRF version"),col=c("blue","red","darkgreen"),lwd=2,bty="n",cex=1.5)
# dev.off()
#}

medW=apply(warmC,c(1,2),median,na.rm=T)
medW[medW>=0.69]=0.69
medW[medW<=(-0.69)]=-0.69
#medW[abs(medW)<=0.27]=NaN
medC=apply(coolC,c(1,2),median,na.rm=T)
medC[medC>=0.69]=0.69
medC[medC<=(-0.69)]=-0.69
#medC[abs(medC)<=0.27]=NaN

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
# tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_CMIP3_median_cool.tiff",sep=""), height=500, width=600,pointsize=20)
# image.plot(lon,lat,medC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_CMIP3_median_warm.tiff",sep=""), height=500, width=600,pointsize=20)
# image.plot(lon,lat,medW,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()

# medW[abs(medW)<=0.27]=NaN
# medC[abs(medC)<=0.27]=NaN
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


