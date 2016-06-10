
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
  fig1=paste("Uwind_",name[n],"_WRF50_Rcomp.tiff",sep="")
  
  f1=open.nc(C_ufile)
  u=var.get.nc(f1,"ua")
  CMIP=apply(u[61,a:(a+4),],2,mean)
  timeC=timeREF[(timeREF[,3]>=date1a[n] & timeREF[,3]<=date2a[n]),1:2]
  
  Ucomp=matrix(0,12,6)
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
    W_ufile=paste("/srv/ccrc/data34/z3478332/CMIP3/WRF_R",r,"_",name2[n],"_uwnd.nc",sep="")
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
  
  mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  yl=c(-10,10)
  tiff(file=fig1, height=500, width=800)
  plot(seq(1:12),Ucomp[,2],ty="l",lwd=4,col="blue",xlab="",ylab="",axes=F,ylim=yl)
  axis(1, at = seq(1:12), labels = c(mnames),cex.axis=1.5) 
  axis(2, at = seq(yl[1],yl[2]),cex.axis=1.5)
  title(main="Average Uwind (m/s)",cex.main=2)
  lines(seq(1:12),Ucomp[,3],ty="l",lwd=4,col="red")
  lines(seq(1:12),Ucomp[,4],ty="l",lwd=4,col="magenta3")
  lines(seq(1:12),Ucomp[,5],ty="l",lwd=4,col="purple")
  lines(seq(1:12),Ucomp[,6],ty="l",lwd=4,col="purple4")
  abline(h=0,col="gray",lty=2)
  legend("topleft",legend=c("NCEP",name[n],"WRF R1","WRF R2","WRF R3"),
         col=c("blue","red","magenta3","purple","purple4"),lwd=4,bty="n",cex=1.5,ncol=2)
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

Ucomp2=matrix(0,12,13)
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
    W_ufile=paste("/srv/ccrc/data34/z3478332/CMIP3/WRF_R",r,"_",name2[n],"_uwnd.nc",sep="")
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
    count=count+1
  }
}

Ucomp1[,7]=apply(Ucomp1[,3:6],1,median)
Ucomp2[,13]=apply(Ucomp2[,1:12],1,median)
  
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
yl=c(-10,10)
tiff(file="Uwind_CMIP3_WRF_all.tiff", height=500, width=800)
  plot(seq(1:12),Ucomp1[,2],ty="l",lwd=4,col="blue",xlab="",ylab="",axes=F,ylim=yl)
  axis(1, at = seq(1:12), labels = c(mnames),cex.axis=1.5) 
  axis(2, at = seq(yl[1],yl[2]),cex.axis=1.5)
  title(main="Average Uwind (m/s)",cex.main=2)
abline(h=0,col="gray",lty=2)
for(i in 3:6) lines(seq(1:12),Ucomp1[,i],ty="l",col="red",lty=2)
for(i in 1:12) lines(seq(1:12),Ucomp2[,i],ty="l",col="purple",lty=2)
lines(seq(1:12),Ucomp1[,2],ty="l",lwd=4,col="blue")
lines(seq(1:12),Ucomp1[,7],ty="l",lwd=4,col="red")
lines(seq(1:12),Ucomp2[,13],ty="l",lwd=4,col="purple")
legend("topleft",legend=c("NCEP","CMIP3 median","WRF50 median"),
         col=c("blue","red","purple"),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()

  
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

corrlist=data.frame(Names=rep("aaa",12),Avecorr_SW=rep(NaN,12),Avecorr_NE=rep(NaN,12),
                    Corrave_SW=rep(NaN,12),Corrave_NE=rep(NaN,12),
                    Point_SW=rep(NaN,12),Point_NE=rep(NaN,12), stringsAsFactors=FALSE)
coolA<-warmA<-array(NaN,dim=c(151,101,12))

count=1
for(r in 1:3)
for(n in 1:4)
{
W_ufile=paste("/srv/ccrc/data34/z3478332/CMIP3/WRF_R",r,"_",name2[n],"_uwnd.nc",sep="")
W_pfile=paste("/srv/ccrc/data34/z3478332/CMIP3/WRF_R",r,"_",name2[n],"_precip.nc",sep="")
#corrlist[count,1]=paste(name2[n]," R",r,sep="")

fig2a=paste("Rain_corrUwind_sig_",name[n],"_WRFR",r,"_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_",name[n],"_WRFR",r,"_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_",name[n],"_WRFR",r,"_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_",name[n],"_WRFR",r,"_warm.tiff",sep="")
print(name[n])

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
#  coolM[i,1]=mean(coolR[,,i]*maskS,na.rm=T)
#  coolM[i,2]=mean(coolR[,,i]*maskN,na.rm=T)
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

### Plot various medians

Cnames=data.frame(CMIP=rep("aaa",12),R=rep("aaa",12), stringsAsFactors=FALSE)
for(i in 1:12){
  x=strsplit(corrlist[i,1]," ")
  Cnames[i,]=x[[1]]
}

cool2=apply(coolA,c(1,2),median)
warm2=apply(warmA,c(1,2),median)
fig2a=paste("Rain_corrUwind_sig_CMIP3_WRF_median_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_CMIP3_WRF_median_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_CMIP3_WRF_median_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_CMIP3_WRF_median_warm.tiff",sep="")

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
fig2a=paste("Rain_corrUwind_sig_",name[n],"_WRF_median_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_",name[n],"_WRF_median_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_",name[n],"_WRF_median_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_",name[n],"_WRF_median_warm.tiff",sep="")

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
  fig2a=paste("Rain_corrUwind_sig_CMIP3_",rlist[r],"_median_cool.tiff",sep="")
  fig3a=paste("Rain_corrUwind_sig_CMIP3_",rlist[r],"_median_warm.tiff",sep="")
  fig2=paste("Rain_corrUwind_CMIP3_",rlist[r],"_median_cool.tiff",sep="")
  fig3=paste("Rain_corrUwind_CMIP3_",rlist[r],"_median_warm.tiff",sep="")
  
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
##Let's stick with 2.5 for now
data=data[-c(1,2,3),]
cmip=apply(data[2:38,5:10],2,median)

data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/20c3m/newcomp_wrf50_v3.csv")
##Use the 50 km AWAP version
data=data[-c(1,2,3),]
tiff(file="~/Documents/Data/CMIP5/CMIP3/Corrs_scatter_avecorr_ref_v4.tiff",height=600,width=600,pointsize=20)
a=5
b=4
plot(data[2:5,a],data[2:5,b],type="p",col="gray50",pch=19,xlab="Coastal correlation",ylab="Inland correlation",xlim=c(-0.7,0.2),ylim=c(-0.7,0.7))
abline(h=0,col="grey")
abline(v=0,col="grey")
points(data[6:17,a],data[6:17,b],col="grey",pch=17)
points(data[1,a],data[1,b],col="black",pch=4,cex=2,lwd=4)
points(median(data[2:5,a]),median(data[2:5,b]),col="gray50",pch=1,cex=2,lwd=4)
points(median(data[6:17,a]),median(data[6:17,b]),col="grey",pch=2,cex=2,lwd=4)
points(cmip[2],cmip[1],col="black",pch=1,cex=2,lwd=4)
legend("bottomright",c("AWAP","CMIP5","CMIP3","WRF"),pch=c(4,1,1,2),col=c("black","black","gray50","grey"),pt.lwd=2)
dev.off()


#########Baby WRF10 play - using WRF50 for wind




#############
##Wind comps, for reference

library(RNetCDF)
f1=open.nc("/srv/ccrc/data34/z3478332/CMIP3/WRF_R2_miroc3.2_uwnd.nc")
time=var.get.nc(f1,"time")
date=cbind(floor(time/100),time%%100)
Uave=matrix(NaN,length(time),5)

uwind=var.get.nc(f1,"uwnd250")
a=23
Uave[,1]=apply(uwind[61,a:(a+4),],2,mean)

uwind=var.get.nc(f1,"uwnd50")
lat=var.get.nc(f1,"lat50")
lon=var.get.nc(f1,"lon50")
I=which(lat>=-36.25 & lat<=-23.75)
J=which(lon>=148.75 & lon<=151.25)
Uave[,2]=apply(uwind[J,I,],3,mean)

f1=open.nc("/srv/ccrc/data34/z3478332/CMIP3/WRF10_R2_miroc3.2_uwnd.nc")
uwind=var.get.nc(f1,"uwnd250")
lat=var.get.nc(f1,"lat")
lon=var.get.nc(f1,"lon")
I=which(lat>=-36.25 & lat<=-23.75)
J=which(lon>=148.75 & lon<=151.25)
Uave[,3]=apply(uwind[J,I,],2,mean)

uwind=var.get.nc(f1,"uwnd50")
lat=var.get.nc(f1,"lat50")
lon=var.get.nc(f1,"lon50")
I=which(lat>=-36.25 & lat<=-23.75)
J=which(lon>=148.75 & lon<=151.25)
Uave[,4]=apply(uwind[J,I,],3,mean)

uwind=var.get.nc(f1,"uwnd10")
lat=var.get.nc(f1,"lat10")
lon=var.get.nc(f1,"lon10")
I=which(lat>=-36.25 & lat<=-23.75)
J=which(lon>=148.75 & lon<=151.25)
Uave[,5]=apply(uwind[J,I,],3,mean)

Ucomp=matrix(0,13,6)
Ucomp[,1]=seq(0,12)
Ucomp[1,2:6]=apply(Uave,2,mean)

for(i in 1:12)
{
  I=which(date[,2]==i)
  Ucomp[i+1,2:6]=apply(Uave[I,],2,mean)
}

Ucorr=matrix(0,5,5)
for(i in 1:5)
  for(j in 1:5)
    Ucorr[i,j]=cor(Uave[,i],Uave[,j])

a=23
b=47
f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")
timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
NCEP=apply(var.get.nc(f1,'uwnd',c(61,b,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)
N=rep(0,12)
for(i in 1:12)
{
  I=which(timeN[,1]<=2009 & timeN[,1]>=1990 & timeN[,2]==i)
  N[i]=mean(NCEP[I])
}

plot(Ucomp[2:13,1],N,type="l",col="black",ylim=c(-10,15),lwd=2,xlab="Month",ylab="Uwind")
for(i in 2:6) lines(Ucomp[2:13,1],Ucomp[2:13,i],type="l",col=i,lwd=2)
legend("bottomright",c("NCEP1","WRF50 2.5 degree","WRF50 0.5 degree","WRF10 2.5 degree","WRF10 0.5 degree","WRF10 0.1 degree"),
       col=c("black",seq(2,6)),lwd=2,ncol=2)

