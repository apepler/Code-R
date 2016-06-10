
##Currently edited to produce figs for WRF1010, but quick to change refs.
rm(list=ls())
setwd("~/Documents/Data/CMIP5/CMIP3")
library(RNetCDF)
library(fields)
library(akima)
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

###50 km mask??

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

f1=open.nc("/srv/ccrc/data13/z3393020/Analyses/share/geo_em_files/geo_em.d02.narclim.nc")
lat1=var.get.nc(f1,"XLAT_M")
lon1=var.get.nc(f1,"XLONG_M")
latt=as.vector(lat1)
lont=as.vector(lon1)
maskS1<-maskN1<-matrix(NaN,325,200)
I=which(lat1>(-33.75) & lat1<(-23.75) & lon1>151.25 & lon1<153.75)
maskN1[I]=1
I=which(lat1>(-38.75) & lat1<(-36.25) & lon1<146.75 & lon1>(131.25-(10/7.5)*(31.25+lat1)))
maskS1[I]=1
I=which(lat1>=(-36.25) & lat1<(-31.25) & lon1<(138.75-(8/5)*(31.25+lat1)) & lon1>(131.25-(10/7.5)*(31.25+lat1)))
maskS1[I]=1

corrlist<-corrlist1<-data.frame(Names=rep("aaa",12),Avecorr_SW=rep(NaN,12),Avecorr_NE=rep(NaN,12),
                    Corrave_SW=rep(NaN,12),Corrave_NE=rep(NaN,12),
                    Point_SW=rep(NaN,12),Point_NE=rep(NaN,12), stringsAsFactors=FALSE)
coolA1<-warmA1<-array(NaN,dim=c(dim(lat1),12))
coolA<-warmA<-array(NaN,dim=c(301,251,12))

count=1
for(r in 1:3)
for(n in 1:4)
{
W_ufile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF10_R",r,"_",name2[n],"_uwnd.nc",sep="")
W_pfile=paste("/srv/ccrc/data36/z3478332/CMIP3/WRF10_R",r,"_",name2[n],"_precip_forcdo.nc",sep="")
corrlist[count,1]=paste(name2[n]," R",r,sep="")

fig2a=paste("Rain_corrUwind_sig_",name[n],"_WRF10R",r,"_cool_postgrid.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_",name[n],"_WRF10R",r,"_warm_postgrid.tiff",sep="")
fig2=paste("Rain_corrUwind_",name[n],"_WRF10R",r,"_cool_postgrid.tiff",sep="")
fig3=paste("Rain_corrUwind_",name[n],"_WRF10R",r,"_warm_postgrid.tiff",sep="")
fig4a=paste("Rain_corrUwind_sig_",name[n],"_WRF10R",r,"_cool_postgrid_small.tiff",sep="")
fig5a=paste("Rain_corrUwind_sig_",name[n],"_WRF10R",r,"_warm_postgrid_small.tiff",sep="")
fig4=paste("Rain_corrUwind_",name[n],"_WRF10R",r,"_cool_postgrid_small.tiff",sep="")
fig5=paste("Rain_corrUwind_",name[n],"_WRF10R",r,"_warm_postgrid_small.tiff",sep="")
print(name[n])

f1=open.nc(W_ufile)
time=var.get.nc(f1,"time")
year=floor(time/100)
month=time%%100
time=cbind(year,month)
Uave=apply(var.get.nc(f1,'uwnd250',c(9,5,1),c(1,5,length(time[,1])),unpack=T),2,mean,na.rm=T)

f1=open.nc(W_pfile)
rain=var.get.nc(f1,"prcp")

#Finally, let's do the two correlation plots
years=seq(min(time[,1]),max(time[,1]))
years=cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
coolM1<-coolM<-matrix(0,length(years[,1]),2)
coolR<-warmR<-array(0,c(dim(lat1),length(years[,1])))
for(i in 1:length(years[,1]))
{
  I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
  years[i,2]=mean(Uave[I])
  coolR[,,i]=apply(rain[,,I],c(1,2),sum)
#  coolM1[i,1]=mean(coolR[,,i]*maskS1,na.rm=T)
#  coolM1[i,2]=mean(coolR[,,i]*maskN1,na.rm=T)
#  rr=as.vector(coolR[,,i])
#  rr[which(is.na(rr))]=0
#  rr2=interp(lont,latt,rr,lon,lat)
#  coolM[i,1]=mean(rr2$z*maskS,na.rm=T)
#  coolM[i,2]=mean(rr2$z*maskN,na.rm=T)
  I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
  years[i,3]=mean(Uave[I])
  warmR[,,i]=apply(rain[,,I],c(1,2),sum)
}

coolC<-warmC<-matrix(0,325,200)
for(i in 1:325)
  for(j in 1:200)
  {
    coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
    warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
  }
coolA1[,,count]=coolC
warmA1[,,count]=warmC

rr=as.vector(coolC)
rr[which(is.na(rr))]=0
rr2=interp(lont,latt,rr,lon,lat)
coolA[,,count]=rr2$z
rr=as.vector(warmC)
rr[which(is.na(rr))]=0
rr2=interp(lont,latt,rr,lon,lat)
warmA[,,count]=rr2$z

# corrlist[count,2]=mean(coolA[,,count]*maskS,na.rm=T)
# corrlist[count,3]=mean(coolA[,,count]*maskN,na.rm=T)
# corrlist[count,4]=cor(coolM[,1],years[,2],use="pairwise.complete.obs")
# corrlist[count,5]=cor(coolM[,2],years[,2],use="pairwise.complete.obs")
# i=which(lat==-35)
# j=which(lon==145)
# corrlist[count,6]=coolA[j,i,count]
# i=which(lat==-30)
# j=which(lon==152.5)
# corrlist[count,7]=coolA[j,i,count]
# 
# ##Version 2 - using original data (not regridded)
# corrlist1[count,2]=mean(coolA1[,,count]*maskS1,na.rm=T)
# corrlist1[count,3]=mean(coolA1[,,count]*maskN1,na.rm=T)
# corrlist1[count,4]=cor(coolM1[,1],years[,2],use="pairwise.complete.obs")
# corrlist1[count,5]=cor(coolM1[,2],years[,2],use="pairwise.complete.obs")
# coolC=coolA[,,count]
# warmC=warmA[,,count]
count=count+1

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
# 
# tiff(file=fig4, height=500, width=600)
# image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# tiff(file=fig5, height=500, width=600)
# image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
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
# tiff(file=fig4a, height=500, width=600)
# image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# tiff(file=fig5a, height=500, width=600)
# image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# 

}

save(coolA,coolA1,warmA,warmA1,file="WRFR10corrs.RData")
load("WRFR10corrs.RData")

### Plot various medians

Cnames=data.frame(CMIP=rep("aaa",12),R=rep("aaa",12), stringsAsFactors=FALSE)
for(i in 1:12){
  x=strsplit(corrlist[i,1]," ")
  Cnames[i,]=x[[1]]
}

cool2=apply(coolA,c(1,2),median)
warm2=apply(warmA,c(1,2),median)
fig2a=paste("Rain_corrUwind_sig_CMIP3_WRF10_median_cool_postregrid.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_CMIP3_WRF10_median_warm_postregrid.tiff",sep="")
fig2=paste("Rain_corrUwind_CMIP3_WRF10_median_cool_postregrid.tiff",sep="")
fig3=paste("Rain_corrUwind_CMIP3_WRF10_median_warm_postregrid.tiff",sep="")
fig4a=paste("Rain_corrUwind_sig_CMIP3_WRF10_median_cool_postregrid_small.tiff",sep="")
fig5a=paste("Rain_corrUwind_sig_CMIP3_WRF10_median_warm_postregrid_small.tiff",sep="")
fig4=paste("Rain_corrUwind_CMIP3_WRF10_median_cool_postregrid_small.tiff",sep="")
fig5=paste("Rain_corrUwind_CMIP3_WRF10_median_warm_postregrid_small.tiff",sep="")

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
tiff(file=fig4, height=500, width=600, pointsize=20)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig5, height=500, width=600, pointsize=20)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

fig4=paste("Rain_corrUwind_sigdot_CMIP3_WRF10_median_cool_postregrid_small_v2.tiff",sep="")
fig5=paste("Rain_corrUwind_sigdot_CMIP3_WRF10_median_warm_postregrid_small_v2.tiff",sep="")

sigmaskC=which(abs(cool2)>=0.45,arr.ind=T)
I=which(lon[sigmaskC[,1]]>155 | lon[sigmaskC[,1]]<135)
if(length(I)>1) sigmaskC=sigmaskC[-I,]
I=which(lat[sigmaskC[,2]]>(-23) | lat[sigmaskC[,2]]<(-39))
if(length(I)>1) sigmaskC=sigmaskC[-I,]
sigmaskW=which(abs(warm2)>=0.45,arr.ind=T)
I=which(lon[sigmaskW[,1]]>155 | lon[sigmaskW[,1]]<135)
if(length(I)>1) sigmaskW=sigmaskW[-I,]
I=which(lat[sigmaskW[,2]]>(-23) | lat[sigmaskW[,2]]<(-39))
if(length(I)>1) sigmaskW=sigmaskW[-I,]

I=which(lon[sigmaskC[,1]] %% 0.5 == 0 & lat[sigmaskC[,2]] %% 0.5 == 0)
if(length(I)>1) sigmaskC=sigmaskC[I,]
I=which(lon[sigmaskW[,1]] %% 0.5 == 0 & lat[sigmaskW[,2]] %% 0.5 == 0)
if(length(I)>1) sigmaskW=sigmaskW[I,]

tiff(file=fig4, height=1000, width=1100, pointsize=20)
image(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(135,155),ylim=c(-39,-23),cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon[sigmaskC[,1]],lat[sigmaskC[,2]],col="black",pch=19,cex=0.75,xlim=c(135,155),ylim=c(-39,-23))
dev.off()
tiff(file=fig5, height=1000, width=1100, pointsize=20)
image(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7) ,xlim=c(135,155),ylim=c(-39,-23),cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon[sigmaskW[,1]],lat[sigmaskW[,2]],col="black",pch=19,cex=0.75,xlim=c(135,155),ylim=c(-39,-23))
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
tiff(file=fig4a, height=500, width=600, pointsize=20)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig5a, height=500, width=600, pointsize=20)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

for(n in 1:4)
{
  I=which(Cnames[,1]==name2[n])
  
  cool2=apply(coolA[,,I],c(1,2),median)
  warm2=apply(warmA[,,I],c(1,2),median)
  fig2=paste("Rain_corrUwind_sigdot_",name[n],"_WRF10_median_cool_postregrid.tiff",sep="")
  fig3=paste("Rain_corrUwind_sigdot_",name[n],"_WRF10_median_warm_postregrid.tiff",sep="")
  fig4=paste("Rain_corrUwind_sigdot_",name[n],"_WRF10_median_cool_postregrid_small.tiff",sep="")
  fig5=paste("Rain_corrUwind_sigdot_",name[n],"_WRF10_median_warm_postregrid_small.tiff",sep="")
  
  cool2[cool2>=0.69]=0.69
  cool2[cool2<=(-0.69)]=-0.69
  warm2[warm2>=0.69]=0.69
  warm2[warm2<=(-0.69)]=-0.69
  
  sigmaskC=which(abs(cool2)>=0.45,arr.ind=T)
  I=which(lon[sigmaskC[,1]]>160 | lon[sigmaskC[,1]]<110)
  if(length(I)>1) sigmaskC=sigmaskC[-I,]
  I=which(lat[sigmaskC[,2]]>(-10) | lat[sigmaskC[,2]]<(-50))
  if(length(I)>1) sigmaskC=sigmaskC[-I,]
  sigmaskW=which(abs(warm2)>=0.45,arr.ind=T)
  I=which(lon[sigmaskW[,1]]>160 | lon[sigmaskW[,1]]<110)
  if(length(I)>1) sigmaskW=sigmaskW[-I,]
  I=which(lat[sigmaskW[,2]]>(-10) | lat[sigmaskW[,2]]<(-50))
  if(length(I)>1) sigmaskW=sigmaskW[-I,]
  
  tiff(file=fig2, height=1000, width=1100, pointsize=20)
  image(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(110,160),ylim=c(-50,-10),cex.axis=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  points(lon[sigmaskC[,1]],lat[sigmaskC[,2]],col="black",pch=19,cex=0.15,xlim=c(110,160),ylim=c(-50,-10))
  dev.off()
  tiff(file=fig3, height=1000, width=1100, pointsize=20)
  image(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(110,160),ylim=c(-50,-10),cex.axis=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  points(lon[sigmaskW[,1]],lat[sigmaskW[,2]],col="black",pch=19,cex=0.15,xlim=c(110,160),ylim=c(-50,-10))
  dev.off()
  
  I=which(lon[sigmaskC[,1]]>155 | lon[sigmaskC[,1]]<135)
  if(length(I)>1) sigmaskC=sigmaskC[-I,]
  I=which(lat[sigmaskC[,2]]>(-23) | lat[sigmaskC[,2]]<(-39))
  if(length(I)>1) sigmaskC=sigmaskC[-I,]
  I=which(lon[sigmaskW[,1]]>155 | lon[sigmaskW[,1]]<135)
  if(length(I)>1) sigmaskW=sigmaskW[-I,]
  I=which(lat[sigmaskW[,2]]>(-23) | lat[sigmaskW[,2]]<(-39))
  if(length(I)>1) sigmaskW=sigmaskW[-I,]
  
  tiff(file=fig4, height=1000, width=1100, pointsize=20)
  image(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(135,155),ylim=c(-39,-23),cex.axis=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  points(lon[sigmaskC[,1]],lat[sigmaskC[,2]],col="black",pch=19,cex=0.15,xlim=c(135,155),ylim=c(-39,-23))
  dev.off()
  tiff(file=fig5, height=1000, width=1100, pointsize=20)
  image(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7) ,xlim=c(135,155),ylim=c(-39,-23),cex.axis=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  points(lon[sigmaskW[,1]],lat[sigmaskW[,2]],col="black",pch=19,cex=0.15,xlim=c(135,155),ylim=c(-39,-23))
  dev.off()
}

rlist=c("R1","R2","R3")
for(r in 1:3)
{
  I=which(Cnames[,2]==rlist[r])
  
  cool2=apply(coolA[,,I],c(1,2),median)
  warm2=apply(warmA[,,I],c(1,2),median)
  fig2=paste("Rain_corrUwind_sigdot_CMIP3_",rlist[r],"_median_cool_postregrid.tiff",sep="")
  fig3=paste("Rain_corrUwind_sigdot_CMIP3_",rlist[r],"_median_warm_postregrid.tiff",sep="")
  fig4=paste("Rain_corrUwind_sigdot_CMIP3_",rlist[r],"_median_cool_postregrid_small.tiff",sep="")
  fig5=paste("Rain_corrUwind_sigdot_CMIP3_",rlist[r],"_median_warm_postregrid_small.tiff",sep="")
  
  cool2[cool2>=0.69]=0.69
  cool2[cool2<=(-0.69)]=-0.69
  warm2[warm2>=0.69]=0.69
  warm2[warm2<=(-0.69)]=-0.69
  
  sigmaskC=which(abs(cool2)>=0.45,arr.ind=T)
  I=which(lon[sigmaskC[,1]]>160 | lon[sigmaskC[,1]]<110)
  if(length(I)>1) sigmaskC=sigmaskC[-I,]
  I=which(lat[sigmaskC[,2]]>(-10) | lat[sigmaskC[,2]]<(-50))
  if(length(I)>1) sigmaskC=sigmaskC[-I,]
  sigmaskW=which(abs(warm2)>=0.45,arr.ind=T)
  I=which(lon[sigmaskW[,1]]>160 | lon[sigmaskW[,1]]<110)
  if(length(I)>1) sigmaskW=sigmaskW[-I,]
  I=which(lat[sigmaskW[,2]]>(-10) | lat[sigmaskW[,2]]<(-50))
  if(length(I)>1) sigmaskW=sigmaskW[-I,]
  
  tiff(file=fig2, height=1000, width=1100, pointsize=20)
  image(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(110,160),ylim=c(-50,-10),cex.axis=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  points(lon[sigmaskC[,1]],lat[sigmaskC[,2]],col="black",pch=19,cex=0.15,xlim=c(110,160),ylim=c(-50,-10))
  dev.off()
  tiff(file=fig3, height=1000, width=1100, pointsize=20)
  image(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(110,160),ylim=c(-50,-10),cex.axis=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  points(lon[sigmaskW[,1]],lat[sigmaskW[,2]],col="black",pch=19,cex=0.15,xlim=c(110,160),ylim=c(-50,-10))
  dev.off()
  
  
  I=which(lon[sigmaskC[,1]]>155 | lon[sigmaskC[,1]]<135)
  if(length(I)>1) sigmaskC=sigmaskC[-I,]
  I=which(lat[sigmaskC[,2]]>(-23) | lat[sigmaskC[,2]]<(-39))
  if(length(I)>1) sigmaskC=sigmaskC[-I,]
  I=which(lon[sigmaskW[,1]]>155 | lon[sigmaskW[,1]]<135)
  if(length(I)>1) sigmaskW=sigmaskW[-I,]
  I=which(lat[sigmaskW[,2]]>(-23) | lat[sigmaskW[,2]]<(-39))
  if(length(I)>1) sigmaskW=sigmaskW[-I,]
  
  tiff(file=fig4, height=1000, width=1100, pointsize=20)
  image(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(135,155),ylim=c(-39,-23),cex.axis=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  points(lon[sigmaskC[,1]],lat[sigmaskC[,2]],col="black",pch=19,cex=0.15,xlim=c(135,155),ylim=c(-39,-23))
  dev.off()
  tiff(file=fig5, height=1000, width=1100, pointsize=20)
  image(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7) ,xlim=c(135,155),ylim=c(-39,-23),cex.axis=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  points(lon[sigmaskW[,1]],lat[sigmaskW[,2]],col="black",pch=19,cex=0.15,xlim=c(135,155),ylim=c(-39,-23))
  dev.off()
}





for(n in 1:4)
{
  I=which(Cnames[,1]==name2[n])

cool2=apply(coolA[,,I],c(1,2),median)
warm2=apply(warmA[,,I],c(1,2),median)
fig2a=paste("Rain_corrUwind_sig_",name[n],"_WRF10_median_cool_postregrid.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_",name[n],"_WRF10_median_warm_postregrid.tiff",sep="")
fig2=paste("Rain_corrUwind_",name[n],"_WRF10_median_cool_postregrid.tiff",sep="")
fig3=paste("Rain_corrUwind_",name[n],"_WRF10_median_warm_postregrid.tiff",sep="")
fig4a=paste("Rain_corrUwind_sig_",name[n],"_WRF10_median_cool_postregrid_small.tiff",sep="")
fig5a=paste("Rain_corrUwind_sig_",name[n],"_WRF10_median_warm_postregrid_small.tiff",sep="")
fig4=paste("Rain_corrUwind_",name[n],"_WRF10_median_cool_postregrid_small.tiff",sep="")
fig5=paste("Rain_corrUwind_",name[n],"_WRF10_median_warm_postregrid_small.tiff",sep="")

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
tiff(file=fig4, height=500, width=600)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig5, height=500, width=600)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

cool2[abs(cool2)<=0.27]=NaN
warm2[abs(warm2)<=0.27]=NaN

tiff(file=fig2a, height=500, width=600)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3a, height=500, width=600)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig4a, height=500, width=600)
image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig5a, height=500, width=600)
image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
}

rlist=c("R1","R2","R3")
for(r in 1:3)
{
  I=which(Cnames[,2]==rlist[r])
  
  cool2=apply(coolA[,,I],c(1,2),median)
  warm2=apply(warmA[,,I],c(1,2),median)
  fig2a=paste("Rain_corrUwind_sig_CMIP3_",rlist[r],"_median_cool_postregrid.tiff",sep="")
  fig3a=paste("Rain_corrUwind_sig_CMIP3_",rlist[r],"_median_warm_postregrid.tiff",sep="")
  fig2=paste("Rain_corrUwind_CMIP3_",rlist[r],"_median_cool_postregrid.tiff",sep="")
  fig3=paste("Rain_corrUwind_CMIP3_",rlist[r],"_median_warm_postregrid.tiff",sep="")
  fig4a=paste("Rain_corrUwind_sig_CMIP3_",rlist[r],"_median_cool_postregrid_small.tiff",sep="")
  fig5a=paste("Rain_corrUwind_sig_CMIP3_",rlist[r],"_median_warm_postregrid_small.tiff",sep="")
  fig4=paste("Rain_corrUwind_CMIP3_",rlist[r],"_median_cool_postregrid_small.tiff",sep="")
  fig5=paste("Rain_corrUwind_CMIP3_",rlist[r],"_median_warm_postregrid_small.tiff",sep="")
  
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
  tiff(file=fig4, height=500, width=600)
  image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=fig5, height=500, width=600)
  image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  
  cool2[abs(cool2)<=0.27]=NaN
  warm2[abs(warm2)<=0.27]=NaN
  
  tiff(file=fig2a, height=500, width=800)
  image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=fig3a, height=500, width=800)
  image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=fig4a, height=500, width=600)
  image.plot(lon,lat,cool2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=fig5a, height=500, width=600)
  image.plot(lon,lat,warm2,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(135,155),ylim=c(-39,-23))
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

data1=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/20c3m/newcomp_wrf50_postregrid.csv")
data1=data1[-c(1,2,3),]
data2=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/20c3m/newcomp_wrf10_postregrid.csv")
data2=data2[-c(1,2,3),]
a=seq(5,13,2)
b=seq(4,12,2)
name=c("avecorr_reg","corr_regave","corr_point","avecorr_regORIG","corr_regaveORIG")

for(i in 1:5)
{
  fname=paste("~/Documents/Data/CMIP5/CMIP3/20c3m/WRF10/akima regrid/Corrs_scatter_",name[i],"_v5.tiff",sep="")
  tiff(file=fname,height=600,width=600,pointsize=20)
  plot(data2[2:5,a[i]],data2[2:5,b[i]],type="p",col="gray50",pch=19,xlab="Coastal correlation",ylab="Inland correlation",xlim=c(-0.7,0.2),ylim=c(-0.7,0.7))
  abline(h=0,col="grey")
  abline(v=0,col="grey")
  points(data2[6:17,a[i]],data2[6:17,b[i]],col="grey",pch=17)
  points(data2[1,a[i]],data2[1,b[i]],col="black",pch=4,cex=2,lwd=4)
  points(median(data2[2:5,a[i]]),median(data2[2:5,b[i]]),col="gray50",pch=1,cex=2,lwd=4)
  points(median(data2[6:17,a[i]]),median(data2[6:17,b[i]]),col="grey",pch=2,cex=2,lwd=4)
  points(cmip[2],cmip[1],col="black",pch=1,cex=2,lwd=4)
  legend("bottomright",c("AWAP","CMIP5","CMIP3","WRF10"),pch=c(4,1,1,2),col=c("black","black","gray50","grey"),pt.lwd=2)
  dev.off()
  
  fname=paste("~/Documents/Data/CMIP5/CMIP3/20c3m/WRF10/akima regrid/Corrs_scatter_",name[i],"_v5_both.tiff",sep="")
  tiff(file=fname,height=600,width=600,pointsize=20)
  plot(data2[2:5,a[i]],data2[2:5,b[i]],type="p",col="gray50",pch=19,xlab="East Coast",ylab="South Coast",xlim=c(-0.7,0.2),ylim=c(-0.7,0.7))
  abline(h=0,col="grey")
  abline(v=0,col="grey")
  
  points(data1[6:17,a[i]],data1[6:17,b[i]],col="gray65",pch=17)
  points(data2[6:17,a[i]],data2[6:17,b[i]],col="gray80",pch=15)
  points(data1[1,a[i]],data1[1,b[i]],col="black",pch=4,cex=2,lwd=4)
  points(median(data1[2:5,a[i]]),median(data1[2:5,b[i]]),col="gray50",pch=1,cex=2,lwd=4)
  points(median(data1[6:17,a[i]]),median(data1[6:17,b[i]]),col="gray65",pch=2,cex=2,lwd=4)
  points(median(data2[6:17,a[i]]),median(data2[6:17,b[i]]),col="gray80",pch=0,cex=2,lwd=4)
  points(cmip[2],cmip[1],col="black",pch=1,cex=2,lwd=4)
  legend("bottomright",c("AWAP","CMIP5","CMIP3","WRF50","WRF10"),
         pch=c(4,1,1,2,0),col=c("black","black","gray50","gray65","gray80"),pt.lwd=2)
  dev.off()
}