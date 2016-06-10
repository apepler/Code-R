library("RNetCDF")
setwd('~/Documents/Data/JE_WRF')
f1=open.nc('P.198411.day.nc')
Den=var.get.nc(f1,'P',c(19,21,2,1),c(1,1,1,61))/100
Gay=var.get.nc(f1,'P',c(32,42,2,1),c(1,1,1,61))/100
GDI=(Den-Gay)
##Hmmm, Deniliquin averages ~ 22 hPA higher. 

##Use Alejandro's output instead for SLP
f2=open.nc('/srv/ccrc/data23/z3444417/studies/Data/WRF50/2D/WRF50_Lon.nc')
lon=var.get.nc(f2,'Lon')
f2=open.nc('/srv/ccrc/data23/z3444417/studies/Data/WRF50/2D/WRF50_Lat.nc')
lat=var.get.nc(f2,'Lat')
f1=open.nc('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/WRF50_nnrp_R2_3h_psl_1979.nc')
Den=var.get.nc(f1,'SLP',c(103,47,1),c(1,1,2920))
Gay=var.get.nc(f1,'SLP',c(116,69,1),c(1,1,2920))
GDI=(Den-Gay)


##Run for all years/months 1979-2009
years=seq(1979,2009)
months=seq(1,12)
times=seq(0,21,3)
GDI<-matrix(0,length(years)*length(months),10)

n=1
for(i in 1:length(years))
  for(j in 1:length(months))
  {
    if(months[j]<10) fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-0',months[j],'.nc',sep="")
    else fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-',months[j],'.nc',sep="")
    f1=open.nc(fname)
    nt=dim.inq.nc(f1,'Time')
    nt<-nt$length
    Den=var.get.nc(f1,'SLP',c(103,47,1),c(1,1,nt))
    Gay=var.get.nc(f1,'SLP',c(116,69,1),c(1,1,nt))
    data1=(Den-Gay)
    time=var.get.nc(f1,'Times')
    time=as.integer(substr(time,12,13))
    close.nc(f1)
    GDI[n,1]=years[i]
    GDI[n,2]=months[j]
    for(m in 1:8)
    {
      I=which(time==times[m])
      GDI[n,2+m]=mean(data1[I])
    }
    n=n+1    
  }

GDI2=GDI #Anomalies
for(j in 1:length(months))
{
  I=which(GDI[,2]==months[j])
  m=colMeans(GDI[I,3:10])
  GDI2[I,3:10]=GDI[I,3:10]-m
}

write.csv(GDI,"GDI_WRF.csv")
write.csv(GDI2,"GDI_WRF_anom.csv")

##Compare to normal GDI
read.table('~/Documents/Timeseries/gdi.txt',header=T,sep="")->gdi
gdi=gdi[86:116,2:13]
comp=matrix(0,12,8)
for(j in 1:length(months))
{
  I=which(GDI[,2]==months[j])
  for(k in 1:8) comp[j,k]=cor(GDI[I,k+2],gdi[,j])
}
write.csv(comp,"GDI_WRF_cor.csv")

##WRF v2 - get ave slp for a 1 unit box around points

years=seq(1979,2009)
months=seq(1,12)
times=seq(0,21,3)
GDI<-matrix(0,length(years)*length(months),10)

n=1
for(i in 1:length(years))
  for(j in 1:length(months))
  {
    if(months[j]<10) fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-0',months[j],'.nc',sep="")
    else fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-',months[j],'.nc',sep="")
    f1=open.nc(fname)
    nt=dim.inq.nc(f1,'Time')
    nt<-nt$length
    Den1=var.get.nc(f1,'SLP',c(102,46,1),c(3,3,nt))
    Den=apply(Den1,3,mean)
    Gay1=var.get.nc(f1,'SLP',c(115,68,1),c(3,3,nt))
    Gay=apply(Gay1,3,mean)
    data1=(Den-Gay)
    time=var.get.nc(f1,'Times')
    time=as.integer(substr(time,12,13))
    close.nc(f1)
    GDI[n,1]=years[i]
    GDI[n,2]=months[j]
    for(m in 1:8)
    {
      I=which(time==times[m])
      GDI[n,2+m]=mean(data1[I])
    }
    n=n+1    
  }

GDI2=GDI #Anomalies
for(j in 1:length(months))
{
  I=which(GDI[,2]==months[j])
  m=colMeans(GDI[I,3:10])
  GDI2[I,3:10]=GDI[I,3:10]-m
}

write.csv(GDI,"GDI_WRF1.csv")
write.csv(GDI2,"GDI_WRF1_anom.csv")

##Compare to normal GDI
read.table('~/Documents/Timeseries/gdi.txt',header=T,sep="")->gdi
gdi=gdi[86:116,2:13]
comp=matrix(0,12,8)
for(j in 1:length(months))
{
  I=which(GDI[,2]==months[j])
  for(k in 1:8) comp[j,k]=cor(GDI[I,k+2],gdi[,j])
}
write.csv(comp,"GDI_WRF1_cor.csv")

##Comparison: What if I use NCEP reanalyses?
##Obviously this will be a bit less close as 2.5 degree resolution
library("RNetCDF")
setwd('~/Documents/Data/JE_WRF')
years=seq(1979,2009)
months=seq(1,12)
times=seq(0,18,6)
GDI<-matrix(0,length(years)*length(months),6)

n=1
for(i in 1:length(years))
{
  fname=paste('/srv/ccrc/data23/z3478332/NCEP/slp.',years[i],'.nc',sep="")
  f1=open.nc(fname)
  nt=dim.inq.nc(f1,'time')
  nt<-nt$length
  time=var.get.nc(f1,'time')
  time2=as.Date((time/24-2),origin="0001-01-01")
  date=matrix(0,length(time2),4)
  date[,1]=as.numeric(format(time2,"%Y"))
  date[,2]=as.numeric(format(time2,"%m"))
  date[,3]=as.numeric(format(time2,"%d"))
  date[,4]=(time %% 24)
  
  Den=var.get.nc(f1,'slp',c(59,51,1),c(1,1,nt),unpack=T)/100
  Gay=var.get.nc(f1,'slp',c(61,47,1),c(1,1,nt),unpack=T)/100
  data1=(Den-Gay)
  close.nc(f1)
  
  for(j in 1:length(months))
  {
    GDI[n,1]=years[i]
    GDI[n,2]=months[j]
    for(m in 1:4)
    {
      I=which(date[,2]==months[j] & date[,4]==times[m])
      GDI[n,2+m]=mean(data1[I])
    }
    n=n+1    
  }
}

GDI2=GDI #Anomalies
for(j in 1:length(months))
{
  I=which(GDI[,2]==months[j])
  m=colMeans(GDI[I,3:6])
  GDI2[I,3:6]=GDI[I,3:6]-m
}

write.csv(GDI,"GDI_NCEP.csv")
write.csv(GDI2,"GDI_NCEP_anom.csv")

##Compare to normal GDI
read.table('~/Documents/Timeseries/gdi.txt',header=T,sep="")->gdi
gdi=gdi[86:116,2:13]
comp=matrix(0,12,4)
for(j in 1:length(months))
{
  I=which(GDI[,2]==months[j])
  for(k in 1:4) comp[j,k]=cor(GDI[I,k+2],gdi[,j])
}
write.csv(comp,"GDI_NCEP_cor.csv")

##NCEP v2 - days +/- GDI
library("RNetCDF")
setwd('~/Documents/Data/JE_WRF')
years=seq(1979,2009)
months=seq(1,12)
times=seq(0,18,6)
GDI<-matrix(0,length(years),length(months))

n=1
for(i in 1:length(years))
{
  fname=paste('/srv/ccrc/data23/z3478332/NCEP/slp.',years[i],'.nc',sep="")
  f1=open.nc(fname)
  nt=dim.inq.nc(f1,'time')
  nt<-nt$length
  time=var.get.nc(f1,'time')
  time2=as.Date((time/24-2),origin="0001-01-01")
  date=matrix(0,length(time2),4)
  date[,1]=as.numeric(format(time2,"%Y"))
  date[,2]=as.numeric(format(time2,"%m"))
  date[,3]=as.numeric(format(time2,"%d"))
  date[,4]=(time %% 24)
  Den=var.get.nc(f1,'slp',c(59,51,1),c(1,1,nt),unpack=T)/100
  Gay=var.get.nc(f1,'slp',c(61,47,1),c(1,1,nt),unpack=T)/100
  data1=(Den-Gay)
  close.nc(f1)
  
  I=which(date[,4]==0)
  date=date[I,]
  data1=data1[I]
  
  for(j in 1:length(months))
  {
    I=which(date[,2]==months[j])
    J=which(data1[I]>0)
    GDI[i,j]=length(J)/length(I) 
  }
}
write.csv(GDI,"GDI_NCEP_prop.csv")

##All years of MSLP averaged across selected months
##Obviously doesn't do cross-year seasons. Yet

setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
library("abind")
library("akima")
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(12)

years=seq(1979,2009)
for(m in 1:12)
{
  NCEP<-array(0,dim=c(144,73,length(years)))
  for(i in 1:length(years))
  {
    fname=paste('/srv/ccrc/data23/z3478332/NCEP/slp.',years[i],'.nc',sep="")
    f1=open.nc(fname)
    nt=dim.inq.nc(f1,'time')
    nt<-nt$length
    time=var.get.nc(f1,'time')
    lat2=var.get.nc(f1,'lat')
    lon2=var.get.nc(f1,'lon')
    time2=as.Date((time/24-2),origin="0001-01-01")
    mm=as.numeric(format(time2,"%m"))
    hh=(time %% 24)
    slp=var.get.nc(f1,'slp',c(1,1,1),c(144,73,nt),unpack=T)/100
    NCEP[,,i]=apply(slp[,,(hh==0 & mm==m)],c(1,2),mean)
    close.nc(f1)
  }
  latn=lat2[lat2>(-60) & lat2<10]
  lonn=lon2[lon2>100 & lon2<180]
  NCEP2=NCEP[(lon2>100 & lon2<180),(lat2>(-60) & lat2<10),]
  
  f2=open.nc('/srv/ccrc/data23/z3444417/studies/Data/WRF50/2D/WRF50_Lon.nc')
  lon=var.get.nc(f2,'Lon')
  f2=open.nc('/srv/ccrc/data23/z3444417/studies/Data/WRF50/2D/WRF50_Lat.nc')
  lat=var.get.nc(f2,'Lat')
  WRF<-array(0,dim=c(215,144,length(years)))
  for(i in 1:length(years))
  {
    if(m<10) fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-0',m,'.nc',sep="")
    else fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-',m,'.nc',sep="")
    f1=open.nc(fname)
    nt=dim.inq.nc(f1,'Time')
    nt<-nt$length
    slp=var.get.nc(f1,'SLP',c(1,1,1),c(215,144,nt))
    time=var.get.nc(f1,'Times')
    time=as.integer(substr(time,12,13))
    close.nc(f1)
    WRF[,,i]=apply(slp[,,time==0],c(1,2),mean,na.rm=T)
  }
  
  latt=as.vector(lat)
  lont=as.vector(lon)
  lat3=seq(7.5,-57.5,-0.5)
  lon3=seq(102.5,177.5,0.5)
  WRF2=array(0,dim=c(length(lonn),length(latn),length(years)))
  WRF3=array(0,dim=c(length(lon3),length(lat3),length(years)))
  for(i in 1:length(years))
  {
    a=as.vector(WRF[,,i])
    b=interp(lont,latt,a,lonn,latn)
    WRF2[,,i]=b$z
    b=interp(lont,latt,a,lon3,lat3)
    WRF3[,,i]=b$z
  }
  
  fout1=paste('Plots/slp/MSLP_NCEP_',m,'.tiff',sep="")
  fout2=paste('Plots/slp/MSLP_WRF_',m,'.tiff',sep="")
  fout3=paste('Plots/slp/MSLP_diff_',m,'.tiff',sep="")
  
  tiff(file=fout1, height=600, width=600)
  contour(lonn,latn[27:1],apply(NCEP2[,27:1,],c(1,2),mean),levels=seq(970,1030,2))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=fout2, height=600, width=600)
  contour(lon3,lat3[131:1],apply(WRF3[,131:1,],c(1,2),mean),levels=seq(970,1030,2))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  ##And now, plot the difference!
  tiff(file=fout3, height=450, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(lonn,latn[27:1],apply(WRF2[,27:1,],c(1,2),mean)-apply(NCEP2[,27:1,],c(1,2),mean),levels=seq(-3,3,0.5),col=cm)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(lonn,latn[27:1],levels=seq(-3,3,0.5),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")
  dev.off()
}

##Correlations
corr<-matrix(0,31,27)
for(i in 1:31) for(j in 1:27) corr[i,j]<-cor(NCEP2[i,j,],WRF2[i,j,])
contour(lonn,latn[27:1],corr)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

plot.new()
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(lonn,latn[27:1],corr,levels=seq(0,1,0.1),col=cm)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(lonn,latn[27:1],levels=seq(0,1,0.1),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")


##Compare WRF and AWAP mean rainfall
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')##Plotty stuff   
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
library("akima")
library("abind")
library("RNetCDF")
##Run for all years/months 1979-2009
years=seq(1979,2009)
months=seq(1,12)
load("DailyGDI.RData")
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(20)
cm[9:12]="white"
source('~/Documents/R/corr.plot.R')
readMat('~/Documents/Data/Useful_ECL.mat')->Useful2
mask2<-t(Useful2$mask)
mask2[is.na(mask2)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

load("Rain_monthly.RData")
load("AWAP_monthly.RData")
latt=as.vector(lat)
lont=as.vector(lon)

WRFm<-AWAPm<-array(0,dim=c(886,691,(length(years)-1)))
  for(i in 1:(length(years)-1))
  {
    I=which(time2[,2]==years[i] & time2[,3]>4 & time2[,3]<=10)
    AWAPm[,,i]=apply(AWAP[,,I],c(1,2),sum)    
    a=as.vector(apply(monthR[,,I],c(1,2),sum))
    b=interp(lont,latt,a,Useful$x,Useful$y)
    WRFm[,,i]=b$z
  }
WRFg=apply(WRFm,c(1,2),mean)
AWAPg=apply(AWAPm,c(1,2),mean)

WRFm<-AWAPm<-array(0,dim=c(886,691,(length(years)-1)))
for(i in 1:(length(years)-1))
{
  I=which((time2[,2]==years[i] & time2[,3]>=11) | (time2[,2]==years[i+1] & time2[,3]<=3))
  AWAPm[,,i]=apply(AWAP[,,I],c(1,2),sum)    
  a=as.vector(apply(monthR[,,I],c(1,2),sum))
  b=interp(lont,latt,a,Useful$x,Useful$y)
  WRFm[,,i]=b$z
}
WRFn=apply(WRFm,c(1,2),mean)
AWAPn=apply(AWAPm,c(1,2),mean)
rm(AWAP,monthR,WRFm,AWAPm)

tiff(file="WRFrain_propAprOct_sea.tiff", height=300, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,(WRFg/(WRFn+WRFg))*t(Useful$mask),levels=seq(0,1,0.1),col=rainbow(10),xlim=c(130,160),ylim=c(-40,-25))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,WRFg,levels=seq(0,1,0.1),col=rainbow(20),xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
dev.off()
tiff(file="AWAPrain_propAprOct_sea.tiff", height=300, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,(AWAPg/(AWAPn+AWAPg))*t(Useful$mask),levels=seq(0,1,0.1),col=rainbow(10),xlim=c(130,160),ylim=c(-40,-25))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,WRFg,levels=seq(0,1,0.1),col=rainbow(10),xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
dev.off()
tiff(file="WRFrain_AprOct_propAWAP_sea.tiff", height=300, width=600)
cc=rainbow(20)
cc2=c(cc[1:5],cc[7:19])
d2=((WRFg-AWAPg)/AWAPg)
I=which(d2>1)
d2[I]=1.05
I=which(d2<(-0.45))
d2[I]=-0.45
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,d2*t(Useful$mask),levels=seq(-0.5,1.1,0.1),col=cc2,xlim=c(130,160),ylim=c(-40,-25))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,d2,levels=seq(-0.5,1.1,0.1),col=cc2,xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "",key.axes = axis(4,labels=c("-40%","-20%","0%","+20%","+40%","+60%","+80%","+100%"),at=seq(-0.4,1,0.2)))
dev.off()

