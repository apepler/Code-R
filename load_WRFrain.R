rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
##Run for all years/months 1979-2009
years=seq(1979,2009)
months=seq(1,12)
times=seq(0,21,3)

##Daily GDI
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
    hh=as.integer(substr(time,12,13))
    close.nc(f1)
    I=which(hh==0)
    if(i==1 & j==1) GDI=data.frame(Year=as.integer(substr(time[I],1,4)),Month=as.integer(substr(time[I],6,7)),Day=as.integer(substr(time[I],9,10)),GDI=data1[I])
    else
    {
      data2=data.frame(Year=as.integer(substr(time[I],1,4)),Month=as.integer(substr(time[I],6,7)),Day=as.integer(substr(time[I],9,10)),GDI=data1[I])
      GDI=rbind(GDI,data2)
    }
  }

save(GDI,file="DailyGDI.RData")

##Now, getting the 00Z rainfall
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
years=seq(1981,2009)
months=seq(1,12)
library("abind")

for(i in 1:length(years))
{
  for(j in 1:length(months))
  {
    if(months[j]<10) fin=paste('/srv/ccrc/data18/z3393242/studies/NARCliM/nnrp/R2/out/wrfhrly_d01_',years[i],'-0',months[j],'-01_00:00:00',sep="")
    else fin=paste('/srv/ccrc/data18/z3393242/studies/NARCliM/nnrp/R2/out/wrfhrly_d01_',years[i],'-',months[j],'-01_00:00:00',sep="")
    f1=open.nc(fin)
    nt=dim.inq.nc(f1,'Time')
    nt<-nt$length
    time=var.get.nc(f1,'Times')
    time2=cbind(as.integer(substr(time,1,4)),as.integer(substr(time,6,7)),as.integer(substr(time,9,10)),as.integer(substr(time,12,13)))
    rainn=var.get.nc(f1,'RAINC',c(1,1,1),c(215,144,nt))
    rainc=var.get.nc(f1,'RAINNC',c(1,1,1),c(215,144,nt))
    I=which(time2[,4]==0)
    rain=rainn[,,I]+rainc[,,I]
    time2=time2[I,]
    
    if(j==1) 
    {
      Rain=rain
      Time=time2
      lat=var.get.nc(f1,'XLAT',c(1,1,1),c(215,144,1))
      lon=var.get.nc(f1,'XLONG',c(1,1,1),c(215,144,1))
    }
    else
    {
      Rain=abind(Rain,rain)
      Time=rbind(Time,time2)
    }
    close.nc(f1)
  }
  dd=dim(Rain)
  Rain2=array(0,dim=dd)
  for(k in 2:dd[3]) Rain2[,,k]=Rain[,,k]-Rain[,,k-1]
  fout=paste('Rain_',years[i],'.RData',sep="")
  save(Rain,Rain2,Time,lat,lon,file=fout)
  rm(Rain,Time)
}

rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
##Run for all years/months 1979-2009
years=seq(1979,2009)
library("abind")
for(i in 1:length(years))
{
  fin=paste('Rain_',years[i],'.RData',sep="")
  load(fin)
  rm(Rain2)
  if(i==1)
  {
    rain=Rain
    time=Time
    
  }  
  else 
  {
    rain=abind(rain,Rain)
    time=rbind(time,Time)
  }
}
dd=dim(rain)
rain2=array(0,dim=dd)
for(k in 2:dd[3]) rain2[,,k]=rain[,,k]-rain[,,k-1]
save(rain,rain2,time,lat,lon,file="Rain_daily.RData")

##Monthly rain total
ym=time[,1]*100+time[,2]
ym2=unique(ym)
ym2=cbind(ym2,floor(ym2/100),ym2 %% 100,rep(0,length(ym2)))
for(i in 1:length(ym2[,1]))
{
  I=which(ym==ym2[i,1])
  ym2[i,4]=length(I)
}

monthR=array(0,dim=c(215,144,length(ym2[,1])))
dd1=time[,1]*10000+time[,2]*100+time[,3]
dd2=ym2[,1]*100+ym2[,4]
I=which(dd1==dd2[1])
monthR[,,1]=rain[,,I]-rain[,,1]
for(i in 2:length(dd2))
{
  I1=which(dd1==dd2[i])
  I2=which(dd1==dd2[i-1])
  monthR[,,i]=rain[,,I1]-rain[,,I2]
}

monthM=array(0,dim=c(215,144,12))
for(i in 1:12)
{
  I=which(ym2[,3]==i)
  monthM[,,i]=apply(monthR[,,I],c(1,2),mean)
}
annR<-coolR<-array(0,dim=c(215,144,length(years)))
for(i in 1:length(years))
{
  I=which(ym2[,2]==years[i])
  annR[,,i]=apply(monthR[,,I],c(1,2),sum)
  I=which(ym2[,2]==years[i] & ym2[,3]>=6 & ym2[,3]<=10)
  coolR[,,i]=apply(monthR[,,I],c(1,2),sum)
}
warmR<-array(0,dim=c(215,144,(length(years)-1)))
for(i in 1:(length(years)-1))
{
  I=which((ym2[,2]==years[i] & ym2[,3]>=11) | (ym2[,2]==years[i+1] & ym2[,3]<=3))
  warmR[,,i]=apply(monthR[,,I],c(1,2),sum)
}
time2=ym2
save(monthR,annR,coolR,warmR,monthM,time2,lat,lon,file="Rain_monthly.RData")

##Still want to compare to the AWAP averages 1979-2009
for(i in 1:length(years))
  for(j in 1:12)
  {
    if(j<10) fin=paste('/media/Seagate Expansion Drive/monthly rainfall/',years[i],'0',j,'.grid',sep="")
    else fin=paste('/media/Seagate Expansion Drive/monthly rainfall/',years[i],j,'.grid',sep="")
    read.table(fin, sep="",skip=6,nrows=691)->data
    as.matrix(data)->data
    data[data<0]=NaN
    data<-t(data[nrow(data):1,])
    if(i==1 & j==1) AWAP=data else AWAP=abind(AWAP,data,along=3)
  }

annA<-coolA<-array(0,dim=c(886,691,length(years)))
for(i in 1:length(years))
{
  I=which(ym2[,2]==years[i])
  annA[,,i]=apply(AWAP[,,I],c(1,2),sum)
  I=which(ym2[,2]==years[i] & ym2[,3]>=6 & ym2[,3]<=10)
  coolA[,,i]=apply(AWAP[,,I],c(1,2),sum)
}
warmA<-array(0,dim=c(886,691,(length(years)-1)))
for(i in 1:(length(years)-1))
{
  I=which((ym2[,2]==years[i] & ym2[,3]>=11) | (ym2[,2]==years[i+1] & ym2[,3]<=3))
  warmA[,,i]=apply(AWAP[,,I],c(1,2),sum)
}
save(AWAP,annA,coolA,warmA,time2,Useful,file="AWAP_monthly.RData")


##Monthly MSLP
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
library("abind")
##Run for all years/months 1979-2009
years=seq(1979,2009)
months=seq(1,12)
for(i in 1:length(years))
  for(j in 1:length(months))
  {
    if(months[j]<10) fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-0',months[j],'.nc',sep="")
    else fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-',months[j],'.nc',sep="")
    f1=open.nc(fname)
    nt=dim.inq.nc(f1,'Time')
    nt<-nt$length
    slp=var.get.nc(f1,'SLP',c(1,1,1),c(215,144,nt))
    time=var.get.nc(f1,'Times')
    hh=as.integer(substr(time,12,13))
    close.nc(f1)
    I=which(hh==0)
    if(i==1 & j==1) MSLP=apply(slp[,,I],c(1,2),mean,na.rm=T)
    else MSLP=abind(MSLP,apply(slp[,,I],c(1,2),mean,na.rm=T),along=3)
  }

save(MSLP,file="WRF_MSLP.RData")
load("Rain_monthly.RData")
rm(monthM,monthR,warmR,coolR,annR)
latt=as.vector(lat)
lont=as.vector(lon)
library("akima")
lat3=seq(-57.5,7.5,0.5)
lon3=seq(102.5,177.5,0.5)

MSLP2=array(0,dim=c(151,131,372))
for(i in 1:372)
{
a=as.vector(MSLP[,,i])
a[which(is.na(a))]=0
b=interp(lont,latt,a,lon3,lat3)
MSLP2[,,i]=b$z
}
save(MSLP,MSLP2,file="WRF_MSLP.RData")

I=which(lon3>=145 & lon3<=150)
J=which(lat3>=(-45) & lat3<=(-10))
MSLP3=apply(MSLP2[I,J,],c(2,3),mean,na.rm=T)
lat4=lat3[J]
STRw=cbind(time2,matrix(0,372,2))
for(i in 1:372)
{
  I=which(MSLP3[,i]==max(MSLP3[,i]))
  STRw[i,5]=-lat4[I] ##So positive latitudes. COnsistency!
  STRw[i,6]=MSLP3[I,i]
}

##Now I need the NCEP version
f1=open.nc('ncep.mslp.nc')
nt=dim.inq.nc(f1,'time')
nt<-nt$length
time=var.get.nc(f1,'time')
lat2=var.get.nc(f1,'lat')
lon2=var.get.nc(f1,'lon')
time2=as.Date((time/24-2),origin="0001-01-01")
yy=as.numeric(format(time2,"%Y"))
slp=var.get.nc(f1,'slp',c(1,1,1),c(144,73,nt),unpack=T)
I=which(yy>=1979 & yy<=2009)
time2=time2[I]
slp=slp[,,I]

I=which(lon2>=145 & lon2<=150)
J=which(lat2>=(-45) & lat2<=(-10))
slp2=apply(slp[I,J,],c(2,3),mean,na.rm=T)
STRn=cbind(STRw[,1:4],matrix(0,372,2))
for(i in 1:372)
{
  y=aspline(lat2[J],slp2[,i],lat4)
  I=which(y$y==max(y$y))
  STRn[i,5]=-(y$x[I[1]]) ##So positive latitudes. Consistency! Where equal, choose lower lat
  STRn[i,6]=y$y[I[1]]
}

STRd=read.csv('../../Timeseries/str.csv')
STRd=cbind(STRd[,1:2],STRd[,4],STRd[,3])
STRn=cbind(STRn[,2:3],STRn[,5:6])
STRw=cbind(STRw[,2:3],STRw[,5:6])
save(STRd,STRn,STRw,file="STR.RData")

##For both, do a new GDI based on interpolated data (comparable)

GDIn=slp[61,47,]-slp[59,51,]
STRn=cbind(STRn,GDIn)
GDIw=MSLP2[99,65,]-MSLP2[86,45,]
STRw=cbind(STRw,GDIw)
save(STRd,STRn,STRw,file="STR.RData")

##
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
load("STR.RData")
aves=cbind(seq(1,12),matrix(0,12,6))
for(i in 1:12)
{
  I=which(STRd[,2]==i)
  aves[i,2:3]=colMeans(STRd[I,3:4])
  aves[i,4:5]=colMeans(STRn[I,3:4])
  aves[i,6:7]=colMeans(STRw[I,3:4])
}
aves=rbind(aves,aves[1,])
##28-41, 1012-1024
plot(aves[,3],-aves[,2],type='l',col='black',xlim=c(1012,1024),ylim=c(-41,-28),xlab="Intensity",ylab="Latitude")
lines(aves[,5],-aves[,4],type='l',col='blue')
lines(aves[,7],-aves[,6],type='l',col='red')
legend("topleft",legend=c("Drosd.","NCEP","WRF"),lty=c(1,1,1),col=c("black","blue","red"))