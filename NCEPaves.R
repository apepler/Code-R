rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
years=seq(1950,2012)
months=seq(1,12)
times=seq(0,18,6)

NCEP<-array(0,dim=c(144,73,length(months),length(times),length(years)))
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
  for(j in 1:length(months))
    for(k in 1:length(times))
  {
    I=which(mm==months[j] & hh==times[k])
    NCEP[,,j,k,i]=apply(slp[,,I],c(1,2),mean,na.rm=T)
  }
  close.nc(f1)
}
NCEP2=apply(NCEP,c(1,2,3,4),mean,na.rm=T)
latn=lat2[lat2>(-43) & lat2<(-22)]
lonn=lon2[lon2>147 & lon2<163]
NCEP2=NCEP2[(lon2>147 & lon2<163),(lat2>(-43) & lat2<(-22)),,]
NCEP3<-matrix(0,length(months)*length(times)*length(latn)*length(lonn),5)
n=1
for(i in 1:length(months))
  for(j in 1:length(times))end process
    for(k in 1:length(latn))
      for(l in 1:length(lonn))
      {
        NCEP3[n,1]=months[i]
        NCEP3[n,2]=times[j]
        NCEP3[n,3]=latn[k]
        NCEP3[n,4]=lonn[l]
        NCEP3[n,5]=NCEP2[l,k,i,j]
        n=n+1    
      }

write.csv(NCEP3,"aveMSLP_19502012.csv")

data=read.csv("Data.csv",header=T)
for(i in 1:length(data[,1]))
{
  I=which(NCEP3[,1]==data[i,4] & NCEP3[,2]==data[i,5] & NCEP3[,3]==data[i,12] & NCEP3[,4]==data[i,11])
  data[i,13]=NCEP3[I,5]
  data[i,14]=data[i,9]-data[i,13]
}
write.csv(data,"Data.csv")

##V2 - Using 1990-2009 period
##Remember, latitudinal average over whole ESB domain, 147.5 - 162.5E
##For lats -22.5:-42.5

rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
years=seq(1990,2009)
months=seq(1,12)
times=seq(0,18,6)

NCEP<-array(0,dim=c(7,9,length(months),length(times),length(years)))
for(i in 1:length(years))
{
  fname=paste('/srv/ccrc/data23/z3478332/NCEP/slp.',years[i],'.nc',sep="")
  f1=open.nc(fname)
  nt=dim.inq.nc(f1,'time')
  nt<-nt$length
  time=var.get.nc(f1,'time')
  time2=as.Date((time/24-2),origin="0001-01-01")
  mm=as.numeric(format(time2,"%m"))
  hh=(time %% 24)
  slp=var.get.nc(f1,'slp',c(60,46,1),c(7,9,nt),unpack=T)/100
  for(j in 1:length(months))
    for(k in 1:length(times))
    {
      I=which(mm==months[j] & hh==times[k])
      NCEP[,,j,k,i]=apply(slp[,,I],c(1,2),mean,na.rm=T)
    }
  close.nc(f1)
}
NCEP2=apply(NCEP,c(2,3,4),mean,na.rm=T)
latn=seq(-22.5,-42.5,-2.5)
aveMSLP<-matrix(0,length(months)*length(times)*length(latn),5)
n=1
for(i in 1:length(months))
  for(j in 1:length(times))
for(k in 1:length(latn))
  {
    aveMSLP[n,1]=months[i]
    aveMSLP[n,2]=times[j]
    aveMSLP[n,3]=latn[k]
    aveMSLP[n,4]=NCEP2[k,i,j]
    n=n+1    
  }

##Now, add in NCEP2!
NCEP<-array(0,dim=c(7,9,length(months),length(times),length(years)))
for(i in 1:length(years))
{
  fname=paste('/srv/ccrc/data23/z3444417/studies/Data/nnrp2/global/psl/raw/mslp.',years[i],'.nc',sep="")
  f1=open.nc(fname)
  nt=dim.inq.nc(f1,'time')
  nt<-nt$length
  time=var.get.nc(f1,'time')
  time2=as.Date((time/24-2),origin="1800-01-03")
  mm=as.numeric(format(time2,"%m"))
  hh=(time %% 24)
  slp=var.get.nc(f1,'mslp',c(60,46,1),c(7,9,nt),unpack=T)/100
  for(j in 1:length(months))
    for(k in 1:length(times))
    {
      I=which(mm==months[j] & hh==times[k])
      NCEP[,,j,k,i]=apply(slp[,,I],c(1,2),mean,na.rm=T)
    }
  close.nc(f1)
}
NCEP2=apply(NCEP,c(2,3,4),mean,na.rm=T)
n=1
for(i in 1:length(months))
  for(j in 1:length(times))
    for(k in 1:length(latn))
    {
      aveMSLP[n,5]=NCEP2[k,i,j]
      n=n+1    
    }
write.csv(aveMSLP,"aveMSLP_19902009.csv")

##Finally, ERAI
##Slightly more of a pain due to the weird way of doing files.
rm(list=ls())
years=seq(1990,2009)
months=seq(1,12)
times=seq(0,18,6)
f1=open.nc('/srv/ccrc/data23/z3444417/studies/Data/ERAI/global/psl/raw/ERAI_mslp_1989-01_1999-12.nc')
f2=open.nc('/srv/ccrc/data23/z3444417/studies/Data/ERAI/global/psl/raw/ERAI_mslp_2000-01_2010-12.nc')
lat=var.get.nc(f1,'latitude')
lon=var.get.nc(f1,'longitude')
time1=var.get.nc(f1,'time')
time1a=as.Date((time1/24-2),origin="1900-01-03")
time2=var.get.nc(f2,'time')
time2a=as.Date((time2/24-2),origin="1900-01-03")

yy1=as.numeric(format(time1a,"%Y"))
mm1=as.numeric(format(time1a,"%m"))
hh1=(time1 %% 24)
yy2=as.numeric(format(time2a,"%Y"))
mm2=as.numeric(format(time2a,"%m"))
hh2=(time2 %% 24)

ERAI<-array(0,dim=c(10,14,length(months),length(times),length(years)))
for(i in 1:length(years))
  for(j in 1:length(months))
  {
    if(years[i]<2000)
    {
      I=which(yy1==years[i] & mm1==months[j])
      slp=var.get.nc(f1,'msl',c(100,76,I[1]),c(10,14,length(I)),unpack=T)/100
      h=hh1[I]
    }
    else
    {
      I=which(yy2==years[i] & mm2==months[j])
      slp=var.get.nc(f2,'msl',c(100,76,I[1]),c(10,14,length(I)),unpack=T)/100
      h=hh2[I]
    }
    for(k in 1:length(times))
    {
      I=which(h==times[k])
      ERAI[,,j,k,i]=apply(slp[,,I],c(1,2),mean,na.rm=T)
    }
  }
latn=seq(-22.5,-42,-1.5)
ERAI2=apply(ERAI,c(2,3,4),mean,na.rm=T)
aveMSLP2<-matrix(0,length(months)*length(times)*length(latn),4)
n=1
for(i in 1:length(months))
  for(j in 1:length(times))
    for(k in 1:length(latn))
    {
      aveMSLP2[n,1]=months[i]
      aveMSLP2[n,2]=times[j]
      aveMSLP2[n,3]=latn[k]
      aveMSLP2[n,4]=ERAI2[k,i,j]
      n=n+1    
    }
write.csv(aveMSLP2,"aveMSLP_19902009_erai_fix.csv")

##ERAI Half
##Slightly more of a pain due to the weird way of doing files.
rm(list=ls())
years=seq(1990,2009)
months=seq(1,12)
times=seq(0,18,6)
f1=open.nc('/srv/ccrc/data23/z3478332/ERAI/ERAI_mslp_1989-01_1999-12_half.nc')
f2=open.nc('/srv/ccrc/data23/z3478332/ERAI/ERAI_mslp_2000-01_2010-12_half.nc')
lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
time1=var.get.nc(f1,'time')
time1a=as.Date((time1/24-2),origin="1900-01-03")
time2=var.get.nc(f2,'time')
time2a=as.Date((time2/24-2),origin="1900-01-03")

yy1=as.numeric(format(time1a,"%Y"))
mm1=as.numeric(format(time1a,"%m"))
hh1=(time1 %% 24)
yy2=as.numeric(format(time2a,"%Y"))
mm2=as.numeric(format(time2a,"%m"))
hh2=(time2 %% 24)

ERAI<-array(0,dim=c(5,7,length(months),length(times),length(years)))
for(i in 1:length(years))
  for(j in 1:length(months))
  {
    if(years[i]<2000)
    {
      I=which(yy1==years[i] & mm1==months[j])
      slp=var.get.nc(f1,'msl',c(51,39,I[1]),c(5,7,length(I)),unpack=T)/100
      h=hh1[I]
    }
    else
    {
      I=which(yy2==years[i] & mm2==months[j])
      slp=var.get.nc(f2,'msl',c(51,39,I[1]),c(5,7,length(I)),unpack=T)/100
      h=hh2[I]
    }
    for(k in 1:length(times))
    {
      I=which(h==times[k])
      ERAI[,,j,k,i]=apply(slp[,,I],c(1,2),mean,na.rm=T)
    }
  }
latn=seq(-24,-42,-3)
ERAI2=apply(ERAI,c(2,3,4),mean,na.rm=T)
aveMSLP2<-matrix(0,length(months)*length(times)*length(latn),4)
n=1
for(i in 1:length(months))
  for(j in 1:length(times))
    for(k in 1:length(latn))
    {
      aveMSLP2[n,1]=months[i]
      aveMSLP2[n,2]=times[j]
      aveMSLP2[n,3]=latn[k]
      aveMSLP2[n,4]=ERAI2[k,i,j]
      n=n+1    
    }
write.csv(aveMSLP2,"aveMSLP_19902009_erai_half.csv")

#### JRA25

rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
years=seq(1990,2009)
months=seq(1,12)
times=seq(0,18,6)

JRA<-array(0,dim=c(7,9,length(months),length(times),length(years)))
for(i in 1:length(years))
  for(j in 1:length(months))
{
  fname=paste('/srv/ccrc/data23/z3478332/JRA25/prmslmsl.anlp6h',years[i],sprintf("%02d",months[j]),'.nc',sep="")
  f1=open.nc(fname)
  time=var.get.nc(f1,'initial_time0_encoded')
  nt=length(time)
  hh=(time %% 100)
  slp=var.get.nc(f1,'PRMSL_GDS0_MSL',c(60,46,1),c(7,9,nt))/100
    for(k in 1:length(times))
    {
      I=which(hh==times[k])
      JRA[,,j,k,i]=apply(slp[,,I],c(1,2),mean,na.rm=T)
    }
  close.nc(f1)
}
JRA2=apply(JRA,c(2,3,4),mean,na.rm=T)
latn=seq(-22.5,-42.5,-2.5)
aveMSLP<-matrix(0,length(months)*length(times)*length(latn),5)
n=1
for(i in 1:length(months))
  for(j in 1:length(times))
    for(k in 1:length(latn))
    {
      aveMSLP[n,1]=months[i]
      aveMSLP[n,2]=times[j]
      aveMSLP[n,3]=latn[k]
      aveMSLP[n,4]=JRA2[k,i,j]
      n=n+1    
    }
write.csv(aveMSLP,"aveMSLP_19902009_jra25.csv")
