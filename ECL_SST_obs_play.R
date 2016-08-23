rm(list=ls())
library(RNetCDF)
setwd('/srv/ccrc/data34/z3478332/')

year=1982:2014

Tasman<-array(NaN,c(33,12,3))
names=c("HADISST","ERSST","OISST")

a2=open.nc("HadISST_sst.nc")
b=var.get.nc(a2,"time_bnds")
time2=as.Date(b[1,],origin="1870-1-1")
lon=var.get.nc(a2,"longitude")
lon[lon<0]=lon[lon<0]+360
lat=var.get.nc(a2,"latitude")
I=which(time2>=as.Date("19820101","%Y%m%d") & time2<as.Date("20160101","%Y%m%d") )
sst=var.get.nc(a2,"sst",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I)),unpack=T)
sstM=array(0,c(length(lon),length(lat),12))
for(i in 1:12) sstM[,,i]=apply(sst[,,which(as.integer(format(time2[I],"%m"))==i & as.integer(format(time2[I],"%Y"))%in%year)],c(1,2),mean,na.rm=T)
mask=!is.na(sst[,,1])
mask[mask==0]=NaN
mask2=mask[c(181:360,1:180),]

n=1
for(y in 1:33)
  for(m in 1:12)
  {
    I=which(lat>=-41 & lat<=-24)
    J=which(lon>=148 & lon<=161)
    Tasman[y,m,1]=max(sst[J,I,n]-sstM[J,I,m],na.rm=T)
    n=n+1
  }


a1=open.nc("noaa.oisst.mnmean.nc")
time1=as.Date(var.get.nc(a1,"time"),origin="1800-01-01")
lat1=var.get.nc(a1,"lat")
lon1=var.get.nc(a1,"lon")
I=which(time1>=as.Date("19820101","%Y%m%d") & time1<as.Date("20160101","%Y%m%d") )
sst1=var.get.nc(a1,"sst",start=c(1,1,I[1]),count=c(length(lon1),length(lat1),length(I)),unpack=T)
sstM=array(0,c(length(lon1),length(lat1),12))
for(i in 1:12) sstM[,,i]=apply(sst1[,,which(as.integer(format(time1[I],"%m"))==i & as.integer(format(time1[I],"%Y"))%in%year)],c(1,2),mean)

n=1
for(y in 1:33)
  for(m in 1:12)
  {
    I=which(lat1>=-41 & lat1<=-24)
    J=which(lon1>=148 & lon1<=161)
    Tasman[y,m,3]=max(sst1[J,I,n]*mask2[J,I]-sstM[J,I,m],na.rm=T)
    n=n+1
  }

a1=open.nc("noaa.ersst.mnmean.nc")
time1=as.Date(var.get.nc(a1,"time"),origin="1800-01-01")
lat1=var.get.nc(a1,"lat")
lon1=var.get.nc(a1,"lon")
I=which(time1>=as.Date("19820101","%Y%m%d") & time1<as.Date("20160101","%Y%m%d") )
sst1=var.get.nc(a1,"sst",start=c(1,1,I[1]),count=c(length(lon1),length(lat1),length(I)),unpack=T)
sstM=array(0,c(length(lon1),length(lat1),12))
for(i in 1:12) sstM[,,i]=apply(sst1[,,which(as.integer(format(time1[I],"%m"))==i & as.integer(format(time1[I],"%Y"))%in%year)],c(1,2),mean)

n=1
for(y in 1:33)
  for(m in 1:12)
  {
    I=which(lat1>=-41 & lat1<=-24)
    J=which(lon1>=148 & lon1<=161)
    Tasman[y,m,2]=max(sst1[J,I,n]-sstM[J,I,m],na.rm=T)
    n=n+1
  }

library(abind)
Tasman2=abind(apply(Tasman[,,],c(1,3),max),apply(Tasman[,3:5,],c(1,3),max),apply(Tasman[,6:8,],c(1,3),max),
              apply(Tasman[,9:11,],c(1,3),max),apply(Tasman[,c(12,1:2),],c(1,3),max),
              apply(Tasman[,5:10,],c(1,3),max),apply(Tasman[,5:10,],c(1,3),max),along=3)
dimnames(Tasman2)[[3]]=c("Annual","Autumn","Winter","Spring","Summer","Cool","Warm")

##Now, fix summer & warm
Tasman2[1:32,,5]=apply(abind(Tasman[1:32,12,],Tasman[2:33,1:2,],along=2),c(1,3),max)
Tasman2[1:32,,7]=apply(abind(Tasman[1:32,11:12,],Tasman[2:33,1:4,],along=2),c(1,3),max)
Tasman2[33,,c(5,7)]=NaN
        
### Okay, now I need all my different ECL databases
### Now, load in the mldb ECL database & collect April-August events for 1982-2006

ECLs<-array(NaN,c(33,7,2))

mldb=read.csv("~/Documents/ECLs/mldb_events.csv")
MLDB<-array(NaN,c(33,12,2))
for(i in 1:33)
  for(j in 1:12)
{
  I=which(mldb$Year==year[i] & mldb$Month==j)
  MLDB[i,j,1]=length(I)
}

erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100_update.csv")
erai$Year=floor(erai$Date1/10000)
erai$Month=floor(erai$Date1/100) %% 100
for(i in 1:33)
  for(j in 1:12)
{
  I=which(erai$Year==year[i] & erai$Month==j)
  MLDB[i,j,2]=length(I)
}

ECLs=abind(apply(MLDB[,,],c(1,3),sum),apply(MLDB[,3:5,],c(1,3),sum),apply(MLDB[,6:8,],c(1,3),sum),
              apply(MLDB[,9:11,],c(1,3),sum),apply(MLDB[,c(12,1:2),],c(1,3),sum),
              apply(MLDB[,5:10,],c(1,3),sum),apply(MLDB[,5:10,],c(1,3),sum),along=3)

##Now, fix summer & warm
ECLs[1:32,,5]=apply(abind(MLDB[1:32,12,],MLDB[2:33,1:2,],along=2),c(1,3),sum)
ECLs[1:32,,7]=apply(abind(MLDB[1:32,11:12,],MLDB[2:33,1:4,],along=2),c(1,3),sum)


corrM=array(0,c(3,7,2))
dimnames(corrM)[[1]]=names
dimnames(corrM)[[2]]=c("Annual","Autumn","Winter","Spring","Summer","Cool","Warm")
  for(j in 1:3)
    for(k in 1:7)
    {
    a=cor.test(ECLs[,2,k],Tasman2[,j,k],use="pairwise.complete.obs")
    corrM[j,k,1]=a$estimate
    corrM[j,k,2]=a$p.value
    }

autocorr<-matrix(0,4,7)

for(i in 1:7)
{
  a=acf(ECLs[,2,i])
  autocorr[1,i]=a$acf[2]
  for(j in 1:3)
  {
    a=acf(Tasman2[1:32,j,i])
    autocorr[j+1,i]=a$acf[2]
  }
}


