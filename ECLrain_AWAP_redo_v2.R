rm(list=ls())
setwd('~/Documents/ECLs')
library(RNetCDF)
#load("/home/nfs/z3478332/Documents/GDI/JE_WRF/Rain_daily.RData")

years2=seq(1970,2009)
rain<-array(0,c(886,691,length(years2),4))
dimnames(rain)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ECLrain_mldb<-ECLrain_ncep<-ECLrain_erai<-rain

dates=seq.POSIXt(ISOdatetime(1970,1,1,0,0,0, tz = "GMT"),ISOdatetime(2009,12,31,0,0,0, tz = "GMT"),by="days")
dates2=as.numeric(format.Date(dates,"%Y%m%d"))
years=floor(dates2/10000)
months=floor(dates2/100)%%100
decades=floor(years/10)

### Now I need to get the ECL days in a clear way, & consider a day's rain as an ECL if any ECL between 00Z Day & 00Z day-1 (since AWAP is tp 9am))

NCEP=rbind(read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ncep1_topo_rad2_proj100/ECLfixes_umelb_ncep1_topo_rad2_proj100_4970.csv"),
           read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ncep1_topo_rad2_proj100/ECLfixes_umelb_ncep1_topo_rad2_proj100_7190.csv"),
           read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ncep1_topo_rad2_proj100/ECLfixes_umelb_ncep1_topo_rad2_proj100_9114.csv"))

NCEP=NCEP[NCEP$Date>=19700000 & NCEP$Date<=20100000,]
NCEP$Date2=as.POSIXct(paste(as.character(NCEP$Date),substr(NCEP$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

ERAI=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")
ERAI=ERAI[ERAI$Date>=19700000 & ERAI$Date<=20100000,]
ERAI$Date2=as.POSIXct(paste(as.character(ERAI$Date),substr(ERAI$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

MLDB=read.csv("~/Documents/ECLs/Speer_ECL_fixes.csv")
MLDB=MLDB[MLDB$date>=19700000 & MLDB$date<=20060000,]
MLDB$Date2=as.POSIXct(paste(as.character(MLDB$date),as.character(MLDB$time),sep=""),format="%Y%m%d%H",tz="GMT")

ECLlist<-matrix(0,length(dates),3)
colnames(ECLlist)=c("MLDB","NCEP","ERAI")
for(i in 2:length(dates))
{
  I=which(MLDB$Date2<=dates[i] & MLDB$Date2>=dates[i]-18*60*60)
  if(length(I)>0) ECLlist[i,1]=1
  I=which(NCEP$Date2<=dates[i] & NCEP$Date2>=dates[i]-18*60*60 & NCEP$Location==1)
  if(length(I)>0) ECLlist[i,2]=1
  I=which(ERAI$Date2<=dates[i] & ERAI$Date2>=dates[i]-18*60*60 & ERAI$Location==1)
  if(length(I)>0) ECLlist[i,3]=1

}

for(y in 1:length(years2))
{
  print(years2[y])
  fname=paste("/srv/ccrc/data02/z3236814/data/AWAP/DAILY/netcdf/rainfall_calib/pre.",years2[y],".nc",sep="")
  a=open.nc(fname)
  tmp=var.get.nc(a,"pre")

  I=which(years==years2[y])  # Within the ECL list
  
  rain[,,yy,1]=apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
  rain[,,yy,2]=apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
  rain[,,yy,3]=apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
  rain[,,yy,4]=apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  
  M=which(ECLlist[I,1]==1)
  if(length(M)>0) {
    ECLrain_mldb[,,yy,1]=apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrain_mldb[,,yy,2]=apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrain_mldb[,,yy,3]=apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrain_mldb[,,yy,4]=apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  
  M=which(ECLlist[I,2]==1)
  if(length(M)>0) {
    ECLrain_ncep[,,yy,1]=apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrain_ncep[,,yy,2]=apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrain_ncep[,,yy,3]=apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrain_ncep[,,yy,4]=apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  
  M=which(ECLlist[I,3]==1)
  if(length(M)>0) {
    ECLrain_erai[,,yy,1]=apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrain_erai[,,yy,2]=apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrain_erai[,,yy,3]=apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrain_erai[,,yy,4]=apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  

}

save(rain,ECLrain_ncep,ECLrain_mldb,ECLrain_erai,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_AWAPcalib_proj100_rad2cv1.RData")

rm(list=ls())
setwd('~/Documents/ECLs')
library(RNetCDF)
#load("/home/nfs/z3478332/Documents/GDI/JE_WRF/Rain_daily.RData")
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->UsefulE
maskE<-t(UsefulE$mask)
maskE[is.na(maskE)]=0
readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
library(RNetCDF)
library(akima)
a=open.nc("/srv/ccrc/data45/z3478332/WRF/output/extracted_data/WRF_d01_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
latt=as.vector(lat)
lont=as.vector(lon)

ratio<-array(0,c(886,691,5,4))
dimnames(ratio)[[3]]<-c("MLDB","NCEP","ERAI","NCEP-WRF","CMIP-WRF")
dimnames(ratio)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_AWAPcalib_proj100_rad2cv1.RData")
ratio[,,1,]=apply(ECLrain_mldb[,,21:37,],c(1,2,4),sum)/apply(rain[,,21:37,],c(1,2,4),sum)
ratio[,,2,]=apply(ECLrain_ncep[,,21:40,],c(1,2,4),sum)/apply(rain[,,21:40,],c(1,2,4),sum)
ratio[,,3,]=apply(ECLrain_erai[,,21:40,],c(1,2,4),sum)/apply(rain[,,21:40,],c(1,2,4),sum)
rm(ECLrain_ncep,ECLrain_erai,ECLrain_mldb,rain)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_NCEP-WRF_proj100_rad2cv1.RData")
ratio2=apply(ECLrain37,c(1,2,4,5),sum)/apply(rain,c(1,2,4,5),sum)
ratio3=apply(ratio2,c(1,2,4),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratio[,,4,i]=b$z
}
rm(ECLrain,ECLrain37,rain)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_CMIP-WRF_proj100_rad2cv0.8.RData")
ratio2=apply(ECLrain37,c(1,2,4,5,6),sum)/apply(rain,c(1,2,4,5,6),sum)
ratio3=apply(ratio2,c(1,2,5),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratio[,,5,i]=b$z
}
rm(ECLrain,ECLrain37,rain)

save(ratio,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_ratio_AWAPcalib_proj100_rad2cv1.RData")






