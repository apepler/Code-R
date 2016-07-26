rm(list=ls())
setwd('~/Documents/ECLs')
library(RNetCDF)
#load("/home/nfs/z3478332/Documents/GDI/JE_WRF/Rain_daily.RData")

years2=seq(1970,2009)
rainC<-array(0,c(886,691,length(years2),4))
dimnames(rainC)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ECLrainC_mldb<-ECLrainC_ncep<-ECLrainC_erai<-rainC

rainW<-array(0,c(886,691,length(years2)-1,4))
dimnames(rainW)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ECLrainW_mldb<-ECLrainW_ncep<-ECLrainW_erai<-rainW
  
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
  months2=months[I]
  
  M=which(months2>=5 & months2<=10)
  rainC[,,yy,1]=apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
  rainC[,,yy,2]=apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
  rainC[,,yy,3]=apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
  rainC[,,yy,4]=apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  
  M=which(months2<=4)
  if(yy>1) {
    rainW[,,yy-1,1]=rainW[,,yy-1,1]+apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    rainW[,,yy-1,2]=rainW[,,yy-1,2]+apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    rainW[,,yy-1,3]=rainW[,,yy-1,3]+apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    rainW[,,yy-1,4]=rainW[,,yy-1,4]+apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  M=which(months2>=11)
  if(yy<length(yearlist)) {
    rainW[,,yy,1]=rainW[,,yy,1]+apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    rainW[,,yy,2]=rainW[,,yy,2]+apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    rainW[,,yy,3]=rainW[,,yy,3]+apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    rainW[,,yy,4]=rainW[,,yy,4]+apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  
  M=which(months2>=5 & months2<=10 & ECLlist[I,1]==1)
  if(length(M)>0) {
    ECLrainC_mldb[,,yy,1]=apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainC_mldb[,,yy,2]=apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainC_mldb[,,yy,3]=apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainC_mldb[,,yy,4]=apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  
  M=which(months2<=4 & ECLlist[I,1]==1)
  if(yy>1 & length(M)>0) {
    ECLrainW_mldb[,,yy-1,1]=ECLrainW_mldb[,,yy-1,1]+apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainW_mldb[,,yy-1,2]=ECLrainW_mldb[,,yy-1,2]+apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainW_mldb[,,yy-1,3]=ECLrainW_mldb[,,yy-1,3]+apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainW_mldb[,,yy-1,4]=ECLrainW_mldb[,,yy-1,4]+apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  M=which(months2>=11 & ECLlist[I,1]==1)
  if(yy<length(yearlist) & length(M)>0) {
    ECLrainW_mldb[,,yy,1]=ECLrainW_mldb[,,yy,1]+apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainW_mldb[,,yy,2]=ECLrainW_mldb[,,yy,2]+apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainW_mldb[,,yy,3]=ECLrainW_mldb[,,yy,3]+apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainW_mldb[,,yy,4]=ECLrainW_mldb[,,yy,4]+apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  
  M=which(months2>=5 & months2<=10 & ECLlist[I,2]==1)
  if(length(M)>0) {
    ECLrainC_ncep[,,yy,1]=apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainC_ncep[,,yy,2]=apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainC_ncep[,,yy,3]=apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainC_ncep[,,yy,4]=apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  
  M=which(months2<=4 & ECLlist[I,2]==1)
  if(yy>1 & length(M)>0) {
    ECLrainW_ncep[,,yy-1,1]=ECLrainW_ncep[,,yy-1,1]+apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainW_ncep[,,yy-1,2]=ECLrainW_ncep[,,yy-1,2]+apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainW_ncep[,,yy-1,3]=ECLrainW_ncep[,,yy-1,3]+apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainW_ncep[,,yy-1,4]=ECLrainW_ncep[,,yy-1,4]+apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  M=which(months2>=11 & ECLlist[I,2]==1)
  if(yy<length(yearlist) & length(M)>0) {
    ECLrainW_ncep[,,yy,1]=ECLrainW_ncep[,,yy,1]+apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainW_ncep[,,yy,2]=ECLrainW_ncep[,,yy,2]+apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainW_ncep[,,yy,3]=ECLrainW_ncep[,,yy,3]+apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainW_ncep[,,yy,4]=ECLrainW_ncep[,,yy,4]+apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  
  M=which(months2>=5 & months2<=10 & ECLlist[I,3]==1)
  if(length(M)>0) {
    ECLrainC_erai[,,yy,1]=apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainC_erai[,,yy,2]=apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainC_erai[,,yy,3]=apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainC_erai[,,yy,4]=apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  
  M=which(months2<=4 & ECLlist[I,3]==1)
  if(yy>1 & length(M)>0) {
    ECLrainW_erai[,,yy-1,1]=ECLrainW_erai[,,yy-1,1]+apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainW_erai[,,yy-1,2]=ECLrainW_erai[,,yy-1,2]+apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainW_erai[,,yy-1,3]=ECLrainW_erai[,,yy-1,3]+apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainW_erai[,,yy-1,4]=ECLrainW_erai[,,yy-1,4]+apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }
  M=which(months2>=11 & ECLlist[I,3]==1)
  if(yy<length(yearlist) & length(M)>0) {
    ECLrainW_erai[,,yy,1]=ECLrainW_erai[,,yy,1]+apply(tmp[,,I[M]],c(1,2),sum,na.rm=T)
    ECLrainW_erai[,,yy,2]=ECLrainW_erai[,,yy,2]+apply(tmp[,,I[M]]>=5,c(1,2),sum,na.rm=T)
    ECLrainW_erai[,,yy,3]=ECLrainW_erai[,,yy,3]+apply(tmp[,,I[M]]>=25,c(1,2),sum,na.rm=T)
    ECLrainW_erai[,,yy,4]=ECLrainW_erai[,,yy,4]+apply(tmp[,,I[M]]>=50,c(1,2),sum,na.rm=T)
  }

}

save(rainC,ECLrainC_ncep,ECLrainC_mldb,ECLrainC_erai,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainC_AWAPcalib_proj100_rad2cv1.RData")
save(rainW,ECLrainW_ncep,ECLrainW_mldb,ECLrainW_erai,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainW_AWAPcalib_proj100_rad2cv1.RData")

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

ratioC<-array(0,c(886,691,5,4))
dimnames(ratioC)[[3]]<-c("MLDB","NCEP","ERAI","NCEP-WRF","CMIP-WRF")
dimnames(ratioC)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ratioW<-ratioC

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainC_AWAPcalib_proj100_rad2cv1.RData")
ratioC[,,1,]=apply(ECLrainC_mldb[,,21:37,],c(1,2,4),sum)/apply(rainC[,,21:37,],c(1,2,4),sum)
ratioC[,,2,]=apply(ECLrainC_ncep[,,21:40,],c(1,2,4),sum)/apply(rainC[,,21:40,],c(1,2,4),sum)
ratioC[,,3,]=apply(ECLrainC_erai[,,21:40,],c(1,2,4),sum)/apply(rainC[,,21:40,],c(1,2,4),sum)
rm(ECLrainC_ncep,ECLrainC_erai,ECLrainC_mldb,rainC)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainW_AWAPcalib_proj100_rad2cv1.RData")
ratioW[,,1,]=apply(ECLrainW_mldb[,,21:36,],c(1,2,4),sum)/apply(rainW[,,21:36,],c(1,2,4),sum)
ratioW[,,2,]=apply(ECLrainW_ncep[,,21:39,],c(1,2,4),sum)/apply(rainW[,,21:39,],c(1,2,4),sum)
ratioW[,,3,]=apply(ECLrainW_erai[,,21:39,],c(1,2,4),sum)/apply(rainW[,,21:39,],c(1,2,4),sum)
rm(ECLrainW_ncep,ECLrainW_erai,ECLrainW_mldb,rainW)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainC_NCEP-WRF_proj100_rad2cv1.RData")
ratio2=apply(ECLrain37C,c(1,2,4,5),sum)/apply(rainC,c(1,2,4,5),sum)
ratio3=apply(ratio2,c(1,2,4),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratioC[,,4,i]=t(b$z)
}
rm(ECLrainC,ECLrain37C,rainC)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainW_NCEP-WRF_proj100_rad2cv1.RData")
ratio2=apply(ECLrain37W,c(1,2,4,5),sum)/apply(rainW,c(1,2,4,5),sum)
ratio3=apply(ratio2,c(1,2,4),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratioW[,,4,i]=b$z
}
rm(ECLrainW,ECLrain37W,rainW)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainC_CMIP-WRF_proj100_rad2cv0.8.RData")
ratio2=apply(ECLrain37C,c(1,2,4,5,6),sum)/apply(rainC,c(1,2,4,5,6),sum)
ratio3=apply(ratio2,c(1,2,5),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratioC[,,5,i]=b$z
}
rm(ECLrainC,ECLrain37C,rainC)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainW_CMIP-WRF_proj100_rad2cv0.8.RData")
ratio2=apply(ECLrain37W,c(1,2,4,5,6),sum)/apply(rainW,c(1,2,4,5,6),sum)
ratio3=apply(ratio2,c(1,2,5),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratioW[,,5,i]=b$z
}
rm(ECLrainW,ECLrain37W,rainW)

save(ratioW,ratioC,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainWC_ratio_AWAPcalib_proj100_rad2cv1.RData")






