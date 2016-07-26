rm(list=ls())
setwd('~/Documents/ECLs')
library(RNetCDF)
#load("/home/nfs/z3478332/Documents/GDI/JE_WRF/Rain_daily.RData")

years2=seq(1970,2009)
rainC<-array(0,c(691,886,length(years2),4))
dimnames(rainC)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
eclrainC_mldb<-eclrainC_ncep<-eclrainC_erai<-rainC

rainW<-array(0,c(691,886,length(years2)-1,4))
dimnames(rainW)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
eclrainW_mldb<-eclrainW_ncep<-eclrainW_erai<-rainW
  
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
colnames(ECLlist)=c("NCEP","ERAI","MLDB")
for(i in 2:length(dates))
{
  I=which(NCEP$Date2<=dates[i] & NCEP$Date2>=dates[i]-18*60*60 & NCEP$Location==1)
  if(length(I)>0) ECLlist[i,1]=1
  I=which(ERAI$Date2<=dates[i] & ERAI$Date2>=dates[i]-18*60*60 & ERAI$Location==1)
  if(length(I)>0) ECLlist[i,2]=1
  I=which(MLDB$Date2<=dates[i] & MLDB$Date2>=dates[i]-18*60*60)
  if(length(I)>0) ECLlist[i,3]=1
}

for(i in 1:length(dates2))
{
  if(dates2[i]%%10000==101) print(dates2[i])
  fname<-paste('/media/Seagate Expansion Drive/Data/daily rainfall/rainfall_',decades[i],'0-',decades[i],'9/rainfall-',years[i],'/r',dates2[i],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  
  y=which(years2==years[i])
  if(months[i]<=4) y2=y-1 else y2=y
  
  if(months[i]>=5 & months[i]<=10)
  {
  rainC[,,y,1]=rainC[,,y,1]+data
  rainC[,,y,2]=rainC[,,y,2]+(data>=5)
  rainC[,,y,3]=rainC[,,y,3]+(data>=25)
  rainC[,,y,4]=rainC[,,y,4]+(data>=50)
  } else 
  {
  if(y2>0 & y2<length(years2))
  {
  rainW[,,y2,1]=rainW[,,y2,1]+data
  rainW[,,y2,2]=rainW[,,y2,2]+(data>=5)
  rainW[,,y2,3]=rainW[,,y2,3]+(data>=25)
  rainW[,,y2,4]=rainW[,,y2,4]+(data>=50)
  }
  }
  
  if(ECLlist[i,1]==1)
  {
    if(months[i]>=5 & months[i]<=10)
    {
    eclrainC_ncep[,,y,1]=eclrainC_ncep[,,y,1]+data
    eclrainC_ncep[,,y,2]=eclrainC_ncep[,,y,2]+(data>=5)
    eclrainC_ncep[,,y,3]=eclrainC_ncep[,,y,3]+(data>=25)
    eclrainC_ncep[,,y,4]=eclrainC_ncep[,,y,4]+(data>=50)
    } else 
    {
      if(y2>0 & y2<length(years2))
      {
        eclrainW_ncep[,,y2,1]=eclrainW_ncep[,,y2,1]+data
        eclrainW_ncep[,,y2,2]=eclrainW_ncep[,,y2,2]+(data>=5)
        eclrainW_ncep[,,y2,3]=eclrainW_ncep[,,y2,3]+(data>=25)
        eclrainW_ncep[,,y2,4]=eclrainW_ncep[,,y2,4]+(data>=50)
      }
      }
  }
  if(ECLlist[i,2]==1)
  {
    if(months[i]>=5 & months[i]<=10)
    {
      eclrainC_erai[,,y,1]=eclrainC_erai[,,y,1]+data
      eclrainC_erai[,,y,2]=eclrainC_erai[,,y,2]+(data>=5)
      eclrainC_erai[,,y,3]=eclrainC_erai[,,y,3]+(data>=25)
      eclrainC_erai[,,y,4]=eclrainC_erai[,,y,4]+(data>=50)
    } else 
    {
      if(y2>0 & y2<length(years2))
      {
        eclrainW_erai[,,y2,1]=eclrainW_erai[,,y2,1]+data
        eclrainW_erai[,,y2,2]=eclrainW_erai[,,y2,2]+(data>=5)
        eclrainW_erai[,,y2,3]=eclrainW_erai[,,y2,3]+(data>=25)
        eclrainW_erai[,,y2,4]=eclrainW_erai[,,y2,4]+(data>=50)
      }
    }
  }
  if(ECLlist[i,3]==1)
  {
    if(months[i]>=5 & months[i]<=10)
    {
      eclrainC_mldb[,,y,1]=eclrainC_mldb[,,y,1]+data
      eclrainC_mldb[,,y,2]=eclrainC_mldb[,,y,2]+(data>=5)
      eclrainC_mldb[,,y,3]=eclrainC_mldb[,,y,3]+(data>=25)
      eclrainC_mldb[,,y,4]=eclrainC_mldb[,,y,4]+(data>=50)
    } else 
    {
      if(y2>0 & y2<length(years2))
      {
        eclrainW_mldb[,,y2,1]=eclrainW_mldb[,,y2,1]+data
        eclrainW_mldb[,,y2,2]=eclrainW_mldb[,,y2,2]+(data>=5)
        eclrainW_mldb[,,y2,3]=eclrainW_mldb[,,y2,3]+(data>=25)
        eclrainW_mldb[,,y2,4]=eclrainW_mldb[,,y2,4]+(data>=50)
      }
    }
    
  }
  rm(data)
}

save(rainC,eclrainC_ncep,eclrainC_mldb,eclrainC_erai,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainC_AWAP_proj100_rad2cv1.RData")
save(rainW,eclrainW_ncep,eclrainW_mldb,eclrainW_erai,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainW_AWAP_proj100_rad2cv1.RData")

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

ratioC<-array(0,c(691,886,5,4))
dimnames(ratioC)[[3]]<-c("MLDB","NCEP","ERAI","NCEP-WRF","CMIP-WRF")
dimnames(ratioC)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ratioW<-ratioC

ratioC[,,1,]=apply(eclrainC_mldb[,,21:37,],c(1,2,4),sum)/apply(rainC[,,21:37,],c(1,2,4),sum)
ratioC[,,2,]=apply(eclrainC_ncep[,,21:40,],c(1,2,4),sum)/apply(rainC[,,21:40,],c(1,2,4),sum)
ratioC[,,3,]=apply(eclrainC_erai[,,21:40,],c(1,2,4),sum)/apply(rainC[,,21:40,],c(1,2,4),sum)
ratioW[,,1,]=apply(eclrainW_mldb[,,21:36,],c(1,2,4),sum)/apply(rainW[,,21:36,],c(1,2,4),sum)
ratioW[,,2,]=apply(eclrainW_ncep[,,21:39,],c(1,2,4),sum)/apply(rainW[,,21:39,],c(1,2,4),sum)
ratioW[,,3,]=apply(eclrainW_erai[,,21:39,],c(1,2,4),sum)/apply(rainW[,,21:39,],c(1,2,4),sum)

save(ratioW,ratioC,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainWC_ratio_AWAP_proj100_rad2cv1.RData")






