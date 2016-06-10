rm(list=ls())
setwd('~/Documents/ECLs')
library(RNetCDF)
#load("/home/nfs/z3478332/Documents/GDI/JE_WRF/Rain_daily.RData")

years2=seq(1970,2009)
rain<-array(0,c(691,886,length(years2),4))
dimnames(rain)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
eclrain_mldb<-eclrain_ncep<-eclrain_erai<-rain
  
dates=seq.POSIXt(ISOdatetime(1970,1,1,0,0,0, tz = "GMT"),ISOdatetime(2009,12,31,0,0,0, tz = "GMT"),by="days")
dates2=as.numeric(format.Date(dates,"%Y%m%d"))
years=floor(dates2/10000)
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
  rain[,,y,1]=rain[,,y,1]+data
  rain[,,y,2]=rain[,,y,2]+(data>=5)
  rain[,,y,3]=rain[,,y,3]+(data>=25)
  rain[,,y,4]=rain[,,y,4]+(data>=50)
  
  if(ECLlist[i,1]==1)
  {
    eclrain_ncep[,,y,1]=eclrain_ncep[,,y,1]+data
    eclrain_ncep[,,y,2]=eclrain_ncep[,,y,2]+(data>=5)
    eclrain_ncep[,,y,3]=eclrain_ncep[,,y,3]+(data>=25)
    eclrain_ncep[,,y,4]=eclrain_ncep[,,y,4]+(data>=50)
  }
  if(ECLlist[i,2]==1)
  {
    eclrain_erai[,,y,1]=eclrain_erai[,,y,1]+data
    eclrain_erai[,,y,2]=eclrain_erai[,,y,2]+(data>=5)
    eclrain_erai[,,y,3]=eclrain_erai[,,y,3]+(data>=25)
    eclrain_erai[,,y,4]=eclrain_erai[,,y,4]+(data>=50)
  }
  if(ECLlist[i,3]==1)
  {
    eclrain_mldb[,,y,1]=eclrain_mldb[,,y,1]+data
    eclrain_mldb[,,y,2]=eclrain_mldb[,,y,2]+(data>=5)
    eclrain_mldb[,,y,3]=eclrain_mldb[,,y,3]+(data>=25)
    eclrain_mldb[,,y,4]=eclrain_mldb[,,y,4]+(data>=50)
  }
  rm(data)
}

save(rain,eclrain_ncep,eclrain_mldb,eclrain_erai,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_AWAP_proj100_rad2cv1.RData")

##### Now, need daily ECL rain for NCEP-WRF, 1990-2009

rm(list=ls())
library(RNetCDF)
setwd('~/Documents/ECLs')
#load("/home/nfs/z3478332/Documents/GDI/JE_WRF/Rain_daily.RData")

yearlist=seq(1990,2009)
rain<-array(0,c(215,144,length(yearlist),3,4))
dimnames(rain)[[4]]=c("R1","R2","R3")
dimnames(rain)[[5]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ECLrain<-ECLrain37<-rain

dates=seq.POSIXt(ISOdatetime(1990,1,1,0,0,0, tz = "GMT"),ISOdatetime(2009,12,31,0,0,0, tz = "GMT"),by="days")
dates2=as.numeric(format.Date(dates,"%Y%m%d"))
years=floor(dates2/10000)
decades=floor(years/10)

proj="100"
rad="rad2cv1"
fixes<-events<-events2<-list()
thresh=c(1.32,1.07,1.39)

for(r in 1:3) {
  fixes[[r]]=read.csv(paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/proj",proj,"/outputUM_ncep_WRFR",r,"_50_",rad,"/ECLfixes_umelb_ncep_wrfR",r,"_proj",proj,"_",rad,"_9009.csv",sep=""))
  fixes[[r]]$Date2=as.POSIXct(paste(as.character(fixes[[r]]$Date),substr(fixes[[r]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
}


ECLlist<-matrix(0,length(dates),3)
colnames(ECLlist)=c("R1","R2","R3")
ECLlist2<-ECLlist
for(i in 2:length(dates))
  for(r in 1:3)
{

    ##For consistency - AWAP is 00Z day N-1 to 00Z day N, and we're using all ECLs from 0600 day N-1 to 00Z day N
    ## Then, since WRF totals are 00Z day N to 00Z day N+1, use all ECLs from 0600 day N to 00Z day N+1
    
  I=which(fixes[[r]]$Date2<=dates[i]+24*60*60 & fixes[[r]]$Date2>=dates[i]+6*60*60 & fixes[[r]]$Location==1)
  if(length(I)>0) ECLlist[i,r]=1
  
  I=which(fixes[[r]]$Date2<=dates[i]+24*60*60 & fixes[[r]]$Date2>=dates[i]+6*60*60 & fixes[[r]]$Location==1 & fixes[[r]]$CV>=thresh[r])
  if(length(I)>0) ECLlist2[i,r]=1
  
}
apply(ECLlist2,2,sum)/20

######### Now, loading the WRF daily rain
#### Can I even load a whole file? No, do individually
cmip="NNRP"
runyears="1950-2010"
ystart=c(1990,1995,2000,2005)

for(w in 1:3)
  for(s in 1:4)
  {
    dir=paste("/srv/ccrc/data30/z3393020/NARCliM/filtered/",cmip,"/R",w,"/",runyears,"/d01/",sep="")
    fname=paste(dir,"CCRC_NARCliM_DAY_",ystart[s],"-",ystart[s]+4,"_pracc_fl.nc",sep="")
    a=open.nc(fname)
    tmp=var.get.nc(a,"pracc_fl")
  
    ytmp=seq(ystart[s],ystart[s]+4)
    I=which(years %in% ytmp)
    years2=years[I]
    for(y in 1:5)
    {
      yy=which(yearlist==ytmp[y])
      I=which(years2==ytmp[y]) # Within the WRF data
      J=which(years==ytmp[y])  # Within the ECL list
      rain[,,yy,w,1]=apply(tmp[,,I],c(1,2),sum,na.rm=T)
      rain[,,yy,w,2]=apply(tmp[,,I]>=5,c(1,2),sum,na.rm=T)
      rain[,,yy,w,3]=apply(tmp[,,I]>=25,c(1,2),sum,na.rm=T)
      rain[,,yy,w,4]=apply(tmp[,,I]>=50,c(1,2),sum,na.rm=T)
      
      J=which(years==ytmp[y])  # Within the ECL list
      K=which(ECLlist[J,w]==1)
      ECLrain[,,yy,w,1]=apply(tmp[,,I[K]],c(1,2),sum,na.rm=T)
      ECLrain[,,yy,w,2]=apply(tmp[,,I[K]]>=5,c(1,2),sum,na.rm=T)
      ECLrain[,,yy,w,3]=apply(tmp[,,I[K]]>=25,c(1,2),sum,na.rm=T)
      ECLrain[,,yy,w,4]=apply(tmp[,,I[K]]>=50,c(1,2),sum,na.rm=T)
      K=which(ECLlist2[J,w]==1)
      ECLrain37[,,yy,w,1]=apply(tmp[,,I[K]],c(1,2),sum,na.rm=T)
      ECLrain37[,,yy,w,2]=apply(tmp[,,I[K]]>=5,c(1,2),sum,na.rm=T)
      ECLrain37[,,yy,w,3]=apply(tmp[,,I[K]]>=25,c(1,2),sum,na.rm=T)
      ECLrain37[,,yy,w,4]=apply(tmp[,,I[K]]>=50,c(1,2),sum,na.rm=T)
    }
  }

save(rain,ECLrain,ECLrain37,ECLlist,ECLlist2,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_NCEP-WRF_proj100_rad2cv1.RData")

####### And now, CMIP-WRF, 1990-2009

rm(list=ls())
library(RNetCDF)
setwd('~/Documents/ECLs')
#load("/home/nfs/z3478332/Documents/GDI/JE_WRF/Rain_daily.RData")
cmip=c("cccma","csiromk3","echam5","miroc")
cmip2=c("CCCMA3.1","CSIROMk3.0","ECHAM5","MIROC3.2")   
runyears="1990-2010"
ystart=c(1990,1995,2000,2005)

yearlist=seq(1990,2009)
rain<-array(0,c(215,144,length(yearlist),4,3,4))
dimnames(rain)[[4]]=cmip
dimnames(rain)[[5]]=c("R1","R2","R3")
dimnames(rain)[[6]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ECLrain<-ECLrain37<-rain

dates=seq.POSIXt(ISOdatetime(1990,1,1,0,0,0, tz = "GMT"),ISOdatetime(2009,12,31,0,0,0, tz = "GMT"),by="days")
dates2=as.numeric(format.Date(dates,"%Y%m%d"))
years=floor(dates2/10000)
decades=floor(years/10)

proj="100"
rad="rad2cv06"
rad2="rad2cv0.8"
thresh=rbind(c(1.01,0.92,1.08),c(1.07,0.86,1.07),c(1.08,0.94,1.12),c(1.05,0.91,0.94))
          
for(c in 1:4)
for(w in 1:3)
{
  fixes=read.csv(paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/proj",proj,"/outputUM_",cmip[c],"_WRFR",w,"_50_",rad,"/ECLfixes_umelb_",cmip[c],"_wrfR",w,"_proj",proj,"_",rad2,"_9009.csv",sep=""))
  fixes$Date2=as.POSIXct(paste(as.character(fixes$Date),substr(fixes$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  ECLlist<-matrix(0,length(dates),2)
  for(i in 2:length(dates))
    {
      ##For consistency - AWAP is 00Z day N-1 to 00Z day N, and we're using all ECLs from 0600 day N-1 to 00Z day N
      ## Then, since WRF totals are 00Z day N to 00Z day N+1, use all ECLs from 0600 day N to 00Z day N+1
      
      I=which(fixes$Date2<=dates[i]+24*60*60 & fixes$Date2>=dates[i]+6*60*60 & fixes$Location==1)
      if(length(I)>0) ECLlist[i,1]=1
      
      I=which(fixes$Date2<=dates[i]+24*60*60 & fixes$Date2>=dates[i]+6*60*60 & fixes$Location==1 & fixes$CV>=thresh[c,w])
      if(length(I)>0) ECLlist[i,2]=1
    }

  for(s in 1:4)
  {
    dir=paste("/srv/ccrc/data30/z3393020/NARCliM/filtered/",cmip2[c],"/R",w,"/",runyears,"/d01/",sep="")
    fname=paste(dir,"CCRC_NARCliM_DAY_",ystart[s],"-",ystart[s]+4,"_pracc_fl.nc",sep="")
    a=open.nc(fname)
    tmp=var.get.nc(a,"pracc_fl")
    
    ytmp=seq(ystart[s],ystart[s]+4)
    I=which(years %in% ytmp)
    years2=years[I]
    for(y in 1:5)
    {
      yy=which(yearlist==ytmp[y])
      I=which(years2==ytmp[y]) # Within the WRF data
      J=which(years==ytmp[y])  # Within the ECL list
      rain[,,yy,c,w,1]=apply(tmp[,,I],c(1,2),sum,na.rm=T)
      rain[,,yy,c,w,2]=apply(tmp[,,I]>=5,c(1,2),sum,na.rm=T)
      rain[,,yy,c,w,3]=apply(tmp[,,I]>=25,c(1,2),sum,na.rm=T)
      rain[,,yy,c,w,4]=apply(tmp[,,I]>=50,c(1,2),sum,na.rm=T)
      
      J=which(years==ytmp[y])  # Within the ECL list
      K=which(ECLlist[J,1]==1)
      ECLrain[,,yy,c,w,1]=apply(tmp[,,I[K]],c(1,2),sum,na.rm=T)
      ECLrain[,,yy,c,w,2]=apply(tmp[,,I[K]]>=5,c(1,2),sum,na.rm=T)
      ECLrain[,,yy,c,w,3]=apply(tmp[,,I[K]]>=25,c(1,2),sum,na.rm=T)
      ECLrain[,,yy,c,w,4]=apply(tmp[,,I[K]]>=50,c(1,2),sum,na.rm=T)
      K=which(ECLlist[J,2]==1)
      ECLrain37[,,yy,c,w,1]=apply(tmp[,,I[K]],c(1,2),sum,na.rm=T)
      ECLrain37[,,yy,c,w,2]=apply(tmp[,,I[K]]>=5,c(1,2),sum,na.rm=T)
      ECLrain37[,,yy,c,w,3]=apply(tmp[,,I[K]]>=25,c(1,2),sum,na.rm=T)
      ECLrain37[,,yy,c,w,4]=apply(tmp[,,I[K]]>=50,c(1,2),sum,na.rm=T)
    }
  }
}

save(rain,ECLrain,ECLrain37,ECLlist,dates,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_CMIP-WRF_proj100_rad2cv0.8.RData")



