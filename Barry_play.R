
##Version 2.0

rm(list=ls())
library(RNetCDF)
setwd("/srv/ccrc/data34/z3478332/NCEP/")
f1=open.nc('slp.1948.nc')
lon=var.get.nc(f1,"lon")
lat=var.get.nc(f1,"lat")
timeN=var.get.nc(f1,'time')
hour=timeN%%24
date=as.Date((timeN/24-2),origin="0001-01-01")
date2=(as.numeric(format(date,"%Y%m%d")))

I1=which(lon==150) # Moruya
J1=which(lat==-35)

I2=which(lon==152.5)
J2=which(lat==-30) # Coffs

data=data.frame(Date=date2,Time=hour)
##If using northern region
data$Moruya=var.get.nc(f1,'slp',c(I1,J1,1),c(1,1,length(date)),unpack=T)/100
data$Coffs=var.get.nc(f1,'slp',c(I2,J2,1),c(1,1,length(date)),unpack=T)/100

slp=data

for(year in 1949:2014)
{
  f1=open.nc(paste("slp.",year,".nc",sep=""))
  timeN=var.get.nc(f1,'time')
  hour=timeN%%24
  date=as.Date((timeN/24-2),origin="0001-01-01")
  date2=(as.numeric(format(date,"%Y%m%d")))
  
  data=data.frame(Date=date2,Time=hour)
  ##If using northern region
  data$Moruya=var.get.nc(f1,'slp',c(I1,J1,1),c(1,1,length(date)),unpack=T)/100
  data$Coffs=var.get.nc(f1,'slp',c(I2,J2,1),c(1,1,length(date)),unpack=T)/100
  
  slp=rbind(slp,data)
}

slp$Diff=slp$Coffs-slp$Moruya
write.csv(slp,file='NCEP_slpgrad_Barry.csv')
write.csv(slp[slp$Time==0,],file='NCEP_slpgrad_Barry_9am.csv')

### Now, ERAI


setwd("/srv/ccrc/data34/z3478332/ERAI/")
f1=open.nc('ERAI_mslp_1979-01_1983-12.nc')
lon=var.get.nc(f1,"longitude")
lat=var.get.nc(f1,"latitude")
timeN=var.get.nc(f1,'time')
hour=timeN%%24
date=as.Date((timeN/24),origin="1900-01-01")
date2=(as.numeric(format(date,"%Y%m%d")))

I1=which(lon==150) # Moruya
J1=which(lat==-36)

I2=which(lon==153)
J2=which(lat==-30) # Coffs

data=data.frame(Date=date2,Time=hour)
##If using northern region
data$Moruya=var.get.nc(f1,'msl',c(I1,J1,1),c(1,1,length(date)),unpack=T)/100
data$Coffs=var.get.nc(f1,'msl',c(I2,J2,1),c(1,1,length(date)),unpack=T)/100

slp=data
y1=c(1984,1989,1994,1999,2005,2011)
y2=c(1988,1993,1998,2004,2010,2014)

for(i in 1:6)
{
  f1=open.nc(paste('ERAI_mslp_',y1[i],'-01_',y2[i],'-12.nc',sep=""))
  timeN=var.get.nc(f1,'time')
  hour=timeN%%24
  date=as.Date((timeN/24),origin="1900-01-01")
  date2=(as.numeric(format(date,"%Y%m%d")))
  
  data=data.frame(Date=date2,Time=hour)
  ##If using northern region
  data$Moruya=var.get.nc(f1,'msl',c(I1,J1,1),c(1,1,length(date)),unpack=T)/100
  data$Coffs=var.get.nc(f1,'msl',c(I2,J2,1),c(1,1,length(date)),unpack=T)/100
  
  slp=rbind(slp,data)
}

slp$Diff=slp$Coffs-slp$Moruya
write.csv(slp,file='ERAI_slpgrad_Barry.csv')
write.csv(slp[slp$Time==0,],file='ERAI_slpgrad_Barry_9am.csv')
