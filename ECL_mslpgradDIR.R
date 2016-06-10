rm(list=ls())
setwd('~/Documents/ECLs')
erai<-read.csv('Algorithm Comparison/UM_rad2_p100.csv')

##Testing - only using lows in 1981 & when in the location
erai=erai[erai$Date>=19810000 & erai$Date<=19820000 & erai$Location==1,]
date1=erai$Date+(as.numeric(erai$Time)-1)/4
erai$Ngrad<-erai$Sgrad<-erai$Egrad<-erai$Wgrad<-erai$GRAD<-0
erai2<-erai3<-erai

library(RNetCDF)
a<-open.nc('/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_1979-01_1983-12.nc')
lat<-var.get.nc(a,'latitude')
lon<-var.get.nc(a,'longitude')
time<-var.get.nc(a,'time')
time1=as.Date((time/24),origin="1900-01-01")
date2=as.double(format(time1,"%Y%m%d"))+(time %% 24)/24
I=which(date2>=19810000 & date2<=19820000)
slp<-var.get.nc(a,'msl',start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I)),unpack=T)/100
date2=date2[I]

library(sp)
for(k in 1:length(erai[,1]))
{
  dist=matrix(NaN,length(lon),length(lat))
  for(i in 1:length(lon))
    for(j in 1:length(lat))
      dist[i,j]=spDistsN1(cbind(lon[i],lat[j]),cbind(erai$Lon[k],erai$Lat[k]))
  
  loc=which(dist==min(dist),arr.ind=T) ##I.e. the lon/lat location
  I=which(date2==date1[k])             ##I.e. the time
  
  ##Gradient over 0-2.5 degrees
  erai$Ngrad[k]=slp[loc[1],loc[2],I]-slp[loc[1],loc[2]-1,I]
  erai$Sgrad[k]=slp[loc[1],loc[2],I]-slp[loc[1],loc[2]+1,I]
  erai$Egrad[k]=slp[loc[1],loc[2],I]-slp[loc[1]-1,loc[2],I]
  erai$Wgrad[k]=slp[loc[1],loc[2],I]-slp[loc[1]+1,loc[2],I]
  
  ##Gradient over 2.5 to 5 degrees
  erai2$Ngrad[k]=slp[loc[1],loc[2]-1,I]-slp[loc[1],loc[2]-2,I]
  erai2$Sgrad[k]=slp[loc[1],loc[2]+1,I]-slp[loc[1],loc[2]+2,I]
  erai2$Egrad[k]=slp[loc[1]-1,loc[2],I]-slp[loc[1]-2,loc[2],I]
  erai2$Wgrad[k]=slp[loc[1]+1,loc[2],I]-slp[loc[1]+2,loc[2],I]
  
  ##Big gradient - 10 degrees away (so prevailing winds)
  erai3$Ngrad[k]=slp[loc[1],loc[2]-1,I]-slp[loc[1],loc[2]-4,I]
  erai3$Sgrad[k]=slp[loc[1],loc[2]+1,I]-slp[loc[1],loc[2]+4,I]
  erai3$Egrad[k]=slp[loc[1],loc[2]-1,I]-slp[loc[1]-4,loc[2],I]
  erai3$Wgrad[k]=slp[loc[1],loc[2]+1,I]-slp[loc[1]+4,loc[2],I]
}

########### Of course, maybe it's better to actually do the mean slope rather than just the absolute difference. 

x<-rle(erai3$ID)
events<-data.frame(ID=x$values,Length=x$lengths,Wgrad=rep(0,length(x$values)),Egrad=rep(0,length(x$values)),
                   Sgrad=rep(0,length(x$values)),Ngrad=rep(0,length(x$values)),Match=rep("aaa",length(x$values)),stringsAsFactors=F)
for(i in 1:length(events3[,1])) 
{
  I<-which(erai3$ID==events3[i,1])
  events3[i,3:6]=apply(erai3[I,18:21],2,mean)
}

dirname=c("W","E","S","N")
for(i in 1:length(events3[,1]))
{
I=which(events3[i,5:6]==apply(events3[i,5:6],1,min))
events3$Match[i]=dirname[I+2]
}

test=data.frame(Low=events$Match,Near=events2$Match,Far=events3$Match,stringsAsFactors=F)

