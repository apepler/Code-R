rm(list=ls())
setwd('~/Documents/ECLs')
ECLs<-read.csv('ECLlist.csv',header=T,sep=";")
years=seq(1980,2008)

ECL_all<-ECL_cool<-nECL_all<-nECL_cool<-array(0,dim=c(691,886,29))
ECL_warm<-nECL_warm<-array(0,dim=c(691,886,30)) ## Have to kill first & last
for(i in 1:length(ECLs[,1]))
{
  pos=ECLs[i,3]-1979
  fname<-paste('/media/Seagate Expansion Drive/daily rainfall/rainfall_',ECLs[i,2],'0-',ECLs[i,2],'9/rainfall-',ECLs[i,3],'/r',ECLs[i,1],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  
  if(ECLs[i,5]==1)
{
  ECL_all[,,pos]=ECL_all[,,pos]+data
  if(ECLs[i,4]>=5 & ECLs[i,4]<=10) ECL_cool[,,pos]=ECL_cool[,,pos]+data else 
    if(ECLs[i,4]>=11) ECL_warm[,,pos+1]=ECL_warm[,,pos+1]+data else ECL_warm[,,pos]=ECL_warm[,,pos]+data
} else {
  nECL_all[,,pos]=nECL_all[,,pos]+data
  if(ECLs[i,4]>=5 & ECLs[i,4]<=10) nECL_cool[,,pos]=nECL_cool[,,pos]+data else 
    if(ECLs[i,4]>=11) nECL_warm[,,pos+1]=nECL_warm[,,pos+1]+data else nECL_warm[,,pos]=nECL_warm[,,pos]+data
}
}
save(ECL_all,ECL_cool,ECL_warm,nECL_all,nECL_cool,nECL_warm,file="ECLrain_NCEP_8008_fix.RData")

####Now, WRF version
##Haha, have Rain_daily for 1979-2009

rm(list=ls())
setwd('~/Documents/ECLs')
ECLs<-read.csv('ECLlist.csv',header=T,sep=";")
load('/home/nfs/z3478332/Documents/Data/JE_WRF/Rain_daily.RData')
I=which(time[,1]>=1980 & time[,1]<=2008)
rm(rain)
rain2=rain2[,,I]
years=seq(1980,2008)
time=time[I,]

ECL_all<-nECL_all<-array(0,dim=c(215,144,29))
ECL_cool<-nECL_cool<-array(0,dim=c(215,144,29))
ECL_warm<-nECL_warm<-array(0,dim=c(215,144,28))
for(i in 1:29)
{
  I=which(ECLs[,6]==1 & ECLs[,3]==years[i])
  ECL_all[,,i]=apply(rain2[,,I],c(1,2),sum)
  I=which(ECLs[,6]==0 & ECLs[,3]==years[i])
  nECL_all[,,i]=apply(rain2[,,I],c(1,2),sum)
  I=which(ECLs[,6]==1 & ECLs[,3]==years[i] & ECLs[,4]>=5 & ECLs[,4]<=10)
  ECL_cool[,,i]=apply(rain2[,,I],c(1,2),sum)
  I=which(ECLs[,6]==0 & ECLs[,3]==years[i] & ECLs[,4]>=5 & ECLs[,4]<=10)
  nECL_cool[,,i]=apply(rain2[,,I],c(1,2),sum)  
}
for(i in 1:28)
{
  I=which(ECLs[,6]==1 & ((ECLs[,3]==years[i] & ECLs[,4]>=11) | (ECLs[,3]==years[i+1] & ECLs[,4]<=4)))
  ECL_warm[,,i]=apply(rain2[,,I],c(1,2),sum)
  I=which(ECLs[,6]==0 & ((ECLs[,3]==years[i] & ECLs[,4]>=11) | (ECLs[,3]==years[i+1] & ECLs[,4]<=4)))
  nECL_warm[,,i]=apply(rain2[,,I],c(1,2),sum)  
}

save(ECL_all,ECL_cool,ECL_warm,nECL_all,nECL_cool,nECL_warm,file="ECLrain_WRF_8008.RData")

##Finally, need to still have the data for the MLDB for comparison. Same period. & We were counting 00Z ECL with the rain to 9am that day. 

rm(list=ls())
setwd('~/Documents/ECLs')
ECLs<-read.csv('ECLlist.csv',header=T,sep=";")
years=seq(1980,2006)
ECLs=ECLs[ECLs[,3]<=2006,]

ECL_all<-ECL_cool<-nECL_all<-nECL_cool<-array(0,dim=c(691,886,27))
ECL_warm<-nECL_warm<-array(0,dim=c(691,886,28)) ## Have to kill first & last
for(i in 1:length(ECLs[,1]))
{
  pos=ECLs[i,3]-1979
  fname<-paste('/media/Seagate Expansion Drive/daily rainfall/rainfall_',ECLs[i,2],'0-',ECLs[i,2],'9/rainfall-',ECLs[i,3],'/r',ECLs[i,1],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  
  if(ECLs[i,7]==1)
  {
    ECL_all[,,pos]=ECL_all[,,pos]+data
    if(ECLs[i,4]>=5 & ECLs[i,4]<=10) ECL_cool[,,pos]=ECL_cool[,,pos]+data else 
      if(ECLs[i,4]>=11) ECL_warm[,,pos+1]=ECL_warm[,,pos+1]+data else ECL_warm[,,pos]=ECL_warm[,,pos]+data
  } else {
    nECL_all[,,pos]=nECL_all[,,pos]+data
    if(ECLs[i,4]>=5 & ECLs[i,4]<=10) nECL_cool[,,pos]=nECL_cool[,,pos]+data else 
      if(ECLs[i,4]>=11) nECL_warm[,,pos+1]=nECL_warm[,,pos+1]+data else nECL_warm[,,pos]=nECL_warm[,,pos]+data
  }
}
save(ECL_all,ECL_cool,ECL_warm,nECL_all,nECL_cool,nECL_warm,file="ECLrain_MLDB_8006_fix.RData")

##Now, days > 25mm
rm(list=ls())
setwd('~/Documents/ECLs')
ECLs<-read.csv('ECLlist.csv',header=T,sep=";")
years=seq(1980,2008)

ECL_all<-ECL_cool<-nECL_all<-nECL_cool<-array(0,dim=c(691,886,29))
ECL_warm<-nECL_warm<-array(0,dim=c(691,886,30)) ## Have to kill first & last
for(i in 1:length(ECLs[,1]))
{
  pos=ECLs[i,3]-1979
  fname<-paste('/media/Seagate Expansion Drive/daily rainfall/rainfall_',ECLs[i,2],'0-',ECLs[i,2],'9/rainfall-',ECLs[i,3],'/r',ECLs[i,1],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  I=which(data>=25)
  data2=matrix(0,691,886)
  data2[I]=1
  
  if(ECLs[i,5]==1)
  {
    ECL_all[,,pos]=ECL_all[,,pos]+data2
    if(ECLs[i,4]>=5 & ECLs[i,4]<=10) ECL_cool[,,pos]=ECL_cool[,,pos]+data2 else 
      if(ECLs[i,4]>=11) ECL_warm[,,pos+1]=ECL_warm[,,pos+1]+data2 else ECL_warm[,,pos]=ECL_warm[,,pos]+data2
  } else {
    nECL_all[,,pos]=nECL_all[,,pos]+data2
    if(ECLs[i,4]>=5 & ECLs[i,4]<=10) nECL_cool[,,pos]=nECL_cool[,,pos]+data2 else 
      if(ECLs[i,4]>=11) nECL_warm[,,pos+1]=nECL_warm[,,pos+1]+data2 else nECL_warm[,,pos]=nECL_warm[,,pos]+data2
  }
}
save(ECL_all,ECL_cool,ECL_warm,nECL_all,nECL_cool,nECL_warm,file="ECLrain_NCEP_8008_25mm.RData")

##Finally, need to still have the data for the MLDB for comparison. Same period. & We were counting 00Z ECL with the rain to 9am that day. 

rm(list=ls())
setwd('~/Documents/ECLs')
ECLs<-read.csv('ECLlist.csv',header=T,sep=";")
years=seq(1980,2006)
ECLs=ECLs[ECLs[,3]<=2006,]

ECL_all<-ECL_cool<-nECL_all<-nECL_cool<-array(0,dim=c(691,886,27))
ECL_warm<-nECL_warm<-array(0,dim=c(691,886,28)) ## Have to kill first & last
for(i in 1:length(ECLs[,1]))
{
  pos=ECLs[i,3]-1979
  fname<-paste('/media/Seagate Expansion Drive/daily rainfall/rainfall_',ECLs[i,2],'0-',ECLs[i,2],'9/rainfall-',ECLs[i,3],'/r',ECLs[i,1],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  I=which(data>=25)
  data2=matrix(0,691,886)
  data2[I]=1
  
  if(ECLs[i,7]==1)
  {
    ECL_all[,,pos]=ECL_all[,,pos]+data2
    if(ECLs[i,4]>=5 & ECLs[i,4]<=10) ECL_cool[,,pos]=ECL_cool[,,pos]+data2 else 
      if(ECLs[i,4]>=11) ECL_warm[,,pos+1]=ECL_warm[,,pos+1]+data2 else ECL_warm[,,pos]=ECL_warm[,,pos]+data2
  } else {
    nECL_all[,,pos]=nECL_all[,,pos]+data2
    if(ECLs[i,4]>=5 & ECLs[i,4]<=10) nECL_cool[,,pos]=nECL_cool[,,pos]+data2 else 
      if(ECLs[i,4]>=11) nECL_warm[,,pos+1]=nECL_warm[,,pos+1]+data2 else nECL_warm[,,pos]=nECL_warm[,,pos]+data2
  }
}
save(ECL_all,ECL_cool,ECL_warm,nECL_all,nECL_cool,nECL_warm,file="ECLrain_MLDB_8006_25mm.RData")