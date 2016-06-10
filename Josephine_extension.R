###########
## Step 1 - COnvert Power databse of events into a dataset of ECL/event by day
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/20CR/ECLs/")

data=read.csv("Power_ECLs.csv")
dates=seq.Date(as.Date("1860/1/1"), as.Date("2012/12/31"),"days")

dates2=cbind(as.numeric(format.Date(dates,"%Y%m%d")),matrix(0,length(dates),5))

for(i in 1:length(data[,1]))
{
  I=which(dates2[,1]==data$Date[i])
  for(j in 1:data$Days[i])
    dates2[I+j-1,2]=as.numeric(data$Type2[i])
}

colnames(dates2)=c("Date","Event","Ens Mems","Ens mean","MLDB","ERAI")

## Step 3 - Add columns looking at proportion of members that identify a given day

years=seq(1871,2012)
decades=seq(1871,2001,10)
annfreq=matrix(0,length(years),57)
decadalfreq=matrix(0,length(decades),57)
data$ERAI<-data$MLDB<-data$EnsMean<-data$EnsMems<-0

#rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/20CR/ECLs/")
data=read.csv("/srv/ccrc/data34/z3478332/20CR/ECLs/Power_ECLs.csv")
data$ERAI<-data$MLDB<-data$EnsMean<-data$EnsMems<-0

##EnsMem
data$EnsMem=0
decades=seq(1871,2001,10)
Y=which(data$Year>=1871 & data$Year<=2012)

dates=seq.Date(as.Date("1871/1/1"), as.Date("2012/12/31"),"days")
dates=cbind(as.numeric(format.Date(dates,"%Y%m%d")),matrix(0,length(dates)))

for(m in 1:56)
{
  print(m)
  dates=seq.Date(as.Date("1871/1/1"), as.Date("2012/12/31"),"days")
  dates=cbind(as.numeric(format.Date(dates,"%Y%m%d")),rep(0,length(dates)))
  for(y in 1:length(decades))
  {
    tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_",decades[y],"-",decades[y]+9,"_",m,".csv",sep=""))
    if(y==1) tmp=rbind(tmp,read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_1871_",m,".csv",sep=""))
    )
    decadalfreq[y,m]=length(unique(tmp$ID))
    I=which(tmp$Location==1)
    tmpday=unique(tmp$Date[I])
    I=which(dates[,1]%in%tmpday)
    dates[I,2]=1
    I=which(dates2[,1]%in%tmpday)
    dates2[I,3]=dates2[I,3]+1
    
    if(y==1) tmp2=tmp else tmp2=rbind(tmp2,tmp)
  }
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_2011-2012_",m,".csv",sep=""))
  I=which(tmp$Location==1)
  tmpday=unique(tmp$Date[I])
  I=which(dates[,1]%in%tmpday)
  dates[I,2]=1
  I=which(dates2[,1]%in%tmpday)
  dates2[I,3]=dates2[I,3]+1
  tmp2=rbind(tmp2,tmp)
  
  for(i in Y[1]:Y[length(Y)])
  {
    I=which(dates[,1]==data$Date[i])
    I2=seq((I-1),(I+data$Days[i]))
    data$EnsMem[i]=data$EnsMem[i]+max(dates[I2,2])
  }
  
  yy=floor(tmp2$Date/10000)
  for(kk in 1:length(years))
  {
    I=which(kk==years[kk] & tmp2$Fix==1)
    annfreq[kk,m]=length(unique(tmp2$ID[I]))
  }
}


### EnsMean

decades=seq(1871,2001,10)
dates=seq.Date(as.Date("1871/1/1"), as.Date("2012/12/31"),"days")
dates=cbind(as.numeric(format.Date(dates,"%Y%m%d")),rep(0,length(dates)))
for(y in 1:length(decades))
{
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_mean/tracks/ECLfixes_",decades[y],"-",decades[y]+9,"_mean.csv",sep=""))
  decadalfreq[y,57]=length(unique(tmp$ID))
  I=which(tmp$Location==1)
  tmpday=unique(tmp$Date[I])
  I=which(dates[,1]%in%tmpday)
  dates[I,2]=1
  I=which(dates2[,1]%in%tmpday)
  dates2[I,4]=1
  
  if(y==1) tmp2=tmp else tmp2=rbind(tmp2,tmp)
}
tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_mean/tracks/ECLfixes_2011-2012_mean.csv",sep=""))
I=which(tmp$Location==1)
tmpday=unique(tmp$Date[I])
I=which(dates[,1]%in%tmpday)
dates[I,2]=1
I=which(dates2[,1]%in%tmpday)
dates2[I,4]=1
tmp2=rbind(tmp2,tmp)

Y=which(data$Year>=1871 & data$Year<=2012)
for(i in Y[1]:Y[length(Y)])
{
  I=which(dates[,1]==data$Date[i])
  I2=seq((I-1),(I+data$Days[i]))
  data$EnsMean[i]=max(dates[I2,2])
}

yy=floor(tmp2$Date/10000)
for(kk in 1:length(years))
{
  I=which(yy==years[kk] & tmp2$Fix==1)
  annfreq[kk,57]=length(unique(tmp2$ID[I]))
}

### MLDB

Y=which(data$Year>=1970 & data$Year<=2006)
mldb=read.csv("MLDB.csv",sep=";")
dates=seq.Date(as.Date("1970/1/1"), as.Date("2006/12/31"),"days")
dates=cbind(as.numeric(format.Date(dates,"%Y%m%d")),rep(0,length(dates)))
I=which(dates[,1]%in%mldb$date)
dates[I,2]=1
I=which(dates2[,1]%in%mldb$date)
dates2[I,5]=1

for(i in Y[1]:Y[length(Y)])
{
  I=which(dates[,1]==data$Date[i])
  I2=seq((I-1),(I+data$Days[i]))
  data$MLDB[i]=max(dates[I2,2])
}


### ERAI

Y=which(data$Year>=1980 & data$Year<=2012)
um=read.csv("UM_rad2_p100.csv")
dates=seq.Date(as.Date("1980/1/1"), as.Date("2012/12/31"),"days")
dates=cbind(as.numeric(format.Date(dates,"%Y%m%d")),rep(0,length(dates)))
I=which(um$Location==1)
I=which(dates[,1]%in%um$Date[I])
dates[I,2]=1
I=which(um$Location==1)
I=which(dates2[,1]%in%um$Date[I])
dates2[I,6]=1

for(i in Y[1]:Y[length(Y)])
{
  I=which(dates[,1]==data$Date[i])
  I2=seq((I-1),(I+data$Days[i]))
  data$ERAI[i]=max(dates[I2,2])
}


save(decadalfreq,annfreq,dates2,data,file="20CR_ECLdates_v3.RData")
