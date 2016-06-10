###########
## Step 1 - COnvert Power databse of events into a dataset of ECL/event by day
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/20CR/ECLs/")

## Step 3 - Add columns looking at proportion of members that identify a given day

years=seq(1871,2012)
decades=seq(1871,2001,10)
#rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/20CR/ECLs/")
decades=seq(1871,2001,10)

dates=seq.Date(as.Date("1871/1/1"), as.Date("2012/12/31"),"days")
dates<-cv<-cbind(as.numeric(format.Date(dates,"%Y%m%d")),matrix(0,length(dates),56))

for(m in 1:56)
{
  print(m)
  for(y in 1:length(decades))
  {
    tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_",decades[y],"-",decades[y]+9,"_",m,".csv",sep=""))
    if(y==1) tmp=rbind(tmp,read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_1871_",m,".csv",sep="")))
    I=which(tmp$Location==1)
    tmpday=unique(tmp$Date[I])
    I=which(dates[,1]%in%tmpday)
    dates[I,1+m]=1
    
    for(k in 1:length(tmpday))
    {
      I=which(dates[,1]%in%tmpday[k])
      J=which(tmp$Date==tmpday[k] & tmp$Location==1)
      cv[I,1+m]=max(tmp$CV[J],na.rm=T)
    }
  
  }
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_2011-2012_",m,".csv",sep=""))
  I=which(tmp$Location==1)
  tmpday=unique(tmp$Date[I])
  I=which(dates[,1]%in%tmpday)
  dates[I,1+m]=1
  for(k in 1:length(tmpday))
  {
    I=which(dates[,1]%in%tmpday[k])
    J=which(tmp$Date==tmpday[k] & tmp$Location==1)
    cv[I,1+m]=max(tmp$CV[J],na.rm=T)
  }
}

save(dates,cv,file="20CR_ECLdates_ens.RData")
