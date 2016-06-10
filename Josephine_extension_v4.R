###########
## Step 1 - COnvert Power databse of events into a dataset of ECL/event by day
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/20CR/ECLs/")

## Step 3 - Add columns looking at proportion of members that identify a given day

years=seq(1871,2012)
decades=seq(1871,2001,10)
annfreq3=array(0,c(length(years),57,4))
decadalfreq3=array(0,c(length(decades),57,4))
thresh=c(1,1.25,1.5,2)

decades=seq(1871,2001,10)
dates=seq.Date(as.Date("1871/1/1"), as.Date("2012/12/31"),"days")
dates=cbind(as.numeric(format.Date(dates,"%Y%m%d")),matrix(0,length(dates)))

for(m in 1:56)
{
  print(m)
  for(y in 1:length(decades))
  {
    tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLevents_",decades[y],"-",decades[y]+9,"_",m,".csv",sep=""))
    if(y==1) tmp=rbind(tmp,read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLevents_1871_",m,".csv",sep=""))
    )
    if(y==1) tmp2=tmp else tmp2=rbind(tmp2,tmp)
  }
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLevents_2011-2012_",m,".csv",sep=""))
  tmp2=rbind(tmp2,tmp)
  
  yy=floor(tmp2$Date1/10000)
  for(kk in 1:length(years))
    for(t in 1:4)
  {
    I=which(yy==years[kk] & tmp2$CV2>=thresh[t])
    annfreq3[kk,m,t]=length(I)
  }
}


### EnsMean

for(y in 1:length(decades))
{
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_mean/tracks/ECLevents_",decades[y],"-",decades[y]+9,"_mean.csv",sep=""))
  if(y==1) tmp2=tmp else tmp2=rbind(tmp2,tmp)
}
tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_mean/tracks/ECLevents_2011-2012_mean.csv",sep=""))
tmp2=rbind(tmp2,tmp)

yy=floor(tmp2$Date1/10000)
for(kk in 1:length(years))
  for(t in 1:4)
  {
    I=which(yy==years[kk] & tmp2$CV2>=thresh[t])
    annfreq3[kk,57,t]=length(I)
  }

for(i in 1:length(decades))
{
  I=which(years>=decades[i] & years<(decades[i]+10))
  decadalfreq3[i,,]=apply(annfreq3[I,,],c(2,3),sum)
}

save(decadalfreq3,annfreq3,file="20CR_ECLdates_thresh.RData")
