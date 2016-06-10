###########
## Step 1 - COnvert Power databse of events into a dataset of ECL/event by day
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/20CR/ECLs/")

## Step 3 - Add columns looking at proportion of members that identify a given day

years=seq(1851,2014)
starts=seq(1851,2001,10)
ends=c(seq(1860,2000,10),2014)

freq_ann=matrix(NaN,length(years),60)
freq_mon=matrix(NaN,length(years)*12,60)

colnames(freq_ann)<-colnames(freq_mon)<-c(paste("Member",1:56),"Mean of Members","Ensemble Mean","ERAI","NCEP1")
rownames(freq_ann)<-years
rownames(freq_mon)<-rep("aaa",length(years)*12)

for(y in 1:length(starts))
{
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_mean/tracks_v2c/ECLevents_",starts[y],"-",ends[y],"_mean.csv",sep=""))
  if(y==1) tmp2=tmp else tmp2=rbind(tmp2,tmp)
}

yy=floor(tmp2$Date1/10000)
mm=floor(tmp2$Date1/100)%%100
n=1
for(kk in 1:length(years))
  {
    I=which(yy==years[kk])
    freq_ann[kk,58]=length(I)
    for(m in 1:12)
    {
      I=which(yy==years[kk] & mm==m)
      freq_mon[n,58]=length(I)
      rownames(freq_mon)[n]=paste(years[kk],m)
      n=n+1
    }
  }

for(mem in 1:56)
{
  for(y in 1:length(starts))
  {
    tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks_v2c/ECLevents_",starts[y],"-",ends[y],"_",mem,".csv",sep=""))
    if(y==1) tmp2=tmp else tmp2=rbind(tmp2,tmp)
  }
  
  yy=floor(tmp2$Date1/10000)
  mm=floor(tmp2$Date1/100)%%100
  n=1
  for(kk in 1:length(years))
  {
    I=which(yy==years[kk])
    freq_ann[kk,mem]=length(I)
    for(m in 1:12)
    {
      I=which(yy==years[kk] & mm==m)
      freq_mon[n,mem]=length(I)
      n=n+1
    }
  }
}

freq_ann[,57]=apply(freq_ann[,1:56],1,mean)
freq_mon[,57]=apply(freq_mon[,1:56],1,mean)

tmp2=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100_update.csv")
yy=floor(tmp2$Date1/10000)
mm=floor(tmp2$Date1/100)%%100
J=which(years==1980)
n=(J-1)*12+1
for(kk in J:length(years))
  {
    I=which(yy==years[kk])
    freq_ann[kk,59]=length(I)
    for(m in 1:12)
    {
      I=which(yy==years[kk] & mm==m)
      freq_mon[n,59]=length(I)
      n=n+1
    }
  }

tmp2=rbind(read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ncep1_topo_rad2_proj100/ECLevents_umelb_ncep1_topo_rad2_proj100_4970.csv"),
read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ncep1_topo_rad2_proj100/ECLevents_umelb_ncep1_topo_rad2_proj100_7190.csv"),
read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ncep1_topo_rad2_proj100/ECLevents_umelb_ncep1_topo_rad2_proj100_9114.csv"))
yy=floor(tmp2$Date1/10000)
mm=floor(tmp2$Date1/100)%%100
J=which(years==1949)
n=(J-1)*12+1
for(kk in J:length(years))
  {
    I=which(yy==years[kk])
    freq_ann[kk,60]=length(I)
    for(m in 1:12)
    {
      I=which(yy==years[kk] & mm==m)
      freq_mon[n,60]=length(I)
      n=n+1
    }
  }

write.csv(freq_ann,file="LAP_AnnualFreq_20CR_v2c.csv")
write.csv(freq_mon,file="LAP_MonthlyFreq_20CR_v2c.csv")








