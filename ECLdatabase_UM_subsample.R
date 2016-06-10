##This is code to take the output from the low tracking software 
##and turn it into a .csv file of ECL track data.

##Note that this script performs only the initial filtering 
##to restrict ECLs based on location, track length and at least one closed fix
##To obtain a better database of ECLs with impacts you can filter for 
##a fix curvature of >= 0.25 & the fix being in the ECL region

##This version of the code requires both the annual version of the data
##as well as that run for the previous summer, to collate the tracks for events
##over January 1.

##Note also that R has issues with large datasets - 
##I tend to separate into sets of 20 years or fewer. Pretty quick. 

##yearS and yearE are start & end years. Must have DJF data for yearS-1
##output is a string to help label the fix & event files (so use "name")

##You will need to set the location to where the track files are stored.
#setwd('~/Charles_bom/output')

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
period=c("9009","6079")
proj=c(240,100)

setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/")
#thresh=1.3

i=3
j=2
k=2
l=1
for(t in seq(0.6,1.6,0.1))
{
  indir=paste("proj",proj[k],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/",sep="")
  infile=paste("ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[k],"_rad2cv06_",period[l],".csv",sep="")
  outdir=paste("proj",proj[k],"/test/",sep="")
  outfile=paste("umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[k],"_rad2cv",t,"_",period[l],".csv",sep="")
  ECLdatabase(indir,infile,outdir,outfile,t)
}

count=data.frame(Thresh=seq(0.6,1.6,0.1),Count=rep(0,11))
for(t in 1:11)
{
  outfile=paste("umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[k],"_rad2cv",count[t,1],"_",period[l],".csv",sep="")  
  a=read.csv(paste("proj",proj[k],"/test/ECLevents_",outfile,".csv",sep=""))
  count[t,2]=length(a[,1])/20
}

thresh=0.8
for(j in 1:3)
  for(k in 2)
    for(l in 1)
    {
      indir=paste("proj",proj[k],"/outputUM_ncep_WRF",wrf[j],"_50_rad2cv06/",sep="")
      infile=paste("ECLfixes_umelb_ncep_wrf",wrf[j],"_proj",proj[k],"_rad2cv06_",period[l],".csv",sep="")
      outdir=indir
      outfile=paste("umelb_ncep_wrf",wrf[j],"_proj",proj[k],"_rad2cv",thresh,"_",period[l],sep="")
      ECLdatabase(indir,infile,outdir,outfile,thresh)
    }

for(i in 1:4)
  for(j in 1:3)
    for(k in 2)
      for(l in 1:2)
      {
        indir=paste("proj",proj[k],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/",sep="")
        infile=paste("ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[k],"_rad2cv06_",period[l],".csv",sep="")
        outdir=indir
        outfile=paste("umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[k],"_rad2cv",thresh,"_",period[l],sep="")
        ECLdatabase(indir,infile,outdir,outfile,thresh)
      }


ECLdatabase<-function(indir,infile,outdir,outfile,thresh)
{
  print(outfile)
  data=read.csv(paste(indir,infile,sep=""))
  data=data[,2:13]
  
  data3=data[data$CV>=thresh,]
  data3$Fix2<-data3$ID2<-rep(0,length(data3$ID))
  data3$ID2[1]<-data3$Fix2[1]<-1
  
  for(i in 2:length(data3$ID))
    if(data3$ID[i]==data3$ID[i-1] & data3$Fix[i]==data3$Fix[i-1]+1) {
      data3$ID2[i]=data3$ID2[i-1]
      data3$Fix2[i]=data3$Fix2[i-1]+1
    } else {
      data3$ID2[i]=data3$ID2[i-1]+1
      data3$Fix2[i]=1
    }
  
  data3[,1:2]=data3[,13:14]
  data2=data3[,1:12]    
  
  ##Collates information for events
  x<-rle(data2[,1])
  events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Date2=rep(0,length(x$values)),
                     MSLP=rep(0,length(x$values)),CV=rep(0,length(x$values)),Loc=rep(0,length(x$values)),Open=rep(0,length(x$values)))
  
  ##Single-fix events
  
  I<-which(events[,2]==1)
  if(length(I)>0)
  {
    y<-match(events[I,1],data2[,1])
    events=events[-I,]
    fixes=data2[-y,]
  } else fixes=data2
  
  for(i in 1:length(events[,1])) 
  {
    I<-which(fixes[,1]==events[i,1])
    events[i,3]=min(fixes[I,3]) ##Date1
    events[i,4]=max(fixes[I,3]) ##Date2
    events[i,5]=min(fixes[I,8]) ##Min central pressure
    events[i,6]=max(fixes[I,9]) ##Max PG
    if(sum(fixes[I,10])>0) events[i,7]=1 ##Ever in ECL region
    events[i,8]=min(fixes[I,5]) ##Ever closed  
  }
  
  I<-which(events[,7]>0 & events[,8]==0) ## At least one fix closed & one in ECL region
  events2<-events[I,1:6]
  include<-match(fixes[,1],events2[,1])
  J<-which(is.na(include)==0)
  fixes2=fixes[J,]
  
  ##Fix numbering and save final data
  
  eventnum=events2[,1]  
  for(i in 1:length(eventnum))
  {
    I<-which(fixes2[,1]==eventnum[i])
    fixes2[I,1]=i
    events2[i,1]=i
  }
  
  ##Add extra columns - count,mslp,cv,rad in loc
  events2$PG2<-events2$MSLP2<-events2$Length2<-rep(0,length(events2$ID))
  for(i in 1:length(events2[,1]))
  {
    I<-which(fixes2$ID==events2[i,1] & fixes2$Location==1)
    events2$Length2[i]=length(I)
    events2$MSLP2[i]=min(fixes2$MSLP[I])
    events2$PG2[i]=max(fixes2$CV[I])
  }
  
  outF=paste(outdir,'ECLfixes_',outfile,'.csv',sep="")
  outE=paste(outdir,'ECLevents_',outfile,'.csv',sep="")
  
  write.csv(fixes2,outF)
  write.csv(events2,outE)
}


###### Add CVmean

######### Need to re-edit the events/fixes so we have CV2 (not PG2) and MEAN CV as well
rm(list=ls())
cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
proj=c(100,100,240,240)
thresh1=c("06","06","1","1")
thresh2=c(0.8,0.8,1.35,1.35)
time=c(9009,6079,9009,6079)

for(i in 1:5)
  for(j in 1:3)
    for(k in 1:4)
    {
      dir=paste("outputUM/proj",proj[k],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv",thresh1[k],"/",sep="")
      file1=paste("ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[k],"_rad2cv",thresh2[k],"_",time[k],".csv",sep="")
      file2=paste("ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[k],"_rad2cv",thresh2[k],"_",time[k],".csv",sep="")
      
      data1=read.csv(paste(dir,file1,sep=""))
      data1=data1[,-1]
      data2=read.csv(paste(dir,file2,sep=""))
      data2=data2[,-1]
      
      colnames(data1)[9]="CVmax"
      
      data1$Rmean<-data1$CVmean<-0      
      for(n in 1:length(data1[,1]))
      {
        I=which(data2$ID==data1$ID[n] & data2$Location==1)
        data1[n,10]=mean(data2$CV[I])
        data1[n,11]=mean(data2$Radius[I])
      }
      colnames(data1)[10]="CVmean"  
      colnames(data1)[11]="Rmean"      
      write.csv(data1,file=paste(dir,file1,sep=""))
    }

