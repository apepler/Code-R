rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)

now<-array(0,c(16,12,5))
dimnames(now)[[1]]<-rep("ncep",16)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now[1:12,,]
intcolE<-intcolF<-10

fixesN<-eventsN<-fixesF<-eventsF<-list()

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")

#### Version 1 - all above threshold
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
    filelistF=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""))
    filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""))

    ff=1
    eventsN[[n]]=read.csv(filelistE[ff])
    fixesN[[n]]=read.csv(filelistF[ff])
    
    eventsN[[n]]$Year=floor(eventsN[[n]]$Date1/10000)
    eventsN[[n]]$Month=floor(eventsN[[n]]$Date1/100)%%100
    fixesN[[n]]$Year=floor(fixesN[[n]]$Date/10000)
    fixesN[[n]]$Month=floor(fixesN[[n]]$Date/100)%%100
    
    fixesN[[n]]$Location2<-0
    I<-which(fixesN[[n]][,7]>=149 & fixesN[[n]][,7]<=154 & fixesN[[n]][,8]<(-37) & fixesN[[n]][,8]>=-41)
    fixesN[[n]]$Location2[I]<-1
    I<-which(fixesN[[n]][,7]>=(149+(37+fixesN[[n]][,8])/2) & fixesN[[n]][,7]<=(154+(37+fixesN[[n]][,8])/2) & fixesN[[n]][,8]<(-31) & fixesN[[n]][,8]>=-37)
    fixesN[[n]]$Location2[I]<-1
    I<-which(fixesN[[n]][,7]>=152 & fixesN[[n]][,7]<=157 & fixesN[[n]][,8]<=(-24) & fixesN[[n]][,8]>=-31)
    fixesN[[n]]$Location2[I]<-1
    
    eventsN[[n]]$CV_loc2<-eventsN[[n]]$Location2<-0
    for(k in 1:length(eventsN[[n]]$ID))
    {
      I=which(fixesN[[n]]$ID==eventsN[[n]]$ID[k] & fixesN[[n]]$Location2==1)
      eventsN[[n]]$Location2[k]=length(I)
      if(length(I)>0) eventsN[[n]]$CV_loc2[k]=max(fixesN[[n]]$CV[k])
    }
    
    if(i<5)
    {
      ff=2
      eventsF[[n]]=read.csv(filelistE[ff])
      fixesF[[n]]=read.csv(filelistF[ff])
      eventsF[[n]]$Year=floor(eventsF[[n]]$Date1/10000)
      eventsF[[n]]$Month=floor(eventsF[[n]]$Date1/100)%%100
      fixesF[[n]]$Year=floor(fixesF[[n]]$Date/10000)
      fixesF[[n]]$Month=floor(fixesF[[n]]$Date/100)%%100
      
      fixesF[[n]]$Location2<-0
      I<-which(fixesF[[n]][,7]>=149 & fixesF[[n]][,7]<=154 & fixesF[[n]][,8]<(-37) & fixesF[[n]][,8]>=-41)
      fixesF[[n]]$Location2[I]<-1
      I<-which(fixesF[[n]][,7]>=(149+(37+fixesF[[n]][,8])/2) & fixesF[[n]][,7]<=(154+(37+fixesF[[n]][,8])/2) & fixesF[[n]][,8]<(-31) & fixesF[[n]][,8]>=-37)
      fixesF[[n]]$Location2[I]<-1
      I<-which(fixesF[[n]][,7]>=152 & fixesF[[n]][,7]<=157 & fixesF[[n]][,8]<=(-24) & fixesF[[n]][,8]>=-31)
      fixesF[[n]]$Location2[I]<-1
      
      eventsF[[n]]$CV_loc2<-eventsF[[n]]$Location2<-0
      for(k in 1:length(eventsF[[n]]$ID))
      {
        I=which(fixesF[[n]]$ID==eventsF[[n]]$ID[k] & fixesF[[n]]$Location2==1)
        eventsF[[n]]$Location2[k]=length(I)
        if(length(I)>0) eventsF[[n]]$CV_loc2[k]=max(fixesF[[n]]$CV[k])
      }
    }
    
    n=n+1  
  }

now<-array(NaN,c(5,3,4,2))
dimnames(now)[[3]]<-c("Events","Days Loc1","Events 2Loc2","Days Loc2")
dimnames(now)[[4]]<-c("Warm","Cold")
future=now
thresh=22

mlist=cbind(c(1:4,11:12),5:10)

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    b=order(eventsN[[n]]$CVmax,decreasing=T)
    thresh2=eventsN[[n]]$CVmax[b[20*thresh]]
    if(is.na(thresh2)) thresh2=min(eventsN[[n]]$CVmax,na.rm=T)
    
    for(m in 1:2)
    {
    now[i,j,1,m]=length(which(eventsN[[n]]$CVmax>=thresh2 & eventsN[[n]]$Month%in%mlist[,m]))
    I=which(fixesN[[n]]$Location==1 & fixesN[[n]]$CV>=thresh2 & fixesN[[n]]$Month%in%mlist[,m])
    now[i,j,2,m]=length(unique(fixesN[[n]]$Date[I]))
    now[i,j,3,m]=length(which(eventsN[[n]]$CVmax>=thresh2 & eventsN[[n]]$Location2>1 & eventsN[[n]]$Month%in%mlist[,m]))
    I=which(fixesN[[n]]$Location2==1 & fixesN[[n]]$CV>=thresh2 & fixesN[[n]]$Month%in%mlist[,m])
    now[i,j,4,m]=length(unique(fixesN[[n]]$Date[I]))
    
    if(i<5)
    {
      future[i,j,1,m]=length(which(eventsF[[n]]$CVmax>=thresh2 & eventsF[[n]]$Month%in%mlist[,m]))
      I=which(fixesF[[n]]$Location==1 & fixesF[[n]]$CV>=thresh2 & fixesF[[n]]$Month%in%mlist[,m])
      future[i,j,2,m]=length(unique(fixesF[[n]]$Date[I]))
      future[i,j,3,m]=length(which(eventsF[[n]]$CVmax>=thresh2 & eventsF[[n]]$Location2>1 & eventsF[[n]]$Month%in%mlist[,m]))
      I=which(fixesF[[n]]$Location2==1 & fixesF[[n]]$CV>=thresh2 & fixesF[[n]]$Month%in%mlist[,m])
      future[i,j,4,m]=length(unique(fixesF[[n]]$Date[I]))
    }
    }
    
  n=n+1
  }

a=future/now
apply(a,c(3,4),mean,na.rm=T)
