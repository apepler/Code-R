rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")
names=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
time=c(9009,6079)
thresh=22
n=1

bombprop=array(0,c(15,4,2))

for(t in 1)
{
  n=1
  
  for(i in 5)
    for(j in 1:3)
      
    {    
      filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_",time[t],".csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_",time[t],".csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_",time[t],"_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_",time[t],"_pg1.3.csv",sep=""))
      
      filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_",time[t],".csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_",time[t],".csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_",time[t],"_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_",time[t],"_pg1.3.csv",sep=""))
      
      filelistE2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_",time[t],"_v2.csv",sep=""),
                   paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_",time[t],"_v2.csv",sep=""),
                   paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_",time[t],"_pg0.8_v2.csv",sep=""),
                   paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_",time[t],"_pg1.3_v2.csv",sep=""))
      
      for(ff in 1:4)
      {
        data=read.csv(filelist1[ff])      
        events=read.csv(filelistE[ff])     
        
        data$Location2=0
        I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
        data$Location2[I]<-1
        I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
        data$Location2[I]<-1
        I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
        data$Location2[I]<-1
        
        data$NDR=NaN
        I=which(data$ID[2:length(data$ID)]==data$ID[1:length(data$ID)-1])+1
        data$NDR[I]= (data$MSLP[I]-data$MSLP[I-1])*sin(60*pi/180)/(6*sin(data$Lat[I]*pi/180))
        
        events$NDR<-events$Length3<-NaN
        if(ff<3) events$CV_loc2=NaN else events$PG_loc2=NaN
        
        for(ee in 1:length(events[,1]))
        {
          I=which(data$ID==events$ID[ee] & data$Location==1 & !is.na(data$NDR))
          if(length(I)>0) events$NDR[ee]=max(data$NDR[I],na.rm=T)
          
          I=which(data$ID==events$ID[ee] & data$Location2==1)
          events$Length3[ee]=length(I)
          if(length(I)>0)
            if(ff<3) events$CV_loc2[ee]=max(data$CV[I]) else events$PG_loc2[ee]=max(data$PG200[I])
        }
        
        write.csv(events,file=filelistE2[ff])
      }
      n=n+1
    }
}


######## Okay, now do properly, w/ v2

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

intcol=c(11,11,11,11)

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")
thresh=22

now<-array(0,c(15,12,4))
dimnames(now)[[1]]<-rep("ncep",15)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
bombN<-bombN_na<-lingerN<-linger2N<-now
future<-bombF<-bombF_na<-now[1:12,,]
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_v2.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8_v2.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3_v2.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079_v2.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8_v2.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3_v2.csv",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelist1[ff])
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      data2=data[data[,intcol[ff]]>=thresh2,]
      mm=floor((data2$Date1%%10000)/100)
      for(m in 1:12)
      {
        I=which(mm==m)
        now[n,m,ff]=length(I)
        J=which(data2$NDR[I]>=1)
        bombN[n,m,ff]=length(J)
        J=which(!is.na(data2$NDR[I]))
        bombN_na[n,m,ff]=length(J)
      }
      
      if(i<5)
      {
        data=read.csv(filelist2[ff])
        data2=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          future[n,m,ff]=length(I)
          J=which(data2$NDR[I]>=1)
          bombF[n,m,ff]=length(J)
          J=which(!is.na(data2$NDR[I]))
          bombF_na[n,m,ff]=length(J)
        } 
      }
    }     
    n=n+1  
  }

########## Look at % with duration 5 in Loc1 & Loc2

intcol=c(11,11,11,11)

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")
thresh=22

now<-array(0,c(15,12,4))
dimnames(now)[[1]]<-rep("ncep",15)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
lingerN<-linger2N<-now
future<-lingerF<-linger2F<-now[1:12,,]
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_v2.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8_v2.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3_v2.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079_v2.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8_v2.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3_v2.csv",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelist1[ff])
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      data2=data[data[,intcol[ff]]>=thresh2,]
      mm=floor((data2$Date1%%10000)/100)
      for(m in 1:12)
      {
        I=which(mm==m)
        now[n,m,ff]=length(I)
        J=which(data2$Length2[I]>=5)
        lingerN[n,m,ff]=length(J)
        J=which(data2$Length3[I]>=5)
        linger2N[n,m,ff]=length(J)
      }
      
      if(i<5)
      {
        data=read.csv(filelist2[ff])
        data2=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          future[n,m,ff]=length(I)
          J=which(data2$Length2[I]>=5)
          lingerF[n,m,ff]=length(J)
          J=which(data2$Length3[I]>=5)
          linger2F[n,m,ff]=length(J)
        } 
      }
    }     
    n=n+1  
  }


prop=apply(lingerN,c(1,3),sum)/apply(now,c(1,3),sum)
prop2=apply(linger2N,c(1,3),sum)/apply(now,c(1,3),sum)