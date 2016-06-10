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
    filelistC=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_compositerain_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.nc",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_compositerain_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.nc",sep=""))
    filelistW=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_compositewind_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.nc",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_compositewind_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.nc",sep=""))
    
    
    ff=1
    eventsN[[n]]=read.csv(filelistE[ff])
    data=read.csv(filelistF[ff])
    a=open.nc(filelistC[ff])
    tmp=var.get.nc(a,"ECLrain")
    tmp[tmp>=600]=NaN
    data$MeanRain=apply(tmp,3,mean,na.rm=T)
    data$MaxRain=apply(tmp,3,max,na.rm=T)
    close.nc(a)
    a=open.nc(filelistW[ff])
    tmp=sqrt(var.get.nc(a,"ECL_U10")^2+var.get.nc(a,"ECL_V10")^2)
    if(max(tmp,na.rm=T)>=100) tmp=sqrt((var.get.nc(a,"ECL_U10")/6)^2+(var.get.nc(a,"ECL_V10")/6)^2)
    data$MeanWind=apply(tmp,3,mean,na.rm=T)
    data$MaxWind=apply(tmp,3,max,na.rm=T)
    close.nc(a)
    fixesN[[n]]=data
    
    eventsN[[n]]$MaxPointRain<-eventsN[[n]]$MaxMeanRain<-eventsN[[n]]$TotalRain<-0
    for(k in 1:length(eventsN[[n]]$ID))
    {
      I=which(fixesN[[n]]$ID==eventsN[[n]]$ID[k] & fixesN[[n]]$Location==1)
      eventsN[[n]]$TotalRain[k]=sum(fixesN[[n]]$MeanRain[I],na.rm=T)
      eventsN[[n]]$MaxMeanRain[k]=max(fixesN[[n]]$MeanRain[I],na.rm=T)
      eventsN[[n]]$MaxPointRain[k]=max(fixesN[[n]]$MaxRain[I],na.rm=T)
      eventsN[[n]]$MaxMeanWind[k]=max(fixesN[[n]]$MeanWind[I],na.rm=T)
      eventsN[[n]]$MaxPointWind[k]=max(fixesN[[n]]$MaxWind[I],na.rm=T)
    }
    
    if(i<5)
    {
      ff=2
      eventsF[[n]]=read.csv(filelistE[ff])
      data=read.csv(filelistF[ff])
      a=open.nc(filelistC[ff])
      tmp=var.get.nc(a,"ECLrain")
      tmp[tmp>=600]=NaN
      data$MeanRain=apply(tmp,3,mean,na.rm=T)
      data$MaxRain=apply(tmp,3,max,na.rm=T)
      close.nc(a)

      a=open.nc(filelistW[ff])
      tmp=sqrt(var.get.nc(a,"ECL_U10")^2+var.get.nc(a,"ECL_V10")^2)
      if(max(tmp,na.rm=T)>=100) tmp=sqrt((var.get.nc(a,"ECL_U10")/6)^2+(var.get.nc(a,"ECL_V10")/6)^2)
      data$MeanWind=apply(tmp,3,mean,na.rm=T)
      data$MaxWind=apply(tmp,3,max,na.rm=T)
      close.nc(a)
            
      fixesF[[n]]=data
      eventsF[[n]]$MaxPointRain<-eventsF[[n]]$MaxMeanRain<-eventsF[[n]]$TotalRain<-0
      for(k in 1:length(eventsF[[n]]$ID))
      {
        I=which(fixesF[[n]]$ID==eventsF[[n]]$ID[k] & fixesF[[n]]$Location==1)
        eventsF[[n]]$TotalRain[k]=sum(fixesF[[n]]$MeanRain[I],na.rm=T)
        eventsF[[n]]$MaxMeanRain[k]=max(fixesF[[n]]$MeanRain[I],na.rm=T)
        eventsF[[n]]$MaxPointRain[k]=max(fixesF[[n]]$MaxRain[I],na.rm=T)
        eventsF[[n]]$MaxMeanWind[k]=max(fixesF[[n]]$MeanWind[I],na.rm=T)
        eventsF[[n]]$MaxPointWind[k]=max(fixesF[[n]]$MaxWind[I],na.rm=T)
      }    
    }
    
    n=n+1  
  }

for(thresh in c(22,15,10,5,2))
{
  countW=array(0,c(12,8,2))
  clist=c(9:11,13:17)
  dimnames(countW)[[1]]<-dimnames(future)[[1]]
  dimnames(countW)[[2]]<-c("Min MSLP","Max Lap","Mean Lap","TotalRain","MaxMeanRain","MaxPointRain","MaxMeanWind","MaxPointWind")
  dimnames(countW)[[3]]<-c("1990-2009","2060-2079")
  countC=countW
  
  for(i in 1:12)
  {
    data=eventsN[[i]]
    mm=floor(data$Date1/100)%%100
    data2=eventsF[[i]]
    mm2=floor(data2$Date1/100)%%100
    for(j in 1)
    {
      b=order(data[,clist[j]],decreasing=F)
      thresh2=data[b[20*thresh],clist[j]]
      if(is.na(thresh2)) thresh2=max(data[,intcolE])
      I=which(data[,clist[j]]<=thresh2 & mm>=5 & mm<=10)
      countC[i,j,1]=length(I)
      I=which(data[,clist[j]]<=thresh2 & (mm>=11 | mm<=4))
      countW[i,j,1]=length(I)    
      I=which(data2[,clist[j]]<=thresh2 & mm2>=5 & mm2<=10)
      countC[i,j,2]=length(I)
      I=which(data2[,clist[j]]<=thresh2 & (mm2>=11 | mm2<=4))
      countW[i,j,2]=length(I)      
    }
    for(j in 2:8)
    {
      b=order(data[,clist[j]],decreasing=T)
      thresh2=data[b[20*thresh],clist[j]]
      if(is.na(thresh2)) thresh2=min(data[,intcolE])
      I=which(data[,clist[j]]>=thresh2 & mm>=5 & mm<=10)
      countC[i,j,1]=length(I)
      I=which(data[,clist[j]]>=thresh2 & (mm>=11 | mm<=4))
      countW[i,j,1]=length(I)    
      I=which(data2[,clist[j]]>=thresh2 & mm2>=5 & mm2<=10)
      countC[i,j,2]=length(I)
      I=which(data2[,clist[j]]>=thresh2 & (mm2>=11 | mm2<=4))
      countW[i,j,2]=length(I)    
    }
  }
  
  change=abind(100*((countC[,,2]/countC[,,1])-1),100*((countW[,,2]/countW[,,1])-1),along=3)
  dimnames(change)[[3]]<-c("Cool","Warm")
  names(dimnames(change))<-c("Source","Metric","Season")
  
  pdf(file=paste("LAP_ECLevents_change_byseason_vsintensitymetric_thresh",thresh,".pdf",sep=""),width=10,height=4,pointsize=10)
  tmp=melt(change)
  print(ggplot(tmp, aes(x = Metric, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-500, 500, 25)) +
          theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))
  dev.off()
  pdf(file=paste("LAP_ECLevents_change_byseason_vsintensitymetric_thresh",thresh,"_trunc.pdf",sep=""),width=10,height=4,pointsize=10)
  tmp=melt(change)
  print(ggplot(tmp, aes(x = Metric, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 100, 25),limits = c(-50, 100)) +
          theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))
  
  dev.off()
}

for(thresh in c(22,15,10,5,2))
{
  meanW=array(0,c(12,5,2))
  clist=c(13:17)
  dimnames(meanW)[[1]]<-dimnames(future)[[1]]
  dimnames(meanW)[[2]]<-colnames(eventsN[[1]])[clist]
  dimnames(meanW)[[3]]<-c("1990-2009","2060-2079")
  meanC=meanW
  
  for(i in 1:12)
  {
    data=eventsN[[i]]
    mm=floor(data$Date1/100)%%100
    data2=eventsF[[i]]
    mm2=floor(data2$Date1/100)%%100
    
    for(j in 1:5)
    {
      b=order(data$CVmax,decreasing=T)
      thresh2=data$CVmax[b[20*thresh]]
      if(is.na(thresh2)) thresh2=min(data$CVmax)
      I=which(data$CVmax>=thresh2 & mm>=5 & mm<=10)
      meanC[i,j,1]=mean(data[I,clist[j]])
      I=which(data$CVmax>=thresh2 & (mm>=11 | mm<=4))
      meanW[i,j,1]=mean(data[I,clist[j]])
      I=which(data2$CVmax>=thresh2 & mm2>=5 & mm2<=10)
      meanC[i,j,2]=mean(data2[I,clist[j]])
      I=which(data2$CVmax>=thresh2 & (mm2>=11 | mm2<=4))
      meanW[i,j,2]=mean(data2[I,clist[j]])  
    }
  }
  
  change=abind(100*((meanC[,,2]/meanC[,,1])-1),100*((meanW[,,2]/meanW[,,1])-1),along=3)
  dimnames(change)[[3]]<-c("Cool","Warm")
  names(dimnames(change))<-c("Source","Metric","Season")
  
  pdf(file=paste("LAP_ECLevents_thresh",thresh,"_change_byseason_meanimpacts.pdf",sep=""),width=8,height=4,pointsize=10)
  tmp=melt(change)
  print(ggplot(tmp, aes(x = Metric, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-500, 500, 25)) +
          theme_bw() + ylab("Percentage change in ECL impacts") + xlab("") +  geom_hline(yintercept = 0))
  dev.off()
  pdf(file=paste("LAP_ECLevents_thresh",thresh,"_change_byseason_meanimpacts_trunc.pdf",sep=""),width=8,height=4,pointsize=10)
  tmp=melt(change)
  print(ggplot(tmp, aes(x = Metric, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 100, 25),limits = c(-50, 100)) +
          theme_bw() + ylab("Percentage change in ECL impacts") + xlab("") +  geom_hline(yintercept = 0))
  
  dev.off()
}

thresh=22
xthresh=c(5,50,13.9,22.2)
countW=array(0,c(12,5,2))
clist=c(14:17)
dimnames(countW)[[1]]<-dimnames(future)[[1]]
dimnames(countW)[[2]]<-c("All","MeanR>=5mm","MaxR>=50mm","MeanW>=50km/h","MaxW>=80km/h")
dimnames(countW)[[3]]<-c("1990-2009","2060-2079")
countC=countW
  
  for(i in 1:12)
  {
    data=eventsN[[i]]
    b=order(data$CVmax,decreasing=T)
    thresh2=data$CVmax[b[20*thresh]]
    if(is.na(thresh2)) thresh2=min(data$CVmax)
    data=data[data$CVmax>=thresh2,]
    mm=floor(data$Date1/100)%%100
    data2=eventsF[[i]]
    data2=data2[data2$CVmax>=thresh2,]
    mm2=floor(data2$Date1/100)%%100
    
    I=which(mm>=5 & mm<=10)
    countC[i,1,1]=length(I)
    I=which(mm2>=5 & mm2<=10)
    countC[i,1,2]=length(I)
    I=which(mm>=11 | mm<=4)
    countW[i,1,1]=length(I)    
    I=which(mm2>=11 | mm2<=4)
    countW[i,1,2]=length(I)    
    
    for(j in 1:4)
    {
      I=which(data[,clist[j]]>=xthresh[j] & mm>=5 & mm<=10)
      countC[i,j+1,1]=length(I)
      I=which(data2[,clist[j]]>=xthresh[j] & mm2>=5 & mm2<=10)
      countC[i,j+1,2]=length(I)
      I=which(data[,clist[j]]>=xthresh[j] & (mm>=11 | mm<=4))
      countW[i,j+1,1]=length(I)
      I=which(data2[,clist[j]]>=xthresh[j] & (mm2>=11 | mm2<=4))
      countW[i,j+1,2]=length(I)
    }
  }
  
  change=abind(100*((countC[,,2]/countC[,,1])-1),100*((countW[,,2]/countW[,,1])-1),along=3)
  dimnames(change)[[3]]<-c("Cool","Warm")
  names(dimnames(change))<-c("Source","Metric","Season")
  
  pdf(file=paste("LAP_ECLevents_thresh",thresh,"_change_byseason_impactthresh.pdf",sep=""),width=8,height=4,pointsize=10)
  tmp=melt(change)
  print(ggplot(tmp, aes(x = Metric, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-500, 500, 25)) +
          theme_bw() + ylab("Percentage change in 90th percentile") + xlab("") +  geom_hline(yintercept = 0))
  dev.off()
  pdf(file=paste("LAP_ECLevents_thresh",thresh,"_change_byseason_impactthresh_trunc.pdf",sep=""),width=8,height=4,pointsize=10)
  tmp=melt(change)
  print(ggplot(tmp, aes(x = Metric, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 100, 25),limits = c(-50, 100)) +
          theme_bw() + ylab("Percentage change in 90th percentile") + xlab("") +  geom_hline(yintercept = 0))
  
  dev.off()





countE<-countF<-matrix(0,15,2)

for(n in 1:15)
{
  countE[n,1]=length(which(eventsN[[n]]$MaxMeanRain>=10))
  countF[n,1]=length(which(fixesN[[n]]$MeanRain>=10))
  
  if(n<=12)
  {
    countE[n,2]=length(which(eventsF[[n]]$MaxMeanRain>=10))
    countF[n,2]=length(which(fixesF[[n]]$MeanRain>=10))
  }
}

count2=cbind(countE[1:12,2]/countE[1:12,1],countF[1:12,2]/countF[1:12,1])
colnames(count2)=c("Events","Fixes")
boxplot(count2-1)

meanE<-meanF<-matrix(0,15,2)

for(n in 1:15)
{
  meanE[n,1]=mean(eventsN[[n]]$TotalRain)
  meanF[n,1]=mean(fixesN[[n]]$MeanRain)
  
  if(n<=12)
  {
    meanE[n,2]=mean(eventsF[[n]]$TotalRain)
    meanF[n,2]=mean(fixesF[[n]]$MeanRain)
  }
}

colnames(meanF)<-c("1990-2009","2060-2079")
boxplot(meanF[1:12,])




######## Version 2 - change in rain properties for events over the 22pa threshold

thresh=22
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
    filelistF=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""))
    filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""))
    filelistC=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_composite_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.nc",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_composite_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.nc",sep=""))
    
    
    ff=1
    data=read.csv(filelistE[ff])
    b=order(data[,intcolE],decreasing=T)
    thresh2=data[b[20*thresh],intcolE]
    if(is.na(thresh2)) thresh2=min(data[,intcolE])
    eventsN[[n]]=data[data[,intcolE]>=thresh2,]
    
    data=read.csv(filelistF[ff])
    a=open.nc(filelistC[ff])
    tmp=var.get.nc(a,"ECLrain")
    tmp[tmp>=600]=NaN
    data$MeanRain=apply(tmp,3,mean,na.rm=T)
    data$MaxRain=apply(tmp,3,max,na.rm=T)
    fixesN[[n]]=data[data[,intcolF]>=thresh2,]
    
    eventsN[[n]]$MaxPointRain<-eventsN[[n]]$MaxMeanRain<-eventsN[[n]]$TotalRain<-0
    for(k in 1:length(eventsN[[n]]$ID))
    {
      I=which(fixesN[[n]]$ID==eventsN[[n]]$ID[k])
      eventsN[[n]]$TotalRain[k]=sum(fixesN[[n]]$MeanRain[I],na.rm=T)
      eventsN[[n]]$MaxMeanRain[k]=max(fixesN[[n]]$MeanRain[I],na.rm=T)
      eventsN[[n]]$MaxPointRain[k]=max(fixesN[[n]]$MaxRain[I],na.rm=T)
    }
    
    if(i<5)
    {
      ff=2
      data=read.csv(filelistE[ff])
      eventsF[[n]]=data[data[,intcolE]>=thresh2,]
      data=read.csv(filelistF[ff])
      a=open.nc(filelistC[ff])
      tmp=var.get.nc(a,"ECLrain")
      tmp[tmp>=600]=NaN
      data$MeanRain=apply(tmp,3,mean,na.rm=T)
      data$MaxRain=apply(tmp,3,max,na.rm=T)
      fixesF[[n]]=data[data[,intcolF]>=thresh2,]
      
      eventsF[[n]]$MaxPointRain<-eventsF[[n]]$MaxMeanRain<-eventsF[[n]]$TotalRain<-0
      for(k in 1:length(eventsF[[n]]$ID))
      {
        I=which(fixesF[[n]]$ID==eventsF[[n]]$ID[k])
        eventsF[[n]]$TotalRain[k]=sum(fixesF[[n]]$MeanRain[I],na.rm=T)
        eventsF[[n]]$MaxMeanRain[k]=max(fixesF[[n]]$MeanRain[I],na.rm=T)
        eventsF[[n]]$MaxPointRain[k]=max(fixesF[[n]]$MaxRain[I],na.rm=T)
      }    
    }
    close.nc(a)
    
    n=n+1  
  }

plot(density(eventsN[[13]]$MaxMeanRain))

countE<-countF<-matrix(0,15,2)

for(n in 1:15)
{
  countE[n,1]=length(which(eventsN[[n]]$MaxMeanRain>=15))
  countF[n,1]=length(which(fixesN[[n]]$MeanRain>=15))
  
  if(n<=12)
  {
    countE[n,2]=length(which(eventsF[[n]]$MaxMeanRain>=15))
    countF[n,2]=length(which(fixesF[[n]]$MeanRain>=15))
  }
}

meanE<-meanF<-matrix(0,15,2)

for(n in 1:15)
{
  meanE[n,1]=mean(eventsN[[n]]$TotalRain)
  meanF[n,1]=mean(fixesN[[n]]$MeanRain)
  
  if(n<=12)
  {
    meanE[n,2]=mean(eventsF[[n]]$TotalRain)
    meanF[n,2]=mean(fixesF[[n]]$MeanRain)
  }
}

colnames(meanF)<-c("1990-2009","2060-2079")
boxplot(meanF[1:12,2]-meanF[1:12,1])


lat=seq(-500,500,50)
lon=seq(-500,500,50)
library(sp)

dist=matrix(0,21,21)
for(i in 1:21)
  for(j in 1:21)
    dist[i,j]=sqrt(lat[i]^2 + lon[j]^2)

dist2<-dist3<-matrix(NaN,21,21)
dist2[dist<=500]=1
dist3[dist<=250]=1

proj=c(100,240)
cv1=c("rad2cv06","rad2cv1")
cv2=c("rad2cv0.8","rad2cv1.35")
n=1

for(yr in c("9009"))
for(i in 5)
  for(j in 1:3)
    for(l in 1)
  {
    filelistF=paste("outputUM/proj",proj[l],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_",cv1[l],"/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[l],"_",cv2[l],"_",yr,".csv",sep="")
    filelistE=paste("outputUM/proj",proj[l],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_",cv1[l],"/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[l],"_",cv2[l],"_",yr,".csv",sep="")
    filelistC=paste("outputUM/proj",proj[l],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_",cv1[l],"/ECLfixes_compositerain_",cmip[i],"_wrf",wrf[j],"_proj",proj[l],"_",cv2[l],"_",yr,".nc",sep="")
    filelistW=paste("outputUM/proj",proj[l],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_",cv1[l],"/ECLfixes_compositewind_",cmip[i],"_wrf",wrf[j],"_proj",proj[l],"_",cv2[l],"_",yr,".nc",sep="")
                
    events=read.csv(filelistE)
    data=read.csv(filelistF)
    a=open.nc(filelistC)
    tmp=var.get.nc(a,"ECLrain")
    tmp[tmp>=500]=NaN
    close.nc(a)
    a=dim(tmp)
    for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist2
    data$MeanRain500=apply(tmp,3,mean,na.rm=T)
    data$MaxRain500=apply(tmp,3,max,na.rm=T)
    for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist3
    data$MeanRain250=apply(tmp,3,mean,na.rm=T)
    data$MaxRain250=apply(tmp,3,max,na.rm=T)
    
    a=open.nc(filelistW)
    tmp=var.get.nc(a,"ECL_WS10")
    close.nc(a)
    a=dim(tmp)
    for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist2
    data$MeanWind500=apply(tmp,3,mean,na.rm=T)
    data$MaxWind500=apply(tmp,3,max,na.rm=T)
    for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist3
    data$MeanWind250=apply(tmp,3,mean,na.rm=T)
    data$MaxWind250=apply(tmp,3,max,na.rm=T)
    
    events$MaxPointWind500<-events$MaxMeanWind500<-events$MaxPointRain500<-events$MaxMeanRain500<-events$TotalRain500<-0
    events$MaxPointWind250<-events$MaxMeanWind250<-events$MaxPointRain250<-events$MaxMeanRain250<-events$TotalRain250<-0
    for(k in 1:length(events$ID))
    {
      I=which(data$ID==events$ID[k] & data$Location==1)
      events$TotalRain500[k]=sum(data$MeanRain500[I],na.rm=T)
      events$MaxMeanRain500[k]=max(data$MeanRain500[I],na.rm=T)
      events$MaxPointRain500[k]=max(data$MaxRain500[I],na.rm=T)
      events$MaxMeanWind500[k]=max(data$MeanWind500[I],na.rm=T)
      events$MaxPointWind500[k]=max(data$MaxWind500[I],na.rm=T)
      events$TotalRain250[k]=sum(data$MeanRain250[I],na.rm=T)
      events$MaxMeanRain250[k]=max(data$MeanRain250[I],na.rm=T)
      events$MaxPointRain250[k]=max(data$MaxRain250[I],na.rm=T)
      events$MaxMeanWind250[k]=max(data$MeanWind250[I],na.rm=T)
      events$MaxPointWind250[k]=max(data$MaxWind250[I],na.rm=T)
    }
    filelistF2=paste("outputUM/proj",proj[l],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_",cv1[l],"/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[l],"_",cv2[l],"_",yr,"_impacts_v2.csv",sep="")
    filelistE2=paste("outputUM/proj",proj[l],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_",cv1[l],"/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj",proj[l],"_",cv2[l],"_",yr,"_impacts_v2.csv",sep="")
    
    write.csv(data[,-1],file=filelistF2)
    write.csv(events[,-1],file=filelistE2)
    }


res=c(50,150)
pg=c("1.2","0.7")
n=1
for(yr in c("9009","6079"))
  for(i in 1:4)
    for(j in 1:3)
      for(l in 1:2)
      {
        filelistF=paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res",res[l],"_",yr,"_pg",pg[l],"_final.csv",sep="")
        filelistE=paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res",res[l],"_",yr,"_pg",pg[l],"_final.csv",sep="")
        filelistC=paste("Alejandro/final/ECLfixes_compositerain_Alejandro_",cmip[i],"_wrf",wrf[j],"_res",res[l],"_",yr,"_pg",pg[l],".nc",sep="")
        filelistW=paste("Alejandro/final/ECLfixes_compositewind_Alejandro_",cmip[i],"_wrf",wrf[j],"_res",res[l],"_",yr,"_pg",pg[l],".nc",sep="")
        
        events=read.csv(filelistE)
        data=read.csv(filelistF)
        data=data[,1:15]
        events=events[,1:13]
        a=open.nc(filelistC)
        tmp=var.get.nc(a,"ECLrain")
        tmp[tmp>=500]=NaN
        close.nc(a)
        a=dim(tmp)
        for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist2
        data$MeanRain500=apply(tmp,3,mean,na.rm=T)
        data$MaxRain500=apply(tmp,3,max,na.rm=T)
        for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist3
        data$MeanRain250=apply(tmp,3,mean,na.rm=T)
        data$MaxRain250=apply(tmp,3,max,na.rm=T)

        a=open.nc(filelistW)
        tmp=var.get.nc(a,"ECL_WS10")
        close.nc(a)
        a=dim(tmp)
        for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist2
        data$MeanWind500=apply(tmp,3,mean,na.rm=T)
        data$MaxWind500=apply(tmp,3,max,na.rm=T)
        for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist3
        data$MeanWind250=apply(tmp,3,mean,na.rm=T)
        data$MaxWind250=apply(tmp,3,max,na.rm=T)
        
        events$MaxPointWind500<-events$MaxMeanWind500<-events$MaxPointRain500<-events$MaxMeanRain500<-events$TotalRain500<-0
        events$MaxPointWind250<-events$MaxMeanWind250<-events$MaxPointRain250<-events$MaxMeanRain250<-events$TotalRain250<-0
        for(k in 1:length(events$ID))
        {
          I=which(data$ID==events$ID[k] & data$Location==1)
          events$TotalRain500[k]=sum(data$MeanRain500[I],na.rm=T)
          events$MaxMeanRain500[k]=max(data$MeanRain500[I],na.rm=T)
          events$MaxPointRain500[k]=max(data$MaxRain500[I],na.rm=T)
          events$MaxMeanWind500[k]=max(data$MeanWind500[I],na.rm=T)
          events$MaxPointWind500[k]=max(data$MaxWind500[I],na.rm=T)
          events$TotalRain250[k]=sum(data$MeanRain250[I],na.rm=T)
          events$MaxMeanRain250[k]=max(data$MeanRain250[I],na.rm=T)
          events$MaxPointRain250[k]=max(data$MaxRain250[I],na.rm=T)
          events$MaxMeanWind250[k]=max(data$MeanWind250[I],na.rm=T)
          events$MaxPointWind250[k]=max(data$MaxWind250[I],na.rm=T)
        }
        filelistF2=paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res",res[l],"_",yr,"_pg",pg[l],"_impacts_v2.csv",sep="")
        filelistE2=paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res",res[l],"_",yr,"_pg",pg[l],"_impacts_v2.csv",sep="")
        
        write.csv(data[,-1],file=filelistF2)
        write.csv(events[,-1],file=filelistE2)
      }



