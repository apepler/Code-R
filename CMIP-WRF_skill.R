rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)
cmip=c("ncep","cccma","csiromk3","echam5","miroc")
wrf=c("R1","R2","R3")

count<-thresh2<-array(NaN,c(5,4))

fixes<-events<-fixesRCM<-eventsRCM<-list()

#### Four methods
#### 22 ECLs p.a, by season
#### Count, Mean(meanR), mean(maxR), MeanR>=5, MaxR>=50


thresh=22

n=1
for(i in 1:5)
{
  tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_rad2cv06_v2/ECLevents_umelb_",cmip[i],"_proj100_rad2cv06_9009.csv",sep=""))
    count[i,1]=length(which(tmp$CV2>=1))
    if(!("CV2" %in% colnames(tmp))) tmp$CV2=tmp$CVmax
    
    b=order(tmp$CV2,decreasing=T)
    thresh2[i,1]=tmp$CV2[b[20*thresh]]
    if(is.na(thresh2[i,1])) thresh2[i,1]=min(tmp$CV2,na.rm=T)
    events[[i]]=tmp[tmp$CV2>=thresh2[i,1],]
    events[[i]]$Year=floor(events[[i]]$Date1/10000)
    events[[i]]$Month=floor(events[[i]]$Date1/100)%%100
    
    tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_rad2cv06_v2/ECLfixes_umelb_",cmip[i],"_proj100_rad2cv06_9009.csv",sep=""))

    fixes[[i]]=tmp[tmp$ID%in%events[[i]]$ID & tmp$CV>=thresh2[i,1],]
    fixes[[i]]$Year=floor(fixes[[i]]$Date/10000)
    fixes[[i]]$Month=floor(fixes[[i]]$Date/100)%%100
    fixes[[i]]$Date2=as.POSIXct(paste(as.character(fixes[[i]]$Date),substr(fixes[[i]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  
  
  for(j in 1:3)
  {
    tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""))
    if(!("CV2" %in% colnames(tmp))) tmp$CV2=tmp$CVmax
        count[i,j+1]=length(which(tmp$CV2>=1))
    b=order(tmp$CV2,decreasing=T)
    thresh2[i,j+1]=tmp$CV2[b[20*thresh]]
    if(is.na(thresh2[i,j+1])) thresh2[i,j+1]=min(tmp$CV2,na.rm=T)
    eventsRCM[[n]]=tmp[tmp$CV2>=thresh2[i,j+1],]
    eventsRCM[[n]]$Year=floor(eventsRCM[[n]]$Date1/10000)
    eventsRCM[[n]]$Month=floor(eventsRCM[[n]]$Date1/100)%%100    
    
    tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""))

    fixesRCM[[n]]=tmp[tmp$ID%in%eventsRCM[[n]]$ID & tmp$CV>=thresh2[i,j+1],]
    fixesRCM[[n]]$Year=floor(fixesRCM[[n]]$Date/10000)
    fixesRCM[[n]]$Month=floor(fixesRCM[[n]]$Date/100)%%100
    fixesRCM[[n]]$Date2=as.POSIXct(paste(as.character(fixesRCM[[n]]$Date),substr(fixesRCM[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    n=n+1  
  }
}


monthly<-array(NaN,c(12,5,4))

n=1
for(i in 1:5)
{
  for(m in 1:12) monthly[m,i,1]=length(which(events[[i]]$Month==m))
  
  for(j in 1:3)
  {
    for(m in 1:12) monthly[m,i,j+1]=length(which(eventsRCM[[n]]$Month==m))
    n=n+1
  }
}

monthly2=cbind(monthly[,1,1],apply(monthly[,2:5,1],1,mean,na.rm=T),apply(monthly[,1,2:4],1,mean),apply(monthly[,2:5,2:4],1,mean))/20

colnames(monthly2)=c("NCEP","CMIP","NCEP-WRF","CMIP-WRF")
clist=c("black","grey","black","grey")
llist=c(1,1,2,2)

plot(1:12,monthly2[,1],type="l",lwd=3,col="black",ylim=c(0,max(monthly2)),
     xlab="Month",ylab="Average ECL frequency")
for(i in 2:4) lines(1:12,monthly2[,i],col=clist[i],lty=llist[i],lwd=3)
legend("bottomright",legend=colnames(monthly2),lty=llist,lwd=3,col=clist,bty="n")
  

############ Locations

fixesALL<-fixesALLR<-list()
lat=seq(-50,0,2.5)
lon=seq(100,180,2.5)

loc=array(0,c(21,33,5,4))

n=1
for(i in 1:5)
{
  tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_rad2cv06_v2/ECLfixes_umelb_",cmip[i],"_proj100_rad2cv06_9009_ALL.csv",sep=""))
  fixesALL[[i]]=tmp[tmp$CV>=thresh2[i,1],]

  for(y in 1:length(lat))
    for(x in 1:length(lon))
      {
        I=which(fixesALL[[i]]$Lon>=lon[x]-1.25 & fixesALL[[i]]$Lon<lon[x]+1.25 & fixesALL[[i]]$Lat>=lat[y]-1.25 & fixesALL[[i]]$Lat<lat[y]+1.25)
        loc[y,x,i,1]=length(I)/20
      }
  
  for(j in 1:3)
  {
    tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_rad2cv06_9009_ALL.csv",sep=""))
    fixesALLR[[n]]=tmp[tmp$CV>=thresh2[i,j+1],]
    for(y in 1:length(lat))
      for(x in 1:length(lon))
      {
        I=which(fixesALLR[[n]]$Lon>=lon[x]-1.25 & fixesALLR[[n]]$Lon<lon[x]+1.25 & fixesALLR[[n]]$Lat>=lat[y]-1.25 & fixesALLR[[n]]$Lat<lat[y]+1.25)
        loc[y,x,i,1+j]=length(I)/20
      }

    n=n+1  
  }
}


library(abind)
loc2=abind(loc[,,1,1],apply(loc[,,2:5,1],c(1,2),mean,na.rm=T),apply(loc[,,1,2:4],c(1,2),mean),apply(loc[,,2:5,2:4],c(1,2),mean),along=3)
names=c("NCEP","CMIP","NCEP-WRF","CMIP-WRF")

library(maps)
ColorBar <- function(brks,cols,labels=NA)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  if(is.na(labels[1])) labels=brks[seq(2, length(brks)-1)]
  
  axis(4, at = seq(1.5, length(brks) - 1.5), tick = TRUE, 
       labels = labels)
}
cm2=gray(seq(1,0.1,-0.15))
bb2=c(-0.5,0,1,2,3,4,5,100)

for(i in 1:4)
{
pdf(file=paste("outputUM/ECL_locations_CORDEX_",names[i],"_proj100_d01_v3.pdf",sep=""),width=8,height=5,pointsize=12)
layout(cbind(1,2),c(1,0.2))
par(mar=c(3,3,3,2))
image(lon,lat,t(loc2[,,i]),xlab="",ylab="",breaks=bb2,
      col=cm2,zlim=c(-Inf,Inf),main=names[i],cex.axis=1.5,cex.main=1.5,
      xlim=c(110,172.5),ylim=c(-45,-10))
map(add=T,lwd=2)
ColorBar(bb2,cm2)
dev.off()
}

### Panels
pdf(file=paste("outputUM/ECL_locations_CORDEX_NCEPvWRF_panel_proj100_d01_v3.pdf",sep=""),width=17,height=5,pointsize=14)
layout(cbind(1,2,3),c(1,1,0.2))
par(mar=c(3,3,3,2))
for(i in c(1,3))
{
image(lon,lat,t(loc2[,,i]),xlab="",ylab="",breaks=bb2,
      col=cm2,zlim=c(-Inf,Inf),main=names[i],cex.axis=1.5,cex.main=2,
      xlim=c(110,172.5),ylim=c(-45,-10))
map(add=T,lwd=2)
}
ColorBar(bb2,cm2)
dev.off()

pdf(file=paste("outputUM/ECL_locations_CORDEX_CMIPvWRF_panel_proj100_d01_v3.pdf",sep=""),width=17,height=5,pointsize=14)
layout(cbind(1,2,3),c(1,1,0.2))
par(mar=c(3,3,3,2))
for(i in c(2,4))
{
  image(lon,lat,t(loc2[,,i]),xlab="",ylab="",breaks=bb2,
        col=cm2,zlim=c(-Inf,Inf),main=names[i],cex.axis=1.5,cex.main=2,
        xlim=c(110,172.5),ylim=c(-45,-10))
  map(add=T,lwd=2)
}
ColorBar(bb2,cm2)
dev.off()


############
## Need the one w/ the impacts

thresh=22

eventsI<-fixesI<-list()

n=1
for(i in 1:5)
{  for(j in 1:3)
  {
    tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""))
    if(!("CV2" %in% colnames(tmp))) tmp$CV2=tmp$CVmax
    count[i,j+1]=length(which(tmp$CV2>=1))
    b=order(tmp$CV2,decreasing=T)
    thresh2[i,j+1]=tmp$CV2[b[20*thresh]]
    if(is.na(thresh2[i,j+1])) thresh2[i,j+1]=min(tmp$CV2,na.rm=T)
    eventsI[[n]]=tmp[tmp$CV2>=thresh2[i,j+1],]
    eventsI[[n]]$Year=floor(eventsI[[n]]$Date1/10000)
    eventsI[[n]]$Month=floor(eventsI[[n]]$Date1/100)%%100    
    
    tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""))
    
    fixesI[[n]]=tmp[tmp$ID%in%eventsI[[n]]$ID & tmp$CV>=thresh2[i,j+1],]
    fixesI[[n]]$Year=floor(fixesI[[n]]$Date/10000)
    fixesI[[n]]$Month=floor(fixesI[[n]]$Date/100)%%100
    fixesI[[n]]$Date2=as.POSIXct(paste(as.character(fixesI[[n]]$Date),substr(fixesI[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    n=n+1  
  }
}

impcol=cbind(c(14:17),c(19:22)) # MaxMeanR, MaxMaxR, MaxMeanW, MaxMaxW
xthresh=c(6,50,13.9,22.2,6,50,13.9,22.2)

stats<-array(0,c(5,3,10))
dimnames(stats)[[1]]=cmip
dimnames(stats)[[2]]=wrf
dimnames(stats)[[3]]=c("Mean TotalRain500","Mean MaxMeanRain500",
                       "MaxMeanR500>6mm","MaxMaxR500>50mm","MaxMeanW500>50km/h","MaxMaxW500>80 km/h",
                       "MaxMeanR250>6mm","MaxMaxR250>50mm","MaxMeanW250>50km/h","MaxMaxW250>80 km/h")

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    stats[i,j,1]=mean(eventsI[[n]]$TotalRain500)
    stats[i,j,2]=mean(eventsI[[n]]$MaxMeanRain500)
    
    for(k in 1:8) stats[i,j,2+k]=length(which(eventsI[[n]][,impcol[k]]>=xthresh[k]))
    n=n+1
  }

apply(stats,c(1,3),mean)

##### What about average rain pattern?

composite<-array(0,c(21,21,15))
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)
cmip=c("ncep","cccma","csiromk3","echam5","miroc")
wrf=c("R1","R2","R3")

thresh=rbind(c(1.32,1.07,1.39),c(1.01,0.92,1.08),c(1.07,0.86,1.07),c(1.08,0.94,1.12),c(1.05,0.91,0.94))

n=1
for(i in 1:5)
{  for(j in 1:3)
{
  filelistC=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_compositewind_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.nc",sep="")
  
  a=open.nc(filelistC)
  tmp=var.get.nc(a,"ECL_WS10")
  tmp[tmp>=500]=NaN
  
  tmp2=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""))
  I=which(tmp2$CV>=thresh[i,j] & tmp2$Location==1)
  composite[,,n]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
  close.nc(a)
  
  n=n+1  
}
}

source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
library(fields)
ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = brks[seq(1, length(brks), subsampleg)])
}

names2=c("NCEP-WRF","CCCMA-WRF","CSIROMK3-WRF","ECHAM5-WRF","MIROC-WRF")

pdf(file="outputUM/ECL_compositewind_mean_vcmip.pdf",width=6.5,height=9)
layout(cbind(c(2,3,1),c(4,5,7),c(6,6,6)),width=c(1,1,0.4))
for(i in 1:5)
{
j=(i*3)-2
  par(mar=c(1,1,3,1))
image(apply(composite[,,j:(j+2)],c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main=names2[i],axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=3,cex=3)
}
ColorBar(c(seq(0,14,1),100),tim.colors(15),subsampleg=2)
dev.off()

pdf(file="outputUM/ECL_compositerain_mean_vwrf_nceponly.pdf",width=12,height=4)
layout(cbind(1,2,3,4),width=c(1,1,1,0.4))
for(i in 1:3)
{
  par(mar=c(1,1,3,1))
  image(composite[,,i],breaks=c(seq(0,14,1),100),col=tim.colors(15),main=paste("R",i,sep=""),axes=F,cex.main=2)
  points(0.5,0.5,pch=4,col="black",lwd=3,cex=3)
}
ColorBar(c(seq(0,14,1),100),tim.colors(15))
dev.off()

stat=matrix(0,5,3)
n=1
for(i in 1:5)
  for(j in 1:3)
  {
  stat[i,j]=mean(composite[,,n])
  n=n+1  
}


source('~/Documents/R/color.palette.R')
pal <- color.palette(c("white","skyblue1","blue","purple4"),c(10,20,10))
pdf(file="outputUM/ECL_compositewind_mean_v3.pdf",width=8.5,height=4)
layout(cbind(1,2,3),width=c(1,1,0.3))
par(mar=c(1,1,3,1))
image(apply(composite[,,1:3],c(1,2),mean),breaks=c(seq(6,12,0.5),100),col=pal(13),main="NCEP-WRF",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=3,cex=3)
image(apply(composite[,,4:15],c(1,2),mean),breaks=c(seq(6,12,0.5),100),col=pal(13),main="CMIP-WRF",axes=F,cex.main=2)
points(0.5,0.5,pch=4,col="black",lwd=3,cex=3)
ColorBar(c(seq(6,12,0.5),100),pal(13),subsampleg=2)
dev.off()



stat=matrix(0,5,3)
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    stat[i,j]=mean(composite[,,n])
    n=n+1  
  }
