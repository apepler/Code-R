######## CORDEX figure for change in ECL frequency
### Do using the 22 p.a. threshold & just cv0.8 threshold

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)

fixesN<-fixesF<-list()
Ethresh=matrix(0,5,3)

lat=seq(-50,0,2.5)
lon=seq(100,180,2.5)
locN22<-locF22<-locF<-locN<-array(0,c(21,33,5,3,12))

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")

#### Version 1 - all above threshold
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    tmp=read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""))
    
    b=order(tmp$CVmax,decreasing=T)
    Ethresh[i,j]=tmp$CVmax[b[20*22]]
    
    filelistF=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_rad2cv06_9009_ALL.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_rad2cv06_6079_ALL.csv",sep=""))
    
    ff=1
    fixesN[[n]]=read.csv(filelistF[ff])
    fixesN[[n]]$Year=floor(fixesN[[n]]$Date/10000)
    fixesN[[n]]$Month=floor(fixesN[[n]]$Date/100)%%100
    
    ##### Remove single-fix events
    x<-rle(fixesN[[n]]$ID)
    events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=8))
    I<-which(events[,2]==1)
    if(length(I)>0)
    {
      y<-match(events[I,1],fixesN[[n]]$ID)
      events=events[-I,]
      fixesN[[n]]=fixesN[[n]][-y,]
    }
    
    for(y in 1:length(lat))
      for(x in 1:length(lon))
        for(m in 1:12)
      {
        I=which(fixesN[[n]]$Lon>=lon[x]-1.25 & fixesN[[n]]$Lon<lon[x]+1.25 & fixesN[[n]]$Lat>=lat[y]-1.25 & fixesN[[n]]$Lat<lat[y]+1.25 & fixesN[[n]]$Month==m)
        locN[y,x,i,j,m]=length(I)/20
        I=which(fixesN[[n]]$Lon>=lon[x]-1.25 & fixesN[[n]]$Lon<lon[x]+1.25 & fixesN[[n]]$Lat>=lat[y]-1.25 & fixesN[[n]]$Lat<lat[y]+1.25 & fixesN[[n]]$Month==m & fixesN[[n]]$CV>=Ethresh[i,j])
        locN22[y,x,i,j,m]=length(I)/20
      }
    
    
    if(i<5)
    {
      ff=2
      fixesF[[n]]=read.csv(filelistF[ff])
      fixesF[[n]]$Year=floor(fixesF[[n]]$Date/10000)
      fixesF[[n]]$Month=floor(fixesF[[n]]$Date/100)%%100
      
      ##### Remove single-fix events
      x<-rle(fixesF[[n]]$ID)
      events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=8))
      I<-which(events[,2]==1)
      if(length(I)>0)
      {
        y<-match(events[I,1],fixesF[[n]]$ID)
        events=events[-I,]
        fixesF[[n]]=fixesF[[n]][-y,]
      }
      
      for(y in 1:length(lat))
        for(x in 1:length(lon))
          for(m in 1:12)
          {
            I=which(fixesF[[n]]$Lon>=lon[x]-1.25 & fixesF[[n]]$Lon<lon[x]+1.25 & fixesF[[n]]$Lat>=lat[y]-1.25 & fixesF[[n]]$Lat<lat[y]+1.25 & fixesF[[n]]$Month==m)
            locF[y,x,i,j,m]=length(I)/20
            I=which(fixesF[[n]]$Lon>=lon[x]-1.25 & fixesF[[n]]$Lon<lon[x]+1.25 & fixesF[[n]]$Lat>=lat[y]-1.25 & fixesF[[n]]$Lat<lat[y]+1.25 & fixesF[[n]]$Month==m & fixesF[[n]]$CV>=Ethresh[i,j])
            locF22[y,x,i,j,m]=length(I)/20
          }
    }
    
    n=n+1  
  }

changeC=apply(locF22[,,1:4,1:3,5:10],c(1,2,3,4),sum)/apply(locN22[,,1:4,1:3,5:10],c(1,2,3,4),sum)
changeW=apply(locF22[,,1:4,1:3,c(1:4,11:12)],c(1,2,3,4),sum)/apply(locN22[,,1:4,1:3,c(1:4,11:12)],c(1,2,3,4),sum)

image(apply(locN,c(1,2),mean))

library(maps)
ColorBar <- function(brks,cols,labels=NA)
{
  par(mar = c(3, 1, 3, 4), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  if(is.na(labels[1])) labels=brks[seq(2, length(brks)-1)]
  
  axis(4, at = seq(1.5, length(brks) - 1.5), tick = TRUE, 
       labels = labels)
}

bb2=c(-10000,c(-50,-25,-10,-5,0,5,10,25,50),10000)
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
cm=pal(10)

pdf(file=paste("ECL_locations_CORDEX_thresh22_change_season_points2.pdf",sep=""),width=8,height=9,pointsize=12)
layout(cbind(c(1,2),c(3,3)),width=c(1,0.2))
par(mar=c(3,3,3,2))
image(lon,lat,t(100*(apply(changeC,c(1,2),mean)-1)),xlab="",ylab="",breaks=bb2,
        col=cm,zlim=c(-Inf,Inf),cex.axis=1.5,cex.main=1.5,main="May-October",
        xlim=c(110,172.5),ylim=c(-45,-10))
tmp=apply(changeC>1,c(1,2),mean) ## Proportion with a positive change
sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=1,lwd=2) ## Add some points
map(add=T,lwd=2)

image(lon,lat,t(100*(apply(changeW,c(1,2),mean)-1)),xlab="",ylab="",breaks=bb2,
      col=cm,zlim=c(-Inf,Inf),cex.axis=1.5,cex.main=1.5,main="November-April",
      xlim=c(110,172.5),ylim=c(-45,-10))
tmp=apply(changeW>1,c(1,2),mean) ## Proportion with a positive change
sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=1,lwd=2) ## Add some points
map(add=T,lwd=2)

ColorBar(bb2,cm,labels=c("-50%","-25%","-10%","-5%","","+5%","+10%","+25%","+50%"))
dev.off()


