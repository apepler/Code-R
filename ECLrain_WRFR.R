##Want to find the proportion of WRF rain that can be attributed to ECLs
##For each version of WRF, need only the ECL days above certain CV thresholds, then attribute subsequent day's rain
##Subsequent day, of course, being the rain at 00Z n+1 - 00Z n for RAINC+RAINNC

###So, I need a nice csv file with day1, day2, and csv for each dataset - I'll prep this in excel

rm(list=ls())
library(RNetCDF)
setwd('~/Documents/ECLs')
load("ECLprops_ann.RData") ##Has WRF lat/lon
library("R.matlab")
library(akima)
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

name=c("R1","R2","R3")
dir=c('/srv/ccrc/data31/z3393020/NARCliM/reanalysis/R1/out/','/srv/ccrc/data27/z3393020/WRF_NNRP/R2/1950-2010/out/','/srv/ccrc/data29/z3236814/NARCliM/reanalysis/R3/out/')

for(n in 1:3)
{
  ###Load monthly data
  years=seq(1980,2009)
  rainA<-rainC<-array(0,dim=c(215,144,30))
  rainW<-array(0,dim=c(215,144,29))
  for(i in 1:length(years))
  {
    ##Get the corresponding rain - in a monthly file? 
    if(i<30)
    {
      f1=open.nc(paste(dir[n],'wrfhrly_d01_',years[i],'-01-01_00:00:00',sep=""))
      f2=open.nc(paste(dir[n],'wrfhrly_d01_',years[i]+1,'-01-01_00:00:00',sep=""))  
      rainA[,,i]=var.get.nc(f2,"RAINC",c(1,1,1),c(215,144,1))+
        var.get.nc(f2,"RAINNC",c(1,1,1),c(215,144,1))-
        var.get.nc(f1,"RAINC",c(1,1,1),c(215,144,1))-
        var.get.nc(f1,"RAINNC",c(1,1,1),c(215,144,1))
      close.nc(f1)
      close.nc(f2)
    } else {
      f1=open.nc(paste(dir[n],'wrfhrly_d01_',years[i],'-01-01_00:00:00',sep=""))
      times=var.get.nc(f1,"Times")
      rainA[,,i]=var.get.nc(f1,"RAINC",c(1,1,dim(times)),c(215,144,1))+
        var.get.nc(f1,"RAINNC",c(1,1,dim(times)),c(215,144,1))-
        var.get.nc(f1,"RAINC",c(1,1,1),c(215,144,1))-
        var.get.nc(f1,"RAINNC",c(1,1,1),c(215,144,1))
      close.nc(f1)
    }

    
    f1=open.nc(paste(dir[n],'wrfhrly_d01_',years[i],'-05-01_00:00:00',sep=""))
    f2=open.nc(paste(dir[n],'wrfhrly_d01_',years[i],'-11-01_00:00:00',sep=""))  
    rainC[,,i]=var.get.nc(f2,"RAINC",c(1,1,1),c(215,144,1))+
      var.get.nc(f2,"RAINNC",c(1,1,1),c(215,144,1))-
      var.get.nc(f1,"RAINC",c(1,1,1),c(215,144,1))-
      var.get.nc(f1,"RAINNC",c(1,1,1),c(215,144,1))
    close.nc(f1)
    close.nc(f2)
    
    if(i<30)
    {
      f1=open.nc(paste(dir[n],'wrfhrly_d01_',years[i],'-11-01_00:00:00',sep=""))
      f2=open.nc(paste(dir[n],'wrfhrly_d01_',years[i]+1,'-05-01_00:00:00',sep=""))  
      rainW[,,i]=var.get.nc(f2,"RAINC",c(1,1,1),c(215,144,1))+
        var.get.nc(f2,"RAINNC",c(1,1,1),c(215,144,1))-
        var.get.nc(f1,"RAINC",c(1,1,1),c(215,144,1))-
        var.get.nc(f1,"RAINNC",c(1,1,1),c(215,144,1))
      close.nc(f1)
      close.nc(f2)
    }    
  }
  
  ##Load ECL list
  Efile=paste('CSV/ECLdays_umelb_wrf',name[n],'_default.csv',sep="")
  ECLs<-read.csv(Efile,header=T,sep=",")
  
  ##V2 - all days with local csv >= 0.25
  ECLs2=ECLs[ECLs[,3]>=0.25,]
  
  date1=cbind(floor(ECLs2[,1]/10000),floor((ECLs2[,1] %% 10000)/100), ECLs2[,1]%%100)
  date2=cbind(floor(ECLs2[,2]/10000),floor((ECLs2[,2] %% 10000)/100), ECLs2[,2]%%100)
  
  ECLrain<-ECLcool<-array(0,dim=c(215,144,30))
  ECLwarm<-array(0,dim=c(215,144,29))
  
  for(i in 1:length(date1[,1]))
  {
    ##Get the corresponding rain - in a monthly file? 
    f1=open.nc(paste(dir[n],'wrfhrly_d01_',date1[i,1],'-',sprintf("%2.2i",date1[i,2]),'-01_00:00:00',sep=""))
    f2=open.nc(paste(dir[n],'wrfhrly_d01_',date2[i,1],'-',sprintf("%2.2i",date2[i,2]),'-01_00:00:00',sep=""))
    
    rain=var.get.nc(f2,"RAINC",c(1,1,date2[i,3]*24-23),c(215,144,1))+
      var.get.nc(f2,"RAINNC",c(1,1,date2[i,3]*24-23),c(215,144,1))-
      var.get.nc(f1,"RAINC",c(1,1,date1[i,3]*24-23),c(215,144,1))-
      var.get.nc(f1,"RAINNC",c(1,1,date1[i,3]*24-23),c(215,144,1))
    
    rain[rain<0]=0 ##Fix for weird oddities
    
    ##Add it to year, index is yy-1979
    ECLrain[,,date2[i,1]-1979]=ECLrain[,,date2[i,1]-1979]+rain
    
    ##Add it to warm/cool depending on month
    if(date2[i,2]>=5 & date2[i,2]<=10) { 
      ECLcool[,,date2[i,1]-1979]=ECLcool[,,date2[i,1]-1979]+rain 
    } else if(date2[i,2]>=11 & date2[i,1]<2009) {
      ECLwarm[,,date2[i,1]-1979]=ECLwarm[,,date2[i,1]-1979]+rain 
    } else if(date2[i,2]<=4 & date2[i,1]>1980) {
      ECLwarm[,,date2[i,1]-1980]=ECLwarm[,,date2[i,1]-1980]+rain }  ##i.e. add to previous year, no initial year
    
    close.nc(f1)
    close.nc(f2)
  }
  
  propA=apply(ECLrain,c(1,2),sum)/apply(rainA,c(1,2),sum)
  propW=apply(ECLwarm,c(1,2),sum)/apply(rainW,c(1,2),sum)
  propC=apply(ECLcool,c(1,2),sum)/apply(rainC,c(1,2),sum)
  
  latt=as.vector(lat)
  lont=as.vector(lon)
  a=as.vector(propA)
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  propA2=b$z
  a=as.vector(propC)
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  propC2=b$z
  a=as.vector(propW)
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  propW2=b$z
  
  ##Best ratio = 450 by 500
  tiff(filename=paste('ECLprop_',name[n],'_ann_cv25.tiff',sep=""),width=450,height=500)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,propA2*t(Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,propA2*t(Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "",    key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%","60%","70%"),at=seq(0,0.7,0.1))) 
  dev.off()
  
  ##Best ratio = 800 by 500
  tiff(filename=paste('ECLprop_',name[n],'_coolwarm_cv25.tiff',sep=""),width=800,height=500)
  par(plt = c(0.06,0.47,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,propC2*t(Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  text(115,-40,"a)",cex=2)
  par(new=TRUE, plt = c(0.47,0.88,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,propW2*t(Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),axes=F,frame.plot=F)
  axis(side = 1, at = seq(146,154,2))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  text(115,-40,"a)",cex=2)
  par(new = "TRUE",plt = c(0.88,0.93,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,propC2*t(Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "", key.axes = axis(4,labels=c("0%","20%","40%","60%","80%","100%"),at=seq(0,1,0.2)))
  dev.off()
  
  ##Redo for top37
  Efile=paste('CSV/ECLdays_umelb_wrf',name[n],'_default.csv',sep="")
  ECLs<-read.csv(Efile,header=T,sep=",")
  a=order(ECLs[,3],decreasing=TRUE)
  ECLs2=ECLs[sort(a[1:1110]),]
  
  date1=cbind(floor(ECLs2[,1]/10000),floor((ECLs2[,1] %% 10000)/100), ECLs2[,1]%%100)
  date2=cbind(floor(ECLs2[,2]/10000),floor((ECLs2[,2] %% 10000)/100), ECLs2[,2]%%100)
  
  ECLrain<-ECLcool<-array(0,dim=c(215,144,30))
  ECLwarm<-array(0,dim=c(215,144,29))
  
  for(i in 1:length(date1[,1]))
  {
    ##Get the corresponding rain - in a monthly file? 
    f1=open.nc(paste(dir[n],'wrfhrly_d01_',date1[i,1],'-',sprintf("%2.2i",date1[i,2]),'-01_00:00:00',sep=""))
    f2=open.nc(paste(dir[n],'wrfhrly_d01_',date2[i,1],'-',sprintf("%2.2i",date2[i,2]),'-01_00:00:00',sep=""))
    
    rain=var.get.nc(f2,"RAINC",c(1,1,date2[i,3]*24-23),c(215,144,1))+
      var.get.nc(f2,"RAINNC",c(1,1,date2[i,3]*24-23),c(215,144,1))-
      var.get.nc(f1,"RAINC",c(1,1,date1[i,3]*24-23),c(215,144,1))-
      var.get.nc(f1,"RAINNC",c(1,1,date1[i,3]*24-23),c(215,144,1))
    
    rain[rain<0]=0 ##Fix for weird oddities
    
    ##Add it to year, index is yy-1979
    ECLrain[,,date2[i,1]-1979]=ECLrain[,,date2[i,1]-1979]+rain
    
    ##Add it to warm/cool depending on month
    if(date2[i,2]>=5 & date2[i,2]<=10) { 
      ECLcool[,,date2[i,1]-1979]=ECLcool[,,date2[i,1]-1979]+rain 
    } else if(date2[i,2]>=11 & date2[i,1]<2009) {
      ECLwarm[,,date2[i,1]-1979]=ECLwarm[,,date2[i,1]-1979]+rain 
    } else if(date2[i,2]<=4 & date2[i,1]>1980) {
      ECLwarm[,,date2[i,1]-1980]=ECLwarm[,,date2[i,1]-1980]+rain }  ##i.e. add to previous year, no initial year
    
    close.nc(f1)
    close.nc(f2)
  }
  
  propA=apply(ECLrain,c(1,2),sum)/apply(rainA,c(1,2),sum)
  propW=apply(ECLwarm,c(1,2),sum)/apply(rainW,c(1,2),sum)
  propC=apply(ECLcool,c(1,2),sum)/apply(rainC,c(1,2),sum)
  
  latt=as.vector(lat)
  lont=as.vector(lon)
  a=as.vector(propA)
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  propA2=b$z
  a=as.vector(propC)
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  propC2=b$z
  a=as.vector(propW)
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  propW2=b$z
  
  ##Best ratio = 250 by 500
  tiff(filename=paste('ECLprop_',name[n],'_ann_37.tiff',sep=""),width=450,height=500)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,propA2*t(Useful$mask),lev=seq(0,0.5,0.05),col=rich.colors(10),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,propA2*t(Useful$mask),lev=seq(0,0.5,0.05),col=rich.colors(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "",    key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%","60%","70%"),at=seq(0,0.7,0.1))) 
  dev.off()
  
  ##Best ratio = 800 by 500
  tiff(filename=paste('ECLprop_',name[n],'_coolwarm_37.tiff',sep=""),width=800,height=500)
  par(plt = c(0.06,0.47,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,propC2*t(Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  text(115,-40,"a)",cex=2)
  par(new=TRUE, plt = c(0.47,0.88,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,propW2*t(Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),axes=F,frame.plot=F)
  axis(side = 1, at = seq(146,154,2))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  text(115,-40,"a)",cex=2)
  par(new = "TRUE",plt = c(0.88,0.93,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,propC2*t(Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "", key.axes = axis(4,labels=c("0%","20%","40%","60%","80%","100%"),at=seq(0,1,0.2)))
  dev.off()
  
}





