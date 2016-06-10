rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

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
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))

#lat=seq(-42.5,-22.5,5)
#lon=seq(147.5,157.5,5)

lat=seq(-40,-25,2.5)
lon=seq(145,160,2.5)

locNC<-locFC<-locFW<-locNW<-array(0,dim=c(7,7,12,4)) #Lon, Lat, Model, Method
intcol=c(10,10,9,9)

thresh=c(10,25,50,100,150)
cm=pal(12)
bb=c(-10000,-100,-50,-25,-10,-5,0,5,10,25,50,100,10000)
blab=c("-100%","-50%","-25%","-10%","-5%","","+5%","+10%","+25%","+50%","+100%")

cm2=gray(seq(1,0.1,-0.15))
#bb2=cbind(c(-0.5,0,0.2,0.5,1,1.5,2,100),c(-0.5,0,0.5,1,1.5,2,3,100),c(-0.5,0,1,2,3,4,5,100),c(-0.5,0,1,2.5,5,7.5,10,100),c(-0.5,0,1,2.5,5,7.5,10,100))
bb2=cbind(c(-0.5,0,0.1,0.2,0.3,0.4,0.5,100),c(-0.5,0,0.2,0.4,0.6,0.8,1,100),c(-0.5,0,0.2,0.5,1,1.5,2,100),c(-0.5,0,0.5,1,1.5,2,3,100),c(-0.5,0,1,2,3,4,5,100))

cm3=pal(10)
bb3=c(-10000,0.1,0.25,0.35,0.45,0.5,0.55,0.65,0.75,0.9,10000)
blab3=c("10%","25%","35%","45%","","55%","65%","75%","90%")

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
names=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
for(t in 1:4)
{
n=1

for(i in 1:4)
  for(j in 1:3)
  {    
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelist1[ff])
      dataa=data[data$Location==1,]
      b=order(dataa[,intcol[ff]],decreasing=T)
      thresh2=dataa[b[20*thresh[t]],intcol[ff]] 
      if(is.na(thresh2)) thresh2=min(dataa[,intcol[ff]])
      dataN=data[data[,intcol[ff]]>=thresh2,]      
      data=read.csv(filelist2[ff])
      dataF=data[data[,intcol[ff]]>=thresh2,]
      
      for(xx in 1:7)
        for(yy in 1:7)
        {     
          I=which(dataN$Lat>=(lat[yy]-1.25) & dataN$Lat<=(lat[yy]+1.25) & dataN$Lon<=(lon[xx]+1.25) & dataN$Lon>=(lon[xx]-1.25))
          mm=floor((dataN[I,]$Date%%10000)/100)
          I=which(mm>=5 & mm<=10)
          locNC[yy,xx,n,ff]=length(I)
          locNW[yy,xx,n,ff]=length(mm)-length(I)     
          
          I=which(dataF$Lat>=lat[yy]-1.25 & dataF$Lat<=lat[yy]+1.25 & dataF$Lon<=lon[xx]+1.25 & dataF$Lon>=lon[xx]-1.25)
          mm=floor((dataF[I,]$Date%%10000)/100)
          I=which(mm>=5 & mm<=10)
          locFC[yy,xx,n,ff]=length(I)
          locFW[yy,xx,n,ff]=length(mm)-length(I)   
        }    
    }      
    n=n+1  
  }

changeC=((locFC/locNC)-1)*100
changeC[locNC==0]=NaN
changeW=((locFW/locNW)-1)*100
changeW[locNW==0]=NaN
changeCa=apply(changeC,c(1,2),median,na.rm=T)
changeWa=apply(changeW,c(1,2),median,na.rm=T)


# pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_change_location_v2_median.pdf",sep=""),width=9,height=4)
# layout(cbind(1,2,3),c(1,1,0.35))
# image(lon,lat,t(changeCa),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# image(lon,lat,t(changeWa),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# ColorBar(bb,cm,blab)
# dev.off()
# 
# pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_location_v2_median.pdf",sep=""),width=9,height=8)
# layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
# image(lon,lat,t(apply(locNC,c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October 1990-09",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# image(lon,lat,t(apply(locFC,c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October 2060-79",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# image(lon,lat,t(apply(locNW,c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April 1990-09",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# image(lon,lat,t(apply(locFW,c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April 2060-79",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# ColorBar(bb2[,t],cm2)
# dev.off()
# print(range(t(apply(locNC,c(1,2),median,na.rm=T))/20,t(apply(locFW,c(1,2),median,na.rm=T))/20))
# 
# 
# pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_change_location_v2_cool_vsmethod.pdf",sep=""),width=9,height=8)
# layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
# for(ff in 1:4)
# {
#   changeCa=apply(changeC[,,,ff],c(1,2),median,na.rm=T)
#   image(lon,lat,t(changeCa),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main=paste("May-October",names[ff]),cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# }
# ColorBar(bb,cm,blab)
# dev.off()
# 
# pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_change_location_v2_warm_vsmethod.pdf",sep=""),width=9,height=8)
# layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
# for(ff in 1:4)
# {
#   changeCa=apply(changeW[,,,ff],c(1,2),median,na.rm=T)
#   image(lon,lat,t(changeCa),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main=paste("May-October",names[ff]),cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# }
# ColorBar(bb,cm,blab)
# dev.off()

names=c("LAP 150 km","LAP 50km","PG 150km","PG 50km")
pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_location_cool_bymethod.pdf",sep=""),width=9,height=8)
layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
for(i in 1:4)
{
  image(lon,lat,t(apply(locNC[,,,i],c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main=names[i],cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb2[,t],cm2)
dev.off()

pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_location_warm_bymethod.pdf",sep=""),width=9,height=8)
layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
for(i in 1:4)
{
  image(lon,lat,t(apply(locNW[,,,i],c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main=names[i],cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb2[,t],cm2)
dev.off()

  names=c("LAP 150 km","LAP 50km","PG 150km","PG 50km")
  pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_location_trend_cool_bymethod.pdf",sep=""),width=9,height=8)
  layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
  for(i in 1:4)
  {
    image(lon,lat,t(apply(changeC[,,,i],c(1,2),median,na.rm=T)),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main=names[i],cex.axis=1.5,cex.main=2)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  }
  ColorBar(bb,cm,blab)
  dev.off()
  
  pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_location_trend_warm_bymethod.pdf",sep=""),width=9,height=8)
  layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
  for(i in 1:4)
  {
    image(lon,lat,t(apply(changeW[,,,i],c(1,2),median,na.rm=T)),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main=names[i],cex.axis=1.5,cex.main=2)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  }
  ColorBar(bb,cm,blab)
  dev.off()

# pdf(file=paste("Figures/ECLfixes_thresh",thresh[t],"_change_location_medianPC.pdf",sep=""),width=9,height=8)
# layout(cbind(c(1,3),c(2,4),c(5,6)),width=c(1,1,0.35),height=c(1,1))
# image(lon,lat,t(apply(changeC,c(1,2),median,na.rm=T)),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October median change",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# image(lon,lat,t(apply(changeW,c(1,2),median,na.rm=T)),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April median change",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# 
# image(lon,lat,t(apply(changeC>0,c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb3,col=cm3,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October % with increase",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# image(lon,lat,t(apply(changeW>0,c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb3,col=cm3,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April % with increase",cex.axis=1.5,cex.main=2)
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# ColorBar(bb,cm,blab)
# ColorBar(bb3,cm3,blab3)
# dev.off()

}


####### Okay, what are the overall projections for ECL fixes in the "near" region, ~500km from coast

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
intcol=c(10,10,9,9)
nowC<-array(NaN,dim=c(12,4,5))
dimnames(nowC)[[1]]<-rep("aaa",12)
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(nowC)[[1]][n]<-paste(cmip[i],wrf[j])
    n=n+1
  }
dimnames(nowC)[[2]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
dimnames(nowC)[[3]]=c(22,15,10,5,2)
nowW<-futureC<-futureW<-nowC

sigens<-matrix(NaN,4,2)
colnames(sigens)=c("Cool","Warm")
rownames(sigens)=c(22,15,10,5,2)

##What about for different duration thresholds?
thresh=c(22,15,10,5,2)
for(t in 1:5)
{
  n=1
  
  for(i in 1:4)
    for(j in 1:3)
    {    
      filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
      filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
      
      for(ff in 1:4)
      {
        data=read.csv(filelist1[ff])
        data$Location2=0
        I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
        data$Location2[I]<-1
        I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
        data$Location2[I]<-1
        I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
        data$Location2[I]<-1
        data=data[data$Location2==1,]
        
        b=order(data[,intcol[ff]],decreasing=T)
        thresh2=data[b[40*thresh[t]],intcol[ff]] ##Assume ~ 3 fixes per event.
        if(is.na(thresh2)) thresh2=min(data[,intcol[ff]])
        dataN=data[data[,intcol[ff]]>=thresh2,]            
        mm=floor((dataN$Date%%10000)/100)
        I=which(mm>=5 & mm<=10)
        nowC[n,ff,t]=length(I)
        I=which((mm>=11 | mm<=4))
        nowW[n,ff,t]=length(I)
        
        data=read.csv(filelist2[ff])
        data$Location2=0
        I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
        data$Location2[I]<-1
        I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
        data$Location2[I]<-1
        I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
        data$Location2[I]<-1
        data=data[data$Location2==1,]        
        dataF=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((dataF$Date%%10000)/100)
        I=which(mm>=5 & mm<=10)
        futureC[n,ff,t]=length(I)
        I=which((mm>=11 | mm<=4))
        futureW[n,ff,t]=length(I)
      }      
      n=n+1  
    }
}

changeC=((futureC/nowC)-1)*100
changeW=((futureW/nowW)-1)*100

for(t in 1:5)
{
  bysource=matrix(NaN,4,3)
  colnames(bysource)=c("GCM","RCM","Method")
  bycmip=matrix(NaN,9,4)
  colnames(bycmip)=cmip
  bywrf<-bymethod<-matrix(NaN,12,3)
  colnames(bywrf)=wrf
  colnames(bymethod)=c("GV","LAP 150km","PG 150km")
  
  tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_bymethod_boxplot_cool.tiff",sep=""),width=500,height=400)
  boxplot(changeC[,,t],main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
  abline(h=0,col="red",lwd=2)
  dev.off()
  tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_bymethod_boxplot_warm.tiff",sep=""),width=500,height=400)
  boxplot(changeW[,,t],main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
  abline(h=0,col="red",lwd=2)
  dev.off()
    
  change4=changeC[,c(1,2,4),t]
    
    for(i in 1:4)
      for(j in 1:3)
        for(k in 1:3)
          bycmip[3*(k-1)+j,i]<-bywrf[3*(i-1)+k,j]<-bymethod[3*(i-1)+j,k]<-change4[3*(i-1)+j,k]
    
    for(i in 1:4)
    {
      a=change4[seq(3*(i-1)+1,3*i),]
      bysource[i,1]=mean(a)
    }
    for(j in 1:3)
    {
      a=change4[seq(j,12,3),]
      bysource[j,2]=mean(a)
    }
    for(k in 1:3)
    {
      a=change4[,k]
      bysource[k,3]=mean(a)
    }
    
    tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_boxplot_bysource_cool.tiff",sep=""),width=500,height=400)
    boxplot(bysource,main="Range of trends by source of uncertainty",ylab="% change in ECLs",ylim=c(-30,30))
    abline(h=0,col="red",lwd=2)
    dev.off()
    tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_boxplot_byGCM_cool.tiff",sep=""),width=500,height=400)
    boxplot(bycmip,main="Range of trends by GCM",ylab="% change in ECLs")
    abline(h=0,col="red",lwd=2)
    dev.off()
    tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_boxplot_byRCM_cool.tiff",sep=""),width=500,height=400)
    boxplot(bywrf,main="Range of trends by RCM",ylab="% change in ECLs")
    abline(h=0,col="red",lwd=2)
    dev.off()
    tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_boxplot_bymethod_cool.tiff",sep=""),width=500,height=400)
    boxplot(bymethod,main="Range of trends by ECL detection method",ylab="% change in ECLs")
    abline(h=0,col="red",lwd=2)
    dev.off()
  
  change4=changeW[,c(1,2,4),t]
  
  for(i in 1:4)
    for(j in 1:3)
      for(k in 1:3)
        bycmip[3*(k-1)+j,i]<-bywrf[3*(i-1)+k,j]<-bymethod[3*(i-1)+j,k]<-change4[3*(i-1)+j,k]
  
  for(i in 1:4)
  {
    a=change4[seq(3*(i-1)+1,3*i),]
    bysource[i,1]=mean(a)
  }
  for(j in 1:3)
  {
    a=change4[seq(j,12,3),]
    bysource[j,2]=mean(a)
  }
  for(k in 1:3)
  {
    a=change4[,k]
    bysource[k,3]=mean(a)
  }
  
  tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_boxplot_bysource_warm.tiff",sep=""),width=500,height=400)
  boxplot(bysource,main="Range of trends by source of uncertainty",ylab="% change in ECLs",ylim=c(-30,30))
  abline(h=0,col="red",lwd=2)
  dev.off()
  tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_boxplot_byGCM_warm.tiff",sep=""),width=500,height=400)
  boxplot(bycmip,main="Range of trends by GCM",ylab="% change in ECLs")
  abline(h=0,col="red",lwd=2)
  dev.off()
  tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_boxplot_byRCM_warm.tiff",sep=""),width=500,height=400)
  boxplot(bywrf,main="Range of trends by RCM",ylab="% change in ECLs")
  abline(h=0,col="red",lwd=2)
  dev.off()
  tiff(file=paste("Figures/ECLfixCLOSE_thresh",thresh[t],"_change_boxplot_bymethod_warm.tiff",sep=""),width=500,height=400)
  boxplot(bymethod,main="Range of trends by ECL detection method",ylab="% change in ECLs")
  abline(h=0,col="red",lwd=2)
  dev.off()
  
}




######Other threshes

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
intcol=c(10,10,9,9)
nowC<-array(NaN,dim=c(12,4,4))
dimnames(nowC)[[1]]<-rep("aaa",12)
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(nowC)[[1]][n]<-paste(cmip[i],wrf[j])
    n=n+1
  }
dimnames(nowC)[[2]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
dimnames(nowC)[[3]]=c("20-15","15-10","10-5","<5")
nowW<-futureC<-futureW<-nowC


##What about for different duration thresholds?
thresh=c(20,15,10,5,0)
for(t in 1:4)
{
  n=1
  
  for(i in 1:4)
    for(j in 1:3)
    {    
      filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
      filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
      
      for(ff in 1:4)
      {
        data=read.csv(filelist1[ff])
        data$Location2=0
        I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
        data$Location2[I]<-1
        I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
        data$Location2[I]<-1
        I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
        data$Location2[I]<-1
        data=data[data$Location2==1,]
        
        b=order(data[,intcol[ff]],decreasing=T)
        thresh1=data[b[40*thresh[t]],intcol[ff]] ##Assume ~ 3 fixes per event.
        if(is.na(thresh1)) thresh1=min(data[,intcol[ff]])
        if(t<4) thresh2=data[b[40*thresh[t+1]],intcol[ff]] else thresh2=9999999

        dataN=data[(data[,intcol[ff]]>=thresh1 & data[,intcol[ff]]<thresh2),]       
        mm=floor((dataN$Date%%10000)/100)
        I=which(mm>=5 & mm<=10)
        nowC[n,ff,t]=length(I)
        I=which((mm>=11 | mm<=4))
        nowW[n,ff,t]=length(I)
        
        data=read.csv(filelist2[ff])
        data$Location2=0
        I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
        data$Location2[I]<-1
        I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
        data$Location2[I]<-1
        I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
        data$Location2[I]<-1
        data=data[data$Location2==1,]        
        dataF=data[(data[,intcol[ff]]>=thresh1 & data[,intcol[ff]]<thresh2),]
        mm=floor((dataF$Date%%10000)/100)
        I=which(mm>=5 & mm<=10)
        futureC[n,ff,t]=length(I)
        I=which((mm>=11 | mm<=4))
        futureW[n,ff,t]=length(I)
      }      
      n=n+1  
    }
}

changeC=((futureC/nowC)-1)*100
changeW=((futureW/nowW)-1)*100

bythreshC<-bythreshW<-matrix(NaN,48,4)
colnames(bythreshC)<-colnames(bythreshW)<-c("20-15","15-10","10-5","<5")

for(i in 1:4)
  for(j in 1:4)
  {
    bythreshC[((12*(i-1)+1):(12*i)),j]<-changeC[,i,j]
    bythreshW[((12*(i-1)+1):(12*i)),j]<-changeW[,i,j]
  }

tiff(file=paste("ECLfixCLOSE_change_boxplot_bythresh_cool.tiff",sep=""),width=500,height=400)
boxplot(bythreshC,main="Range of trends by threshold (ECLs p.a.)",ylab="% change in ECLs",ylim=c(-100,100))
abline(h=0,col="red",lwd=2)
dev.off()
tiff(file=paste("ECLfixCLOSE_change_boxplot_bythresh_warm.tiff",sep=""),width=500,height=400)
boxplot(bythreshW,main="Range of trends by threshold (ECLs p.a.)",ylab="% change in ECLs",ylim=c(-100,100))
abline(h=0,col="red",lwd=2)
dev.off()




rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

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
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))

#lat=seq(-42.5,-22.5,5)
#lon=seq(147.5,157.5,5)

lat=seq(-40,-25,2.5)
lon=seq(145,160,2.5)

locNC<-locFC<-locFW<-locNW<-array(0,dim=c(7,7,12,4)) #Lon, Lat, Model, Method
intcol=c(10,10,9,9)
intcolE=rep(10,4)
cm=pal(12)
bb=c(-10000,-100,-50,-25,-10,-5,0,5,10,25,50,100,10000)
blab=c("-100%","-50%","-25%","-10%","-5%","","+5%","+10%","+25%","+50%","+100%")

cm2=gray(seq(1,0.1,-0.15))
#bb2=cbind(c(-0.5,0,0.2,0.5,1,1.5,2,100),c(-0.5,0,0.5,1,1.5,2,3,100),c(-0.5,0,1,2,3,4,5,100),c(-0.5,0,1,2.5,5,7.5,10,100),c(-0.5,0,1,2.5,5,7.5,10,100))
bb2=cbind(c(-0.5,0,0.1,0.2,0.3,0.4,0.5,100),c(-0.5,0,0.2,0.4,0.6,0.8,1,100),c(-0.5,0,0.2,0.5,1,1.5,2,100),c(-0.5,0,0.5,1,1.5,2,3,100),c(-0.5,0,1,2,3,4,5,100))

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
names=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
thresh=22
  n=1
  
  for(i in 1:4)
    for(j in 1:3)
    {    
      filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
      
      filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
      filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
      
      for(ff in 1:4)
      {
        data=read.csv(filelistE[ff])
        b=order(data[,intcol[ff]],decreasing=T)
        thresh2=data[b[20*thresh],intcolE[ff]]
        if(is.na(thresh2)) thresh2=min(data[,intcolE[ff]])
        
        data=read.csv(filelist1[ff])
        dataN=data[data[,intcol[ff]]>=thresh2,]      
        data=read.csv(filelist2[ff])
        dataF=data[data[,intcol[ff]]>=thresh2,]
        
        for(xx in 1:7)
          for(yy in 1:7)
          {     
            I=which(dataN$Lat>=(lat[yy]-1.25) & dataN$Lat<=(lat[yy]+1.25) & dataN$Lon<=(lon[xx]+1.25) & dataN$Lon>=(lon[xx]-1.25))
            mm=floor((dataN[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            locNC[yy,xx,n,ff]=length(I)
            locNW[yy,xx,n,ff]=length(mm)-length(I)     
            
            I=which(dataF$Lat>=lat[yy]-1.25 & dataF$Lat<=lat[yy]+1.25 & dataF$Lon<=lon[xx]+1.25 & dataF$Lon>=lon[xx]-1.25)
            mm=floor((dataF[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            locFC[yy,xx,n,ff]=length(I)
            locFW[yy,xx,n,ff]=length(mm)-length(I)   
          }    
      }      
      n=n+1  
    }
  
  changeC=((locFC/locNC)-1)*100
  changeC[locNC==0]=NaN
  changeW=((locFW/locNW)-1)*100
  changeW[locNW==0]=NaN
  changeCa=apply(changeC,c(1,2),median,na.rm=T)
  changeWa=apply(changeW,c(1,2),median,na.rm=T)
  

  bb2=c(-0.5,0,0.25,0.5,0.75,1,1.25,100)
  pdf(file=paste("Figures/ECLfixes_Eventthresh",thresh,"_location_v2_mean.pdf",sep=""),width=9,height=8)
  layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
  image(lon,lat,t(apply(locNC,c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2,col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October 1990-09",cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  image(lon,lat,t(apply(locFC,c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2,col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October 2060-79",cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  image(lon,lat,t(apply(locNW,c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2,col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April 1990-09",cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  image(lon,lat,t(apply(locFW,c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2,col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April 2060-79",cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  ColorBar(bb2,cm2)
  dev.off()
  print(range(t(apply(locNC,c(1,2),median,na.rm=T))/20,t(apply(locFW,c(1,2),median,na.rm=T))/20))



########## Fixes above the event threshold
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

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
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))

#lat=seq(-42.5,-22.5,5)
#lon=seq(147.5,157.5,5)

lat=seq(-40,-25,2.5)
lon=seq(145,160,2.5)


intcol=c(10,10,9,9)
intcolE=rep(10,4)
cm=pal(12)
bb=c(-10000,-100,-50,-25,-10,-5,0,5,10,25,50,100,10000)
blab=c("-100%","-50%","-25%","-10%","-5%","","+5%","+10%","+25%","+50%","+100%")

cm3=pal(10)
bb3=c(-10000,0.1,0.25,0.35,0.45,0.5,0.55,0.65,0.75,0.9,10000)
blab3=c("10%","25%","35%","45%","","55%","65%","75%","90%")

cm2=gray(seq(1,0.1,-0.15))
#bb2=cbind(c(-0.5,0,0.2,0.5,1,1.5,2,100),c(-0.5,0,0.5,1,1.5,2,3,100),c(-0.5,0,1,2,3,4,5,100),c(-0.5,0,1,2.5,5,7.5,10,100),c(-0.5,0,1,2.5,5,7.5,10,100))
bb2=cbind(c(-0.5,0,0.1,0.2,0.3,0.4,0.5,100),c(-0.5,0,0.2,0.4,0.6,0.8,1,100),c(-0.5,0,0.25,0.5,0.75,1,1.5,100),c(-0.5,0,0.25,0.5,0.75,1,1.5,100),c(-0.5,0,1,2,3,4,5,100))

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
names=c("LAP 150km","LAP 50km","PG 150km","PG 50km")


thresh=c(2,5,10,22)
t=4
for(t in 4)
{
  n=1
  locNC<-locFC<-locFW<-locNW<-array(0,dim=c(7,7,12,4)) #Lon, Lat, Model, Method
  for(i in 1:4)
    for(j in 1:3)
    {    
      filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
      
      filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
      filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
      
      for(ff in 1:4)
      {
        data=read.csv(filelistE[ff])
        b=order(data[,intcol[ff]],decreasing=T)
        thresh2=data[b[20*thresh[t]],intcolE[ff]]
        if(is.na(thresh2)) thresh2=min(data[,intcolE[ff]])
        
        data=read.csv(filelist1[ff])
        dataN=data[data[,intcol[ff]]>=thresh2,]      
        data=read.csv(filelist2[ff])
        dataF=data[data[,intcol[ff]]>=thresh2,]
        
        for(xx in 1:7)
          for(yy in 1:7)
          {     
            I=which(dataN$Lat>=(lat[yy]-1.25) & dataN$Lat<=(lat[yy]+1.25) & dataN$Lon<=(lon[xx]+1.25) & dataN$Lon>=(lon[xx]-1.25))
            mm=floor((dataN[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            locNC[yy,xx,n,ff]=length(I)
            locNW[yy,xx,n,ff]=length(mm)-length(I)     
            
            I=which(dataF$Lat>=lat[yy]-1.25 & dataF$Lat<=lat[yy]+1.25 & dataF$Lon<=lon[xx]+1.25 & dataF$Lon>=lon[xx]-1.25)
            mm=floor((dataF[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            locFC[yy,xx,n,ff]=length(I)
            locFW[yy,xx,n,ff]=length(mm)-length(I)   
          }    
      }      
      n=n+1  
    }
  
  
  changeC=((locFC/locNC)-1)*100
  changeC[locNC==0]=NaN
  changeW=((locFW/locNW)-1)*100
  changeW[locNW==0]=NaN
  changeCa=apply(changeC,c(1,2),median,na.rm=T)
  changeWa=apply(changeW,c(1,2),median,na.rm=T)
  
  
#   pdf(file=paste("Figures/ECLfixes_Eventthresh",thresh[t],"_location_v2_mean.pdf",sep=""),width=9,height=8)
#   layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
#   image(lon,lat,t(apply(locNC,c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October 1990-09",cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   image(lon,lat,t(apply(locFC,c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October 2060-79",cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   image(lon,lat,t(apply(locNW,c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April 1990-09",cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   image(lon,lat,t(apply(locFW,c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April 2060-79",cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   ColorBar(bb2[,t],cm2)
#   dev.off()
#   pdf(file=paste("Figures/ECLfixes_Eventthresh",thresh[t],"_location_v2_median.pdf",sep=""),width=9,height=8)
#   layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
#   image(lon,lat,t(apply(locNC,c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October 1990-09",cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   image(lon,lat,t(apply(locFC,c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October 2060-79",cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   image(lon,lat,t(apply(locNW,c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April 1990-09",cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   image(lon,lat,t(apply(locFW,c(1,2),median,na.rm=T))/20,xlab="",ylab="",breaks=bb2[,t],col=cm2,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April 2060-79",cex.axis=1.5,cex.main=2)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   ColorBar(bb2[,t],cm2)
#   dev.off()
#   
#   names=c("LAP 150 km","LAP 50km","PG 150km","PG 50km")
#   pdf(file=paste("Figures/ECLfixes_Eventthresh",thresh[t],"_location_trend_mean_cool_bymethod.pdf",sep=""),width=9,height=8)
#   layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
#   for(i in 1:4)
#   {
#     image(lon,lat,t(apply(changeC[,,,i],c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main=names[i],cex.axis=1.5,cex.main=2)
#     contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   }
#   ColorBar(bb,cm,blab)
#   dev.off()
#   
#   pdf(file=paste("Figures/ECLfixes_Eventthresh",thresh[t],"_location_trend_mean_warm_bymethod.pdf",sep=""),width=9,height=8)
#   layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.35),height=c(1,1))
#   for(i in 1:4)
#   {
#     image(lon,lat,t(apply(changeW[,,,i],c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main=names[i],cex.axis=1.5,cex.main=2)
#     contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   }
#   ColorBar(bb,cm,blab)
#   dev.off()
  
    pdf(file=paste("Figures/ECLfixes_Eventthresh",thresh[t],"_trend_location_medianPC.pdf",sep=""),width=9,height=8)
    layout(cbind(c(1,3),c(2,4),c(5,6)),width=c(1,1,0.35),height=c(1,1))
    image(lon,lat,t(apply(changeC,c(1,2),median,na.rm=T)),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October median change",cex.axis=1.5,cex.main=2)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    image(lon,lat,t(apply(changeW,c(1,2),median,na.rm=T)),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April median change",cex.axis=1.5,cex.main=2)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    
    image(lon,lat,t(apply(changeC>0,c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb3,col=cm3,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="May-October % with increase",cex.axis=1.5,cex.main=2)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    image(lon,lat,t(apply(changeW>0,c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb3,col=cm3,zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="November-April % with increase",cex.axis=1.5,cex.main=2)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    ColorBar(bb,cm,blab)
    ColorBar(bb3,cm3,blab3)
    dev.off()
}


###### Save colorbars for reference with large text

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
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))

cm=pal(12)
bb=c(-10000,-100,-50,-25,-10,-5,0,5,10,25,50,100,10000)
blab=c("-100%","-50%","-25%","-10%","-5%","","+5%","+10%","+25%","+50%","+100%")

cm2=gray(seq(1,0.1,-0.15))
bb2=c(-0.5,0,0.5,1,1.5,2,3,100)

pdf(file="Location_colorbars.pdf",width=3,height=5)
layout(cbind(1,2),height=1,width=c(1,1))
ColorBar(bb2,cm2)
ColorBar(bb,cm,blab)
dev.off()

############ Whole cordex frequencies

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

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
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))

#lat=seq(-42.5,-22.5,5)
#lon=seq(147.5,157.5,5)

lat=seq(-45,-10,2.5)
lon=seq(110,175,2.5)
intcol<-intcolE<-10

cm=pal(12)
bb=c(-10000,-100,-50,-25,-10,-5,0,5,10,25,50,100,10000)
blab=c("-100%","-50%","-25%","-10%","-5%","","+5%","+10%","+25%","+50%","+100%")

cm2=gray(seq(1,0.1,-0.15))
#bb2=cbind(c(-0.5,0,0.2,0.5,1,1.5,2,100),c(-0.5,0,0.5,1,1.5,2,3,100),c(-0.5,0,1,2,3,4,5,100),c(-0.5,0,1,2.5,5,7.5,10,100),c(-0.5,0,1,2.5,5,7.5,10,100))
bb2=c(-0.5,0,1,2,3,4,5,100)

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
names=c("LAP 150km","LAP 50km","PG 150km","PG 50km")

thresh=c(2,5,10,22)
t=4
n=1

locC<-locW<-array(0,dim=c(15,27))
filelistE=paste("outputUM/proj100/outputUM_ncep_rad2cv06/ECLevents_umelb_ncep_proj100_rad2cv06_9009.csv",sep="")
filelist1=paste("outputUM/proj100/outputUM_ncep_rad2cv06/ECLfixes_umelb_ncep_proj100_rad2cv06_9009_ALL.csv",sep="")

data=read.csv(filelistE)
b=order(data[,intcol],decreasing=T)
thresh2=data[b[20*thresh[t]],intcolE]
if(is.na(thresh2)) thresh2=min(data[,intcolE])

data=read.csv(filelist1)
dataN=data[data[,intcol]>=thresh2,]      

for(xx in 1:27)
  for(yy in 1:15)
  {     
    I=which(dataN$Lat>=(lat[yy]-1.25) & dataN$Lat<=(lat[yy]+1.25) & dataN$Lon<=(lon[xx]+1.25) & dataN$Lon>=(lon[xx]-1.25))
    mm=floor((dataN[I,]$Date%%10000)/100)
    I=which(mm>=5 & mm<=10)
    locC[yy,xx]=length(I)
    locW[yy,xx]=length(mm)-length(I)     
  }


  n=1
  locNC<-locFC<-locFW<-locNW<-array(0,dim=c(15,27,15)) #Lon, Lat, Model, Method
  for(i in 1:5)
    for(j in 1:3)
    {    
      filelistE=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_9009.csv",sep="")
      filelist1=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_rad2cv06_9009_ALL.csv",sep="")
      filelist2=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_rad2cv06_6079_ALL.csv",sep="")
      
      data=read.csv(filelistE)
      b=order(data[,intcol],decreasing=T)
      thresh2=data[b[20*thresh[t]],intcolE]
        if(is.na(thresh2)) thresh2=min(data[,intcolE])
        
        data=read.csv(filelist1)
        dataN=data[data[,intcol]>=thresh2,]      
        if(i>1)
          {data=read.csv(filelist2)
        dataF=data[data[,intcol]>=thresh2,]
        }
        
        for(xx in 1:27)
          for(yy in 1:15)
          {     
            I=which(dataN$Lat>=(lat[yy]-1.25) & dataN$Lat<=(lat[yy]+1.25) & dataN$Lon<=(lon[xx]+1.25) & dataN$Lon>=(lon[xx]-1.25))
            mm=floor((dataN[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            locNC[yy,xx,n]=length(I)
            locNW[yy,xx,n]=length(mm)-length(I)     
            
            if(i>1)
            {
            I=which(dataF$Lat>=lat[yy]-1.25 & dataF$Lat<=lat[yy]+1.25 & dataF$Lon<=lon[xx]+1.25 & dataF$Lon>=lon[xx]-1.25)
            mm=floor((dataF[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            locFC[yy,xx,n]=length(I)
            locFW[yy,xx,n]=length(mm)-length(I)   
            }
          }    
      n=n+1  
    }
  

  
  pdf(file=paste("Figures/ECLfixes_Eventthresh",thresh[t],"_location_v2_mean_CORDEX_LAP.pdf",sep=""),width=5,height=4)
  layout(c(1,2,3))
  
  image(lon,lat,t(locC+locW)/20,xlab="",ylab="",breaks=bb2,col=cm2,zlim=c(-1000,1000),xlim=c(110,175),ylim=c(-45,-10),main="NCEP",cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  image(lon,lat,t(apply(locNC[,,1:3]+locNW[,,1:3],c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2,col=cm2,
        zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="NCEP-WRF",cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  image(lon,lat,t(apply(locNW[,,4:15]+locNC[,,4:15],c(1,2),mean,na.rm=T))/20,xlab="",ylab="",breaks=bb2,col=cm2,
        zlim=c(-1000,1000),xlim=c(145,160),ylim=c(-40,-25),main="CMIP-WRF",cex.axis=1.5,cex.main=2)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  ColorBar(bb2,cm2)
  dev.off()
  
  
}


