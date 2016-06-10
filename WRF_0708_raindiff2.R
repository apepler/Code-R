rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3
dir='~/output/outputUM_wrf_2007_all/'
events<-eventsNT<-fixes<-fixesNT<-comp<-compNT<-list()
compa<-compNTa<-array(0,c(21,21,3))
compb<-compNTb<-array(0,c(101,101,3))

n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    events[[n]]=rbind(read.csv(paste(dir,"ECLevents_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                      read.csv(paste(dir,"ECLevents_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    
    eventsNT[[n]]=rbind(read.csv(paste(dir,"ECLevents_",dom,"_2007_R",r,"_notopo_",cat[c],".csv",sep="")),
                        read.csv(paste(dir,"ECLevents_",dom,"_2008_R",r,"_notopo_",cat[c],".csv",sep="")))
    eventsNT[[n]]$Year=floor(eventsNT[[n]]$Date1/10000)
    eventsNT[[n]]$Month=floor(eventsNT[[n]]$Date1/100)%%100
    
    fixes[[n]]=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                     read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    
    fixesNT[[n]]=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,"_notopo_",cat[c],".csv",sep="")),
                       read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,"_notopo_",cat[c],".csv",sep="")))
    fixesNT[[n]]$Year=floor(fixesNT[[n]]$Date/10000)
    fixesNT[[n]]$Month=floor(fixesNT[[n]]$Date/100)%%100
    
    fixes[[n]]$Location2<-0
    I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    
    fixesNT[[n]]$Location2<-0
    I<-which(fixesNT[[n]][,7]>=149 & fixesNT[[n]][,7]<=154 & fixesNT[[n]][,8]<(-37) & fixesNT[[n]][,8]>=-41)
    fixesNT[[n]]$Location2[I]<-1
    I<-which(fixesNT[[n]][,7]>=(149+(37+fixesNT[[n]][,8])/2) & fixesNT[[n]][,7]<=(154+(37+fixesNT[[n]][,8])/2) & fixesNT[[n]][,8]<(-31) & fixesNT[[n]][,8]>=-37)
    fixesNT[[n]]$Location2[I]<-1
    I<-which(fixesNT[[n]][,7]>=152 & fixesNT[[n]][,7]<=157 & fixesNT[[n]][,8]<=(-24) & fixesNT[[n]][,8]>=-31)
    fixesNT[[n]]$Location2[I]<-1
    
    
    ### Rain stuff

    if(dom=="d01") a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLrain_0708_",cat[c],".nc",sep="")) else
      a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLrain_d02_0708_",cat[c],".nc",sep=""))
    tmp=var.get.nc(a,"ECLrain")
    fixes[[n]]$MeanRain=apply(tmp,3,mean,na.rm=T)
    fixes[[n]]$MaxRain=apply(tmp,3,max,na.rm=T)
    I=which(fixes[[n]]$Location2==1)
    if(dom=="d01") {
      fixes[[n]]$MeanRainS=apply(tmp[,1:10,],3,mean,na.rm=T)
      fixes[[n]]$MeanRainN=apply(tmp[,12:21,],3,mean,na.rm=T)
      compa[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T) 
      } else {
        fixes[[n]]$MeanRainS=apply(tmp[,1:50,],3,mean,na.rm=T)
        fixes[[n]]$MeanRainN=apply(tmp[,52:101,],3,mean,na.rm=T)
        compb[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
      }
    comp[[n]]=tmp
    
    if(dom=="d01") a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_notopo/out/ECLrain_0708_",cat[c],".nc",sep="")) else
      a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_notopo/out/ECLrain_d02_0708_",cat[c],".nc",sep=""))
    tmp=var.get.nc(a,"ECLrain")
    fixesNT[[n]]$MeanRain=apply(tmp,3,mean,na.rm=T)
    fixesNT[[n]]$MaxRain=apply(tmp,3,max,na.rm=T)
    I=which(fixesNT[[n]]$Location2==1)
    if(dom=="d01") {
      compNTa[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T) 
      fixesNT[[n]]$MeanRainS=apply(tmp[,1:10,],3,mean,na.rm=T)
      fixesNT[[n]]$MeanRainN=apply(tmp[,12:21,],3,mean,na.rm=T)
      } else {
        compNTb[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
        fixesNT[[n]]$MeanRainS=apply(tmp[,1:50,],3,mean,na.rm=T)
        fixesNT[[n]]$MeanRainN=apply(tmp[,52:101,],3,mean,na.rm=T)
      }
    compNT[[n]]=tmp
    
    n=n+1     
  }

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
#bb2=c(-100,-4,-3,-2,-1,0,1,2,3,4,100)
bb2=c(-100,-5,-2,-1,0,1,2,5,100)
layout(cbind(1,2,3))
par(mar=c(1,1,3,1))
for(i in 1:3) image(compNTb[,,i]-compb[,,i],col=pal(8),breaks=bb2,main=paste("RCM",i),axes=F)

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = brks[seq(1, length(brks), subsampleg)])
}
layout(cbind(1,2,3),width=c(1,1,0.3))
image(apply(compb,c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main="Control",axes=F,cex.main=2)
image(apply(compNTb,c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main="NoTopo",axes=F,cex.main=2)
ColorBar(c(seq(0,14,1),100),tim.colors(15),subsampleg=2)

layout(cbind(1,2,3),width=c(1,1,0.3))
image(apply(compNTa-compa,c(1,2),mean),col=pal(10),breaks=bb2,main="50 km",axes=F,cex.main=2)
image(apply(compNTb-compb,c(1,2),mean),col=pal(10),breaks=bb2,main="10 km",axes=F,cex.main=2)
ColorBar(col=pal(10),brks=bb2)

## Playing with N/S rain dist
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))
par(mar=c(1,1,3,1))
n=1
I=which(fixes[[n]]$MeanRainN>=1.2*fixes[[n]]$MeanRainS & fixes[[n]]$Location==1)
image(apply(comp[[n]][,,I],c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main="N",axes=F,cex.main=2)
I=which(fixes[[n]]$MeanRainN<1.2*fixes[[n]]$MeanRainS & fixes[[n]]$MeanRainS<1.2*fixes[[n]]$MeanRainN & fixes[[n]]$Location==1)
image(apply(comp[[n]][,,I],c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main="C",axes=F,cex.main=2)
I=which(fixes[[n]]$MeanRainS>=1.2*fixes[[n]]$MeanRainN & fixes[[n]]$Location==1)
image(apply(comp[[n]][,,I],c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main="S",axes=F,cex.main=2)
ColorBar(c(seq(0,14,1),100),tim.colors(15),subsampleg=2)

n=1
plot(fixes[[n]]$MeanRainS/fixes[[n]]$MeanRainN,fixes[[n]]$MeanRain,xlab="S/N Ratio",ylab="Mean rain",
     pch=4,lwd=2,cex=1.5,log="x")


layout(1)
par(mar=c(5,4,4,2))
meanR<-matrix(0,512,12)
for(i in 1:3)
{
  a=density(fixes[[i+3]]$MeanRain,from=0,to=15)
  meanR[,i]=a$y
  a=density(fixesNT[[i+3]]$MeanRain,from=0,to=15)
  meanR[,i+3]=a$y
  a=density(fixes[[i+3]]$MeanRain[fixes[[i]]$Location2==1],from=0,to=15)
  meanR[,i+6]=a$y
  a=density(fixesNT[[i+3]]$MeanRain[fixesNT[[i]]$Location2==1],from=0,to=15)
  meanR[,i+9]=a$y
}

plot(a$x,apply(meanR[,1:3],1,mean),type="l",col="blue",lwd=3,
     xlab="Mean rainfall (mm)",ylab="Frequency",ylim=c(0,0.3))
lines(a$x,apply(meanR[,4:6],1,mean),col="red",lwd=3)
lines(a$x,apply(meanR[,7:9],1,mean),col="blue",lwd=3,lty=2)
lines(a$x,apply(meanR[,10:12],1,mean),col="red",lwd=3,lty=2)
legend("topright",legend=c("All lows","All lows, NoTopo","Coastal lows","Coastal lows, NoTopo"),
       lwd=3,lty=c(1,1,2,2),col=c("blue","red","blue","red"))



###### What about SST changes for R2?

rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3
dir='~/output/outputUM_wrf_2007_all/'
events<-eventsBRAN<-fixes<-fixesBRAN<-comp<-compBRAN<-events_noeac<-fixes_noeac<-comp_noeac<-list()
compa<-compBRANa<-comp_noeaca<-array(0,c(21,21,3))
compb<-compBRANb<-comp_noeacb<-array(0,c(101,101,3))
n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {

    fixes[[n]]=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                     read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    
    fixesBRAN[[n]]=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,"_BRAN_",cat[c],".csv",sep="")),
                       read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,"_BRAN_",cat[c],".csv",sep="")))
    fixesBRAN[[n]]$Year=floor(fixesBRAN[[n]]$Date/10000)
    fixesBRAN[[n]]$Month=floor(fixesBRAN[[n]]$Date/100)%%100
    
    fixes_noeac[[n]]=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")),
                       read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")))
    fixes_noeac[[n]]$Year=floor(fixes_noeac[[n]]$Date/10000)
    fixes_noeac[[n]]$Month=floor(fixes_noeac[[n]]$Date/100)%%100
    
    fixes_noeac[[n]]$Location2<-0
    I<-which(fixes_noeac[[n]][,7]>=149 & fixes_noeac[[n]][,7]<=154 & fixes_noeac[[n]][,8]<(-37) & fixes_noeac[[n]][,8]>=-41)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=(149+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,7]<=(154+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,8]<(-31) & fixes_noeac[[n]][,8]>=-37)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=152 & fixes_noeac[[n]][,7]<=157 & fixes_noeac[[n]][,8]<=(-24) & fixes_noeac[[n]][,8]>=-31)
    fixes_noeac[[n]]$Location2[I]<-1
    
    fixes[[n]]$Location2<-0
    I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    
    fixesBRAN[[n]]$Location2<-0
    I<-which(fixesBRAN[[n]][,7]>=149 & fixesBRAN[[n]][,7]<=154 & fixesBRAN[[n]][,8]<(-37) & fixesBRAN[[n]][,8]>=-41)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=(149+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,7]<=(154+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,8]<(-31) & fixesBRAN[[n]][,8]>=-37)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=152 & fixesBRAN[[n]][,7]<=157 & fixesBRAN[[n]][,8]<=(-24) & fixesBRAN[[n]][,8]>=-31)
    fixesBRAN[[n]]$Location2[I]<-1
    
    
    ### Rain stuff
    
    if(dom=="d01") a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLrain_0708_",cat[c],".nc",sep="")) else
      a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLrain_d02_0708_",cat[c],".nc",sep=""))
    tmp=var.get.nc(a,"ECLrain")
    fixes[[n]]$MeanRain=apply(tmp,3,mean,na.rm=T)
    fixes[[n]]$MaxRain=apply(tmp,3,max,na.rm=T)
    I=which(fixes[[n]]$Location2==1)
    if(dom=="d01") compa[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T) else compb[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
    comp[[n]]=tmp
    
    if(dom=="d01") a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN/out/ECLrain_0708_",cat[c],".nc",sep="")) else
      a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN/out/ECLrain_d02_0708_",cat[c],".nc",sep=""))
    tmp=var.get.nc(a,"ECLrain")
    fixesBRAN[[n]]$MeanRain=apply(tmp,3,mean,na.rm=T)
    fixesBRAN[[n]]$MaxRain=apply(tmp,3,max,na.rm=T)
    I=which(fixesBRAN[[n]]$Location2==1)
    if(dom=="d01") compBRANa[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T) else compBRANb[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
    compBRAN[[n]]=tmp
    
    if(dom=="d01") a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN_noeac/out/ECLrain_0708_",cat[c],".nc",sep="")) else
      a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN_noeac/out/ECLrain_d02_0708_",cat[c],".nc",sep=""))
    tmp=var.get.nc(a,"ECLrain")
    fixes_noeac[[n]]$MeanRain=apply(tmp,3,mean,na.rm=T)
    fixes_noeac[[n]]$MaxRain=apply(tmp,3,max,na.rm=T)
    I=which(fixes_noeac[[n]]$Location2==1)
    if(dom=="d01") comp_noeaca[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T) else comp_noeacb[,,r]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
    comp_noeac[[n]]=tmp
    
    n=n+1     
  }

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}
layout(cbind(1,2,3),width=c(1,1,0.3))
image(apply(compBRANa,c(1,2),mean),breaks=c(seq(0,11,1),100),col=tim.colors(12),main="BRAN",axes=F)
image(apply(comp_noeaca,c(1,2),mean),breaks=c(seq(0,11,1),100),col=tim.colors(12),main="NoEAC",axes=F)
ColorBar(c(seq(0,11,1),100),tim.colors(12),subsampleg=2)

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
bb2=c(-100,-4,-3,-2,-1,0,1,2,3,4,100)
layout(cbind(1,2,3),width=c(1,1,0.3))
image(apply(comp_noeaca-compBRANa,c(1,2),mean),col=pal(10),breaks=bb2,main="50 km",axes=F,cex.main=2)
image(apply(comp_noeacb-compBRANb,c(1,2),mean),col=pal(10),breaks=bb2,main="10 km",axes=F,cex.main=2)
ColorBar(col=pal(10),brks=bb2)

layout(cbind(1,2,3),width=c(1,1,1))
for(i in 1:3) image(comp_noeacb[,,i]-compBRANb[,,i],col=pal(10),breaks=bb2,main=paste("RCM",i),axes=F,cex.main=2)


meanR<-matrix(0,512,12)
for(i in 1)
{
  a=density(fixesBRAN[[i]]$MeanRain,from=0,to=15)
  meanR[,i]=a$y
  a=density(fixes_noeac[[i]]$MeanRain,from=0,to=15)
  meanR[,i+3]=a$y
  a=density(fixesBRAN[[i]]$MeanRain[fixesBRAN[[i]]$Location2==1],from=0,to=15)
  meanR[,i+6]=a$y
  a=density(fixes_noeac[[i]]$MeanRain[fixes_noeac[[i]]$Location2==1],from=0,to=15)
  meanR[,i+9]=a$y
}


plot(a$x,meanR[,1],type="l",col="blue",lwd=3,
     xlab="Mean rainfall (mm)",ylab="Frequency",ylim=c(0,0.3))
lines(a$x,meanR[,4],col="red",lwd=3)
lines(a$x,meanR[,7],col="blue",lwd=3,lty=2)
lines(a$x,meanR[,10],col="red",lwd=3,lty=2)
legend("topright",legend=c("All lows","All lows, NoEAC","Coastal lows","Coastal lows, NoEAC"),
       lwd=3,lty=c(1,1,2,2),col=c("blue","red","blue","red"))


###### More assessment of rain pdf stuff

rm(list=ls())
library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3
dirs=c("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007/out/",
       "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007_notopo/out/",
       "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN/out/",
       "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_noeac/out/")

dayrain<-array(NaN,c(731,4,4))
for(i in 1:4) dayrain[,,i]=as.matrix(read.table(paste(dirs[i],"dayrain_d01.txt",sep="")))

Thresh=cbind(c(1,2,5,10,25,Inf),c(1,5,25,50,100,Inf),c(1,9,22.5,45,81,Inf),c(1,2,4,9,30,Inf))
Rcount=array(0,c(5,4,4))
dimnames(Rcount)[[2]]=c("Mean","Max","Cells>5","Cells>25")
dimnames(Rcount)[[3]]=c("Control","NoTopo","BRAN","NoEAC")
for(r in 1:4)
  for(s in 1:4)
  for(t in 1:5)
  {
    I=which(dayrain[,s,r]>=Thresh[t,s] & dayrain[,s,r]<Thresh[t+1,s] )
    Rcount[t,s,r]=length(I)
  }

dayECL<-array(0,c(731,4,4))
dimnames(dayECL)[[2]]=c("Control","NoTopo","BRAN","NoEAC")
dimnames(dayECL)[[3]]=c("Loc1","Loc1 +1","Loc2","Loc2+1")
daylist=seq(ISOdate(2007,1,1),ISOdate(2008,12,31),by="1 day")
date=as.numeric(format.Date(daylist,"%Y%m%d"))
type=c("","_notopo","_BRAN","_BRAN_noeac")
dom="d01"
r=2
dir="~/output/outputUM_wrf_2007_all/"
for(t in 1:4)
{
  fixes=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,type[t],"_",cat[c],".csv",sep="")),
              read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,type[t],"_",cat[c],".csv",sep="")))
  fixes$Year=floor(fixes$Date/10000)
  fixes$Month=floor(fixes$Date/100)%%100
  fixes$Location2<-0
  I<-which(fixes[,7]>=149 & fixes[,7]<=154 & fixes[,8]<(-37) & fixes[,8]>=-41)
  fixes$Location2[I]<-1
  I<-which(fixes[,7]>=(149+(37+fixes[,8])/2) & fixes[,7]<=(154+(37+fixes[,8])/2) & fixes[,8]<(-31) & fixes[,8]>=-37)
  fixes$Location2[I]<-1
  I<-which(fixes[,7]>=152 & fixes[,7]<=157 & fixes[,8]<=(-24) & fixes[,8]>=-31)
  fixes$Location2[I]<-1
  
  for(i in 1:length(date))
  {
    I=which(fixes$Date==date[i] & fixes$Location==1)
    if(length(I)>0) dayECL[i,t,1]<-1
    if(length(I)>0) dayECL[(i-1):(i+1),t,2]<-1
    I=which(fixes$Date==date[i] & fixes$Location2==1)
    if(length(I)>0) dayECL[i,t,3]<-1
    if(length(I)>0) dayECL[(i-1):(i+1),t,4]<-1
  }
}

totrain=matrix(0,5,4)
rownames(totrain)=c("All","Loc1","Loc1 +1","Loc2","Loc2+1")
colnames(totrain)=c("Control","NoTopo","BRAN","NoEAC")
totrain[1,]=apply(dayrain[,1,],2,sum)
for(i in 1:4) totrain[i+1,]=apply(dayECL[,,i]*dayrain[,1,],2,sum)

totrain[1,]=apply(dayrain[,4,]>0,2,sum)
for(i in 1:4) totrain[i+1,]=apply(dayECL[,,i]*(dayrain[,4,]>0),2,sum)