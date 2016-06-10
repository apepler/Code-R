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

dir=c(36,34,34,34,34,36,36,36)
name1=c("default","default_SST-2","default_SST-1","default_SST+1","default_SST+2","BRAN","BRAN_noeac","BRAN_coarseSST")
name2=c("default","SST-2","SST-1","SST+1","SST+2","BRAN","BRAN_noeac","BRAN_coarseSST")

meanR=matrix(0,8,4)
rownames(meanR)=c("default","SST-2","SST-1","SST+1","SST+2","BRAN","BRAN_noeac","BRAN_coarseSST")
colnames(meanR)=c("Monthly rain","Wettest monthly total","Average hourly extreme","Wettest hour")
for(n in 1:8)
{
  dir1=paste("/srv/ccrc/data",dir[n],"/z3478332/WRF/output/ERAI_R2_nudging_",name1[n],"/out/",sep="")
  fname="wrfhrly_d01_2007-06-01_00:00:00"
  
  a=open.nc(paste(dir1,fname,sep=""))
  latW=var.get.nc(a,"XLAT",c(1,1,1),c(215,144,1))
  lonW=var.get.nc(a,"XLONG",c(1,1,1),c(215,144,1))
  rain1=apply(var.get.nc(a,"PREC_ACC_C"),c(1,2),sum)+apply(var.get.nc(a,"PREC_ACC_NC"),c(1,2),sum)  
  rain2=apply(var.get.nc(a,"PREC_ACC_C"),c(1,2),sum)+apply(var.get.nc(a,"PREC_ACC_NC"),c(1,2),max)  
  lont=as.vector(lonW)
  latt=as.vector(latW)
  rr=as.vector(rain1)
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  meanR[n,1]=mean(rr2$z*t(Useful$mask),na.rm=T)
  meanR[n,2]=max(rr2$z*t(Useful$mask),na.rm=T)
  
  rr=as.vector(rain2)
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  meanR[n,3]=mean(rr2$z*t(Useful$mask),na.rm=T)
  meanR[n,4]=max(rr2$z*t(Useful$mask),na.rm=T)
}

meanR3=matrix(0,30,8)
rownames(meanR)=c("default","SST-2","SST-1","SST+1","SST+2","BRAN","BRAN_noeac","BRAN_coarseSST")
colnames(meanR)=c("Monthly rain","Wettest monthly total","Average hourly extreme","Wettest hour")
meanR3b<-meanR3
for(n in 1:8)
{
  dir1=paste("/srv/ccrc/data",dir[n],"/z3478332/WRF/output/ERAI_R2_nudging_",name1[n],"/out/",sep="")
  fname="wrfhrly_d01_2007-06-01_00:00:00"
  a=open.nc(paste(dir1,fname,sep=""))
  time=var.get.nc(a,"Times")
  day=as.numeric(substr(time,9,10))
  latW=var.get.nc(a,"XLAT",c(1,1,1),c(215,144,1))
  lonW=var.get.nc(a,"XLONG",c(1,1,1),c(215,144,1))
  rain1=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")  
  lont1=as.vector(lonW)
  latt1=as.vector(latW)
  
  fname="wrfhrly_d02_2007-06-01_00:00:00"
  a=open.nc(paste(dir1,fname,sep=""))
  latW=var.get.nc(a,"XLAT",c(1,1,1),c(325,200,1))
  lonW=var.get.nc(a,"XLONG",c(1,1,1),c(325,200,1))
  rain2=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")  
  lont2=as.vector(lonW)
  latt2=as.vector(latW)
  
  for(t in 1:30)
  {
    I=which(day==t)
    rr=as.vector(apply(rain1[,,I],c(1,2),sum))
    rr[which(is.na(rr))]=0
    rr2=interp(lont1,latt1,rr,Useful$x,Useful$y)
    meanR3[t,n]=mean(rr2$z*t(Useful$mask),na.rm=T)
    
    rr=as.vector(apply(rain2[,,I],c(1,2),sum))
    rr[which(is.na(rr))]=0
    rr2=interp(lont2,latt2,rr,Useful$x,Useful$y)
    meanR3b[t,n]=mean(rr2$z*t(Useful$mask),na.rm=T)
  }
  
}

#setwd('/home/nfs/z3478332/Documents/ECLs/')
#save(meanR,meanR2,meanR3,meanR3b,file="sst_test.RData")

# Lisa - load in ACTUAL AWAP rain
# Remember lag - day rain is ~ equiv to AWAP for following day

rainA<-array(0,dim=c(691,886,30))
meanA<-rep(0,30)
month=c(rep(6,29),7)
day=c(2:30,1)
Adir="/media/Seagate Expansion Drive/Data/daily rainfall/rainfall_2000-2009/rainfall-2007/"

for(i in 1:30)
{
  fname<-paste(Adir,"r2007",sprintf("%2.2i",month[i]),sprintf("%2.2i",day[i]),'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->rain
  as.matrix(rain)->rain
  rain[rain<0]=0
  rainA[,,i]<-rain[nrow(rain):1,]
  meanA[i]=mean(rainA[,,i]*Useful$mask,na.rm=T)
}

meanW1<-array(0,c(215,144,30,8))
meanW2<-array(0,c(325,200,30,8))

for(n in 1:8)
{
  dir1=paste("/srv/ccrc/data",dir[n],"/z3478332/WRF/output/ERAI_R2_nudging_",name1[n],"/out/",sep="")
  fname="wrfhrly_d01_2007-06-01_00:00:00"
  a=open.nc(paste(dir1,fname,sep=""))
  time=var.get.nc(a,"Times")
  day=as.numeric(substr(time,9,10))
  latW=var.get.nc(a,"XLAT",c(1,1,1),c(215,144,1))
  lonW=var.get.nc(a,"XLONG",c(1,1,1),c(215,144,1))
  rain1=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")  
  lont1=as.vector(lonW)
  latt1=as.vector(latW)
  
  fname="wrfhrly_d02_2007-06-01_00:00:00"
  a=open.nc(paste(dir1,fname,sep=""))
  latW=var.get.nc(a,"XLAT",c(1,1,1),c(325,200,1))
  lonW=var.get.nc(a,"XLONG",c(1,1,1),c(325,200,1))
  rain2=var.get.nc(a,"PREC_ACC_C")+var.get.nc(a,"PREC_ACC_NC")  
  lont2=as.vector(lonW)
  latt2=as.vector(latW)
  
  for(t in 1:30)
  {
    I=which(day==t)
    meanW1[,,t,n]=apply(rain1[,,I],c(1,2),sum)
    meanW2[,,t,n]=apply(rain2[,,I],c(1,2),sum)
  }
  
}

save(meanW1,meanW2,meanA,file="sst_test_grid.RData")

##### How to convert AWAP to WRF grid instead of vv?
n=1
dir1=paste("/srv/ccrc/data",dir[n],"/z3478332/WRF/output/ERAI_R2_nudging_",name1[n],"/out/",sep="")
fname="wrfhrly_d01_2007-06-01_00:00:00"
a=open.nc(paste(dir1,fname,sep=""))
latW1=var.get.nc(a,"XLAT",c(1,1,1),c(215,144,1))
lonW1=var.get.nc(a,"XLONG",c(1,1,1),c(215,144,1))
fname="wrfhrly_d02_2007-06-01_00:00:00"
a=open.nc(paste(dir1,fname,sep=""))
latW2=var.get.nc(a,"XLAT",c(1,1,1),c(325,200,1))
lonW2=var.get.nc(a,"XLONG",c(1,1,1),c(325,200,1))

rainA1<-array(NaN,c(215,144,30))
rainA2<-array(NaN,c(325,200,30))

for(i in 2:214)
  for(j in 2:143)
  {
    I=which(Useful$x>=mean(lonW1[(i-1):i,j]) & Useful$x<=mean(lonW1[i:(i+1),j]))
    J=which(Useful$y>=mean(latW1[i,(j-1):j]) & Useful$y<=mean(latW1[i,j:(j+1)]))  
    if(length(I)>1 & length(J)>1) rainA1[i,j,]=apply(rainA[J,I,],3,mean,na.rm=T)
  }

for(i in 2:324)
  for(j in 2:199)
  {
    I=which(Useful$x>=mean(lonW2[(i-1):i,j]) & Useful$x<=mean(lonW2[i:(i+1),j]))
    J=which(Useful$y>=mean(latW2[i,(j-1):j]) & Useful$y<=mean(latW2[i,j:(j+1)]))  
    if(length(I)>1 & length(J)>1) rainA2[i,j,]=apply(rainA[J,I,],3,mean,na.rm=T)
  }

save(rainW1,rainW2,rainA,rainA1,rainA2,latW1,lonW1,latW2,lonW2,Useful,file="sst_test_grid.RData")

##########

setwd('/home/nfs/z3478332/Documents/ECLs/')
load("sst_test_grid.RData")
rms1<-array(NaN,c(215,144,8))
rms2<-array(NaN,c(325,200,8))

for(i in 1:215)
  for(j in 1:144)
    for(k in 1:8)
      rms1[i,j,k]=sqrt( mean( (rainW1[i,j,,k]-rainA1[i,j,])^2 , na.rm = TRUE ) )

for(i in 1:325)
  for(j in 1:200)
    for(k in 1:8)
      rms2[i,j,k]=sqrt( mean( (rainW2[i,j,,k]-rainA2[i,j,])^2 , na.rm = TRUE ) )

dayR1<-dayR2<-matrix(0,30,8)

mask1<-array(NaN,c(215,144))
mask2<-array(NaN,c(325,200))

for(i in 2:214)
  for(j in 2:143)
  {
    I=which(Useful$x>=mean(lonW1[(i-1):i,j]) & Useful$x<=mean(lonW1[i:(i+1),j]))
    J=which(Useful$y>=mean(latW1[i,(j-1):j]) & Useful$y<=mean(latW1[i,j:(j+1)]))  
    if(length(I)>1 & length(J)>1) mask1[i,j]=mean(mask[I,J])
  }

for(i in 2:324)
  for(j in 2:199)
  {
    I=which(Useful$x>=mean(lonW2[(i-1):i,j]) & Useful$x<=mean(lonW2[i:(i+1),j]))
    J=which(Useful$y>=mean(latW2[i,(j-1):j]) & Useful$y<=mean(latW2[i,j:(j+1)]))  
    if(length(I)>1 & length(J)>1) mask2[i,j]=mean(mask[I,J])
  }

mask1[mask1>=0.5]=1
mask1[mask1<0.5]=NaN
mask2[mask2>=0.5]=1
mask2[mask2<0.5]=NaN


for(i in 1:30)
  for(j in 1:8)
  {
    data=(rainW1[,,i,j]-rainA1[,,i])^2
    dayR1[i,j]=sqrt( mean( as.vector(data*mask1) , na.rm = TRUE ) )
    data=(rainW2[,,i,j]-rainA2[,,i])^2
    dayR2[i,j]=sqrt( mean( as.vector(data*mask2) , na.rm = TRUE ) )
  }

plot(1:30,dayR2[,1],type="l",col=1,ylim=c(0,60),lwd=2)
for(i in 2:4) lines(1:30,dayR2[,i+4],col=i,lwd=2)
legend("topleft",legend=c("default","BRAN","BRAN_noeac","BRAN_coarseSST"),ncol=2,col=1:4,lwd=2)


############## 
## Comparing ECL stats
## Default = p100, rad 2, cv1

rm(list=ls())
setwd('/home/nfs/z3478332/Documents/ECLs/')
dir="/home/nfs/z3478332/output/outputUM_wrf_sst/"
name1=c("R2","R2_-2","R2_-1","R2_+1","R2_+2","R2_BRAN","R2_noeac","R2_coarseSST")
name2=c("default","SST-2","SST-1","SST+1","SST+2","BRAN","BRAN_noeac","BRAN_coarseSST")


for(p in c("p100","p240"))
  for(r in c("rad2","rad5"))
  {
    events1<-events2<-fixes1<-fixes2<-list()
    
    for(n in 1:8)
    {
      events1[[n]]=read.csv(paste(dir,"ECLevents_d01__",name1[n],"_",r,"_",p,".csv",sep=""))
      events2[[n]]=read.csv(paste(dir,"ECLevents_d02__",name1[n],"_",r,"_",p,".csv",sep=""))
      fixes1[[n]]=read.csv(paste(dir,"ECLfixes_d01__",name1[n],"_",r,"_",p,".csv",sep=""))
      fixes2[[n]]=read.csv(paste(dir,"ECLfixes_d02__",name1[n],"_",r,"_",p,".csv",sep=""))
      
      fixes1[[n]]$Date2=fixes1[[n]]$Date+(as.numeric(fixes1[[n]]$Time)-1)/4
      fixes2[[n]]$Date2=fixes2[[n]]$Date+(as.numeric(fixes2[[n]]$Time)-1)/4
    }
    
    
    dates1<-dates2<-cbind(seq(20070601,20070630.75,0.25),matrix(NaN,120,9))
    ##Needs a control, like ERAI w/ same criteria.
    erai=read.csv(paste(dir,"ECLfixes_erai_",r,"_",p,".csv",sep=""))
    erai$Date2=erai$Date+(as.numeric(erai$Time)-1)/4
    
    for(i in 1:120)
    {
      I=which(erai$Date2==dates1[i,1] & erai$Location==1)
      if(length(I)>0) dates1[i,2]<-dates2[i,2]<-max(erai$CV[I])
      
      for(n in 1:8)
      {
        I=which(fixes1[[n]]$Date2==dates1[i,1] & fixes1[[n]]$Location==1)
        if(length(I)>0) dates1[i,n+2]=max(fixes1[[n]]$CV[I])
        I=which(fixes2[[n]]$Date2==dates2[i,1] & fixes2[[n]]$Location==1)
        if(length(I)>0) dates2[i,n+2]=max(fixes2[[n]]$CV[I])
      }
    }
    
    ymax=ceiling(max(max(dates1[,2:9],na.rm=T),max(dates2[,2:9],na.rm=T)))
    
    tiff(file=paste("SSTcomp","_",p,"_",r,"_eac.tiff",sep=""),width=600,height=400)
    plot(dates1[,1],dates1[,3],type="l",col=1,lwd=2,
         xlab="Date",ylab="Maximum CV",ylim=c(0,ymax),main=paste("June 2007, d01,",p,r))
    for(i in 2:4) lines(dates1[,1],dates1[,6+i],col=i,lwd=2)
    legend("topleft",legend=c("default","BRAN","No EAC","Coarse SST"),
           lwd=2,col=1:4,ncol=2,bty='n')
    dev.off()
    
    tiff(file=paste("SSTcomp","_",p,"_",r,"_add.tiff",sep=""),width=600,height=400)
    plot(dates1[,1],dates1[,3],type="l",col=1,lwd=2,
         xlab="Date",ylab="Maximum CV",ylim=c(0,ymax),main=paste("June 2007, d01,",p,r))
    for(i in 2:5) lines(dates1[,1],dates1[,2+i],col=i,lwd=2)
    legend("topleft",legend=c("default","SST-2","SST-1","SST+1","SST+2"),
           lwd=2,col=1:5,ncol=2,bty='n')
    dev.off()
    
    #     tiff(file=paste("SSTcomp","_",p,"_",r,".tiff",sep=""),width=600,height=400)
    #     plot(dates1[,1],dates1[,2],type="l",col="black",lwd=2,
    #          xlab="Date",ylab="Maximum CV",ylim=c(0,ymax),main=paste("June 2007",p,r))
    #     lines(dates1[,1],dates1[,3],col="blue",lwd=2)
    #     lines(dates1[,1],dates1[,8],col="red",lwd=2)
    #     lines(dates1[,1],dates2[,3],col="blue",lwd=2,lty=2)
    #     lines(dates1[,1],dates2[,8],col="red",lwd=2,lty=2)
    #     legend("topleft",legend=c("ERAI","d01 default","d01 BRAN","d02 default","d02 BRAN"),
    #            lwd=2,lty=c(1,1,1,2,2),col=c("black","blue","red","blue","red"),ncol=2,bty='n')
    #     dev.off()
  }

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
setwd('/home/nfs/z3478332/Documents/ECLs/')
dir="/home/nfs/z3478332/output/outputUM_wrf_sst/"
mask<-t(Useful$mask)
mask[is.na(mask)]=0
f1<-open.nc("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_2005-01_2010-12.nc")
lonP=var.get.nc(f1,'longitude')
latP=var.get.nc(f1,'latitude')
time=var.get.nc(f1,'time')
hh=time%%24
time=as.Date(time/24,origin="1900-01-01")
date=as.numeric(substr(time,1,4))*10000+as.numeric(substr(time,6,7))*100+
  as.numeric(substr(time,9,10))

name1=c("R2","R2_BRAN","R2_noeac","R2_coarseSST")
name2=c("default","BRAN","BRAN_noeac","BRAN_coarseSST")
keydates=c(20070626,20070627)

for(p in c("p100","p240"))
  for(r in c("rad2","rad5"))
  {
    fixes<-list()
    fixes[[1]]=read.csv(paste(dir,"ECLfixes_erai_",r,"_",p,".csv",sep=""))
    fixes[[1]]$Date2=fixes[[1]]$Date+(as.numeric(fixes[[1]]$Time)-1)/4
    for(n in 1:4)
    {
      fixes[[n+1]]=read.csv(paste(dir,"ECLfixes_d01__",name1[n],"_",r,"_",p,".csv",sep=""))    
      fixes[[n+1]]$Date2=fixes[[n+1]]$Date+(as.numeric(fixes[[n+1]]$Time)-1)/4
    }
    
    for(t in 1:2)
    {
      pdf(file=paste("ECL_track_SST_",keydates[t],"_d01_",p,"_",r,".pdf",sep=""),width=5,height=5,pointsize=12)
      par(mar=c(4,4,2,2))
      I=which(date==keydates[t] & hh==0)
      MSLP=var.get.nc(f1,"msl",c(1,1,I),c(length(lonP),length(latP),1),unpack=T)/100
      #pdf(file=paste("Fig8_forposter_",keydates[t],".pdf",sep=""),width=9,height=10,pointsize=36)
      filled.contour3(Useful$x,Useful$y,mask,col=c("white","lightgreen"),xlim=c(145,170),ylim=c(-45,-20))
      contour(lonP,latP[length(latP):1],MSLP[,length(latP):1],levels=seq(992,1050,4),add=T,lwd=2,labcex=1) #19980808 = 29 (1 + 7*4)
      
      dates=data.frame(Date=seq(keydates[t]-3,keydates[t]+3,0.25),
                       lat_erai=rep(NaN,25),lon_erai=rep(NaN,25),
                       lat_default=rep(NaN,25),lon_default=rep(NaN,25),
                       lat_BRAN=rep(NaN,25),lon_BRAN=rep(NaN,25),
                       lat_noeac=rep(NaN,25),lon_noeac=rep(NaN,25),
                       lat_coarseSST=rep(NaN,25),lon_coarseSST=rep(NaN,25))
      
      for(i in 1:5)
      {
        b=fixes[[i]]
        I=which(b$Date==keydates[t] & b$Time=="00:00")
        a=unique(b$ID[I])
        for(j in 1:length(a)) 
        {
          c=b[b$ID==a[j],]
          lines(c$Lon,c$Lat,col=i+1,lwd=3)
        }
      }
      
      legend("topright",legend=c("ERAI",name2),
             col=2:6,lwd=3) 
      dev.off()
    }
  }
    

