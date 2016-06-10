setwd("~/Documents/ECLs/WRFruns/Ensemble")

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful2$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

###############
##
##  Need a figure with PDFs of monthly rain difference for all 3 versions & both resolutions
##
###############

raindiff<-matrix(NaN,20,6)
colnames(raindiff)=rep("aaa",6)
n=1
for(d in c("D01","D02"))
  for(w in c("R1","R2","R3"))
  {
    a=read.csv(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_",w,"_ensemble_notopo/out/",d,"_totalESBrain.csv",sep=""))
    raindiff[,n]=((a[,3]/a[,2])-1)*100
    colnames(raindiff)[n]=paste(d,w)
    n=n+1
  }

a=density(raindiff[,1:3],na.rm=T,from=-60,to=10)
dens=data.frame(X=a$x,d01=a$y,d02=rep(0,length(a$y)),
                R1=rep(0,length(a$y)),R2=rep(0,length(a$y)),R3=rep(0,length(a$y)),
                R1a=rep(0,length(a$y)),R2a=rep(0,length(a$y)),R3a=rep(0,length(a$y)))
a=density(raindiff[,4:6],na.rm=T,from=-60,to=10)
dens[,3]=a$y
for(i in 1:6){
  a=density(raindiff[,i],na.rm=T,from=-60,to=10)
  dens[,i+3]=a$y
} 
clist=c(1,1,2:4,2:4)
tlist=c(1,2,1,1,1,2,2,2)
wlist=c(4,4,2,2,2,2,2,2)
plot(NA,xlim=c(-60,10),ylim=range(dens[,2:9]),xlab="% change in ESB mean rainfall",
     ylab="Frequency",main="June 2007 ESB rainfall change for no topography",cex.main=1.2)
abline(v=0,col="gray",lwd=2,lty=3)
for(i in 8:1) lines(dens[,1],dens[,i+1],col=clist[i],lty=tlist[i],lwd=wlist[i])
legend("topleft",legend=c("All","R1","R2","R3"),col=1:4,lwd=c(4,2,2,2),lty=1,bty='n')
legend("topright",legend=c("50km","10km"),col=1,lwd=2,lty=1:2,bty='n')

######
##
## Now, something about distribution of intensity across all R versions
##
######

makePDF = function(data1,data2,xlabel,mlab="June 2007 ECL statistics for 20 runs",labloc="topleft") {
  a=density(data1,na.rm=T)
  b=density(data2,na.rm=T)
  
  lims=range(data1,data2,na.rm=T)
  if((lims[2]-lims[1])<10)
  {
    lims[1]=floor(lims[1])
    lims[2]=ceiling(lims[2])
  } else {
    lims[1]=floor(lims[1]/5)*5
    lims[2]=ceiling(lims[2]/5)*5
  }
  
  plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
       xlab=xlabel,ylab="Frequency",cex.main=1.2,
       main=mlab)
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=c("Control","NoTopo"),
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
}


cat=c("d01_p100_rad2cv05","d01_p100_rad5cv015","d02_p100_rad2cv05","d02_p100_rad5cv015")
count<-countNT<-matrix(0,60,4)

for(c in 1:4)
{
  n=1
  for(r in 1:3)
    for(day in 27:31)
      for(hour in c("00","06","12","18"))
      {
        dir1=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep="")
        dir2=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_notopo/",cat[c],"/",sep="")
        
        if(n==1)
        {
          events=read.csv(paste(dir1,"ECLevents_200705",day,hour,".csv",sep=""))
          fixes=read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep=""))
          events_notopo=read.csv(paste(dir2,"ECLevents_200705",day,hour,".csv",sep=""))
          fixes_notopo=read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep=""))
        } else {
          events=rbind(events,read.csv(paste(dir1,"ECLevents_200705",day,hour,".csv",sep="")))
          fixes=rbind(fixes,read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep="")))
          events_notopo=rbind(events_notopo,read.csv(paste(dir2,"ECLevents_200705",day,hour,".csv",sep="")))
          fixes_notopo=rbind(fixes_notopo,read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep="")))
        }
        
        count[n,c]=length(read.csv(paste(dir1,"ECLevents_200705",day,hour,".csv",sep=""))[,1])
        countNT[n,c]=length(read.csv(paste(dir2,"ECLevents_200705",day,hour,".csv",sep=""))[,1])
        n=n+1    
      }
  
  pdf(file=paste("FixCV_Jun07_all_",cat[c],".pdf",sep=""),width=6,height=4)
  makePDF(fixes$CV[fixes$Location==1],fixes_notopo$CV[fixes_notopo$Location==1],xlab="Intensity",mlab="Strength of low - all systems in June 2007","topright")
  dev.off()
  pdf(file=paste("EventCV_Jun07_all_",cat[c],".pdf",sep=""),width=6,height=4)
  makePDF(events$CV2,events_notopo$CV2,xlab="Intensity",mlab="Maximum strength of low - all systems in June 2007","topright")
  dev.off()
  pdf(file=paste("EventMSLP_Jun07_all_",cat[c],".pdf",sep=""),width=6,height=4)
  makePDF(events$MSLP2,events_notopo$MSLP2,xlab="MSLP",mlab="Maximum strength of low - all systems in June 2007","topleft")
  dev.off()
  pdf(file=paste("EventL_Jun07_all_",cat[c],".pdf",sep=""),width=6,height=4)
  makePDF(events$Length2,events_notopo$Length2,xlab="Intensity",mlab="Duration of low in ESB - all systems in June 2007","topright")
  dev.off()
  
  fixes$Date2=fixes$Date+(as.numeric(fixes$Time)-1)/4
  fixes_notopo$Date2=fixes_notopo$Date+(as.numeric(fixes_notopo$Time)-1)/4
  dates=data.frame(Date=seq(20070601,20070630.75,0.25),
                   CV=rep(NaN,120),MSLP=rep(NaN,120),CV_NT=rep(NaN,120),MSLP_NT=rep(NaN,120))
  for(i in 1:120)
  {
    I=which(fixes$Date2==dates$Date[i] & fixes$Location==1)
    if(length(I)>0){
      dates[i,2]=mean(fixes$CV[I])
      dates[i,3]=mean(fixes$MSLP[I])
    } 
    
    I=which(fixes_notopo$Date2==dates$Date[i] & fixes_notopo$Location==1)
    if(length(I)>0){
      dates[i,4]=mean(fixes_notopo$CV[I])
      dates[i,5]=mean(fixes_notopo$MSLP[I])
    } 
  }
  
  pdf(file=paste("DailyCV_Jun07_all_",cat[c],".pdf",sep=""),width=7,height=4)
  plot(dates$Date,dates[,2],type="l",lwd=3,col="blue",
       ylim=c(0,ceiling(max(dates[,c(2,4)],na.rm=T))),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
  lines(dates$Date,dates[,4],lwd=3,col="red")
  legend("topleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3)
  dev.off()

}

for(c in 1)
{
  pdfE<-pdfF<-pdfEnt<-pdfFnt<-data.frame(X=seq(0,5,length.out=512),R1=rep(0,512),R2=rep(0,512),R3=rep(0,512),All=rep(0,512))
  pval<-matrix(0,4,2)
  for(r in 1:3)
  {
    n=1
    for(day in 27:31)
      for(hour in c("00","06","12","18"))
      {
        dir1=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep="")
        dir2=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_notopo/",cat[c],"/",sep="")
        
        if(n==1)
        {
          events=read.csv(paste(dir1,"ECLevents_200705",day,hour,".csv",sep=""))
          events_notopo=read.csv(paste(dir2,"ECLevents_200705",day,hour,".csv",sep=""))
          
          a=read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep=""))
          fixes=a[a$Location==1,]
          a=read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep=""))
          fixes_notopo=a[a$Location==1,]
        } else {
          events=rbind(events,read.csv(paste(dir1,"ECLevents_200705",day,hour,".csv",sep="")))
          events_notopo=rbind(events_notopo,read.csv(paste(dir2,"ECLevents_200705",day,hour,".csv",sep="")))
          a=read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep=""))
          fixes=rbind(fixes,a[a$Location==1,])
          a=read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep=""))
          fixes_notopo=rbind(fixes_notopo,a[a$Location==1,])
        
        
        }
        n=n+1    
      }
  
    a=density(events$CV2,na.rm=T,from=0,to=5)
    pdfE[,r+1]=a$y
    a=density(events_notopo$CV2,na.rm=T,from=0,to=5)
    pdfEnt[,r+1]=a$y
    a=density(fixes$CV,na.rm=T,from=0,to=5)
    pdfF[,r+1]=a$y
    a=density(fixes_notopo$CV,na.rm=T,from=0,to=5)
    pdfFnt[,r+1]=a$y
    
    a=ks.test(events$CV2,events_notopo$CV2)
    pval[r,1]=a$p.value
    a=ks.test(fixes$CV,fixes_notopo$CV)
    pval[r,2]=a$p.value
    
    if(r==1)
    {
      events2=events
      fixes2=fixes
      eventsNT2=events_notopo
      fixesNT2=fixes_notopo
    } else {
      events2=rbind(events2,events)
      eventsNT2=rbind(eventsNT2,events_notopo)
      fixes2=rbind(fixes2,fixes)
      fixesNT2=rbind(fixesNT2,fixes_notopo)
    }
  }
  a=ks.test(events2$CV2,eventsNT2$CV2)
  pval[4,1]=a$p.value
  a=ks.test(fixes2$CV,fixesNT2$CV)
  pval[4,2]=a$p.value
  
  
  a=density(events2$CV2,na.rm=T,from=0,to=5)
  pdfE[,5]=a$y
  a=density(eventsNT2$CV2,na.rm=T,from=0,to=5)
  pdfEnt[,5]=a$y
  a=density(fixes2$CV,na.rm=T,from=0,to=5)
  pdfF[,5]=a$y
  a=density(fixesNT2$CV,na.rm=T,from=0,to=5)
  pdfFnt[,5]=a$y
  
  pdfE2=cbind(pdfE[,1],pdfEnt[,2:5]-pdfE[,2:5])
  pdfF2=cbind(pdfE[,1],pdfFnt[,2:5]-pdfF[,2:5])
  
  
  yy=c(-signif(max(abs(pdfE2[,2:5]))+0.1,digits=1),signif(max(abs(pdfE2[,2:5]))+0.1,digits=1))
  pdf(file=paste("EventCV_Jun07_PDFchange_",cat[c],"_vwrf.pdf",sep=""),width=6,height=4)
  plot(pdfE2[,1],pdfE2[,5],lwd=4,col="black",type="l",
       xlim=c(0,5),ylim=yy,xlab="Intensity",ylab="Change in PDF")
  for(i in 2:4) lines(pdfE2[,1],pdfE2[,i],col=i,lwd=2)
  abline(h=0,col="gray",lwd=2,lty=2)
  legend("topright",legend=c("All","R1","R2","R3"),col=1:4,lwd=c(4,2,2,2),bty="n")
  dev.off()
  
  yy=c(-signif(max(abs(pdfF2[,2:5]))+0.1,digits=1),signif(max(abs(pdfF2[,2:5]))+0.1,digits=1))
  pdf(file=paste("FixCV_Jun07_PDFchange_",cat[c],"_vwrf.pdf",sep=""),width=6,height=4)
  plot(pdfF2[,1],pdfF2[,5],lwd=4,col="black",type="l",
       xlim=c(0,5),ylim=yy,xlab="Intensity",ylab="Change in PDF")
  for(i in 2:4) lines(pdfF2[,1],pdfF2[,i],col=i,lwd=2)
  abline(h=0,col="gray",lwd=2,lty=2)
  legend("topright",legend=c("All","R1","R2","R3"),col=1:4,lwd=c(4,2,2,2),bty="n")
  dev.off()
}


cat=c("d01_p100_rad2cv05","d01_p100_rad5cv015","d02_p100_rad2cv05","d02_p100_rad5cv015")

for(c in 1:4)
{
  n=1
  dates=seq(20070601,20070630.75,0.25)
  CV<-CVnt<-matrix(NaN,length(dates),60)
  for(r in 1:3)
    for(day in 27:31)
      for(hour in c("00","06","12","18"))
      {
        dir1=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep="")
        dir2=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_notopo/",cat[c],"/",sep="")
        
        fixes=read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep=""))
        fixes_notopo=read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep=""))
        fixes$Date2=fixes$Date+(as.numeric(fixes$Time)-1)/4
        fixes_notopo$Date2=fixes_notopo$Date+(as.numeric(fixes_notopo$Time)-1)/4
        
        for(i in 1:120)
        {
          I=which(fixes$Date2==dates[i] & fixes$Location==1)
          if(length(I)>0) CV[i,n]=mean(fixes$CV[I])          
          I=which(fixes_notopo$Date2==dates[i] & fixes_notopo$Location==1)
          if(length(I)>0) CVnt[i,n]=mean(fixes_notopo$CV[I])
        }

        n=n+1    
      }
  
  CV2=cbind(apply(CV[,1:20],1,mean,na.rm=T),apply(CV[,21:40],1,mean,na.rm=T),apply(CV[,41:60],1,mean,na.rm=T))
  CVnt2=cbind(apply(CVnt[,1:20],1,mean,na.rm=T),apply(CVnt[,21:40],1,mean,na.rm=T),apply(CVnt[,41:60],1,mean,na.rm=T))
  yy=range(CV,CVnt,na.rm=T)
    
  pdf(file=paste("DailyCV_scatter_",cat[c],"_vwrf.pdf",sep=""),width=4,height=4)
  plot(0:ceiling(yy[2]),0:ceiling(yy[2]),type="l",col="black",lwd=3,xlab="Default case",ylab="No topography")
  for(i in 1:3) points(CV2[,i],CVnt2[,i],col=i+1,pch=4,lwd=2)
  legend("topleft",legend=c("R1","R2","R3"),col=2:4,pch=4,pt.lwd=2,bty="n")
  dev.off()
     
#   pdf(file=paste("DailyCV_Jun07_all_",cat[c],"_v2.pdf",sep=""),width=7,height=4)
#   plot(NA,xlim=range(dates),ylim=c(0,ceiling(yy[2])),xlab="Date",ylab="Intensity",main="Average intensity - all systems in June 2007")
#   for(i in 1:60) lines(dates,CV[,i],col=rgb(0,0,1,1/4))
#   for(i in 1:60) lines(dates,CVnt[,i],col=rgb(1,0,0,1/4))
#   lines(dates,apply(CV,1,mean,na.rm=T),col="blue",lwd=4)
#   lines(dates,apply(CVnt,1,mean,na.rm=T),col="red",lwd=4)
#   legend("topleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=4,bty="n")
#   dev.off()
#   
#   pdf(file=paste("DailyCV_Jun07_all_",cat[c],"_vWRF.pdf",sep=""),width=7,height=4)
#   plot(NA,xlim=range(dates),ylim=c(0,ceiling(yy[2])),xlab="Date",ylab="Intensity",main="Average intensity - all systems in June 2007")
#   for(i in 1:3) lines(dates,CV2[,i],col=i+1,lwd=4)
#   for(i in 1:3) lines(dates,CVnt2[,i],col=i+1,lty=2,lwd=4)
#   legend("topright",legend=c("Control","NoTopo"),lty=c(1,2),lwd=4,bty="n")
#   legend("topleft",legend=c("R1","R2","R3"),col=2:4,lwd=4,bty="n")
#   dev.off()
#   
#   prop=cbind(apply(1-is.na(CV),1,sum),apply(1-is.na(CVnt),1,sum))*100/60
# 
#   pdf(file=paste("DailyProp_Jun07_all_",cat[c],".pdf",sep=""),width=7,height=4)
#   plot(NA,xlim=range(dates),ylim=c(0,110),xlab="Date",ylab="%",main="% of ensemble members with a low present")
#   lines(dates,prop[,1],col="blue",lwd=4)
#   lines(dates,prop[,2],col="red",lwd=4)
#   legend("topleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=4,bty="n")
#   dev.off()
#   
#   prop2=cbind(apply(1-is.na(CV[,1:20]),1,sum),apply(1-is.na(CV[,21:40]),1,sum),apply(1-is.na(CV[,41:60]),1,sum))*100/20
#   prop2nt=cbind(apply(1-is.na(CVnt[,1:20]),1,sum),apply(1-is.na(CVnt[,21:40]),1,sum),apply(1-is.na(CVnt[,41:60]),1,sum))*100/20
#   
#   pdf(file=paste("DailyProp_Jun07_all_",cat[c],"_vwrf.pdf",sep=""),width=7,height=4)
#   plot(NA,xlim=range(dates),ylim=c(0,130),xlab="Date",ylab="%",main="% of ensemble members with a low present")
#   for(i in 1:3) lines(dates,prop2[,i],col=i+1,lwd=4)
#   for(i in 1:3) lines(dates,prop2nt[,i],col=i+1,lty=2,lwd=4)
#   legend("topright",legend=c("Control","NoTopo"),lty=c(1,2),lwd=4,bty="n")
#   legend("topleft",legend=c("R1","R2","R3"),col=2:4,lwd=4,bty="n")
#   dev.off()
  
}


######### Re-do location stuff for a particular event

## Mean ECL track

##Tracks for ECLs - June 2007
library(RNetCDF)
f1<-open.nc("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_2005-01_2010-12.nc")
lonP=var.get.nc(f1,'longitude')
latP=var.get.nc(f1,'latitude')
time=var.get.nc(f1,'time')
hh=time%%24
time=as.Date(time/24,origin="1900-01-01")
date=as.numeric(substr(time,1,4))*10000+as.numeric(substr(time,6,7))*100+
  as.numeric(substr(time,9,10))

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

cat=c("d01_p100_rad2cv05","d01_p100_rad5cv015","d02_p100_rad2cv05","d02_p100_rad5cv015")
keydates=c(20070608,20070616,20070619,20070626,20070628)

for(cc in 1:4)
{
  n=1
  fixes<-fixes_notopo<-list()
  for(r in 1:3)
    for(day in 27:31)
      for(hour in c("00","06","12","18"))
      {
        dir1=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[cc],"/",sep="")
        dir2=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_notopo/",cat[cc],"/",sep="")
        
        fixes[[n]]=read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep=""))
        fixes_notopo[[n]]=read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep=""))
        fixes[[n]]$Date2=fixes[[n]]$Date+(as.numeric(fixes[[n]]$Time)-1)/4
        fixes_notopo[[n]]$Date2=fixes_notopo[[n]]$Date+(as.numeric(fixes_notopo[[n]]$Time)-1)/4
        n=n+1
      }
  
for(t in 1:5)
{
  pdf(file=paste("ECL_meantrack_",keydates[t],"_",cat[cc],".pdf",sep=""),width=5,height=5,pointsize=12)
  par(mar=c(4,4,2,2))
  I=which(date==keydates[t] & hh==0)
  MSLP=var.get.nc(f1,"msl",c(1,1,I),c(length(lonP),length(latP),1),unpack=T)/100
  #pdf(file=paste("Fig8_forposter_",keydates[t],".pdf",sep=""),width=9,height=10,pointsize=36)
  filled.contour3(Useful$x,Useful$y,mask,col=c("white","lightgreen"),xlim=c(145,170),ylim=c(-45,-20))
  contour(lonP,latP[length(latP):1],MSLP[,length(latP):1],levels=seq(992,1050,4),add=T,lwd=2,labcex=1) #19980808 = 29 (1 + 7*4)
  
  for(i in 1:60)
  {
    b=fixes[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,col=rgb(0,0,1,1/4))
      if(i==1 & j==1) fixesshort=c else fixesshort=rbind(fixesshort,c)
    }
    
    b=fixes_notopo[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,col=rgb(1,0,0,1/4))
      if(i==1 & j==1) fixesshort_NT=c else fixesshort_NT=rbind(fixesshort_NT,c)
    }
  }
  
  fixesshort$Date2=fixesshort$Date+(as.numeric(fixesshort$Time)-1)/4
  fixesshort_NT$Date2=fixesshort_NT$Date+(as.numeric(fixesshort_NT$Time)-1)/4
  dates=data.frame(Date=seq(keydates[t]-3,keydates[t]+3,0.25),count=rep(NaN,25),
                   lat=rep(NaN,25),lon=rep(NaN,25),count_NT=rep(NaN,25),lat_NT=rep(NaN,25),lon_NT=rep(NaN,25))
  
  for(i in 1:length(dates[,1]))
  {
    I=which(fixesshort$Date2==dates$Date[i])
    dates$count[i]=length(I)
    if(length(I)>10)
    {
      dates$lat[i]=mean(fixesshort$Lat[I])
      dates$lon[i]=mean(fixesshort$Lon[I])
    }
    
    I=which(fixesshort_NT$Date2==dates$Date[i])
    dates$count_NT[i]=length(I)
    if(length(I)>10)
    {
      dates$lat_NT[i]=mean(fixesshort_NT$Lat[I])
      dates$lon_NT[i]=mean(fixesshort_NT$Lon[I])
    }
  }
  
  lines(dates$lon,dates$lat,col="blue",lwd=3)
  lines(dates$lon_NT,dates$lat_NT,col="red",lwd=3)
  I=which(dates$Date==keydates[t])
  points(dates$lon[I],dates$lat[I],col="blue",lwd=3,cex=2,pch=4)
  points(dates$lon_NT[I],dates$lat_NT[I],col="red",lwd=3,cex=2,pch=4)
  legend("topright",legend=c("Control","NoTopo"),
         col=c("blue","red"),lwd=3) 
  dev.off()
}
}


####### What about ECL locations?

lat=seq(-40,-25,2.5)
lon=seq(145,160,2.5)

cat=c("d01_p100_rad2cv05","d01_p100_rad5cv015","d02_p100_rad2cv05","d02_p100_rad5cv015")

for(c in 1:4)
{
  n=1
  loc<-locNT<-matrix(0,7,7)
  
  for(r in 1:3)
    for(day in 27:31)
      for(hour in c("00","06","12","18"))
      {
        dir1=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep="")
        dir2=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_notopo/",cat[c],"/",sep="")
        
        if(n==1)
        {
          fixes=read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep=""))
          fixes_notopo=read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep=""))
        } else {
          fixes=rbind(fixes,read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep="")))
          fixes_notopo=rbind(fixes_notopo,read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep="")))
        }
        n=n+1    
      }
  
  for(i in 1:7)
    for(j in 1:7)
    {
      I=which(fixes$Lat>=lat[i]-1.25 & fixes$Lat<lat[i]+1.25 & fixes$Lon>=lon[j]-1.25 & fixes$Lon<lon[j]+1.25)
      loc[i,j]=length(I)/60
      
      I=which(fixes_notopo$Lat>=lat[i]-1.25 & fixes_notopo$Lat<lat[i]+1.25 & fixes_notopo$Lon>=lon[j]-1.25 & fixes_notopo$Lon<lon[j]+1.25)
      locNT[i,j]=length(I)/60
    }
  
  bb=c(-0.5,0,0.5,1,2,3,5,100)
  cols=gray(seq(1,0.1,-0.15))
  ColorBar <- function(brks,cols,vert=T,subsampleg=1)
  {
    par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1.2)
    image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
          xlab = '', ylab = '')
    box()
    axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
         labels = brks[seq(2, length(brks)-1, subsampleg)])
  }
  
  pdf(file=paste("ECL_locations_",cat[c],".pdf",sep=""),width=8,height=4,pointsize=12)
  layout(cbind(1,2,3),c(1,1,0.3))
  image(lons,lats,t(loc),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="Control",cex.axis=1.5,cex.main=1.5)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  image(lons,lats,t(locNT),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="No Topography",cex.axis=1.5,cex.main=1.5)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  ColorBar(bb,cols)
  dev.off()
  
  locdiff=locNT-loc
  bb=c(-100,-2,-1,-0.5,0.5,1,2,100)
  source('~/Documents/R/color.palette.R')
  pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
  cm=pal(7)
  cm[4]="white"
  
  pdf(file=paste("ECL_locations_",cat[c],"_change.pdf",sep=""),width=4.5,height=4,pointsize=12)
  layout(cbind(1,2),c(1,0.3))
  image(lons,lats,t(locdiff),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="Change in ECL frequency")
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  ColorBar(bb,cm)
  dev.off()  
}
  
  
  