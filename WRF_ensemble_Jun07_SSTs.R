rm(list=ls())
library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful2$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

setwd("~/Documents/ECLs/WRFruns/SSTensemble")


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
SSTtype=c("ERAI","BRAN","NoEAC","DoubleEAC")
count<-matrix(0,15,4)

events<-fixes<-list()

c=1
n=1
for(r in 1:3)
  for(day in 27:31)
    for(hour in c("00"))
    {
      dir=c(paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_noeac/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_2eac/",cat[c],"/",sep=""))
      
      
      for(k in 1:4)
      {
        count[n,k]=length(read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep=""))[,1])
        if(n==1)
        {
          
          events[[k]]=read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep=""))
          fixes[[k]]=read.csv(paste(dir[k],"ECLfixes_200705",day,hour,".csv",sep=""))
        } else {
          
          events[[k]]=rbind(events[[k]],read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep="")))
          fixes[[k]]=rbind(fixes[[k]],read.csv(paste(dir[k],"ECLfixes_200705",day,hour,".csv",sep="")))
        }
        
      }
      n=n+1    
    }

ks.test(events[[2]]$CV2,events[[3]]$CV2)
ks.test(fixes[[2]]$CV[fixes[[2]]$Location==1],fixes[[3]]$CV[fixes[[3]]$Location==1])

a=density(fixes[[1]]$CV[fixes[[1]]$Location==1],from=0,to=4)
FixDens=cbind(a$x,a$y,matrix(0,length(a$x),3))
for(i in 2:4)
{
  a=density(fixes[[i]]$CV[fixes[[i]]$Location==1],from=0,to=4)
  FixDens[,i+1]=a$y
}

pdf(file=paste("FixCV_Jun07_all_",cat[c],"_SST.pdf",sep=""),width=6,height=4)
clist=c("black","blue","deepskyblue","blue4")
plot(FixDens[,1],FixDens[,2],col=clist[1],xlim=c(0,4),ylim=c(0,1),type="l",lwd=3,
     xlab="Intensity",ylab="Frequency")
for(i in 2:4) lines(FixDens[,1],FixDens[,i+1],col=clist[i],lwd=3)
lines(FixDens[,1],FixDens[,2],col=clist[1],lwd=3)
legend("topright",legend=SSTtype,col=clist,lwd=3,bty='n')
dev.off()

a=density(events[[1]]$CV2,from=0,to=4)
EvDens=cbind(a$x,a$y,matrix(0,length(a$x),3))
for(i in 2:4)
{
  a=density(events[[i]]$CV2,from=0,to=4)
  EvDens[,i+1]=a$y
}

pdf(file=paste("EventCV_Jun07_all_",cat[c],"_SST.pdf",sep=""),width=6,height=4)
clist=c("black","blue","deepskyblue","blue4")
plot(FixDens[,1],FixDens[,2],col=clist[1],xlim=c(0,4),ylim=c(0,1),type="l",lwd=3,
     xlab="Intensity",ylab="Frequency")
for(i in 2:4) lines(FixDens[,1],FixDens[,i+1],col=clist[i],lwd=3)
lines(FixDens[,1],FixDens[,2],col=clist[1],lwd=3)
legend("topright",legend=SSTtype,col=clist,lwd=3,bty='n')
dev.off()

pval=matrix(0,4,4)
for(i in 1:4)
  for(j in 1:4)
  {
    a=(ks.test(fixes[[i]]$CV[fixes[[i]]$Location==1],fixes[[j]]$CV[fixes[[j]]$Location==1]))
    pval[i,j]=a$p.value
  }


for(k in 1:4)  fixes[[k]]$Date2=fixes[[k]]$Date+(as.numeric(fixes[[k]]$Time)-1)/4
dates=data.frame(Date=seq(20070601,20070630.75,0.25),
                 ERAI=rep(NaN,120),BRAN=rep(NaN,120),BRAN_noeac=rep(NaN,120),BRAN_2eac=rep(NaN,120))
for(i in 1:120)
  for(k in 1:4)
  {
    I=which(fixes[[k]]$Date2==dates$Date[i] & fixes[[k]]$Location==1)
    if(length(I)>0) dates[i,k+1]=mean(fixes[[k]]$CV[I])
  } 

pdf(file=paste("DailyCV_Jun07_all_",cat[c],".pdf",sep=""),width=7,height=4)
plot(dates$Date,dates[,2],type="l",lwd=3,col=clist[1],
     ylim=c(0,ceiling(max(dates[,2:5],na.rm=T))),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
for(k in 2:4) lines(dates$Date,dates[,k+1],lwd=3,col=clist[k])
legend("topleft",legend=SSTtype,col=clist,lwd=3,bty='n')
dev.off()

###### Splitting into different WRF versions

pdfE<-pdfF<-array(0,c(512,3,4))
pval<-array(0,c(4,4,3,2))
c=1
for(r in 1:3)
{
  events<-fixes<-list()
  n=1
  for(day in 27:31)
    for(hour in c("00"))
    {
      dir=c(paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_noeac/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_2eac/",cat[c],"/",sep=""))
      
      
      for(k in 1:4)
      {
        count[n,k]=length(read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep=""))[,1])
        if(n==1)
        {
          
          events[[k]]=read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep=""))
          fixes[[k]]=read.csv(paste(dir[k],"ECLfixes_200705",day,hour,".csv",sep=""))
        } else {
          
          events[[k]]=rbind(events[[k]],read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep="")))
          fixes[[k]]=rbind(fixes[[k]],read.csv(paste(dir[k],"ECLfixes_200705",day,hour,".csv",sep="")))
        }          
      }
      n=n+1    
    }    
  
  for(k in 1:4)
  {
    a=density(events[[k]]$CV2,na.rm=T,from=0,to=4)
    pdfE[,r,k]=a$y
    a=density(fixes[[k]]$CV[fixes[[k]]$Location==1],from=0,to=4)
    pdfF[,r,k]=a$y
    
    for(i in 1:4)
      for(j in 1:4)
      {
        a=(ks.test(fixes[[i]]$CV[fixes[[i]]$Location==1],fixes[[j]]$CV[fixes[[j]]$Location==1]))
        pval[i,j,r,1]=a$p.value
        a=(ks.test(events[[i]]$CV2,events[[j]]$CV2))
        pval[i,j,r,2]=a$p.value
      }
  }
}

pdfE2=cbind(EvDens[,1],EvDens[,5]-EvDens[,4],pdfE[,,4]-pdfE[,,3])
pdfF2=cbind(FixDens[,1],FixDens[,5]-FixDens[,4],pdfF[,,4]-pdfF[,,3])

pdf(file=paste("EventCV_Jun07_PDFchange_",cat[c],"_vwrf_2eac_noeac.pdf",sep=""),width=6,height=4)
plot(pdfE2[,1],pdfE2[,2],lwd=4,col="black",type="l",
     xlim=c(0,4),ylim=c(-0.2,0.2),xlab="Intensity",ylab="Change in PDF")
for(i in 2:4) lines(pdfE2[,1],pdfE2[,i+1],col=i,lwd=2)
abline(h=0,col="gray",lwd=2,lty=2)
legend("topleft",legend=c("All","R1","R2","R3"),col=1:4,lwd=c(4,2,2,2),bty="n",ncol=2)
dev.off()

pdf(file=paste("FixCV_Jun07_PDFchange_",cat[c],"_vwrf_2eac_noeac.pdf",sep=""),width=6,height=4)
plot(pdfF2[,1],pdfF2[,2],lwd=4,col="black",type="l",
     xlim=c(0,4),ylim=c(-0.3,0.3),xlab="Intensity",ylab="Change in PDF")
for(i in 2:4) lines(pdfF2[,1],pdfF2[,i+1],col=i,lwd=2)
abline(h=0,col="gray",lwd=2,lty=2)
legend("bottomright",legend=c("All","R1","R2","R3"),col=1:4,lwd=c(4,2,2,2),bty="n",ncol=2)
dev.off()


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

cc=1
n=1
fixes<-fixesB<-fixes_NE<-fixes_2E<-list()
for(r in 1:3)
  for(day in 27:31)
    for(hour in c("00"))
    {
      dir=c(paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[cc],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN/",cat[cc],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_noeac/",cat[cc],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_2eac/",cat[cc],"/",sep=""))
      
      fixes[[n]]=read.csv(paste(dir[1],"ECLfixes_200705",day,hour,".csv",sep=""))
      fixes[[n]]$Date2=fixes[[n]]$Date+(as.numeric(fixes[[n]]$Time)-1)/4
      fixesB[[n]]=read.csv(paste(dir[2],"ECLfixes_200705",day,hour,".csv",sep=""))
      fixesB[[n]]$Date2=fixesB[[n]]$Date+(as.numeric(fixesB[[n]]$Time)-1)/4
      fixes_NE[[n]]=read.csv(paste(dir[3],"ECLfixes_200705",day,hour,".csv",sep=""))
      fixes_NE[[n]]$Date2=fixes_NE[[n]]$Date+(as.numeric(fixes_NE[[n]]$Time)-1)/4
      fixes_2E[[n]]=read.csv(paste(dir[4],"ECLfixes_200705",day,hour,".csv",sep=""))
      fixes_2E[[n]]$Date2=fixes_2E[[n]]$Date+(as.numeric(fixes_2E[[n]]$Time)-1)/4
      
      n=n+1
    }

for(t in 1:5)
{
  pdf(file=paste("ECL_meantrack_",keydates[t],"_",cat[cc],"_controlBRAN.pdf",sep=""),width=5,height=5,pointsize=12)
  par(mar=c(4,4,2,2))
  I=which(date==keydates[t] & hh==0)
  MSLP=var.get.nc(f1,"msl",c(1,1,I),c(length(lonP),length(latP),1),unpack=T)/100
  #pdf(file=paste("Fig8_forposter_",keydates[t],".pdf",sep=""),width=9,height=10,pointsize=36)
  filled.contour3(Useful$x,Useful$y,mask,col=c("white","lightgreen"),xlim=c(145,170),ylim=c(-45,-20))
  contour(lonP,latP[length(latP):1],MSLP[,length(latP):1],levels=seq(992,1050,4),add=T,lwd=2,labcex=1) #19980808 = 29 (1 + 7*4)
  
  for(i in 1:15)
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
    
    b=fixesB[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,col=rgb(1,0,0,1/4))
      if(i==1 & j==1) fixesshortB=c else fixesshortB=rbind(fixesshortB,c)
    }
  }
  
  fixesshort$Date2=fixesshort$Date+(as.numeric(fixesshort$Time)-1)/4
  fixesshortB$Date2=fixesshortB$Date+(as.numeric(fixesshortB$Time)-1)/4
  dates=data.frame(Date=seq(keydates[t]-3,keydates[t]+3,0.25),count=rep(NaN,25),
                   lat=rep(NaN,25),lon=rep(NaN,25),countB=rep(NaN,25),latB=rep(NaN,25),lonB=rep(NaN,25))
  
  for(i in 1:length(dates[,1]))
  {
    I=which(fixesshort$Date2==dates$Date[i])
    dates$count[i]=length(I)
    if(length(I)>10)
    {
      dates$lat[i]=mean(fixesshort$Lat[I])
      dates$lon[i]=mean(fixesshort$Lon[I])
    }
    
    I=which(fixesshortB$Date2==dates$Date[i])
    dates$countB[i]=length(I)
    if(length(I)>10)
    {
      dates$latB[i]=mean(fixesshortB$Lat[I])
      dates$lonB[i]=mean(fixesshortB$Lon[I])
    }
  }
  
  lines(dates$lon,dates$lat,col="blue",lwd=3)
  lines(dates$lonB,dates$latB,col="red",lwd=3)
  I=which(dates$Date==keydates[t])
  points(dates$lon[I],dates$lat[I],col="blue",lwd=3,cex=2,pch=4)
  points(dates$lonB[I],dates$latB[I],col="red",lwd=3,cex=2,pch=4)
  legend("topright",legend=c("Control","BRAN"),
         col=c("blue","red"),lwd=3) 
  dev.off()
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
  
  
  rm(fix2,fix2N,fix22)
  for(i in 1:15)
  {
    b=fixesB[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      #lines(c$Lon,c$Lat,col="blue")
      if(i==1 & j==1) fix2=c else fix2=rbind(fix2,c)
    }
    
    b=fixes_NE[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      #lines(c$Lon,c$Lat,col="deepskyblue")
      if(i==1 & j==1) fix2N=c else fix2N=rbind(fix2N,c)
    }
    
    b=fixes_2E[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      #lines(c$Lon,c$Lat,col="blue4")
      if(i==1 & j==1) fix22=c else fix22=rbind(fix22,c)
    }
    
  }
  
  dates=data.frame(Date=seq(keydates[t]-3,keydates[t]+3,0.25),
                   count=rep(NaN,25),lat=rep(NaN,25),lon=rep(NaN,25),
                   count_NE=rep(NaN,25),lat_NE=rep(NaN,25),lon_NE=rep(NaN,25),
                   count_2E=rep(NaN,25),lat_2E=rep(NaN,25),lon_2E=rep(NaN,25))
  
  for(i in 1:length(dates[,1]))
  {
    I=which(fix2$Date2==dates$Date[i])
    dates$count[i]=length(I)
    if(length(I)>10)
    {
      dates$lat[i]=mean(fix2$Lat[I])
      dates$lon[i]=mean(fix2$Lon[I])
    }
    
    I=which(fix2N$Date2==dates$Date[i])
    dates$count_NE[i]=length(I)
    if(length(I)>10)
    {
      dates$lat_NE[i]=mean(fix2N$Lat[I])
      dates$lon_NE[i]=mean(fix2N$Lon[I])
    }
    
    I=which(fix22$Date2==dates$Date[i])
    dates$count_2E[i]=length(I)
    if(length(I)>10)
    {
      dates$lat_2E[i]=mean(fix22$Lat[I])
      dates$lon_2E[i]=mean(fix22$Lon[I])
    }}
  
  lines(dates$lon,dates$lat,col="blue",lwd=3)
  lines(dates$lon_NE,dates$lat_NE,col="deepskyblue",lwd=3)
  lines(dates$lon_2E,dates$lat_2E,col="blue4",lwd=3)
  I=which(dates$Date==keydates[t])
  points(dates$lon[I],dates$lat[I],col="blue",lwd=3,cex=2,pch=4)
  points(dates$lon_NE[I],dates$lat_NE[I],col="deepskyblue",lwd=3,cex=2,pch=4)
  points(dates$lon_2E[I],dates$lat_2E[I],col="blue4",lwd=3,cex=2,pch=4)
  legend("topright",legend=c("BRAN","NoEac","DoubleEAC"),
         col=c("blue","deepskyblue","blue4"),lwd=3) 
  dev.off()
}



######## Mean loactions

lat=seq(-40,-25,2.5)
lon=seq(145,160,2.5)

cat=c("d01_p100_rad2cv05","d01_p100_rad5cv015","d02_p100_rad2cv05","d02_p100_rad5cv015")
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


c=1

loc<-array(0,c(7,7,4))
for(k in 1:15)
  for(i in 1:7)
    for(j in 1:7)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-1.25 & fixes[[k]]$Lat<lat[i]+1.25 & fixes[[k]]$Lon>=lon[j]-1.25 & fixes[[k]]$Lon<lon[j]+1.25 & fixes[[k]]$Location==1)
      loc[i,j,1]=loc[i,j,1]+length(I)
      I=which(fixesB[[k]]$Lat>=lat[i]-1.25 & fixesB[[k]]$Lat<lat[i]+1.25 & fixesB[[k]]$Lon>=lon[j]-1.25 & fixesB[[k]]$Lon<lon[j]+1.25 & fixesB[[k]]$Location==1)
      loc[i,j,2]=loc[i,j,2]+length(I)
      I=which(fixes_NE[[k]]$Lat>=lat[i]-1.25 & fixes_NE[[k]]$Lat<lat[i]+1.25 & fixes_NE[[k]]$Lon>=lon[j]-1.25 & fixes_NE[[k]]$Lon<lon[j]+1.25 & fixes_NE[[k]]$Location==1)
      loc[i,j,3]=loc[i,j,3]+length(I)
      I=which(fixes_2E[[k]]$Lat>=lat[i]-1.25 & fixes_2E[[k]]$Lat<lat[i]+1.25 & fixes_2E[[k]]$Lon>=lon[j]-1.25 & fixes_2E[[k]]$Lon<lon[j]+1.25 & fixes_2E[[k]]$Location==1)
      loc[i,j,4]=loc[i,j,4]+length(I)
    }

loc=loc/15

pdf(file=paste("ECL_locations_",cat[c],".pdf",sep=""),width=12,height=4,pointsize=12)
layout(cbind(1,2,3,4),c(1,1,1,0.3))

for(k in c(3,2,4))
{
  image(lon,lat,t(loc[,,k]),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=SSTtype[k],cex.axis=1.5,cex.main=1.5)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb,cols)
dev.off()

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm2=pal(10)
bb2=c(-10000,-50,-25,-10,-5,0,5,10,25,50,10000)
blab=c("-50%","-25%","-10%","-5%","","+5%","+10%","+25%","+50%")
bb3=c(-100,-3,-2,-1,-0.5,0,0.5,1,2,3,100)


pdf(file=paste("ECL_locations_",cat[c],"_change.pdf",sep=""),width=12,height=4,pointsize=12)
layout(cbind(1,2,3,4),c(1,1,1,0.3))
tmp=loc[,,3]-loc[,,2]
tmp[loc[,,3]==0 & loc[,,2]==0]=NaN
image(lon,lat,t(tmp),xlab="",ylab="",breaks=bb3,col=cm2,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="NoEAC",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,t(loc[,,2]),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="BRAN",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
tmp=loc[,,4]-loc[,,2]
tmp[loc[,,4]==0 & loc[,,2]==0]=NaN
image(lon,lat,t(tmp),xlab="",ylab="",breaks=bb3,col=cm2,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="2EAC",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb3,cm2)
dev.off()



#### Divide domain into different regions, and count # ECL days in each region for each version

ECLdays<-array(0,c(15,4,4))
dimnames(ECLdays)[[1]]<-rep("aaa",15)
dimnames(ECLdays)[[2]]<-c("All","S coast","C coast","N coast")
dimnames(ECLdays)[[3]]<-SSTtype

n=1
for(r in 1:3)
  for(day in 27:31)
    for(hour in c("00"))
    {
      dimnames(ECLdays)[[1]][n]<-paste("R",r,": ",day," May",sep="")
      
      dir=c(paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_noeac/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_2eac/",cat[c],"/",sep=""))
      
      for(k in 1:4)
      {
        data=read.csv(paste(dir[k],"ECLfixes_200705",day,hour,".csv",sep=""))    
        I=which(data$Location==1)
        ECLdays[n,1,k]=length(unique(data$Date[I]))
        
        data$Location2=0
        I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-36) & data$Lat>=-40)
        data$Location2[I]<-1
        I=which(data$Location2==1)
        ECLdays[n,2,k]=length(unique(data$Date[I]))
        
        data$Location2=0
        I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-36)
        data$Location2[I]<-1
        I=which(data$Location2==1)
        ECLdays[n,3,k]=length(unique(data$Date[I]))
        
        data$Location2=0
        I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
        data$Location2[I]<-1
        I=which(data$Location2==1)
        ECLdays[n,4,k]=length(unique(data$Date[I]))       
        
      }
      
      n=n+1    
    }

### Daily ECLs
c=1
daylist=20070601:20070630

ECLdays<-array(NaN,c(30,15,4))
dimnames(ECLdays)[[1]]<-daylist
dimnames(ECLdays)[[2]]<-rep("aaa",15)
dimnames(ECLdays)[[3]]<-SSTtype


n=1
for(r in 1:3)
  for(day in 27:31)
    for(hour in c("00"))
    {
      dimnames(ECLdays)[[2]][n]<-paste("R",r,": ",day," May",sep="")
      
      dir=c(paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_noeac/",cat[c],"/",sep=""),
            paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble_BRAN_2eac/",cat[c],"/",sep=""))
      
      for(k in 1:4)
      {
        data=read.csv(paste(dir[k],"ECLfixes_200705",day,hour,".csv",sep="")) 
        for(i in 1:30)
        {
          I=which(data$Date==daylist[i] & data$Location==1)
          if(length(I)>0) ECLdays[i,n,k]=max(data$CV[I])
        }
      }
      
      n=n+1    
    }

ECLdays2=!is.na(ECLdays)
ECLdays3=(ECLdays2==1 & ECLdays>=1)
tmp=apply(ECLdays3,c(1,3),mean)*100

SSTtype=c("Control","BRAN","NoEAC","2EAC")
plot(daylist,tmp[,1],col="black",type="l",lwd=3,xlab="Date",ylab="% of members with an ECL present",ylim=c(0,115))
clist=c("black","grey","blue","red")
for(i in 2:4) lines(daylist,tmp[,i],lwd=3,col=clist[i])
legend("topleft",SSTtype,lwd=3,col=clist,ncol=4,bty="n")

tmp=apply(ECLdays[,,3]-ECLdays[,,2],1,mean,na.rm=T)

meanCV<-matrix(0,15,4)
for(i in 1:15)
{
  meanCV[i,1]=mean(fixes[[i]]$CV[fixes[[i]]$Location==1])
  meanCV[i,2]=mean(fixesB[[i]]$CV[fixesB[[i]]$Location==1])
  meanCV[i,3]=mean(fixes_NE[[i]]$CV[fixes_NE[[i]]$Location==1])
  meanCV[i,4]=mean(fixes_2E[[i]]$CV[fixes_2E[[i]]$Location==1])
}
