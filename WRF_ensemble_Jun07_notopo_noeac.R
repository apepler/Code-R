###### Okay, so compare the ensemble-results for notopo and noeac
###### With the results for doing them both (individual case)
###### Currently just R2

rm(list=ls())
library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->Useful
readMat('~/Documents/GDI/Useful.mat')->Useful2
mask<-t(Useful2$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

setwd("~/Documents/ECLs/WRFruns/")
d1="/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007-06/"

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
dom2=c("d01","d01","d02","d02")
cat2=c("rad2_p100","rad5_p100","rad2_p100","rad5_p100")

SSTtype=c("Control","NoTopo","NoEAC","Both")
count<-matrix(0,5,4)

events<-fixes<-list()

c=1
n=1
for(r in 2)
  for(day in 27:31)
    for(hour in c("00"))
    {
      dir=c(paste(d1,"ERAI_R",r,"_ensemble/",cat[c],"/",sep=""),
            paste(d1,"ERAI_R",r,"_ensemble_notopo/",cat[c],"/",sep=""),
            paste(d1,"ERAI_R",r,"_ensemble_BRAN_noeac/",cat[c],"/",sep=""))

      
      for(k in 1:3)
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

events[[4]]=read.csv(paste(d1,"ECLevents_",dom2[c],"_R",r,"_notopo_noeac_",cat2[c],".csv",sep=""))
fixes[[4]]=read.csv(paste(d1,"ECLfixes_",dom2[c],"_R",r,"_notopo_noeac_",cat2[c],".csv",sep=""))
count[1,4]=length(events[[4]][,1])

a=density(fixes[[1]]$CV[fixes[[1]]$Location==1],from=0,to=4)
FixDens=cbind(a$x,a$y,matrix(0,length(a$x),3))
for(i in 2:4)
{
  a=density(fixes[[i]]$CV[fixes[[i]]$Location==1],from=0,to=4)
  FixDens[,i+1]=a$y
}

clist=c("black","blue","red","purple")
plot(FixDens[,1],FixDens[,2],col=clist[1],xlim=c(0,4),ylim=c(0,1),type="l",lwd=3,
     xlab="Intensity",ylab="Frequency")
for(i in 2:4) lines(FixDens[,1],FixDens[,i+1],col=clist[i],lwd=3)
lines(FixDens[,1],FixDens[,2],col=clist[1],lwd=3)
legend("topright",legend=SSTtype,col=clist,lwd=3,bty='n')


for(k in 1:4)  fixes[[k]]$Date2=fixes[[k]]$Date+(as.numeric(fixes[[k]]$Time)-1)/4
dates=data.frame(Date=seq(20070601,20070630.75,0.25),
                 ERAI=rep(NaN,120),BRAN=rep(NaN,120),BRAN_noeac=rep(NaN,120),BRAN_2eac=rep(NaN,120))
for(i in 1:120)
  for(k in 1:4)
  {
    I=which(fixes[[k]]$Date2==dates$Date[i] & fixes[[k]]$Location==1)
    if(length(I)>0) dates[i,k+1]=mean(fixes[[k]]$CV[I])
  } 

plot(dates$Date,dates[,2],type="l",lwd=3,col=clist[1],
     ylim=c(0,ceiling(max(dates[,2:5],na.rm=T))),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
for(k in 2:4) lines(dates$Date,dates[,k+1],lwd=3,col=clist[k])
legend("topleft",legend=SSTtype,col=clist,lwd=3,bty='n')

##Tracks??
