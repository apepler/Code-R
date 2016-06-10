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

load('~/Documents/Data/Useful_WRF.RData')
dailymaxwind<-dailymaxwind_nt<-array(0,c(325,200,30,20))

n=1
for(day in 27:31)
  for(hour in c("00","06","12","18"))
  {
    dir1=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble/out/200705",day,hour,"/",sep="")
    dir2=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/200705",day,hour,"/",sep="")    
    
    for(i in 1:30)
    {
    fname=paste("wrfxtrm_d02_2007-06-",sprintf("%2.2d",i),"_00:00:00",sep="") 
    a=open.nc(paste(dir1,fname,sep=""))
    latW=var.get.nc(a,"XLAT",c(1,1,1),c(325,200,1))
    lonW=var.get.nc(a,"XLONG",c(1,1,1),c(325,200,1))
    dailymaxwind[,,i,n]=var.get.nc(a,"SPDUV10MAX")
    close.nc(a)
    a=open.nc(paste(dir2,fname,sep=""))
    dailymaxwind_nt[,,i,n]=var.get.nc(a,"SPDUV10MAX")
    close.nc(a)
    }
    n=n+1
  }

thresholds=c(50,60,70,80,90)
threshms=thresholds*1000/3600

windthresh<-windthresh_nt<-array(0,c(30,5,20))
ncells=length(which(is.na(mask10a)==0))

for(g in 1:5)
for(n in 1:20)
  for(i in 1:30)
  {
    I=which(dailymaxwind[,,i,n]*mask10a>=threshms[g])
    if(length(I)>0) windthresh[i,g,n]=(length(I)/ncells)*100
    I=which(dailymaxwind_nt[,,i,n]*mask10a>=threshms[g])
    if(length(I)>0) windthresh_nt[i,g,n]=(length(I)/ncells)*100
  }


pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D02_dailywindthresh_notopo.pdf",
    width=6,height=4)
plot(seq(20070601,20070630),apply(windthresh[,1,],1,mean),type="l",lwd=3,col="blue",
     ylim=c(0,30),xlab="Date",ylab="Area (%)",main="Daily mean area of ESB exceeding hourly wind thresholds")
lines(seq(20070601,20070630),apply(windthresh[,2,],1,mean),lwd=3,col="blue",lty=2)
for(i in 1:2) lines(seq(20070601,20070630),apply(windthresh_nt[,i,],1,mean),lwd=3,col="red",lty=i)
legend("topleft",legend=c("50 km/h","60 km/h"),col="black",lwd=3,bty="n",lty=c(1,2))
legend("topright",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3,bty="n")
dev.off()


a=density(apply(windthresh[,1,],2,max),na.rm=T)
b=density(apply(windthresh_nt[,1,],2,max),na.rm=T)
lims=range(apply(windthresh[,1,],2,max),apply(windthresh_nt[,1,],2,max))

lims[1]=floor(lims[1]/5)*5
lims[2]=ceiling(lims[2]/5)*5

#pdf(file=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D02_Jun",day,"_area100mm_notopo.pdf",sep=""),
pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D02_maxarea50kmh_notopo.pdf",
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Wind area (%)",ylab="Frequency",cex.main=1.2,
     main=paste("Percentage area above 50 km/h"))
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()
#}

mamax=apply(dailymaxwind,c(1,2,4),max)
for(i in 1:20) mamax[,,i]=mamax[,,i]*mask10a
mamax_nt=apply(dailymaxwind_nt,c(1,2,4),max)
for(i in 1:20) mamax_nt[,,i]=mamax_nt[,,i]*mask10a

a=density(as.vector(mamax[!is.na(mamax)])*36/10,na.rm=T)
b=density(as.vector(mamax_nt[!is.na(mamax_nt)])*36/10,na.rm=T)
lims=range(as.vector(mamax[!is.na(mamax)]),as.vector(mamax_nt[!is.na(mamax_nt)]),0)*36/10

lims[1]=floor(lims[1]/5)*5
lims[2]=ceiling(lims[2]/5)*5

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D02_maxdailywind_allpoints_notopo2.pdf",
    width=6,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Wind speed (km/h)",ylab="Frequency",cex.main=1.2,
     main=paste("Spatial distribution of windiest hour"))
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("Control","NoTopo"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")   
dev.off()