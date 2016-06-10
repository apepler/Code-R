setwd('~/output/outputUM_wrf_cases/')

makePDF = function(data1,data2,xlabel,labloc="topleft") {
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
       main="Case study ECL statistics")
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=c("Control","NoTopo"),
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
}

wrfv=c("R1","R2","R3")
fixes <- fixesNT <- list()

for(i in 1:3)
{
  fixes[[i]]=read.csv(paste('ECLfixes_cases_',wrfv[i],'_nudging.csv',sep=""))
  fixesNT[[i]]=read.csv(paste('ECLfixes_cases_',wrfv[i],'_nudging_notopo.csv',sep=""))
}
fixes2=rbind(fixes[[1]],fixes[[2]],fixes[[3]])
fixesNT2=rbind(fixesNT[[1]],fixesNT[[2]],fixesNT[[3]])

pdf(file="D01_ECLs_fixCV.pdf",width=6,height=4)
makePDF(fixes2$CV[fixes2$Location==1],fixesNT2$CV[fixesNT2$Location==1],"Curvature for all lows in the ESB",labloc="topright") 
dev.off()
pdf(file="D01_ECLs_fixMSLP.pdf",width=6,height=4)
makePDF(fixes2$MSLP[fixes2$Location==1],fixesNT2$MSLP[fixesNT2$Location==1],"MSLP for all lows in the ESB")
dev.off()
  
fixes2$Date2=fixes2$Date+(as.numeric(fixes2$Time)-1)/4
fixesNT2$Date2=fixesNT2$Date+(as.numeric(fixesNT2$Time)-1)/4

date1=c(200706,200711,200808)
eday=c(30,30,31)

for(n in 1:3)
{
dates=data.frame(Date=seq((date1[n]*100+1),(date1[n]*100+eday[n]+0.75),0.25),
                 CV=rep(NaN,eday[n]*4),MSLP=rep(NaN,eday[n]*4),CV_NT=rep(NaN,eday[n]*4),MSLP_NT=rep(NaN,eday[n]*4))
for(i in 1:length(dates[,1]))
{
  I=which(fixes2$Date2==dates$Date[i] & fixes2$Location==1)
  if(length(I)>0){
    dates[i,2]=mean(fixes2$CV[I])
    dates[i,3]=mean(fixes2$MSLP[I])
  } 
  
  I=which(fixesNT2$Date2==dates$Date[i] & fixesNT2$Location==1)
  if(length(I)>0){
    dates[i,4]=mean(fixesNT2$CV[I])
    dates[i,5]=mean(fixesNT2$MSLP[I])
  } 
}

pdf(file=paste("D01_ECLs_dailyCV_",date1[n],".pdf",sep=""),width=7,height=4)
plot(dates$Date,dates[,2],type="l",lwd=3,col="blue",
     ylim=c(0,3),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
lines(dates$Date,dates[,4],lwd=3,col="red")
legend("topleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3)
dev.off()
pdf(file=paste("D01_ECLs_dailyMSLP_",date1[n],".pdf",sep=""),width=7,height=4)
plot(dates$Date,dates[,3],type="l",lwd=3,col="blue",
     ylim=c(970,1030),xlab="Date",ylab="MSLP",main="Mean ECL MSLP across all runs")
lines(dates$Date,dates[,5],lwd=3,col="red")
legend("bottomleft",legend=c("Control","NoTopo"),col=c("blue","red"),lwd=3)
dev.off()
}

  

######### Similar tests for all of 2007


  