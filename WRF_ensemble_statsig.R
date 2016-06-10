rm(list=ls())

####
## WRF total rain

for(rm in c("R1","R2","R3"))
  for(res in c("D01","D02"))
  {
    read.csv(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_",rm,"_ensemble_notopo/out/",res,"_totalESBrain.csv",sep=""))->meanR
    bootstats<-bootstats2<-matrix(0,20,4)
    for(thresh in 1:20)
    {
      means<-rep(0,1000)
      for(i in 1:1000)
      {
        a=sample(meanR[,4],thresh,replace=T)
        means[i]=mean(a)
      }
      
      bootstats[thresh,1]=mean(means)
      bootstats[thresh,2]=sd(means)
      bootstats[thresh,3]=quantile(means,0.95)-quantile(means,0.05)
      bootstats[thresh,4]=sum(means<0)/length(means)
      
      means<-rep(0,1000)
      for(i in 1:1000)
      {
        a=sample(meanR[,4],thresh,replace=F)
        means[i]=mean(a)
      }
      
      bootstats2[thresh,1]=mean(means)
      bootstats2[thresh,2]=sd(means)
      bootstats2[thresh,3]=quantile(means,0.95)-quantile(means,0.05)
      bootstats2[thresh,4]=sum(means<0)/length(means)
    }
    pdf(file=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_",rm,"_ensemble_notopo/",res,"_totalESBraindiff_bootvariance.pdf",sep=""),
        width=5,height=5)
    plot(seq(1:20),abs(bootstats[,2]/bootstats[,1]),type="l",lwd=2,col="red",ylim=c(0,0.5),
         xlab="Number of members",ylab="Ratio of standard deviation to mean")
    lines(seq(1:20),abs(bootstats2[,2]/bootstats2[,1]),type="l",lwd=2,col="blue")
    legend("topright",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
    dev.off()
  }




load('/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out/WRFdailyrain_ensemble_d01.RData')
load('~/Documents/Data/Useful_WRF.RData')
dailyR<-dailyR_nt<-matrix(0,30,20)
for(i in 1:30)
  for(n in 1:20)
  {
    dailyR[i,n]=mean(dailyrain[,,i,n]*mask50a,na.rm=T)
    dailyR_nt[i,n]=mean(dailyrain_nt[,,i,n]*mask50a,na.rm=T)
  }

maxR<-matrix(0,20,3)
maxR[,1]=apply(dailyR,2,max)
maxR[,2]=apply(dailyR_nt,2,max)
maxR[,3]=maxR[,2]-maxR[,1]

bootstats<-matrix(NaN,20,2)
for(thresh in 1:20)
{
  pvals<-matrix(0,1000,2)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=T)
    b=sample(1:20,thresh,replace=F)
    a=ks.test(maxR[a,1],maxR[a,2])
    pvals[i,1]=a$p.value
    a=ks.test(maxR[b,1],maxR[b,2])
    pvals[i,2]=a$p.value
  }
  for(i in 1:2) bootstats[thresh,i]=100*sum(pvals[,i]<0.05)/length(pvals[,i])
}
pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_maxESBrain_bootKStest.pdf",
    width=5,height=5)
plot(seq(1:20),bootstats[,1],type="l",lwd=2,col="red",ylim=c(0,100),
     xlab="Number of members",ylab="% of time a ks-test has p<0.05")
lines(seq(1:20),bootstats[,2],type="l",lwd=2,col="blue")
legend("bottomright",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
dev.off()


####### ECL stats

setwd("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/out")
events<-fixes<-events_notopo<-fixes_notopo <- list()

n=1
for(day in 27:31)
  for(hour in c("00","06","12","18"))
  {
    dir1="/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R2_ensemble/d01_p60/"
    dir2="/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R2_ensemble_notopo/d01_p60/"
    events[[n]]=read.csv(paste(dir1,"ECLevents_200705",day,hour,".csv",sep=""))
    fixes[[n]]=read.csv(paste(dir1,"ECLfixes_200705",day,hour,".csv",sep=""))
    events_notopo[[n]]=read.csv(paste(dir2,"ECLevents_200705",day,hour,".csv",sep=""))
    fixes_notopo[[n]]=read.csv(paste(dir2,"ECLfixes_200705",day,hour,".csv",sep=""))
    
    n=n+1
  }


bootstats<-bootstats2<-matrix(0,20,4)

for(thresh in 1:20)
{
  pvals<-pvalsR<-matrix(0,1000,4)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=T)
    b=sample(1:20,thresh,replace=F)
    
    E1=events[[a[1]]]
    E2=events[[b[1]]]
    E1nt=events_notopo[[a[1]]]
    E2nt=events_notopo[[b[1]]]
    
    if(thresh>1)
      for(nn in 2:thresh)
      {
        E1=rbind(E1,events[[a[nn]]])
        E2=rbind(E2,events[[b[nn]]])
        E1nt=rbind(E1nt,events_notopo[[a[nn]]])
        E2nt=rbind(E2nt,events_notopo[[b[nn]]])
      }
    
    a=ks.test(E1$CV2,E1nt$CV2)
    pvalsR[i,1]=a$p.value
    a=ks.test(E1$MSLP2,E1nt$MSLP2)
    pvalsR[i,2]=a$p.value
    a=ks.test(E1$Rad2,E1nt$Rad2)
    pvalsR[i,3]=a$p.value
    a=ks.test(E1$Length2,E1nt$Length2)
    pvalsR[i,4]=a$p.value
    
    a=ks.test(E2$CV2,E2nt$CV2)
    pvals[i,1]=a$p.value
    a=ks.test(E2$MSLP2,E2nt$MSLP2)
    pvals[i,2]=a$p.value
    a=ks.test(E2$Rad2,E2nt$Rad2)
    pvals[i,3]=a$p.value
    a=ks.test(E2$Length2,E2nt$Length2)
    pvals[i,4]=a$p.value
  }
  
  for(i in 1:4)
  {
    bootstats[thresh,i]=sum(pvalsR[,i]<0.05)/length(pvalsR[,i])
    bootstats2[thresh,i]=sum(pvals[,i]<0.05)/length(pvals[,i])
  }
}

bootstats=bootstats*100
bootstats2=bootstats2*100

type=c("CV","MSLP","Rad","Length")
for(i in 1:4)
{
  pdf(file=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_ECL",type[i],"_bootKS.pdf",sep=""),
      width=5,height=5)
  plot(seq(1:20),bootstats[,i],type="l",lwd=2,col="red",ylim=c(0,100),
       xlab="Number of members",ylab="% of time a KS test has p<0.05")
  lines(seq(1:20),bootstats2[,i],type="l",lwd=2,col="blue")
  legend("topleft",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
  dev.off()
}

bootstats<-bootstats2<-matrix(0,20,3)

for(thresh in 1:20)
{
  pvals<-pvalsR<-matrix(0,1000,3)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=T)
    b=sample(1:20,thresh,replace=F)
    
    E1=fixes[[a[1]]]
    E2=fixes[[b[1]]]
    E1nt=fixes_notopo[[a[1]]]
    E2nt=fixes_notopo[[b[1]]]
    
    if(thresh>1)
      for(nn in 2:thresh)
      {
        E1=rbind(E1,fixes[[a[nn]]])
        E2=rbind(E2,fixes[[b[nn]]])
        E1nt=rbind(E1nt,fixes_notopo[[a[nn]]])
        E2nt=rbind(E2nt,fixes_notopo[[b[nn]]])
      }
    
    a=ks.test(E1$CV[E1$Location==1],E1nt$CV[E1$Location==1])
    pvalsR[i,1]=a$p.value
    a=ks.test(E1$MSLP[E1$Location==1],E1nt$MSLP[E1nt$Location==1])
    pvalsR[i,2]=a$p.value
    a=ks.test(E1$Rad[E1$Location==1],E1nt$Rad[E1nt$Location==1])
    pvalsR[i,3]=a$p.value
    
    
    a=ks.test(E2$CV[E2$Location==1],E2nt$CV[E2nt$Location==1])
    pvals[i,1]=a$p.value
    a=ks.test(E2$MSLP[E2$Location==1],E2nt$MSLP[E2nt$Location==1])
    pvals[i,2]=a$p.value
    a=ks.test(E2$Rad[E2$Location==1],E2nt$Rad[E2nt$Location==1])
    pvals[i,3]=a$p.value
  }
  
  for(i in 1:3)
  {
    bootstats[thresh,i]=100*sum(pvalsR[,i]<0.05)/length(pvalsR[,i])
    bootstats2[thresh,i]=100*sum(pvals[,i]<0.05)/length(pvals[,i])
  }
}

type=c("CV","MSLP","Rad")
for(i in 1:3)
{
  pdf(file=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_ECLfix_",type[i],"_bootKS.pdf",sep=""),
      width=5,height=5)
  plot(seq(1:20),bootstats[,i],type="l",lwd=2,col="red",ylim=c(0,100),
       xlab="Number of members",ylab="% of time a KS test has p<0.05")
  lines(seq(1:20),bootstats2[,i],type="l",lwd=2,col="blue")
  legend("topleft",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
  dev.off()
}

CVmax<-MSLPmin<-matrix(0,20,2)
for(i in 1:20)
{
  CVmax[i,1]=max(events[[i]]$CV2)
  CVmax[i,2]=max(events_notopo[[i]]$CV2)
  MSLPmin[i,1]=min(events[[i]]$MSLP2)
  MSLPmin[i,2]=min(events_notopo[[i]]$MSLP2)
}

bootstats<-matrix(NaN,20,2)
for(thresh in 5:20)
{
  pvals<-matrix(0,1000,2)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=T)
    b=sample(1:20,thresh,replace=F)
    a=t.test(CVmax[a,1],CVmax[a,2])
    pvals[i,1]=a$p.value
    a=t.test(CVmax[b,1],CVmax[b,2])
    pvals[i,2]=a$p.value
  }
  for(i in 1:2) bootstats[thresh,i]=100*sum(pvals[,i]<0.05)/length(pvals[,i])
}
pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_ECL_maxCV_bootTtest.pdf",
    width=5,height=5)
plot(seq(1:20),bootstats[,1],type="l",lwd=2,col="red",ylim=c(0,100),
     xlab="Number of members",ylab="% of time a ks-test has p<0.05")
lines(seq(1:20),bootstats[,2],type="l",lwd=2,col="blue")
legend("bottomright",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
dev.off()

CVmax=cbind(CVmax,CVmax[,2]-CVmax[,1])
MSLPmin=cbind(MSLPmin,MSLPmin[,2]-MSLPmin[,1])

bootstats<-bootstats2<-matrix(0,20,4)
for(thresh in 1:20)
{
  means<-rep(0,1000)
  for(i in 1:1000)
  {
    a=sample(MSLPmin[,3],thresh,replace=T)
    means[i]=mean(a)
  }
  
  bootstats[thresh,1]=mean(means)
  bootstats[thresh,2]=sd(means)
  bootstats[thresh,3]=quantile(means,0.95)-quantile(means,0.05)
  bootstats[thresh,4]=sum(means<0)/length(means)
  
  means<-rep(0,1000)
  for(i in 1:1000)
  {
    a=sample(MSLPmin[,3],thresh,replace=F)
    means[i]=mean(a)
  }
  
  bootstats2[thresh,1]=mean(means)
  bootstats2[thresh,2]=sd(means)
  bootstats2[thresh,3]=quantile(means,0.95)-quantile(means,0.05)
  bootstats2[thresh,4]=sum(means<0)/length(means)
}
pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_ECL_minMSLP_bootvariance.pdf",
    width=5,height=5)
plot(seq(1:20),abs(bootstats[,2]/bootstats[,1]),type="l",lwd=2,col="red",ylim=c(0,0.5),
     xlab="Number of members",ylab="Ratio of standard deviation to mean")
lines(seq(1:20),abs(bootstats2[,2]/bootstats2[,1]),type="l",lwd=2,col="blue")
legend("topright",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
dev.off()

########## Individual events?


keydates=c(20070626)
for(tt in 1:2)
{
  CVmax2<-matrix(NaN,20,2)
  for(i in 1:20)
  {
    b=fixes[[i]]
    I=which(b$Date==keydates[tt] & b$Time=="00:00")
    a=unique(b$ID[I])
    if(length(a)>0) CVmax2[i,1]=max(b$CV[b$Location==1 & b$ID%in%a])
    b=fixes_notopo[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    if(length(a)>0) CVmax2[i,2]=max(b$CV[b$Location==1 & b$ID%in%a])  
  }
  
  bootstats<-matrix(NaN,20,2)
  for(thresh in 10:20)
  {
    pvals<-matrix(0,1000,2)
    for(i in 1:1000)
    {
      a=sample(1:20,thresh,replace=T)
      b=sample(1:20,thresh,replace=F)
      a=t.test(CVmax2[a,1],CVmax2[a,2])
      pvals[i,1]=a$p.value
      a=t.test(CVmax2[b,1],CVmax2[b,2])
      pvals[i,2]=a$p.value
    }
    for(i in 1:2) bootstats[thresh,i]=100*sum(pvals[,i]<0.05)/length(pvals[,i])
  }
  pdf(file=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble_notopo/D01_ECL_maxCV_bootTtest_",keydates[tt],".pdf",sep=""),
      width=5,height=5)
  plot(seq(1:20),bootstats[,1],type="l",lwd=2,col="red",ylim=c(0,100),
       xlab="Number of members",ylab="% of time a ks-test has p<0.05")
  lines(seq(1:20),bootstats[,2],type="l",lwd=2,col="blue")
  legend("bottomright",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
  dev.off()
}