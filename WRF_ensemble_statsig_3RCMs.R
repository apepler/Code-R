rm(list=ls())

####
## WRF total rain

rainT<-rainNT<-diff<-matrix(0,20,3)
rm=c("R1","R2","R3")
for(i in 1:3)
{
  read.csv(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_",rm[i],"_ensemble_notopo/out/D01_totalESBrain.csv",sep=""))->meanR
  rainT[,i]=meanR[,2]
  rainNT[,i]=meanR[,3]
  diff[,i]=meanR[,4]
}

bootstats<-bootstats2<-matrix(NaN,20,4)
for(thresh in 1:20)
{
  means<-rep(0,1000)
  pvals<-matrix(0,1000,2)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=T)
    means[i]=mean(diff[a,])
    b=ks.test(rainT[a,],rainNT[a,])
    pvals[i,1]=b$p.value
    if(thresh>5)
      {
      b=t.test(rainT[a,],rainNT[a,])
      pvals[i,2]=b$p.value
    }
  }
  
  bootstats[thresh,1]=mean(means)
  bootstats[thresh,2]=sd(means)
  bootstats[thresh,3]=100*sum(pvals[,1]<0.05)/length(pvals[,1])
  if(thresh>5) bootstats[thresh,4]=100*sum(pvals[,2]<0.05)/length(pvals[,2])
  
  means<-rep(0,1000)
  pvals<-matrix(0,1000,2)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=F)
    means[i]=mean(diff[a,])
    b=ks.test(rainT[a,],rainNT[a,])
    pvals[i,1]=b$p.value
    if(thresh>5)
    {
      b=t.test(rainT[a,],rainNT[a,])
      pvals[i,2]=b$p.value
    }
  }
  
  bootstats2[thresh,1]=mean(means)
  bootstats2[thresh,2]=sd(means)
  bootstats2[thresh,3]=100*sum(pvals[,1]<0.05)/length(pvals[,1])
  if(thresh>5) bootstats2[thresh,4]=100*sum(pvals[,2]<0.05)/length(pvals[,2])
}


pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_ensembleNT_figs/D01all_totalESBraindiff_bootvariance.pdf",
    width=5,height=5)
plot(seq(1:20),abs(bootstats[,2]/bootstats[,1]),type="l",lwd=2,col="red",ylim=c(0,0.5),
     xlab="Number of members",ylab="Ratio of standard deviation to mean")
lines(seq(1:20),abs(bootstats2[,2]/bootstats2[,1]),type="l",lwd=2,col="blue")
legend("topright",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
dev.off()

### All for CV/MSLP changes

setwd("/srv/ccrc/data36/z3478332/WRF/output/ERAI_ensembleNT_figs")
CVmax<-CVmaxNT<-MSLPmin<-MSLPminNT<-matrix(0,20,3)
for(i in 1:3)
{
  n=1
for(day in 27:31)
  for(hour in c("00","06","12","18"))
  {
    dir1=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_",rm[i],"_ensemble/d01_p60/",sep="")
    dir2=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_",rm[i],"_ensemble_notopo/d01_p60/",sep="")
    events=read.csv(paste(dir1,"ECLevents_200705",day,hour,".csv",sep=""))
    events_NT=read.csv(paste(dir2,"ECLevents_200705",day,hour,".csv",sep=""))
    
    CVmax[n,i]=max(events$CV2)
    CVmaxNT[n,i]=max(events_NT$CV2)
    MSLPmin[n,i]=min(events$MSLP2)
    MSLPminNT[n,i]=min(events_NT$MSLP2)
    n=n+1
  }
}

CVdiff=CVmaxNT-CVmax
MSLPdiff=MSLPminNT-MSLPmin

bootstats<-bootstats2<-matrix(NaN,20,4)
for(thresh in 1:20)
{
  means<-rep(0,1000)
  pvals<-matrix(0,1000,2)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=T)
    means[i]=mean(CVdiff[a,])
    b=ks.test(CVmax[a,],CVmaxNT[a,])
    pvals[i,1]=b$p.value
    if(thresh>5)
    {
      b=t.test(CVmax[a,],CVmaxNT[a,])
      pvals[i,2]=b$p.value
    }
  }
  
  bootstats[thresh,1]=mean(means)
  bootstats[thresh,2]=sd(means)
  bootstats[thresh,3]=100*sum(pvals[,1]<0.05)/length(pvals[,1])
  if(thresh>5) bootstats[thresh,4]=100*sum(pvals[,2]<0.05)/length(pvals[,2])
  
  means<-rep(0,1000)
  pvals<-matrix(0,1000,2)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=F)
    means[i]=mean(CVdiff[a,])
    b=ks.test(CVmax[a,],CVmaxNT[a,])
    pvals[i,1]=b$p.value
    if(thresh>5)
    {
      b=t.test(CVmax[a,],CVmaxNT[a,])
      pvals[i,2]=b$p.value
    }
  }
  
  bootstats2[thresh,1]=mean(means)
  bootstats2[thresh,2]=sd(means)
  bootstats2[thresh,3]=100*sum(pvals[,1]<0.05)/length(pvals[,1])
  if(thresh>5) bootstats2[thresh,4]=100*sum(pvals[,2]<0.05)/length(pvals[,2])
}

pdf(file="/srv/ccrc/data36/z3478332/WRF/output/ERAI_ensembleNT_figs/D01all_ECL_maxCV_bootvariance.pdf",
    width=5,height=5)
plot(seq(1:20),abs(bootstats[,2]/bootstats[,1]),type="l",lwd=2,col="red",ylim=c(0,0.5),
     xlab="Number of members",ylab="Ratio of standard deviation to mean")
lines(seq(1:20),abs(bootstats2[,2]/bootstats2[,1]),type="l",lwd=2,col="blue")
legend("topright",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
dev.off()


########## Distribution across all events

fname<-fname_NT<-matrix("aaa",20,3)
n=1
for(i in 1:3)
{
  n=1
for(day in 27:31)
  for(hour in c("00","06","12","18"))
  {
    fname[n,i]=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",i,"_ensemble/d01_p60/ECLfixes_200705",day,hour,".csv",sep="")
    fname_NT[n,i]=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",i,"_ensemble_notopo/d01_p60/ECLfixes_200705",day,hour,".csv",sep="")
    n=n+1
  }
}


bootstats<-bootstats2<-matrix(NaN,20,4)

for(thresh in 1:20)
{
  pvals<-pvalsR<-matrix(NaN,1000,4)
  for(i in 1:1000)
  {
    a=sample(1:20,thresh,replace=T)
    b=sample(1:20,thresh,replace=F)
    
    E1=read.csv(fname[a[1],1])
    E2=read.csv(fname[b[1],1])
    E1nt=read.csv(fname_NT[a[1],1])
    E2nt=read.csv(fname_NT[b[1],1])
    
    for(rr in 2:3)
    {
      E1=rbind(E1,read.csv(fname[a[1],rr]))
      E2=rbind(E2,read.csv(fname[b[1],rr]))
      E1nt=rbind(E1nt,read.csv(fname_NT[a[1],rr]))
      E2nt=rbind(E2nt,read.csv(fname_NT[b[1],rr]))
    }
    
    if(thresh>1)
      for(nn in 2:thresh)
        for(rr in 1:3)
      {
          E1=rbind(E1,read.csv(fname[a[nn],rr]))
          E2=rbind(E2,read.csv(fname[b[nn],rr]))
          E1nt=rbind(E1nt,read.csv(fname_NT[a[nn],rr]))
          E2nt=rbind(E2nt,read.csv(fname_NT[b[nn],rr]))
      }
    
    a=ks.test(E1$CV[E1$Location==1],E1nt$CV[E1nt$Location==1])
    pvalsR[i,1]=a$p.value
    a=ks.test(E1$MSLP[E1$Location==1],E1nt$MSLP[E1nt$Location==1])
    pvalsR[i,2]=a$p.value
    if(thresh>5)
    {
    a=t.test(E1$CV[E1$Location==1],E1nt$CV[E1nt$Location==1])
    pvalsR[i,3]=a$p.value
    a=t.test(E1$MSLP[E1$Location==1],E1nt$MSLP[E1nt$Location==1])
    pvalsR[i,4]=a$p.value
    }
    
    a=ks.test(E2$CV[E2$Location==1],E2nt$CV[E2nt$Location==1])
    pvals[i,1]=a$p.value
    a=ks.test(E2$MSLP[E2$Location==1],E2nt$MSLP[E2nt$Location==1])
    pvals[i,2]=a$p.value
    if(thresh>5)
    {
      a=t.test(E2$CV[E2$Location==1],E2nt$CV[E2nt$Location==1])
      pvals[i,3]=a$p.value
      a=t.test(E2$MSLP[E2$Location==1],E2nt$MSLP[E2nt$Location==1])
      pvals[i,4]=a$p.value
    }
  }
  
  for(i in 1:4)
  {
    bootstats[thresh,i]=100*sum(pvalsR[,i]<0.05)/length(pvalsR[,i])
    bootstats2[thresh,i]=100*sum(pvals[,i]<0.05)/length(pvals[,i])
  }
}

type=c("CV","MSLP")
for(i in 1:2)
{
  pdf(file=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_ensembleNT_figs/D01all_ECLfix_",type[i],"_bootKStest.pdf",sep=""),
      width=5,height=5)
  plot(seq(1:20),bootstats[,i],type="l",lwd=2,col="red",ylim=c(0,100),
       xlab="Number of members",ylab="% of time a KS test has p<0.05")
  lines(seq(1:20),bootstats2[,i],type="l",lwd=2,col="blue")
  legend("topleft",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
  dev.off()
  
  pdf(file=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_ensembleNT_figs/D01all_ECLfix_",type[i],"_bootTtest.pdf",sep=""),
      width=5,height=5)
  plot(seq(1:20),bootstats[,i+2],type="l",lwd=2,col="red",ylim=c(0,100),
       xlab="Number of members",ylab="% of time a t-test has p<0.05")
  lines(seq(1:20),bootstats2[,i+2],type="l",lwd=2,col="blue")
  legend("topleft",lwd=2,col=c("red","blue"),legend=c("With replacement","No replacement"))
  dev.off()
}
