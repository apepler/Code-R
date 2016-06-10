rm(list=ls())
setwd("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R2_ensemble/")
events<-fixes<-list()
a=list.dirs(path=".")

dirlist=c("d01-as-d02_p60_topo","d01-as-d02_p240_topo","d02_p60_topo","d02_p240_topo")

for(i in 1:2)
{
  events[[i]]=read.csv(paste(dirlist[i],"/ECLevents_2007052700_d01.csv",sep=""))
  events[[i+2]]=read.csv(paste(dirlist[i+2],"/ECLevents_2007052700_d02.csv",sep=""))
  fixes[[i]]=read.csv(paste(dirlist[i],"/ECLfixes_2007052700_d01.csv",sep=""))
  fixes[[i+2]]=read.csv(paste(dirlist[i+2],"/ECLfixes_2007052700_d02.csv",sep=""))
}

dateCV<-cbind(seq(20070601,20070630.75,0.25),matrix(NaN,120,4))
for(n in 1:4)
{
  a=fixes[[n]]
  a$Date2=a$Date+(as.numeric(a$Time)-1)/4
  
for(i in 1:120)
{
  I=which(a$Date2==dateCV[i,1] & a$Location==1)
  if(length(I)>0) dateCV[i,n+1]=mean(a$CV[I])
}
}

lcol=c("red","blue","red","blue")
ltype=c(1,1,2,2)
plot(dateCV[,1],dateCV[,2],col=lcol[1],lwd=2,xlab="Date",ylab="Curvature",ylim=c(0,5),type="l")
for(i in c(2,4))
{
  lines(dateCV[,1],dateCV[,i+1],col=lcol[i],lwd=2,lty=ltype[i])
}
legend("topleft",legend=c("d01 p60","d01 p240","d02 p60","d02 p240"),lty=ltype,lwd=2,col=lcol)


