setwd("~/Documents/ECLs/Algorithm Comparison/Alejandro/")
read.csv("ERAI-75-6-global-150-6_v12-1-1.csv")->data
data=data[,-seq(16,19)]
data=data[,-19]
data$date=data[,21]+data[,20]*100+data[,19]*10000
data$ASPloc=0
I<-which(data[,4]>=149 & data[,4]<=161 & data[,3]<(-37) & data[,3]>=-41)
data[I,24]<-1
I<-which(data[,4]>=(149+(37+data[,3])/2) & data[,4]<=161 & data[,3]<(-31) & data[,3]>=-37)
data[I,24]<-1
I<-which(data[,4]>=152 & data[,4]<=161 & data[,3]<=(-24) & data[,3]>=-31)
data[I,24]<-1

##Need to do an events basis too? 

x<-rle(data[,1])
events<-cbind(x$values,matrix(data=0,nrow=length(x$values),ncol=8))

for(i in 1:length(events[,1])) 
{
  I<-which(data[,1]==events[i,1])
  events[i,2]=min(data[I,23]) ##Date1
  events[i,3]=max(data[I,23]) ##Date2
  events[i,4]=min(data[I,5]) ##Min central pressure
  events[i,5]=max(data[I,7]) ##Max PG
  J=which(data[I,24]==1 & data[I,17]==-1)
  if(length(J)>0) events[i,6]=1 ##Ever closed in ASP region
  J=which(data[I,17]==-1)
  events[i,7]=length(J) ##Duration of closed
  J=which(data[I,24]==1)
  if(length(J)>0)
  {
    events[i,8]=max(data[I[J],5]) ##Min P in region
    events[i,9]=max(data[I[J],7]) #Max PG in region
  }
}

events<-data.frame(events[,1],events[,7],events[,6],events[,2:5],events[,8:9])
names(events)=c('ID','Length','Ever closed in region','Date1','Date2','MSLP','PG','MSLP (loc)','PG (loc)')

fixes<-data.frame(ID=data[,1],Length=data[,18],Date=data[,23],Time=data[,22],Closed=data[,17],Lon=data[,4],Lat=data[,3],MSLP=data[,5],Grad=data[,7],Location=data[,24])

write.csv(fixes,file='Alejandro_ERAI_8009_v12.csv')
write.csv(events,file='Alejandro_ERAI_8009_v12_events.csv')

fixes$Fixes=0
fixes[11,1]=1
for(i in 2:length(fixes[,1]))
  if(fixes[i,1]==fixes[i-1,1]) fixes[11,i]=fixes[11,i-1]+1 else fixes[11,i]=1


fixes2=fixes[fixes$Grad>=0.01,]
fixes2$ID2=0
fixes2$fix=0
for(i in 2:length(fixes2[,1]))
{
  if(fixes[1,i]==fixes)
}

x<-rle(data[,1])
events2<-cbind(x$values,matrix(data=0,nrow=length(x$values),ncol=8))

for(i in 1:length(events[,1])) 
{
  I<-which(data[,1]==events[i,1])
  events2[i,2]=min(data[I,23]) ##Date1
  events[i,3]=max(data[I,23]) ##Date2
  events[i,4]=min(data[I,5]) ##Min central pressure
  events[i,5]=max(data[I,7]) ##Max PG
  J=which(data[I,24]==1 & data[I,17]==-1)
  if(length(J)>0) events[i,6]=1 ##Ever closed in ASP region
  J=which(data[I,17]==-1)
  events[i,7]=length(J) ##Duration of closed
  J=which(data[I,24]==1)
  if(length(J)>0)
  {
    events[i,8]=max(data[I[J],5]) ##Min P in region
    events[i,9]=max(data[I[J],7]) #Max PG in region
  }
}

events<-data.frame(events[,1],events[,7],events[,6],events[,2:5],events[,8:9])
names(events)=c('ID','Length','Ever closed in region','Date1','Date2','MSLP','PG','MSLP (loc)','PG (loc)')


I=which(events[,2]>1)
events2=events[I,]
include<-match(data[,1],events2[,1])
J<-which(is.na(include)==0)
data2=data[J,]

write.csv(data2,file='Alejandro_ERAI_8009_6.csv')
write.csv(events2,file='Alejandro_ERAI_8009_6_events.csv')

I=which(data2[,24]==1 & data2[,17]==-1)
write.csv(unique(data2[I,23]),file='Alejandro_ERAI_8009_6_days.csv')
I=which(data[,24]==1 & data[,17]==-1)
write.csv(unique(data[I,23]),file='Alejandro_ERAI_8009_3_days.csv')

#####Comparison of event tracks
setwd("~/Documents/ECLs/")
read.csv('CSV/ECLfixes_unsw_ERAI.csv',sep=";")->mine
read.csv('CSV/ECLfixes_unsw_ERA925.csv',sep=";")->e925
read.csv('CSV/ECLfixes_Fei.csv')->fei
read.csv('CSV/ECLfixes_mldb.csv')->mldb
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

##Boxing day 1998
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(130,170),ylim=c(-50,-20),main="Boxing Day ECL, December 1998")
lines(mine[(mine[,1]==2262),6],mine[(mine[,1]==2262),7],type="b",col="blue",pch=1)
lines(mine[(mine[,1]==2263),6],mine[(mine[,1]==2263),7],type="b",col="blue",pch=3)
lines(mine[(mine[,1]==2264),6],mine[(mine[,1]==2264),7],type="b",col="blue",pch=2)
lines(e925[(e925[,1]==1216),6],e925[(e925[,1]==1216),7],type="b",col="red",pch=1)
lines(e925[(e925[,1]==1217),6],e925[(e925[,1]==1217),7],type="b",col="red",pch=2)
lines(fei[(fei[,1]==444),7],fei[(fei[,1]==444),8],type="b",col="green",pch=2)
points(mldb[(mldb[,2]>=19981224 & mldb[,2]<=19981229),4],mldb[(mldb[,2]>=19981224 & mldb[,2]<=19981229),3],col="black",pch=4,cex=2)
legend("bottomright",legend=c("MLDB","Pepler","Pepler 925hPa","Fei"),pch=4,col=c("black","blue","red","green"))

##Queens Birthday, 8-9 June 2007
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(130,170),ylim=c(-50,-20),main="Pasha Bulker ECL, June 2007")
lines(mine[(mine[,1]==3216),6],mine[(mine[,1]==3216),7],type="b",col="blue",pch=1)
lines(mine[(mine[,1]==3217),6],mine[(mine[,1]==3217),7],type="b",col="blue",pch=3)
lines(mine[(mine[,1]==3218),6],mine[(mine[,1]==3218),7],type="b",col="blue",pch=2)
lines(e925[(e925[,1]==1698),6],e925[(e925[,1]==1698),7],type="b",col="red",pch=1)
lines(e925[(e925[,1]==1699),6],e925[(e925[,1]==1699),7],type="b",col="red",pch=2)
lines(e925[(e925[,1]==1700),6],e925[(e925[,1]==1700),7],type="b",col="red",pch=3)
lines(e925[(e925[,1]==1701),6],e925[(e925[,1]==1701),7],type="b",col="red",pch=4)
lines(fei[(fei[,1]==629),7],fei[(fei[,1]==629),8],type="b",col="green",pch=2)
legend("bottomright",legend=c("Pepler","Pepler 925hPa","Fei"),pch=4,col=c("blue","red","green"))

##Ex-TC, Feb 1990
I=which(mine[,3]>=19900202 & mine[,3]<=19900205)
unique(mine[I,1])

contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(130,170),ylim=c(-50,-20),main="ex-TC, February 1990")
lines(mine[(mine[,1]==1205),6],mine[(mine[,1]==1205),7],type="b",col="blue",pch=1)
lines(mine[(mine[,1]==1206),6],mine[(mine[,1]==1206),7],type="b",col="blue",pch=2)
lines(mine[(mine[,1]==1207),6],mine[(mine[,1]==1207),7],type="b",col="blue",pch=3)
lines(mine[(mine[,1]==1208),6],mine[(mine[,1]==1208),7],type="b",col="blue",pch=4)
lines(e925[(e925[,1]==641),6],e925[(e925[,1]==641),7],type="b",col="red",pch=1)
lines(e925[(e925[,1]==642),6],e925[(e925[,1]==642),7],type="b",col="red",pch=2)
lines(fei[(fei[,1]==252),7],fei[(fei[,1]==252),8],type="b",col="green",pch=2)
lines(mldb[mldb[,1]==3,4],mldb[mldb[,1]==3,3],type="b",col="black",pch=4,cex=2)
legend("bottomright",legend=c("MLDB","Pepler","Pepler 925hPa","Fei"),pch=4,col=c("black","blue","red","green"))

####Try w/ inter-reanalysis comp.
rm(list=ls())
setwd("~/Documents/ECLs/")
read.csv('CSV/ECLfixes_unsw_NCEP1.csv')->ncep1
ncep1=ncep1[,c(-5,-6)]
read.csv('CSV/ECLfixes_unsw_NCEP2.csv')->ncep2
ncep2=ncep2[,c(-5,-6)]
read.csv('CSV/ECLfixes_unsw_JRA25.csv')->jra
jra=jra[,2:11]
read.csv('CSV/ECLfixes_unsw_ERAI.csv',sep=";")->erai
read.csv('CSV/ECLfixes_mldb.csv')->mldb
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

I=which(erai[,3]>=19981225 & erai[,3]<=19981227)
unique(erai[I,1]) 

tiff(file="26Dec98_reanalyses.tiff",height=600,width=500)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(140,170),ylim=c(-50,-20),main="Boxing Day ECL, December 1998",cex.main=2)
lines(ncep1[(ncep1[,1]==2823),6],ncep1[(ncep1[,1]==2823),7],lwd=2,col="blue",pch=1)
lines(ncep1[(ncep1[,1]==2824),6],ncep1[(ncep1[,1]==2824),7],lwd=2,col="blue",pch=2)
lines(ncep1[(ncep1[,1]==2825),6],ncep1[(ncep1[,1]==2825),7],lwd=2,col="blue",pch=3)
lines(ncep2[(ncep2[,1]==548),6],ncep2[(ncep2[,1]==548),7],lwd=2,col="red",pch=1)
lines(ncep2[(ncep2[,1]==549),6],ncep2[(ncep2[,1]==549),7],lwd=2,col="red",pch=2)
lines(jra[(jra[,1]==429),6],jra[(jra[,1]==429),7],lwd=2,col="purple",pch=1)
lines(jra[(jra[,1]==430),6],jra[(jra[,1]==430),7],lwd=2,col="purple",pch=2)
lines(erai[(erai[,1]==2262),6],erai[(erai[,1]==2262),7],lwd=2,col="orange",pch=1)
lines(erai[(erai[,1]==2263),6],erai[(erai[,1]==2263),7],lwd=2,col="orange",pch=2)
lines(erai[(erai[,1]==2264),6],erai[(erai[,1]==2264),7],lwd=2,col="orange",pch=3)
points(mldb[(mldb[,2]>=19981224 & mldb[,2]<=19981229),4],mldb[(mldb[,2]>=19981224 & mldb[,2]<=19981229),3],col="black",pch=4,cex=2)
legend("bottomright",legend=c("MLDB","NCEP1","NCEP2","JRA25","ERAI"),pch=4,col=c("black","blue","red","purple","orange"),cex=2)
dev.off()

##Mine v. Alejandro
rm(list=ls())
setwd("~/Documents/ECLs/")
read.csv('CSV/ECLfixes_unsw_NCEP1.csv')->ncep1
ncep1=ncep1[,c(-5,-6)]
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

read.csv("Algorithm Comparison/Alejandro/Alejandro_NCEP1_8009_6.csv")->ale
ale=ale[,c(1,2,23,22,18,5,4,6,8,24)]
ale[ale[,5]==-1,5]=0

I=which(ale[,3]>=19981225 & ale[,3]<=19981227)
unique(ale[I,2]) 

tiff(file="26Dec98_ale.tiff",height=600,width=500)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main="Sydney-Hobart ECL, 25-27 December 1998",cex.main=2)
lines(ncep1[(ncep1[,1]==2823),6],ncep1[(ncep1[,1]==2823),7],lwd=3,col="blue",pch=1)
lines(ncep1[(ncep1[,1]==2824),6],ncep1[(ncep1[,1]==2824),7],lwd=3,col="blue",pch=2)
lines(ncep1[(ncep1[,1]==2825),6],ncep1[(ncep1[,1]==2825),7],lwd=3,col="blue",pch=3)
lines(ale[(ale[,2]==468),6],ale[(ale[,2]==468),7],lwd=3,col="red",lty=3)
lines(ale[(ale[,2]==469),6],ale[(ale[,2]==469),7],lwd=3,col="red",lty=3)
lines(ale[(ale[,2]==468 & ale[,5]==0),6],ale[(ale[,2]==468 & ale[,5]==0),7],lwd=3,col="red",pch=1)
lines(ale[(ale[,2]==469 & ale[,5]==0),6],ale[(ale[,2]==469 & ale[,5]==0),7],lwd=3,col="red",pch=2)
points(mldb[(mldb[,2]>=19981224 & mldb[,2]<=19981229),4],mldb[(mldb[,2]>=19981224 & mldb[,2]<=19981229),3],col="black",pch=4,cex=3)
legend("bottomright",legend=c("MLDB","Pepler","Browning"),pch=4,col=c("black","blue","red"),cex=2)
dev.off()

tiff(file="7Jun07_ale.tiff",height=600,width=500)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main="Pasha Bulker, 6-8 June 2007",cex.main=2)
lines(ncep1[(ncep1[,1]==3270),6],ncep1[(ncep1[,1]==3270),7],lwd=3,col="blue",pch=1)
lines(ncep1[(ncep1[,1]==3271),6],ncep1[(ncep1[,1]==3271),7],lwd=3,col="blue",pch=1)
lines(ale[(ale[,2]==678),6],ale[(ale[,2]==678),7],lwd=3,col="red",lty=3)
lines(ale[(ale[,2]==679),6],ale[(ale[,2]==679),7],lwd=3,col="red",lty=3)
lines(ale[(ale[,2]==678 & ale[,5]==0),6],ale[(ale[,2]==678 & ale[,5]==0),7],lwd=3,col="red",pch=1)
lines(ale[(ale[,2]==679 & ale[,5]==0),6],ale[(ale[,2]==679 & ale[,5]==0),7],lwd=3,col="red",pch=2)
legend("bottomright",legend=c("Pepler","Browning"),pch=4,col=c("blue","red"),cex=2)
dev.off()

tiff(file="28Jun07_ale.tiff",height=600,width=500)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main="Gippsland low, 25-29 June 2007",cex.main=2)
lines(ncep1[(ncep1[,1]==3276),6],ncep1[(ncep1[,1]==3276),7],lwd=2,col="blue",pch=1)
lines(ncep1[(ncep1[,1]==3277),6],ncep1[(ncep1[,1]==3277),7],lwd=2,col="blue",pch=1)
lines(ale[(ale[,2]==682),6],ale[(ale[,2]==682),7],lwd=3,col="red",lty=3)
lines(ale[(ale[,2]==683),6],ale[(ale[,2]==683),7],lwd=3,col="red",lty=3)
lines(ale[(ale[,2]==684),6],ale[(ale[,2]==684),7],lwd=3,col="red",lty=3)
lines(ale[(ale[,2]==685),6],ale[(ale[,2]==685),7],lwd=3,col="red",lty=3)
lines(ale[(ale[,2]==682 & ale[,5]==0),6],ale[(ale[,2]==682 & ale[,5]==0),7],lwd=3,col="red",pch=1)
lines(ale[(ale[,2]==683 & ale[,5]==0),6],ale[(ale[,2]==683 & ale[,5]==0),7],lwd=3,col="red",pch=2)
lines(ale[(ale[,2]==684 & ale[,5]==0),6],ale[(ale[,2]==684 & ale[,5]==0),7],lwd=3,col="red",pch=1)
lines(ale[(ale[,2]==685 & ale[,5]==0),6],ale[(ale[,2]==685 & ale[,5]==0),7],lwd=3,col="red",pch=2)
legend("bottomright",legend=c("Pepler","Browning"),pch=4,col=c("blue","red"),cex=2)
dev.off()


##Looking at Ale's issues
rm(list=ls())
setwd("~/Documents/ECLs/")
read.csv('CSV/ECLfixes_mldb.csv')->mldb
read.csv('CSV/ECLfixes_unsw_NCEP1.csv')->ncep1
read.csv('CSV/ECLfixes_unsw_ERAI.csv',sep=";")->erai
ncep1=ncep1[,c(-5,-6)]
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

read.csv("Algorithm Comparison/Alejandro/Alejandro_NCEP1_8009_6.csv")->ale
ale=ale[,c(1,2,23,22,18,5,4,6,8)]
ale[ale[,5]==-1,5]=0
read.csv("Algorithm Comparison/Alejandro/ERAI-75-6-global-150-6-3_v8.csv")->aleE1
aleE1$date=aleE1[,21]+aleE1[,20]*100+aleE1[,19]*10000
aleE1=aleE1[,c(1,1,23,22,17,4,3,5,7)]
aleE1[aleE1[,5]==-1,5]=0
read.csv("Algorithm Comparison/Alejandro/ERAI-75-6-global-300-6-3_v8.csv")->aleE3
aleE3$date=aleE3[,21]+aleE3[,20]*100+aleE3[,19]*10000
aleE3=aleE3[,c(1,1,23,22,17,4,3,5,7)]
aleE3[aleE3[,5]==-1,5]=0



trackplot(19981225,19981227,"26Dec98","Sydney-Hobart, 25-27 December 1998")
trackplot(20070606,20070608,"7Jun07","Pasha Bulker, 6-8 June 2007")
trackplot(20070625,20070629,"28Jun07","Gippsland low, 25-29 June 2007")
trackplot(19920218,19920220,"19Feb92","ex-TC, 18-20 Feb 1992")
trackplot<-function(d1,d2,fname,title="")
{
tiff(file=paste(fname,"_ale_ERAIissue.tiff",sep=""),height=600,width=500)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main=title,cex.main=2)
I=which(ncep1[,3]>=d1 & ncep1[,3]<=d2)
a=unique(ncep1[I,1])
for(i in 1:length(a)) lines(ncep1[(ncep1[,1]==a[i]),6],ncep1[(ncep1[,1]==a[i]),7],lwd=3,col="blue",pch=4, type="b")
#I=which(erai[,3]>=d1 & erai[,3]<=d2)
#a=unique(erai[I,1])
#for(i in 1:length(a)) lines(erai[(erai[,1]==a[i]),6],erai[(erai[,1]==a[i]),7],lwd=3,col="purple",pch=4, type="b")
I=which(ale[,3]>=d1 & ale[,3]<=d2)
a=unique(ale[I,2])
for(i in 1:length(a)) lines(ale[(ale[,2]==a[i] & ale[,5]==0),6],ale[(ale[,2]==a[i] & ale[,5]==0),7],lwd=3,col="red",pch=4,type="b")
I=which(aleE1[,3]>=d1 & aleE1[,3]<=d2)
a=unique(aleE1[I,2])
for(i in 1:length(a)) lines(aleE1[(aleE1[,2]==a[i] & aleE1[,5]==0),6],aleE1[(aleE1[,2]==a[i] & aleE1[,5]==0),7],lwd=3,col="orange",pch=4,type="b")
I=which(aleE3[,3]>=d1 & aleE3[,3]<=d2)
a=unique(aleE3[I,2])
for(i in 1:length(a)) lines(aleE3[(aleE3[,2]==a[i] & aleE3[,5]==0),6],aleE3[(aleE3[,2]==a[i] & aleE3[,5]==0),7],lwd=3,col="green",pch=4,type="b")
#points(mldb[(mldb[,2]>=d1 & mldb[,2]<=d2),4],mldb[(mldb[,2]>=d1 & mldb[,2]<=d2),3],col="black",pch=4,cex=3)
legend("bottomright",legend=c("P NCEP1","A NCEP1","A ERAI150","A ERAI300"),pch=4,col=c("blue","red","orange","green"),cex=2)
dev.off()
}

trackplot2(19981225,19981227,"26Dec98","Sydney-Hobart, 25-27 December 1998")
trackplot2(20070606,20070608,"7Jun07","Pasha Bulker, 6-8 June 2007")
trackplot2(20070625,20070629,"28Jun07","Gippsland low, 25-29 June 2007")
trackplot2(19920218,19920220,"19Feb92","ex-TC, 18-20 Feb 1992")
trackplot2<-function(d1,d2,fname,title="")
{
  tiff(file=paste(fname,"_ale_ERAIissue2.tiff",sep=""),height=600,width=500)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main=title,cex.main=2)
  I=which(erai[,3]>=d1 & erai[,3]<=d2)
  a=unique(erai[I,1])
  for(i in 1:length(a)) lines(erai[(erai[,1]==a[i]),6],erai[(erai[,1]==a[i]),7],lwd=3,col="blue",pch=4, type="b")
  I=which(aleE1[,3]>=d1 & aleE1[,3]<=d2)
  a=unique(aleE1[I,2])
#  for(i in 1:length(a)) lines(aleE1[aleE1[,2]==a[i],6],aleE1[aleE1[,2]==a[i],7],lwd=3,col="purple",pch=4,lty=3)
  for(i in 1:length(a)) lines(aleE1[(aleE1[,2]==a[i] & aleE1[,5]==0),6],aleE1[(aleE1[,2]==a[i] & aleE1[,5]==0),7],lwd=3,col="purple",pch=4,type="b")
  I=which(aleE3[,3]>=d1 & aleE3[,3]<=d2)
  a=unique(aleE3[I,2])
#  for(i in 1:length(a)) lines(aleE3[aleE3[,2]==a[i],6],aleE3[aleE3[,2]==a[i],7],lwd=3,col="red",pch=4,lty=3)
  for(i in 1:length(a)) lines(aleE3[(aleE3[,2]==a[i] & aleE3[,5]==0),6],aleE3[(aleE3[,2]==a[i] & aleE3[,5]==0),7],lwd=3,col="red",pch=4,type="b")
  #points(mldb[(mldb[,2]>=d1 & mldb[,2]<=d2),4],mldb[(mldb[,2]>=d1 & mldb[,2]<=d2),3],col="black",pch=4,cex=3)
  legend("bottomright",legend=c("P ERAI150","A ERAI150","A ERAI300"),pch=4,col=c("blue","purple","red"),cex=2)
  dev.off()
}




rm(list=ls())
setwd("~/Documents/ECLs/")
read.csv('CSV/ECLfixes_unsw_ncep1_250.csv')->ncep1
ncep1=ncep1[,c(-5,-6)]
read.csv('CSV/ECLfixes_unsw_wrf_250.csv')->wrf
wrf=wrf[,-1]

read.csv('CSV/ECLfixes_mldb.csv')->mldb
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

d1=20070606
d2=20070608
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),cex.main=2)
I=which(ncep1[,3]>=d1 & ncep1[,3]<=d2)
a=unique(ncep1[I,1])
for(i in 1:length(a)) lines(ncep1[(ncep1[,1]==a[i]),6],ncep1[(ncep1[,1]==a[i]),7],lwd=3,col="blue",pch=4, type="b")
I=which(wrf[,3]>=d1 & wrf[,3]<=d2)
a=unique(wrf[I,1])
for(i in 1:length(a)) lines(wrf[(wrf[,1]==a[i]),6],wrf[(wrf[,1]==a[i]),7],lwd=3,col="red",pch=4, type="b")

d1=20080806
d2=20080808
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),cex.main=2)
I=which(ncep1[,3]>=d1 & ncep1[,3]<=d2)
a=unique(ncep1[I,1])
for(i in 1:length(a)) lines(ncep1[(ncep1[,1]==a[i]),6],ncep1[(ncep1[,1]==a[i]),7],lwd=3,col="blue",pch=4, type="b")
I=which(wrf[,3]>=d1 & wrf[,3]<=d2)
a=unique(wrf[I,1])
for(i in 1:length(a)) lines(wrf[(wrf[,1]==a[i]),6],wrf[(wrf[,1]==a[i]),7],lwd=3,col="red",pch=4, type="b")

trackplotW(19981225,19981227,"26Dec98","Sydney-Hobart, 25-27 December 1998")
trackplotW(20070606,20070608,"7Jun07","Pasha Bulker, 6-8 June 2007")
trackplotW(20070625,20070629,"28Jun07","Gippsland low, 25-29 June 2007")
trackplotW(19920218,19920220,"19Feb92","ex-TC, 18-20 Feb 1992")
trackplotW(19920218,19920220,"19Feb92","ex-TC, 18-20 Feb 1992")
trackplotW(20080805,20080808,"5Aug08","Missing low, 5-8 August 2008")
trackplotW<-function(d1,d2,fname,title="")
{
  tiff(file=paste(fname,"_wrf.tiff",sep=""),height=600,width=700)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,180),ylim=c(-45,-25),main=title,cex.main=2)
  I=which(ncep1$Date>=d1 & ncep1$Date<=d2)
  a=unique(ncep1$ID[I])
  for(i in 1:length(a)) lines(ncep1$Lon[ncep1$ID==a[i]],ncep1$Lat[ncep1$ID==a[i]],lwd=3,col="blue",pch=4, type="b")
  I=which(wrf$Date>=d1 & wrf$Date<=d2)
  a=unique(wrf$ID[I])
  for(i in 1:length(a)) lines(wrf$Lon[wrf$ID==a[i]],wrf$Lat[wrf$ID==a[i]],lwd=3,col="red",pch=4, type="b")
  I=which(mldb[,2]>=d1 & mldb[,2]<=d2)
  legend("bottomright",legend=c("NCEP1","NCEP-WRF"),pch=4,col=c("blue","red"),cex=2)
  dev.off()
}


####ERAI
rm(list=ls())
setwd("~/Documents/ECLs/Algorithm Comparison/")
read.csv('Mine.csv')->mine
read.csv('Mine_cv4.csv')->mine2
read.csv('Ale.csv')->Ale
read.csv('MLDB.csv')->mldb
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

trackplot2(19981225,19981227,"26Dec98","Sydney-Hobart, 25-27 December 1998")
trackplot2(20070606,20070608,"7Jun07","Pasha Bulker, 6-8 June 2007")
trackplot2(20070625,20070629,"28Jun07","Gippsland low, 25-29 June 2007")
trackplot2(19920218,19920220,"19Feb92","ex-TC, 18-20 Feb 1992")
trackplot<-function(d1,d2,fname,title="")
{
  tiff(file=paste(fname,"_ERAI_all.tiff",sep=""),height=600,width=500)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main=title,cex.main=2)
  I=which(mine[,3]>=d1 & mine[,3]<=d2)
  a=unique(mine[I,1])
  for(i in 1:length(a)) lines(mine[(mine[,1]==a[i] & mine[,9]>=0.25),6],mine[(mine[,1]==a[i] & mine[,9]>=0.25),7],lwd=3,col="red",pch=4, type="b")
  I=which(Ale[,3]>=d1 & Ale[,3]<=d2)
  a=unique(Ale[I,1])
  for(i in 1:length(a)) lines(Ale[(Ale[,1]==a[i] & Ale[,5]==-1),6],Ale[(Ale[,1]==a[i] & Ale[,5]==-1),7],lwd=3,col="orange",pch=4, type="b")
  I=which(mldb[,3]>=d1 & mldb[,3]<=d2)
  if(length(I)>0)
  {
    points(mldb[I,6],mldb[I,7],col="black",pch=4,cex=3)
    legend("bottomright",legend=c("Pepler","Browning","Speer"),pch=4,col=c("red","orange","black"),cex=2)
  } else legend("bottomright",legend=c("Pepler","Browning"),pch=4,col=c("red","orange"),cex=2)
  dev.off()
}
trackplot2<-function(d1,d2,fname,title="")
{
  tiff(file=paste(fname,"_ERAI.tiff",sep=""),height=600,width=500)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main=title,cex.main=2)
  I=which(mine2[,3]>=d1 & mine2[,3]<=d2)
  a=unique(mine2[I,1])
  for(i in 1:length(a)) lines(mine2[(mine2[,1]==a[i]),6],mine2[(mine2[,1]==a[i]),7],lwd=3,col="red",pch=4, type="b")
  I=which(Ale[,3]>=d1 & Ale[,3]<=d2)
  a=unique(Ale[I,1])
  for(i in 1:length(a)) lines(Ale[(Ale[,1]==a[i] & Ale[,5]==-1),6],Ale[(Ale[,1]==a[i] & Ale[,5]==-1),7],lwd=3,col="orange",pch=4, type="b")
  I=which(mldb[,3]>=d1 & mldb[,3]<=d2)
  if(length(I)>0)
  {
    points(mldb[I,6],mldb[I,7],col="black",pch=4,cex=3)
    legend("bottomright",legend=c("Pepler","Browning","Speer"),pch=4,col=c("red","orange","black"),cex=2)
  } else legend("bottomright",legend=c("Pepler","Browning"),pch=4,col=c("red","orange"),cex=2)
  dev.off()
}

###Try again with, say, NCEP1 and NCEP1
rm(list=ls())
setwd('~/Documents/ECLs')
ncep1<-read.csv('CSV/ECLfixes_unsw_ncep1_250.csv')
ncep1$Time=sprintf("%2.2i:00",ncep1$Time)
ncep2<-read.csv('CSV/ECLfixes_unsw_ncep2_250.csv')
erai<-read.csv('CSV/ECLfixes_unsw_erai_300.csv')
jra25<-read.csv('CSV/ECLfixes_unsw_jra_250.csv')
merra<-read.csv('CSV/ECLfixes_unsw_merra_250.csv')
cfsr<-read.csv('CSV/ECLfixes_unsw_cfsr_250.csv')
wrf<-read.csv('CSV/ECLfixes_unsw_wrf_250_ext_ext.csv')
datasets=list(ncep1,ncep2,erai,jra25,merra,cfsr)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

trackplot(19981225,19981227,"26Dec98","Sydney-Hobart, 25-27 December 1998")
trackplot(20070606,20070608,"7Jun07","Pasha Bulker, 6-8 June 2007")
trackplot(20070625,20070629,"28Jun07","Gippsland low, 25-29 June 2007")
trackplot(19920218,19920220,"19Feb92","ex-TC, 18-20 Feb 1992")

trackplot<-function(d1,d2,fname,title="")
{
  tiff(file=paste(fname,"_reanalyses.tiff",sep=""),height=600,width=500)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main=title,cex.main=2)
  cols=c("blue","red","orange","green","purple","lightblue")
  
  for(i in 1:6)
  {
    data=datasets[[i]]
    I=which(data$Date>=d1 & data$Date<=d2 & data$CV>=0.25)
    a=unique(data$ID[I])
    for(j in 1:length(a))
    {
      I=which(data$ID==a[j])
      lines(data$Lon[I],data$Lat[I],lwd=3,col=cols[i],pch=4, type="b")
    }
  }
  legend("bottomright",legend=c("NCEP1","NCEP2","ERAI","JRA25","MERRA","CFSR"),pch=4,col=cols,cex=2)
  dev.off()
}



####ERAI
rm(list=ls())
setwd("~/Documents/ECLs/Algorithm Comparison/")
read.csv('Mine_cv4.csv')->mine
read.csv('MLDB.csv')->mldb
read.csv('UM_cv35.csv')->UM
read.csv('Ale_v12.csv')->Ale
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

trackplot2(19981225,19981227,"26Dec98","Sydney-Hobart, 25-27 December 1998")
trackplot2(20070606,20070608,"7Jun07","Pasha Bulker, 6-8 June 2007")
trackplot2(20070625,20070629,"28Jun07","Gippsland low, 25-29 June 2007")
trackplot2(19920218,19920220,"19Feb92","ex-TC, 18-20 Feb 1992")
trackplot2(19860805,19860807,"6Aug86","5-7 August 1986")

trackplot2<-function(d1,d2,fname,title="")
{
  tiff(file=paste(fname,"_new.tiff",sep=""),height=600,width=500)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),main=title,cex.main=2)
  I=which(mine[,3]>=d1 & mine[,3]<=d2)
  a=unique(mine[I,1])
  for(i in 1:length(a)) lines(mine[(mine[,1]==a[i]),6],mine[(mine[,1]==a[i]),7],lwd=3,col="red",pch=4, type="b")
  I=which(UM[,4]>=d1 & UM[,4]<=d2)
  a=unique(UM[I,1])
  for(i in 1:length(a)) lines(UM[(UM[,1]==a[i]),7],UM[(UM[,1]==a[i]),8],lwd=3,col="orange",pch=4, type="b")
  I=which(Ale[,4]>=d1 & Ale[,4]<=d2)
  a=unique(Ale[I,1])
  for(i in 1:length(a)) lines(Ale[(Ale[,1]==a[i]),6],Ale[(Ale[,1]==a[i]),7],lwd=3,col="darkgreen",pch=4, type="b")
  I=which(mldb[,3]>=d1 & mldb[,3]<=d2)
  if(length(I)>0)
  {
    points(mldb[I,6],mldb[I,7],col="black",pch=4,cex=3)
    legend("bottomright",legend=c("MLD","LAPB","LAPM","PG"),pch=4,col=c("black","red","orange","darkgreen"),cex=2)
  } else legend("bottomright",legend=c("LAPB","LAPM","PG"),pch=4,col=c("red","orange","darkgreen"),cex=2)
  dev.off()
}

