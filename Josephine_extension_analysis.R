###########
## Step 1 - COnvert Power databse of events into a dataset of ECL/event by day
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/20CR/ECLs/")
load("20CR_ECLdates_v3.RData")
load("20CR_ECLdates_ens.RData")

## Decadal distribution of # ECLs/day

yy=floor(dates2[,1]/10000)
decadalfreq2=matrix(0,14,6)
rownames(decadalfreq2)<-decades<-seq(1871,2001,10)
threshes<-c(1,14,28,42,50,56)
colnames(decadalfreq2)<-c("Any","25%","50%","75%","90%","100%")

for(i in 1:length(decades))
  for(j in 1:length(threshes))
  {
    if(j<length(threshes)) I=which(yy>=decades[i] & yy<=decades[i]+9 & dates2[,3]>=threshes[j] & dates2[,3]<threshes[j+1]) else
      I=which(yy>=decades[i] & yy<=decades[i]+9 & dates2[,3]>=threshes[j])
    decadalfreq2[i,j]=length(I)
  }

## Annual distribution
yy=floor(dates2[,1]/10000)
years=unique(yy)
annfreq2=matrix(0,length(years),7)
rownames(annfreq2)<-years
colnames(annfreq2)<-threshes<-c(1,14,28,32,42,50,56)

for(i in 1:length(years))
  for(j in 1:length(threshes))
  {
    if(j<length(threshes)) I=which(yy==years[i] & dates2[,3]>=threshes[j] & dates2[,3]<=threshes[j+1]) else
      I=which(yy==years[i] & dates2[,3]>=threshes[j])
    annfreq2[i,j]=length(I)
  }

annfreq=cbind(annfreq,apply(annfreq[,1:56],1,mean))

years2=seq(1871,2012,1)
corr<-matrix(NaN,58,58)

for(i in 1:58)
  for(j in 1:58)
    if(i!=j)
  {
    I=which(years>=1911 & years<=1959) 
    corr[i,j]=cor(annfreq[I,i],annfreq[I,j])
    }

mean(corr[1:56,1:56],na.rm=T)


###### During 1911-1920, how consistent are ensemble & ensemble mean?
I=which(dates2[,1]>=19110000 & dates2[,1]<=19200000)
dates3=dates2[I,]

I=which(dates3[,3]>=50)
mean(dates3[I,4])

I=which(dates3[,4]>=1)
quantile(dates3[I,3],c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))

H=length(which(dates3[,3]>=42 & dates3[,4]==1))
FA=length(which(dates3[,3]>=42 & dates3[,4]==0))
M=length(which(dates3[,3]<42 & dates3[,4]==1))

H/(H+M)
FA/(H+FA)
H/(H+M+FA)

### Playing with a random distribution

test=matrix(0,365,10000)
for(j in 1:10000)
for(i in 1:56) 
  {
  tmp=sample(1:365,33)
  test[tmp,j]=test[tmp,j]+1
}
decadalfreq2=rbind(rep(0,length(threshes)),decadalfreq2)
rownames(decadalfreq2)[1]="Uniform distribution"
for(j in 1:length(threshes))
{
  if(j<length(threshes)) I=which(test>=threshes[j] & test<=threshes[j+1]) else I=which(test>=threshes[j])
  decadalfreq2[1,j]=length(I)/1000
}

## 

##### Replicate Josephine figures

## Figure 2

pdf(file="Boxplot_decadal_vmean.pdf",height=4,width=7)
boxplot(t(decadalfreq[,1:56]),ylim=c(0,225),axes=F,bty="o",col="gray",xlab="Decade",ylab="Count of ECLs")
points(1:14,decadalfreq[,57],pch=4,cex=1,lwd=3)
axis(2,cex.axis=0.7)
axis(1,at=seq(-0.5,14.5,1),labels=c(NA,decades,NA),cex.axis=0.7)
dev.off()


um=read.csv("events_UM_rad2_p100.csv")
umdec=rep(NaN,length(decades))
for(i in 1:length(decades))
{
  I=which(um$Date1>=decades[i]*10000 & um$Date2<(decades[i]+10)*10000)
  if(length(I)>0) umdec[i]=length(unique(um$ID[I]))
}

umann=cbind(seq(1980,2009),rep(0,30))
for(i in 1:30) umann[i,2]=length(which(um$Year==umann[i,1]))

pdf(file="Boxplot_decadal_vmeanERAI_title.pdf",height=4,width=7)
boxplot(t(decadalfreq[,1:56]),ylim=c(0,250),axes=F,bty="o",col="gray",xlab="Decade",ylab="Count of ECLs",main="Figure 2")
points(1:14,decadalfreq[,57],pch=4,cex=1,lwd=3)
axis(2,at=seq(-50,250,50),cex.axis=0.7)
axis(1,at=seq(-0.5,14.5,1),labels=c(NA,decades,2011),cex.axis=0.7)
points(12:14,umdec[12:14],pch=17,col="grey")
dev.off()

I=which(years>=1980 & years<=2009)
ave=apply(annfreq[I,],2,mean)

## Figure 3

pdf(file="Barplot_decadal_title.pdf",height=4,width=7)
barplot(t(decadalfreq2[,6:1]/10),legend.text=T,axisnames=F,xlab="Decade",ylab="Number of ECL days per year",
        ylim=c(0,350),axes=F,main="Figure 3")
axis(1,at=c(-50,50))
axis(1,at=seq(0.05,17.2,1.205),labels=seq(1871,2011,10),cex.axis=0.7)
axis(2,at=seq(0,400,50),cex.axis=0.7)
dev.off()



##### Hit rate/FA stuff

I=which(dates2[,1]>=19800000 & dates2[,1]<=20090000)
dates3=dates2[I,]

thresh=matrix(0,57,5)
thresh[,1]=c(1:56,0)
colnames(thresh)=c("Thresh","HR","FAR","CSI","Count")

for(i in 1:56)
{
  thresh[i,5]=length(which(dates3[,3]>=i))
  H=length(which(dates3[,6]==1 & dates3[,3]>=i))
  M=length(which(dates3[,6]==1 & dates3[,3]<i))
  FA=length(which(dates3[,6]==0 & dates3[,3]>=i))
  
  thresh[i,2]=H/(H+M)
  thresh[i,3]=FA/(H+FA)
  thresh[i,4]=H/(H+M+FA)
}

H=length(which(dates3[,6]==1 & dates3[,4]==1))
M=length(which(dates3[,6]==1 & dates3[,4]==0))
FA=length(which(dates3[,6]==0 & dates3[,4]==1))

thresh[57,2]=H/(H+M)
thresh[57,3]=FA/(H+FA)
thresh[57,4]=H/(H+M+FA)


###### Decadal stuff - amount by all members
apply(decadalfreq2[11:15,],2,mean)/10
apply(decadalfreq2[,2:7],1,sum)/10

I=which(dates2[,1]>=19800000 & dates2[,1]<=20070000)
apply(dates2[I,2:6],2,sum)/27

apply(decadal)

#### Testing the ECL C events

I=which((data$Type=="ECL C" | data$Type=="ECL C L" |  data$Type=="ECL NC") & data$Date>=18710000)
data3=data[I,]


I=which(data3$EnsMean==1 | data3$EnsMem>=50)

###### Doing some locations of events - 1890 & 1990

for(m in 1:56)
{
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_1911-1920_",m,".csv",sep=""))
  I=which(tmp$Date>=19130000 & tmp$Date<=19140000)
  if(m==1) fix1900=tmp[I,] else fix1900=rbind(fix1900,tmp[I,])
  
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_1991-2000_",m,".csv",sep=""))
  I=which(tmp$Date>=19980000 & tmp$Date<=19990000)
  if(m==1) fix2000=tmp[I,] else fix2000=rbind(fix2000,tmp[I,])
}

lat=-50:-20
lon=140:175
locs<-array(0,c(length(lat),length(lon),2))

for(i in 1:length(lat))
  for(j in 1:length(lon))
  {
    I=which(fix1900$Lat>=lat[i]-0.5 & fix1900$Lat<lat[i]+0.5 & fix1900$Lon>=lon[j]-0.5 & fix1900$Lon<lon[j]+0.5)
    locs[i,j,1]=length(I)
    I=which(fix2000$Lat>=lat[i]-0.5 & fix2000$Lat<lat[i]+0.5 & fix2000$Lon>=lon[j]-0.5 & fix2000$Lon<lon[j]+0.5)
    locs[i,j,2]=length(I)
  }

locs=locs/56

library(fields)
library("R.matlab")
readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

ColorBar <- function(brks,cols,labels=NA)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  if(is.na(labels[1])) labels=brks[seq(2, length(brks)-1)]
  
  axis(4, at = seq(1.5, length(brks) - 1.5), tick = TRUE, 
       labels = labels)
}

cm=gray(seq(1,0.1,-0.15))
bb=c(-0.5,0,0.25,0.5,1,1.5,2,100)

layout(cbind(1,2,3),width=c(1,1,0.3))
image(lon,lat,t(locs[,,1]),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(140,175),ylim=c(-50,-20),main="1913",cex.axis=1.5,cex.main=2)
map(xlim=c(140,175),ylim=c(-50,-20),add=T,lwd=2)
image(lon,lat,t(locs[,,2]),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(140,175),ylim=c(-50,-20),main="1998",cex.axis=1.5,cex.main=2)
map(xlim=c(140,175),ylim=c(-50,-20),add=T,lwd=2)
ColorBar(bb,cm)

library(maps)
map(xlim=c(140,175),ylim=c(-50,-20),add=T)

### Alternative to Josephine - look at individual events?

lat=-50:-20
lon=140:175
library(maps)

pdf(file="ECL_location_events_v3_title.pdf",width=12,height=4.5)
layout(rbind(cbind(4,4,4),cbind(1,2,3)),width=c(1,1,1),height=c(0.5,4))
par(mar=c(5,5,4,1)+0.1)
image(lon,lat,matrix(NA,length(lon),length(lat)),
      xlab=expression(paste("Longitude [",degree,"E]")),ylab=expression(paste("Latitude [",degree,"N]")),cex.lab=1.5,
      breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(140,175),ylim=c(-50,-20),main="13-15 May 1913",cex.axis=1.5,cex.main=1.5)
map(xlim=c(140,175),ylim=c(-50,-20),lwd=2,add=T,col="gray")
for(m in 1:56)
{
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_1911-1920_",m,".csv",sep=""))
  I=which(tmp$Date>=19130513 & tmp$Date<=19130515)
  a=unique(tmp$ID[I])
  for(i in 1:length(a))
  {
    I=which(tmp$ID==a[i])
    lines(tmp$Lon[I],tmp$Lat[I])
  }
}

image(lon,lat,matrix(NA,length(lon),length(lat)),xlab=expression(paste("Longitude [",degree,"E]")),cex.lab=1.5,ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(140,175),ylim=c(-50,-20),main="25-27 February 1919",cex.axis=1.5,cex.main=1.5)
map(xlim=c(140,175),ylim=c(-50,-20),lwd=2,add=T,col="gray")
for(m in 1:56)
{
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_1911-1920_",m,".csv",sep=""))
  I=which(tmp$Date>=19190225 & tmp$Date<=19190227)
  a=unique(tmp$ID[I])
  for(i in 1:length(a))
  {
    I=which(tmp$ID==a[i])
    lines(tmp$Lon[I],tmp$Lat[I])
  }
}

image(lon,lat,matrix(NA,length(lon),length(lat)),xlab=expression(paste("Longitude [",degree,"E]")),cex.lab=1.5,ylab="",breaks=bb,col=cm,zlim=c(-1000,1000),xlim=c(140,175),ylim=c(-50,-20),main="7-8 August 1998",cex.axis=1.5,cex.main=1.5)
map(xlim=c(140,175),ylim=c(-50,-20),lwd=2,add=T,col="gray")
for(m in 1:56)
{
  tmp=read.csv(paste("/srv/ccrc/data34/z3478332/20CR/ECLs/output/ensemble_56mem/tracks/ECLfixes_1991-2000_",m,".csv",sep=""))
  I=which(tmp$Date>=19980807 & tmp$Date<=19980808)
  a=unique(tmp$ID[I])
  for(i in 1:length(a))
  {
    I=which(tmp$ID==a[i])
    lines(tmp$Lon[I],tmp$Lat[I])
  }
}

title("Figure 4", line = -3, outer = TRUE,cex.main=3)
dev.off()


#####

### Hit rate/FA for each ens member vs the ensemble mean
I=which(dates2[,1]>=19110000 & dates2[,1]<=20120000)
dates3=dates2[I,]
I=which(dates[,1]>=19110000 & dates[,1]<=20120000)
dates4=dates[I,]

thresh=matrix(0,57,5)
thresh[,1]=c(1:56,0)
colnames(thresh)=c("Thresh","HR","FAR","CSI","Count")

for(i in 1:56)
{
  H=length(which(dates3[,2]==1 & dates4[,i+1]==1))
  M=length(which(dates3[,2]==1 & dates4[,i+1]==0))
  FA=length(which(dates3[,2]==0 & dates4[,i+1]==1))
  
  thresh[i,2]=H/(H+M)
  thresh[i,3]=FA/(H+FA)
  thresh[i,4]=H/(H+M+FA)
}

H=length(which(dates3[,2]==1 & dates3[,4]==1))
M=length(which(dates3[,2]==1 & dates3[,4]==0))
FA=length(which(dates3[,2]==0 & dates3[,4]==1))

thresh[57,2]=H/(H+M)
thresh[57,3]=FA/(H+FA)
thresh[57,4]=H/(H+M+FA)

####Make a new decadalfreq

thresh=c(1,1.25,1.5,2)
decades=seq(1871,2001,10)
decadalfreq2<-array(0,c(length(decades),56,length(thresh)))

for(i in 1:length(decades))
  for(j in 1:56)
    for(k in 1:length(thresh))
      decadalfreq2[i,j,k]=length(which(cv[,j+1]>=thresh[k] & cv[,1]>=decades[i]*10000 & cv[,1]<(decades[i]+10)*10000))

dimnames(decadalfreq2)[[1]]=decades
dimnames(decadalfreq2)[[3]]=thresh

boxplot(t(decadalfreq[,1:56]),ylim=c(0,max(decadalfreq[,1:56])),bty="o",col="gray",xlab="Decade",ylab="Count of ECLs")
for(t in 1:4) 
  {
  png(file=paste("Boxplot_decadal_cv",thresh[t],".png",sep=""),height=400,width=600)
  boxplot(t(decadalfreq2[,,t]),ylim=c(0,max(decadalfreq2[,,t])),bty="o",col="gray",xlab="Decade",ylab="Count of ECLs",main=paste("Thresh",thresh[t]))
dev.off()
}

source("20CR_ECLdates_thresh.RData")
thresh=c(1,1.25,1.5,2)
decades=seq(1871,2001,10)
decadalfreq2<-array(0,c(length(decades),56,length(thresh)))

for(t in 1:4) 
{
  png(file=paste("Boxplot_decadal_cv",thresh[t],".png",sep=""),height=400,width=700)
  boxplot(t(decadalfreq3[,1:56,t]),ylim=c(0,max(decadalfreq3[,1:56,t])),axes=F,bty="o",col="gray",xlab="Decade",ylab="Count of ECLs",main=paste("Intensity >",thresh[t],"hpa/(deg.lat^2)"))
  axis(2)
  axis(1,at=seq(-0.5,14.5,1),labels=c(NA,decades,NA),lwd.ticks=3)
  dev.off()
}