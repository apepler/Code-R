rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
library(maps)

cat="p100_rad2cv0.5"
cat2="rad2_p100_cvthresh1"
events<-fixes<-list()

dirs=c("ERA-nonudge","ERA-nonudge_notopo","ERA-nonudge","ERA-nonudge_notopo")
dom=c("d01","d01","d02","d02")
tmp=c("_2","","_2","")

for(n in 1:4)
{
  events[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLevents_",dirs[n],"_",dom[n],"_",cat2,"_typing.csv",sep=""))
  events[[n]]$Year=floor(events[[n]]$Date1/10000)
  events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
  events[[n]][events[[n]]==-Inf]=NaN
  
  fixes[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLfixes_",dirs[n],"_",dom[n],"_",cat2,"_typing.csv",sep=""))
  fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
  fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
  fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
}

years=1990:2009
count=matrix(NaN,20,6)

for(i in 1:20)
  for(j in 1:4)
  {
    I=which(events[[j]]$Year==years[i])
    count[i,j]=length(I)
  }

apply(count,2,mean,na.rm=T)
t.test(count[,1],count[,2])

count2=array(NaN,c(20,12,4))
for(i in 1:20)
  for(j in 1:12)
    for(k in 1:4)
    {
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j)
      count2[i,j,k]=length(I)
    }

library(abind)
count3=abind(apply(count2[,,],c(1,3),sum),apply(count2[,3:5,],c(1,3),sum),apply(count2[,6:8,],c(1,3),sum),
             apply(count2[,9:11,],c(1,3),sum),apply(count2[,c(12,1:2),],c(1,3),sum),
             apply(count2[,5:10,],c(1,3),sum),apply(count2[,5:10,],c(1,3),sum),along=3)
count3[1:19,,5]=apply(abind(count2[1:19,12,],count2[2:20,1:2,],along=2),c(1,3),sum)
count3[1:19,,7]=apply(abind(count2[1:19,11:12,],count2[2:20,1:4,],along=2),c(1,3),sum)
count3[20,,c(5,7)]=NaN

sigM=matrix(0,2,7)
for(i in 1:7){
  a=t.test(count3[,1,i],count3[,2,i])
  sigM[1,i]=a$p.value
  a=t.test(count3[,3,i],count3[,4,i])
  sigM[2,i]=a$p.value  
  }


#### Change vs location

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

loc<-array(NaN,c(11,17,4,3))

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:4)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1)
      loc[i,j,k,1]=length(I)/20
      loc[i,j,k,2]=mean(fixes[[k]]$CV[I],na.rm=T)
      
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Fix==1)
      loc[i,j,k,3]=length(I)/20
    }

### Plot where ECLs are
pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

pdf(file=paste(figdir,"ECL_location_d01_NoNudgevNoTopo_Genesis_cv0.5.pdf",sep=""),width=7,height=8)
bb2=c(-10000,seq(0,2.5,length.out=11),10000)
cm=pal1(12)
layout(cbind(c(1,2),c(3,3)),width=c(1,0.2))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,2,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 NoNudge",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
image(lon,lat,t(loc[,,1,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 NoNudge NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

####### Significance

loc<-array(NaN,c(20,11,17,4,3))

for(y in 1:20)
for(i in 1:11)
  for(j in 1:17)
    for(k in 1:4)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1 & fixes[[k]]$Year==years[y])
      loc[y,i,j,k,1]=length(I)
      loc[y,i,j,k,2]=mean(fixes[[k]]$CV[I],na.rm=T)
      
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Fix==1 & fixes[[k]]$Year==years[y])
      loc[y,i,j,k,3]=length(I)
    }

type=c("_",NA,"_genesis_")

for(i in 3)
{
  pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_d01_Change",type[i],"cv0.5_sig3.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(-50,50,10),10000)
  cm=pal(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  tmp=t(100*((apply(loc[,,,2,i],c(2,3),sum)/apply(loc[,,,1,i],c(2,3),sum))-1))
  I=which(t(apply(loc[,,,1,i],c(2,3),mean))<0.5)
  tmp[I]=NaN
  image(lon,lat,tmp,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main="% Change",cex.axis=1.5,cex.main=1.5)
  
  pval<-matrix(0,length(lon),length(lat))
  for(y in 1:length(lat))
    for(x in 1:length(lon))
    {
      a=t.test(loc[,y,x,2,i],loc[,y,x,1,i])
      pval[x,y]=a$p.value
    }
  
  sigmask=which(pval<=0.05 & !is.na(tmp),arr.ind=T)
  points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=4,lwd=3,xlim=c(110,175),ylim=c(-45,-10))
  
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}




##### Significant things

count=array(NaN,c(4,5))
dimnames(count)[[1]]=c("NoNudge d01","NoTopo d01","NoNudge d02","NoTopo d02")
dimnames(count)[[2]]=c("All","CV>=2","Bombs","Formed in region","Formed elsewhere")

for(i in 1:4)
{
  
  count[i,1]=length(events[[i]]$CV2)
  count[i,2]=length(which(events[[i]]$CV2>=2))
  count[i,3]=sum(events[[i]]$Bomb)
  count[i,4]=length(which(events[[i]]$EnteredFormed==1))
  count[i,5]=length(which(events[[i]]$EnteredFormed!=1))
}



########## What about locations of all cyclones, and locations of cyclone genesis?
cat2="rad2_p100_cv0.5"
fixesALL<-list()
tmp=c("_2","")
for(n in 1:2) fixesALL[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLfixes_",dirs[n],"_d01_",cat2,"_ALL.csv",sep=""))
for(n in 1:2) fixesALL[[n]]$Year=floor(fixesALL[[n]]$Date/10000)

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

locALL<-array(NaN,c(20,11,17,2,2))

for(y in 1:20)
for(i in 1:11)
  for(j in 1:17)
    for(k in 1:2)
    {
      I=which(fixesALL[[k]]$Lat>=lat[i]-2.5 & fixesALL[[k]]$Lat<lat[i]+2.5 & fixesALL[[k]]$Lon>=lon[j]-2.5 & fixesALL[[k]]$Lon<lon[j]+2.5 & fixesALL[[k]]$Year==years[y])
      J=unique(fixesALL[[k]]$ID[I])
      
      locALL[y,i,j,k,1]=length(I)
      locALL[y,i,j,k,2]=length(which(fixesALL[[k]]$Fix[I]==1))
    }

pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

library(maps)
for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_",dirs[i],"_ALL_cv0.5.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(0,50,length.out=11),10000)
  cm=pal1(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  image(lon,lat,t(apply(locALL[,,,i,1],c(2,3),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main=paste("R2",dirs[i]),cex.axis=1.5,cex.main=1.5)
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_",dirs[i],"_ALL_cv0.5_genesis.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(0,10,length.out=8),10000)
  cm=pal1(9)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  image(lon,lat,t(apply(locALL[,,,i,2],c(2,3),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main=paste("R2",dirs[i]),cex.axis=1.5,cex.main=1.5)
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}


type=c("ALL2","ALL2_genesis")

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_Change_",type[i],"_cv0.5.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(-100,100,20),10000)
  cm=pal(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  image(lon,lat,t(100*((apply(locALL[,,,2,i],c(2,3),sum)/apply(locALL[,,,1,i],c(2,3),sum))-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main="% Change",cex.axis=1.5,cex.main=1.5)
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_Change_",type[i],"_cv0.5_sig3.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(-50,50,10),10000)
  cm=pal(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  tmp=t(100*((apply(locALL[,,,2,i],c(2,3),sum)/apply(locALL[,,,1,i],c(2,3),sum))-1))
  I=which(t(apply(locALL[,,,1,i],c(2,3),mean))<1)
  tmp[I]=NaN
  image(lon,lat,tmp,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main="% Change",cex.axis=1.5,cex.main=1.5)
  
  pval<-matrix(0,length(lon),length(lat))
  for(y in 1:length(lat))
    for(x in 1:length(lon))
    {
      a=t.test(locALL[,y,x,2,i],locALL[,y,x,1,i])
      pval[x,y]=a$p.value
    }
  
  sigmask=which(pval<=0.05 & !is.na(tmp),arr.ind=T)
  points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=4,lwd=3,xlim=c(110,175),ylim=c(-45,-10))
  
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}



