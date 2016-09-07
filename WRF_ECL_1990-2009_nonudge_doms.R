rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
library(maps)

cat="p100_rad2cv1"
cat2="rad2_p100_cv1.0"
events<-fixes<-list()

dirs=c("ERA-nonudge","ERA-nonudge_notopo","ERA-nonudge","ERA-nonudge_notopo")
dom=c("d01","d01","d02","d02")
tmp=c("_2","","_2","")

for(n in 1:4)
{
  events[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLevents_",dirs[n],"_",dom[n],"_",cat2,"_typing_impactsC2.csv",sep=""))
  events[[n]]$Year=floor(events[[n]]$Date1/10000)
  events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
  events[[n]][events[[n]]==-Inf]=NaN
  
  fixes[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLfixes_",dirs[n],"_",dom[n],"_",cat2,"_typing_impactsC2.csv",sep=""))
  fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
  fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
  fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
}

years=1990:2009
count=matrix(NaN,20,4)

for(i in 1:20)
  for(j in 1:4)
  {
    I=which(events[[j]]$Year==years[i])
    count[i,j]=length(I)
  }

apply(count,2,mean,na.rm=T)
cor(count[,1],count[,2])
t.test(count[,2],count[,1])

count2=array(NaN,c(20,12,4))
for(i in 1:20)
  for(j in 1:12)
    for(k in 1:4)
    {
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j)
      count2[i,j,k]=length(I)
    }

tmp=apply(count2,c(2,3),mean)
plot(1:12,tmp[,1],lwd=2,type="l",ylim=c(0,3.5))
for(i in 2:4) lines(1:12,tmp[,i],lwd=2,col=i)

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
ks.test(events[[1]]$CV2,events[[2]]$CV2) # Not significantly different

cvthresh=c(seq(1,4,0.5),NaN)
cvcount=array(0,c(7,4))
for(x in 1:7)
  for(j in 1:4)
    cvcount[x,j]=length(which(events[[j]]$CV2>=cvthresh[x] & events[[j]]$CV2<cvthresh[x+1]))


########## Repeating for CV>=2

count2=array(NaN,c(20,12,4))
for(i in 1:20)
  for(j in 1:12)
    for(k in 1:4)
    {
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j & events[[k]]$CV2>=2)
      count2[i,j,k]=length(I)
    }

tmp=apply(count2,c(2,3),mean)
plot(1:12,tmp[,1],lwd=2,type="l",ylim=c(0,1))
for(i in 2:4) lines(1:12,tmp[,i],lwd=2,col=i)

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

pdf(file=paste(figdir,"ECL_location_d02_NoNudgevNoTopo.pdf",sep=""),width=8.5,height=4)
bb2=c(-10000,seq(0,20,length.out=11),10000)
cm=pal1(12)
layout(cbind(1,2,3),c(1,1,0.3))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,3,1]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="R2 NoNudge",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
image(lon,lat,t(loc[,,4,1]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="R2 NoNudge NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()


pdf(file=paste(figdir,"ECL_location_d02_NoTopoChange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(100*((loc[,,4,1]/loc[,,3,1])-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      cex.axis=1,cex.main=1)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_d02_NoTopo_CVchange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-0.2,0.2,length.out=9),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(loc[,,4,2]-loc[,,3,2]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      cex.axis=1,cex.main=1)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_d02_NoNudgevNoTopo_Genesis.pdf",sep=""),width=7,height=8)
bb2=c(-10000,seq(0,2.5,length.out=11),10000)
cm=pal1(12)
layout(cbind(c(1,2),c(3,3)),width=c(1,0.2))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,3,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 NoNudge",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
image(lon,lat,t(loc[,,4,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 NoNudge NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_d02_NoNudgevNoTopo_Genesis_change.pdf",sep=""),width=7,height=4)
bb2=c(-10000,seq(-1,1,length.out=11),10000)
cm=pal(12)
layout(cbind(1,2),c(1,0.2))
par(mar=c(2,2,2,0))
image(lon,lat,t(loc[,,4,3]-loc[,,3,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      cex.axis=1,cex.main=1)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()



##### Significant things

count=array(NaN,c(4,13))
dimnames(count)[[1]]=c("NoNudge d01","NoTopo d01","NoNudge d02","NoTopo d02")
dimnames(count)[[2]]=c("All","CV>=2","Bombs","Formed in region","Formed elsewhere","EC","SC","Mixed","Mean rain>6 mm/6hr","Mean rain>10 mm/6hr",
                       "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")

for(i in 1:4)
{
  
  count[i,1]=length(events[[i]]$CV2)
  count[i,2]=length(which(events[[i]]$CV2>=2))
  count[i,3]=sum(events[[i]]$Bomb)
  count[i,4]=length(which(events[[i]]$EnteredFormed==1))
  count[i,5]=length(which(events[[i]]$EnteredFormed!=1))
  count[i,6]=length(which(events[[i]]$HartType=="EC"))
  count[i,7]=length(which(events[[i]]$HartType=="SC"))
  count[i,8]=length(which(events[[i]]$HartType=="Mixed"))
  count[i,9]=length(which(events[[i]]$MaxMeanRain500>=6))
  count[i,10]=length(which(events[[i]]$MaxMeanRain500>=10))
  count[i,11]=length(which(events[[i]]$MaxPointRain500>=50))
  count[i,12]=length(which(events[[i]]$MaxMeanWind500>=13.9))
  count[i,13]=length(which(events[[i]]$MaxPointWind500>=22.2))
}



########## What about locations of all cyclones, and locations of cyclone genesis?

fixesALL<-list()
tmp=c("_2","","_2","")
for(n in 1:2) fixesALL[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLfixes_",dirs[n],"_",cat2,"_ALL.csv",sep=""))

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

locALL<-array(NaN,c(11,17,2,2))

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:2)
    {
      I=which(fixesALL[[k]]$Lat>=lat[i]-2.5 & fixesALL[[k]]$Lat<lat[i]+2.5 & fixesALL[[k]]$Lon>=lon[j]-2.5 & fixesALL[[k]]$Lon<lon[j]+2.5)
      J=unique(fixesALL[[k]]$ID[I])
      
      locALL[i,j,k,1]=length(I)/20
      locALL[i,j,k,2]=length(which(fixesALL[[k]]$Fix[I]==1))/20
    }

pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

library(maps)
for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_",dirs[i],"_ALL2.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(0,25,length.out=11),10000)
  cm=pal1(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  image(lon,lat,t(locALL[,,i,1]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main=paste("R2",dirs[i]),cex.axis=1.5,cex.main=1.5)
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_",dirs[i],"_ALL2_genesis.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(0,3.5,length.out=8),10000)
  cm=pal1(9)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  image(lon,lat,t(locALL[,,i,2]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main=paste("R2",dirs[i]),cex.axis=1.5,cex.main=1.5)
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}

type=c("ALL2","ALL2_genesis")

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_Change_",type[i],".pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(-100,100,20),10000)
  cm=pal(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  print(range(locALL[,,2,i]/locALL[,,1,i],na.rm=T))
  image(lon,lat,t(100*((locALL[,,2,i]/locALL[,,1,i])-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main="% Change",cex.axis=1.5,cex.main=1.5)
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}

