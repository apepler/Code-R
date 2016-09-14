rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
library(maps)

cat="p100_rad2cv1"
cat2="rad2_p100_cv1.0"
events<-fixes<-list()

dirs=c("ERA-nudge","ERA-nonudge","ERA-nonudge_notopo","ERA-nudge","ERA-nonudge","ERA-nonudge_notopo")
dom=c("d01","d01","d01","d02","d02","d02")
tmp=c("_2","_2","","_2","_2","")

for(n in 1:6)
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
count=matrix(NaN,20,6)

for(i in 1:20)
  for(j in 1:6)
  {
    I=which(events[[j]]$Year==years[i])
    count[i,j]=length(I)
  }

apply(count,2,mean,na.rm=T)
t.test(count[,1],count[,2])

count2=array(NaN,c(20,12,6))
for(i in 1:20)
  for(j in 1:12)
    for(k in 1:6)
    {
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j)
      count2[i,j,k]=length(I)
    }

tmp=apply(count2,c(2,3),mean)
plot(1:12,tmp[,1],lwd=2,type="l",ylim=c(0,3.5))
for(i in 2:3) lines(1:12,tmp[,i],lwd=2,col=i)
for(i in 1:3) lines(1:12,tmp[,i+3],lwd=2,col=i,lty=2)

library(abind)
count3=abind(apply(count2[,,],c(1,3),sum),apply(count2[,3:5,],c(1,3),sum),apply(count2[,6:8,],c(1,3),sum),
             apply(count2[,9:11,],c(1,3),sum),apply(count2[,c(12,1:2),],c(1,3),sum),
             apply(count2[,5:10,],c(1,3),sum),apply(count2[,5:10,],c(1,3),sum),along=3)
count3[1:19,,5]=apply(abind(count2[1:19,12,],count2[2:20,1:2,],along=2),c(1,3),sum)
count3[1:19,,7]=apply(abind(count2[1:19,11:12,],count2[2:20,1:4,],along=2),c(1,3),sum)
count3[20,,c(5,7)]=NaN

sigM=matrix(0,2,7)
for(i in 1:7){
  a=t.test(count3[,2,i],count3[,3,i])
  sigM[1,i]=a$p.value
  a=t.test(count3[,5,i],count3[,6,i])
  sigM[2,i]=a$p.value  
  }
ks.test(events[[1]]$CV2,events[[2]]$CV2) # Not significantly different

cvthresh=c(seq(1,4,0.5),NaN)
cvcount=array(0,c(7,6))
for(x in 1:7)
  for(j in 1:6)
    cvcount[x,j]=length(which(events[[j]]$CV2>=cvthresh[x] & events[[j]]$CV2<cvthresh[x+1]))


########## Repeating for CV>=2

count2=array(NaN,c(20,12,6))
for(i in 1:20)
  for(j in 1:12)
    for(k in 1:6)
    {
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j & events[[k]]$CV2>=2)
      count2[i,j,k]=length(I)
    }

tmp=apply(count2,c(2,3),mean)

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

loc<-array(NaN,c(11,17,6,3))

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:6)
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
image(lon,lat,t(loc[,,5,1]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="R2 NoNudge",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
image(lon,lat,t(loc[,,6,1]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="R2 NoNudge NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()


pdf(file=paste(figdir,"ECL_location_d02_NoTopoChange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(100*((loc[,,6,1]/loc[,,5,1])-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      cex.axis=1,cex.main=1)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_d02_NoTopo_CVchange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-0.2,0.2,length.out=9),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(loc[,,6,2]-loc[,,5,2]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      cex.axis=1,cex.main=1)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_d02_NoNudgevNoTopo_Genesis.pdf",sep=""),width=7,height=8)
bb2=c(-10000,seq(0,2.5,length.out=11),10000)
cm=pal1(12)
layout(cbind(c(1,2),c(3,3)),width=c(1,0.2))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,5,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 NoNudge",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
image(lon,lat,t(loc[,,6,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 NoNudge NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_d02_NoNudgevNoTopo_Genesis_change.pdf",sep=""),width=7,height=4)
bb2=c(-10000,seq(-1,1,length.out=11),10000)
cm=pal(12)
layout(cbind(1,2),c(1,0.2))
par(mar=c(2,2,2,0))
image(lon,lat,t(loc[,,6,3]-loc[,,5,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      cex.axis=1,cex.main=1)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()



##### Significant things

count=array(NaN,c(6,13))
dimnames(count)[[1]]=c("NoNudge d01","NoTopo d01","NoNudge d02","NoTopo d02")
dimnames(count)[[2]]=c("All","CV>=2","Bombs","Formed in region","Formed elsewhere","EC","SC","Mixed","Mean rain>6 mm/6hr","Mean rain>10 mm/6hr",
                       "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")

for(i in 1:6)
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

count3=array(NaN,c(6,8))
dimnames(count3)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(count3)[[2]]=c("All","CV>=2.5","Bombs","Mean rain>6 mm/6hr","Mean rain>9 mm/6hr",
                        "Max rain>50 mm/6hr","Mean wind> 40 km/hr","Max wind> 80 km/hr")

  for(i in 1:6)
  {
    ev=fixes[[i]][fixes[[i]]$Location2>0,] 
    
    count3[i,1]=length(ev$CV)
    count3[i,2]=length(which(ev$CV>=2))
    count3[i,3]=length(which(ev$NDR>=1))
    count3[i,4]=length(which(ev$MeanRain500>=6))
    count3[i,5]=length(which(ev$MeanRain500>=9))
    count3[i,6]=length(which(ev$MaxRain500>=50))
    count3[i,7]=length(which(ev$MeanWind500>=11.1))
    count3[i,8]=length(which(ev$MaxWind500>=22.2))
  }





########## What about locations of all cyclones, and locations of cyclone genesis?

fixesALL<-list()
tmp=c("_2","_2","")
for(n in 1:3) fixesALL[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLfixes_",dirs[n],"_",cat2,"_ALL.csv",sep=""))
for(n in 1:3) fixesALL[[n]]$Year=floor(fixesALL[[n]]$Date/10000)

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

locALL<-array(NaN,c(20,11,17,3,2))

for(y in 1:20)
for(i in 1:11)
  for(j in 1:17)
    for(k in 1:3)
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
  image(lon,lat,t(100*((apply(locALL[,,,3,i],c(2,3),sum)/apply(locALL[,,,2,i],c(2,3),sum))-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main="% Change",cex.axis=1.5,cex.main=1.5)
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_Change_",type[i],"_sig3.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(-50,50,10),10000)
  cm=pal(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  tmp=t(100*((apply(locALL[,,,3,i],c(2,3),sum)/apply(locALL[,,,2,i],c(2,3),sum))-1))
  I=which(t(apply(locALL[,,,2,i],c(2,3),mean))<1)
  tmp[I]=NaN
  tmp[which(tmp>=50)]=49
  tmp[which(tmp<=-50)]=-49  
  image(lon,lat,tmp,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main="% Change",cex.axis=1.5,cex.main=1.5)
  
  pval<-matrix(0,length(lon),length(lat))
  for(y in 1:length(lat))
    for(x in 1:length(lon))
    {
      a=t.test(locALL[,y,x,3,i],locALL[,y,x,2,i])
      pval[x,y]=a$p.value
    }
  
  sigmask=which(pval<=0.05 & !is.na(tmp),arr.ind=T)
  points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=4,lwd=3,xlim=c(110,175),ylim=c(-45,-10))
  
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}



lat=seq(-60,-10,10)
lon=seq(100,180,10)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

locALL<-array(NaN,c(20,6,9,3,2))

for(y in 1:20)
  for(i in 1:6)
    for(j in 1:9)
      for(k in 1:3)
      {
        I=which(fixesALL[[k]]$Lat>=lat[i]-5 & fixesALL[[k]]$Lat<lat[i]+5 & fixesALL[[k]]$Lon>=lon[j]-5 & fixesALL[[k]]$Lon<lon[j]+5 & fixesALL[[k]]$Year==years[y])
        J=unique(fixesALL[[k]]$ID[I])
        
        locALL[y,i,j,k,1]=length(I)
        locALL[y,i,j,k,2]=length(which(fixesALL[[k]]$Fix[I]==1))
      }

pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

type=c("ALL2","ALL2_genesis")

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_Change_",type[i],"_sig2.pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(-100,100,20),10000)
  cm=pal(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  image(lon,lat,t(100*((apply(locALL[,,,3,i],c(2,3),sum)/apply(locALL[,,,2,i],c(2,3),sum))-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main="% Change",cex.axis=1.5,cex.main=1.5)
  
  pval<-matrix(0,length(lat),length(lon))
  for(x in 1:length(lat))
    for(y in 1:length(lon))
    {
      a=t.test(locALL[,x,y,3,i],locALL[,x,y,2,i])
      pval[x,y]=a$p.value
    }
  
  sigmask=which(pval<=0.05,arr.ind=T)
  points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,lwd=3,xlim=c(110,175),ylim=c(-45,-10))
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}




####### Test - Stuart's typing

types=c("ET","SSL","IT","CL")
count4<-array(0,c(20,2,4))
for(i in 1:2)
  for(j in 1:4)
    for(y in 1:20)
    count4[y,i,j]=length(which(events[[i+1]]$TypeSB==types[j] & events[[i+1]]$Year==years[y]))

for(i in 1:4) print(t.test(count4[,1,i],count4[,2,i]))

#### Repeating that box plot thing

count=array(NaN,c(6,5,4))
dimnames(count)[[1]]=c("ERA-nudge d01","ERA-nonudge d01","ERA-nonudge_notopo d01","ERA-nudge d02","ERA-nonudge d02","ERA-nonudge_notopo d02")
cvthresh=c(1,1.5,2,2.5,3)
dimnames(count)[[2]]=cvthresh
dimnames(count)[[3]]=c("Events","Fixes","Days","CoastDays")

for(i in 1:6)
  for(j in 1:5)
  {
    count[i,j,1]=length(which(events[[i]]$CV2>=cvthresh[j]))
    count[i,j,2]=length(which(fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh[j]))
    count[i,j,3]=length(unique(fixes[[i]]$Date[fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh[j]]))
    count[i,j,4]=length(unique(fixes[[i]]$Date[fixes[[i]]$Location2==1 & fixes[[i]]$CV>=cvthresh[j]]))
  }

