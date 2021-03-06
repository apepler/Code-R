rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
dirs=c("ERA-nudge","ERA-nonudge","ERA-nonudge_notopo")
cat="p100_rad2cv1"
cat2="rad2_p100_cv1.0"
dom="d01"

events<-fixes<-list()

for(n in 1:3)
{
  events[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,"/",dom,"/ECLevents_",dirs[n],"_",cat2,"_typing_impactsC2.csv",sep=""))
  events[[n]]$Year=floor(events[[n]]$Date1/10000)
  events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
  events[[n]][events[[n]]==-Inf]=NaN
  
  fixes[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,"/",dom,"/ECLfixes_",dirs[n],"_",cat2,"_typing_impactsC2.csv",sep=""))
  fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
  fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
  fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
}

years=1990:2009
count=matrix(NaN,20,3)

for(i in 1:20)
  for(j in 1:3)
{
  I=which(events[[j]]$Year==years[i])
  count[i,j]=length(I)
  }
apply(count,2,mean,na.rm=T)
cor(count[,2],count[,3])

count2=array(NaN,c(20,12,3))
for(i in 1:20)
  for(j in 1:12)
  for(k in 1:3)
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

sigM=rep(0,7)
for(i in 1:7){
  a=t.test(count3[,2,i],count3[,3,i])
  sigM[i]=a$p.value }
ks.test(events[[2]]$CV2,events[[3]]$CV2) # Not significantly different

cvthresh=c(seq(1,4,0.5),NaN)
cvcount=array(0,c(7,3))
for(x in 1:7)
  for(j in 1:3)
    cvcount[x,j]=length(which(events[[j]]$CV2>=cvthresh[x] & events[[j]]$CV2<cvthresh[x+1]))


######### Matching

match=eventmatch(events[[2]][events[[2]]$Year>=1990,],fixes[[2]][fixes[[2]]$Year>=1990,],events[[3]],fixes[[3]])
match2=match[,6:8]-match[,1:3]
apply(match2,2,mean,na.rm=T)
plot(match[,1],match2[,1])
abline(h=0,col="red")

###### ERAI comp

erai_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
erai_E$Year=floor(erai_E$Date1/10000)
erai_E$Month=floor(erai_E$Date1/100) %% 100
erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")
erai$Year=floor(erai$Date/10000)
erai$Month=floor(erai$Date/100) %% 100
erai$Date2=as.POSIXct(paste(as.character(erai$Date),substr(erai$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

years=1990:2009
erai_ann=rep(0,20)
for(i in 1:20) erai_ann[i]=length(which(erai_E$Year==years[i]))

eracomp=matrix(0,4,3)
rownames(eracomp)=c("Cor","HR","FAR","CSI")
colnames(eracomp)=c("Nudge","No-nudge","No-nudge notopo")
for(i in 1:3) 
  {
  eracomp[1,i]=cor(count[11:30,i],erai_ann)
  eracomp[2:4,i]=CSI_days(erai[erai$Year>=1990,],fixes[[i]][fixes[[i]]$Year>=1990,])
}

#### Change vs location

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

loc<-array(NaN,c(11,17,3,3))

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:3)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1)
      loc[i,j,k,1]=length(I)/20
      loc[i,j,k,2]=mean(fixes[[k]]$CV[I],na.rm=T)
      
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Fix==1)
      loc[i,j,k,3]=length(I)/20
    }

### Plot where ECLs are
pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo.pdf",sep=""),width=8.5,height=4)
bb2=c(-10000,seq(0,20,length.out=11),10000)
cm=pal1(12)
layout(cbind(1,2,3),c(1,1,0.3))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,2,1]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="R2 NoNudge",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
image(lon,lat,t(loc[,,3,1]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="R2 NoNudge NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()


pdf(file=paste(figdir,"ECL_location_NoTopoChange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(100*((loc[,,3,1]/loc[,,2,1])-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      cex.axis=1,cex.main=1)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_NoTopo_CVchange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-0.2,0.2,length.out=9),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(loc[,,3,2]-loc[,,2,2]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      cex.axis=1,cex.main=1)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_Genesis.pdf",sep=""),width=7,height=8)
bb2=c(-10000,seq(0,2.5,length.out=11),10000)
cm=pal1(12)
layout(cbind(c(1,2),c(3,3)),width=c(1,0.2))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,2,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 NoNudge",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
image(lon,lat,t(loc[,,3,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 NoNudge NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_Genesis_change.pdf",sep=""),width=7,height=4)
bb2=c(-10000,seq(-1,1,length.out=11),10000)
cm=pal(12)
layout(cbind(1,2),c(1,0.2))
par(mar=c(2,2,2,0))
image(lon,lat,t(loc[,,3,3]-loc[,,2,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      cex.axis=1,cex.main=1)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()





##### Significant things

count=array(NaN,c(3,13))
dimnames(count)[[1]]=c("NoNudge","NoTopo")
dimnames(count)[[2]]=c("All","CV>=2","Bombs","Formed in region","Formed elsewhere","EC","SC","Mixed","Mean rain>6 mm/6hr","Mean rain>10 mm/6hr",
                       "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")

for(i in 1:3)
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

for(n in 1:3) fixesALL[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,"/",dom,"/ECLfixes_",dirs[n],"_",cat2,"_ALL.csv",sep=""))

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

locALL<-array(NaN,c(11,17,3,2))

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:3)
    {
      I=which(fixesALL[[k]]$Lat>=lat[i]-2.5 & fixesALL[[k]]$Lat<lat[i]+2.5 & fixesALL[[k]]$Lon>=lon[j]-2.5 & fixesALL[[k]]$Lon<lon[j]+2.5)
      J=unique(fixesALL[[k]]$ID[I])
      
      locALL[i,j,k,1]=length(I)/20
      locALL[i,j,k,2]=length(which(fixesALL[[k]]$Fix[I]==1))/20
    }

pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

library(maps)
for(i in 1:3)
{
pdf(file=paste(figdir,"ECL_location_",dirs[i],"_ALL.pdf",sep=""),width=7,height=4)
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

for(i in 1:3)
{
  pdf(file=paste(figdir,"ECL_location_",dirs[i],"_ALL_genesis.pdf",sep=""),width=7,height=4)
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

type=c("ALL","ALL_genesis")

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_NoNudgevNoTopo_Change_",type[i],".pdf",sep=""),width=7,height=4)
  bb2=c(-10000,seq(-100,100,20),10000)
  cm=pal(12)
  layout(cbind(1,2),width=c(1,0.2))
  par(mar=c(3,3,3,1))
  print(range(locALL[,,3,i]/locALL[,,2,i],na.rm=T))
  image(lon,lat,t(100*((locALL[,,3,i]/locALL[,,2,i])-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
        main="% Change",cex.axis=1.5,cex.main=1.5)
  map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}
