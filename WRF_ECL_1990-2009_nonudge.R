rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
dirs=c("ERA-nudge","ERA-nonudge","ERA-nonudge_notopo")
cat="p100_rad2cv1"
cat2="rad2_p100_cv1.0"

events<-fixes<-list()

for(n in 1:3)
{
  events[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,"/ECLevents_",dirs[n],"_",cat2,".csv",sep=""))
  events[[n]]$Year=floor(events[[n]]$Date1/10000)
  events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
  events[[n]][events[[n]]==-Inf]=NaN
  
  fixes[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,"/ECLfixes_",dirs[n],"_",cat2,".csv",sep=""))
  fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
  fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
  fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
}

years=1980:2009
count=matrix(NaN,30,3)

for(i in 1:30)
  for(j in 1:3)
{
  I=which(events[[j]]$Year==years[i])
  count[i,j]=length(I)
  }
count[1:10,3]=NaN
apply(count[11:20,],2,mean,na.rm=T)
cor(count[11:20,2],count[11:20,3])

count2=array(NaN,c(30,12,3))
for(i in 1:30)
  for(j in 1:12)
  for(k in 1:3)
  {
    I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j)
    count2[i,j,k]=length(I)
  }

count2[c(1:10,21:30),,3]=NaN
a=apply(count2,c(2,3),mean,na.rm=T)
clist=c("red","blue","green")
plot(1:12,a[,1],type="l",lwd=2,col=clist[1],ylim=c(0,4),xlab="Month",ylab="Count")
for(i in 2:3) lines(1:12,a[,i],lwd=2,col=clist[i])
legend("topleft",c("Nudge","NoNudge","NoTopo"),col=clist,lwd=2,bty="n")


ks.test(events[[2]]$CV2,events[[3]]$CV2) # Not significantly different

cvthresh=c(seq(1,4,0.5),NaN)
cvcount=array(0,c(7,3))
for(x in 1:7)
  for(j in 1:3)
    cvcount[x,j]=length(which(events[[j]]$CV2>=cvthresh[x] & events[[j]]$CV2<cvthresh[x+1] &
                              events[[j]]$Year>=1990 & events[[j]]$Year<=1999))

makePDF(events[[2]]$CV2[events[[2]]$Year>=1990],events[[3]]$CV2)

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

loc<-array(NaN,c(11,17,3,2))

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:3)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1 &
                fixes[[k]]$Year>=1990 & fixes[[k]]$Year<=1999)
      loc[i,j,k,1]=length(I)/10
      loc[i,j,k,2]=mean(fixes[[k]]$CV[I],na.rm=T)
    }

### Plot where ECLs are
pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

pdf(file=paste(figdir,"ECL_location_19901999_NoNudgevNoTopo.pdf",sep=""),width=8.5,height=4)
bb2=c(-10000,seq(0,20,length.out=11),10000)
cm=pal1(12)
layout(cbind(1,2,3),c(1,1,0.3))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,2]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="R2 NoNudge",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
image(lon,lat,t(loc[,,3]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="R2 NoNudge NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()


pdf(file=paste(figdir,"ECL_location_19901999_NoTopoChange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(100*((loc[,,3,1]/loc[,,2,1])-1)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      cex.axis=1,cex.main=1)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"ECL_location_19901999_NoTopo_CVchange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-0.2,0.2,length.out=9),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(loc[,,3,2]-loc[,,2,2]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      cex.axis=1,cex.main=1)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()