rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
library(maps)

cat="p100_rad2cv1"
cat2="rad2_p100_cv1.0"
events<-fixes<-list()

dirs=c("ERA-nonudge","ERA-nudge","ERA-nonudge","ERA-nudge","ERA-nonudge_notopo","ERA-nudge_notopo","ERA-nonudge_notopo","ERA-nudge_notopo")
dom=c("d01","d01","d02","d02","d01","d01","d02","d02")
tmp=c("_2","_2","_2","_2","","","","")

years=1990:2009

for(n in 1:7)
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

###d02 doesn't exist for nudge_notopo - make an empty set
fixes[[8]]<-fixes[[7]][1:2,]
I=which(fixes[[8]]$Year==1990)
fixes[[8]]=fixes[[8]][-I,]
events[[8]]<-events[[7]][1:2,]
I=which(events[[8]]$Year==1990)
events[[8]]=events[[8]][-I,]

###Ad ERAI for ref

n=9
events[[n]]<-read.csv("outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100_update.csv")
events[[n]]$Year=floor(events[[n]]$Date1/10000)
events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
events[[n]][events[[n]]==-Inf]=NaN
events[[n]]<-events[[n]][events[[n]]$Year>=1990 & events[[n]]$Year<=2009,]

fixes[[n]]<-read.csv("outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100_update.csv")
fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
fixes[[n]]<-fixes[[n]][fixes[[n]]$Year>=1990 & fixes[[n]]$Year<=2009,]

##Quick check

years=1990:2009
count=matrix(NaN,20,9)
colnames(count)=c(paste(dirs,dom),"ERAI")

cvthresh=1

for(i in 1:20)
  for(j in 1:9)
  {
    I=which(events[[j]]$Year==years[i] & events[[j]]$CV2>=cvthresh)
    count[i,j]=length(I)
  }

apply(count,2,mean,na.rm=T)
apply(count[,5:8],2,mean,na.rm=T)/apply(count[,1:4],2,mean,na.rm=T)

### Now, doing the proper assessment of seasonal change/significance

cvthresh=1
count2=array(NaN,c(20,12,9))
for(i in 1:20)
  for(j in 1:12)
    for(k in 1:9)
    {
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j & events[[k]]$CV2>=cvthresh)
      count2[i,j,k]=length(I)
    }

library(abind)
count3=abind(apply(count2[,,],c(1,3),sum),apply(count2[,3:5,],c(1,3),sum),apply(count2[,6:8,],c(1,3),sum),
             apply(count2[,9:11,],c(1,3),sum),apply(count2[,c(12,1:2),],c(1,3),sum),
             apply(count2[,5:10,],c(1,3),sum),apply(count2[,5:10,],c(1,3),sum),along=3)
count3[1:19,,5]=apply(abind(count2[1:19,12,],count2[2:20,1:2,],along=2),c(1,3),sum)
count3[1:19,,7]=apply(abind(count2[1:19,11:12,],count2[2:20,1:4,],along=2),c(1,3),sum)
count3[20,,c(5,7)]=NaN
dimnames(count3)[[3]]=c("Ann","MAM","JJA","SON","DJF","MJJASO","NDJFMA")
dimnames(count3)[[2]]=c(paste(dirs,dom),"ERAI")


## Significant difference - seasonal

sigM=array(0,c(7,7))
dimnames(sigM)[[2]]=c("Ann","MAM","JJA","SON","DJF","MJJASO","NDJFMA")
dimnames(sigM)[[1]]=c("Count","Change d01","p d01","Change nudge","p nudge","Change d02","p d02")

for(i in 1:7)
  {
  sigM[1,i]=mean(count3[,9,i],na.rm=T)
  n=2
  for(j in 1:3)
  {
  sigM[n,i]=100*((mean(count3[,j+4,i],na.rm=T)/mean(count3[,j,i],na.rm=T))-1)
  a=t.test(count3[,j,i],count3[,j+4,i])
  sigM[n+1,i]=a$p.value
  n=n+2
  }
}

## Near coast etc

count=array(NaN,c(9,5,5))
dimnames(count)[[1]]=c(paste(dirs,dom),"ERAI")
cvthresh=c(1,1.5,2,2.5,3)
dimnames(count)[[2]]=cvthresh
dimnames(count)[[3]]=c("Events","Fixes","CoastFixes","Days","CoastDays")

for(i in 1:9)
  for(j in 1:5)
  {
    count[i,j,1]=length(which(events[[i]]$CV2>=cvthresh[j])) 
    count[i,j,2]=length(which(fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh[j])) 
    count[i,j,3]=length(which(fixes[[i]]$Location2==1 & fixes[[i]]$CV>=cvthresh[j])) 
    count[i,j,4]=length(unique(fixes[[i]]$Date[fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh[j]]))
    count[i,j,5]=length(unique(fixes[[i]]$Date[fixes[[i]]$Location2==1 & fixes[[i]]$CV>=cvthresh[j]]))
  }

change=100*((count[5:7,,]/count[1:3,,])-1)

## Change in intensity

for(j in 1:3) print(ks.test(events[[j]]$CV2,events[[j+4]]$CV2)) # Not significantly different
for(j in 1:3) print(ks.test(events[[j]]$Length2,events[[j+4]]$Length2)) # Not significantly different
for(j in 1:3) print(t.test(events[[j]]$Rad2,events[[j+4]]$Rad2)) # Not significantly different
for(j in 1:3) makePDF(events[[j]]$MSLP2,events[[j+4]]$MSLP2) # Not significantly different

cvthresh=c(seq(1,2.5,0.5),Inf)
cvcount=array(0,c(4,9))
colnames(cvcount)=c(paste(dirs,dom),"ERAI")
for(x in 1:4)
  for(j in 1:9)
    cvcount[x,j]=length(which(events[[j]]$CV2>=cvthresh[x] & events[[j]]$CV2<cvthresh[x+1]))

cvcount[,5:7]/cvcount[,1:3]

mthresh=c(Inf,seq(1015,985,-5),-Inf)
mcount=array(0,c(8,9))
colnames(mcount)=c(paste(dirs,dom),"ERAI")
rownames(mcount)=c(">1015",mthresh[3:8],"<985")
for(x in 1:8)
  for(j in 1:9)
    mcount[x,j]=length(which(events[[j]]$MSLP2<=mthresh[x] & events[[j]]$MSLP2>mthresh[x+1]))

mcount[,5:7]/mcount[,1:3]

#### Change vs location

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region

loc<-array(NaN,c(20,11,17,9))

for(y in 1:20)
for(i in 1:11)
  for(j in 1:17)
    for(k in 1:9)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1 & fixes[[k]]$Year==years[y])
      loc[y,i,j,k,1]=length(I)
      loc[y,i,j,k,2]=mean(fixes[[k]]$CV[I],na.rm=T)
      
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Fix==1 & fixes[[k]]$Year==years[y])
      loc[y,i,j,k,3]=length(I)
    }


change=100*((loc[,,,5:7,]/loc[,,,1:3,])-1)
dimnames(change)[[4]]<-nam<-c("NoNudge_d01","Nudge_d01","NoNudge_d02")
dimnames(change)[[5]]=c("Count","CV","Genesis","MeanRain")

## Figure 5

for(i in 1:3)
{
pdf(file=paste(figdir,"ECL_location_NoTopoChange_",nam[i],".pdf",sep=""),width=4.5,height=4)
bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
layout(cbind(1,2),width=c(1,0.35))
par(mar=c(3,3,3,1))

tmp=t(100*((apply(loc[,,,i+4,1],c(2,3),sum)/apply(loc[,,,i,1],c(2,3),sum))-1))
I=which(t(apply(loc[,,,i,1],c(2,3),mean))<0.5)
tmp[I]=NaN
image(lon,lat,tmp,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(145,161),ylim=c(-42,-23),
      main="ECL Frequency",cex.axis=1,cex.main=1)
pval<-matrix(0,length(lon),length(lat))
for(y in 1:length(lat))
  for(x in 1:length(lon))
  {
    a=t.test(loc[,y,x,i+4,1],loc[,y,x,i,1])
    pval[x,y]=a$p.value
  }
sigmask=which(pval<=0.05 & !is.na(tmp),arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=4,lwd=2,cex=2,xlim=c(145,161),ylim=c(-42,-23))
map(xlim=c(145,161),ylim=c(-42,-23),add=T,lwd=2)
box()
ColorBar(bb2,cm)
dev.off()
}

### Change in rainfall

loc<-array(NaN,c(11,17,2,5))

for(i in 1:11)
    for(j in 1:17)
      for(k in 1:2)
      {
        I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1)
        J=which(fixes[[k+4]]$Lat>=lat[i]-2.5 & fixes[[k+4]]$Lat<lat[i]+2.5 & fixes[[k+4]]$Lon>=lon[j]-2.5 & fixes[[k+4]]$Lon<lon[j]+2.5 & fixes[[k+4]]$Location==1)
        
        loc[i,j,k,1]=length(I)/20
        loc[i,j,k,2]=mean(fixes[[k]]$MeanRain500[I],na.rm=T)
        loc[i,j,k,3]=mean(fixes[[k+4]]$MeanRain500[J],na.rm=T)
        loc[i,j,k,4]=100*((loc[i,j,k,3]/loc[i,j,k,2])-1)
        
        if(length(I)>5 & length(J)>5)
        {
        a=t.test(fixes[[k]]$MeanRain500[I],fixes[[k+4]]$MeanRain500[J],na.rm=T)
        loc[i,j,k,5]=a$p.value
        }
      }

for(i in 1:2)
{
  pdf(file=paste(figdir,"ECL_location_NoTopoChange_",nam[i],"_MeanRain500.pdf",sep=""),width=4.5,height=4)
  bb2=c(-10000,seq(-40,40,10),10000)
  cm=pal(10)
  layout(cbind(1,2),width=c(1,0.35))
  par(mar=c(3,3,3,1))
  
  tmp=loc[,,i,4]
  I=which(loc[,,i,1]<0.5)
  tmp[I]=NaN
  image(lon,lat,t(tmp),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(145,161),ylim=c(-42,-23),
        main="Change in MeanRain500 (%)",cex.axis=1,cex.main=1)
  sigmask=which(loc[,,i,5]<=0.05 & !is.na(tmp),arr.ind=T)
  points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,lwd=2,cex=2,xlim=c(145,161),ylim=c(-42,-23))
  map(xlim=c(145,161),ylim=c(-42,-23),add=T,lwd=2)
  box()
  ColorBar(bb2,cm)
  dev.off()
}





#### Change by location

count4=array(NaN,c(20,9,5))
dimnames(count4)[[2]]=c(paste(dirs,dom),"ERAI")
dimnames(count4)[[3]]=c("All","Lat>=33","Lat>=-35","Lat< -35","Lat < -37")

for(y in 1:20)
  for(i in 1:9)
  {
    ev=fixes[[i]][fixes[[i]]$Location>0 & fixes[[i]]$Year==years[y],] 
    
    count4[y,i,1]=length(ev$CV)
    count4[y,i,2]=length(which(ev$Lat>=-33))
    count4[y,i,3]=length(which(ev$Lat>=-35))
    count4[y,i,4]=length(which(ev$Lat<(-35)))
    count4[y,i,5]=length(which(ev$Lat<(-37)))
  }
apply(count4[,5:7,],c(2,3),mean)/apply(count4[,1:3,],c(2,3),mean)

##### Change in extremes or subgroups of ECLs, e.g. formed w/in v w/out domain

count=array(NaN,c(20,9,13))
dimnames(count)[[2]]=c(paste(dirs,dom),"ERAI")
dimnames(count)[[3]]=c("All","CV>=2","Bombs","Formed in region","Formed elsewhere","EC","SC","Mixed","Mean rain>6 mm/6hr","Mean rain>10 mm/6hr",
                       "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")

for(y in 1:20)
for(i in 1:9)
{
  Y=which(events[[i]]$Year==years[y])
  count[y,i,1]=length(events[[i]]$CV2[Y])
  count[y,i,2]=length(which(events[[i]]$CV2[Y]>=2))
  count[y,i,3]=sum(events[[i]]$Bomb[Y])
  count[y,i,4]=length(which(events[[i]]$EnteredFormed[Y]==1))
  count[y,i,5]=length(which(events[[i]]$EnteredFormed[Y]!=1))
  count[y,i,6]=length(which(events[[i]]$HartType[Y]=="EC"))
  count[y,i,7]=length(which(events[[i]]$HartType[Y]=="SC"))
  count[y,i,8]=length(which(events[[i]]$HartType[Y]=="Mixed"))
  count[y,i,9]=length(which(events[[i]]$MaxMeanRain500[Y]>=6))
  count[y,i,10]=length(which(events[[i]]$MaxMeanRain500[Y]>=10))
  count[y,i,11]=length(which(events[[i]]$MaxPointRain500[Y]>=50))
  count[y,i,12]=length(which(events[[i]]$MaxMeanWind500[Y]>=13.9))
  count[y,i,13]=length(which(events[[i]]$MaxPointWind500[Y]>=22.2))
}

apply(count[,5:7,],c(2,3),mean)/apply(count[,1:3,],c(2,3),mean)

### Changes in rainfall

cvthresh1=2
cvthresh2=Inf
count3=array(NaN,c(7,9))
dimnames(count3)[[1]]=c(paste(dirs[1:7],dom[1:7]))
dimnames(count3)[[2]]=c("Count","Bombs","Mean Rain","Max Rain","Mean rain>6 mm/6hr","Mean rain>9 mm/6hr",
                        "Max rain>50 mm/6hr","Mean wind> 40 km/hr","Max wind> 80 km/hr")

  for(i in 1:7)
  {
    ev=fixes[[i]][fixes[[i]]$Location2>0 & fixes[[i]]$CV>=cvthresh1 & fixes[[i]]$CV<cvthresh2,] 
    
    count3[i,1]=length(ev$CV)
    count3[i,2]=length(which(ev$NDR>=1))
    count3[i,3]=mean(ev$MeanRain500,na.rm=T)
    count3[i,4]=mean(ev$MaxRain500,na.rm=T)
    count3[i,5]=length(which(ev$MeanRain500>=6))
    count3[i,6]=length(which(ev$MeanRain500>=9))
    count3[i,7]=length(which(ev$MaxRain500>=50))
    count3[i,8]=length(which(ev$MeanWind500>=11.1))
    count3[i,9]=length(which(ev$MaxWind500>=22.2))
  }
count3[5:7,]/count3[1:3,]

t.test(fixes[[2]]$MaxRain500[fixes[[2]]$Location2>0],fixes[[6]]$MaxRain500[fixes[[6]]$Location2>0])





### Genesis locations - for weaker precursor lows

dirs=c("ERA-nonudge","ERA-nudge","ERA-nonudge_notopo","ERA-nudge_notopo")
dom=c("d01","d01","d01","d01")
tmp=c("_2","_2","","")

cat="p100_rad2cv0.5"
cat2="rad2_p100_cvthresh1"

eventsW<-fixesW<-list()
for(n in 1:4)
{
  eventsW[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLevents_",dirs[n],"_",dom[n],"_",cat2,"_typing.csv",sep=""))
  eventsW[[n]]$Year=floor(eventsW[[n]]$Date1/10000)
  eventsW[[n]]$Month=floor(eventsW[[n]]$Date1/100)%%100
  eventsW[[n]][eventsW[[n]]==-Inf]=NaN
  
  fixesW[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,tmp[n],"/",dom[n],"/ECLfixes_",dirs[n],"_",dom[n],"_",cat2,"_typing.csv",sep=""))
  fixesW[[n]]$Year=floor(fixesW[[n]]$Date/10000)
  fixesW[[n]]$Month=floor(fixesW[[n]]$Date/100)%%100
  fixesW[[n]]$Date2=as.POSIXct(paste(as.character(fixesW[[n]]$Date),substr(fixesW[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
}

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

loc<-array(NaN,c(20,11,17,4,3))

for(y in 1:20)
  for(i in 1:11)
    for(j in 1:17)
      for(k in 1:4)
      {
        I=which(fixesW[[k]]$Lat>=lat[i]-2.5 & fixesW[[k]]$Lat<lat[i]+2.5 & fixesW[[k]]$Lon>=lon[j]-2.5 & fixesW[[k]]$Lon<lon[j]+2.5 & fixesW[[k]]$Location==1 & fixesW[[k]]$Year==years[y])
        loc[y,i,j,k,1]=length(I)
        loc[y,i,j,k,2]=mean(fixesW[[k]]$CV[I],na.rm=T)
        
        I=which(fixesW[[k]]$Lat>=lat[i]-2.5 & fixesW[[k]]$Lat<lat[i]+2.5 & fixesW[[k]]$Lon>=lon[j]-2.5 & fixesW[[k]]$Lon<lon[j]+2.5 & fixesW[[k]]$Fix==1 & fixesW[[k]]$Year==years[y])
        loc[y,i,j,k,3]=length(I) ## Genesis locations
      }

## Figure 6

pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

pdf(file=paste(figdir,"ECL_location_d01_NudgevNoTopo_Genesis_cv0.5_panel.pdf",sep=""),width=14,height=4)
bb1=c(-10000,seq(0,2.5,length.out=11),10000)
cm1=pal1(12)
layout(cbind(1,3,2,4),width=c(1,0.2,1,0.2))
par(mar=c(3,3,3,1))
image(lon,lat,t(apply(loc[,,,2,3],c(2,3),mean)),xlab="",ylab="",breaks=bb1,col=cm1,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="R2 Nudge Cyclogenesis",cex.axis=1.5,cex.main=1.5)
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)

bb2=c(-10000,seq(-0.5,0.5,length.out=11),10000)
cm2=pal(12)
tmp=t(apply(loc[,,,4,3],c(2,3),mean)-apply(loc[,,,2,3],c(2,3),mean))
I=which(t(apply(loc[,,,2,3],c(2,3),mean))<0.25)
tmp[I]=NaN
image(lon,lat,tmp,xlab="",ylab="",breaks=bb2,col=cm2,zlim=c(-Inf,Inf),xlim=c(110,175),ylim=c(-45,-10),
      main="Change - NoTopography",cex.axis=1.5,cex.main=1.5)
pval<-matrix(0,length(lon),length(lat))
for(y in 1:length(lat))
  for(x in 1:length(lon))
  {
    a=t.test(loc[,y,x,2,1],loc[,y,x,1,1])
    pval[x,y]=a$p.value
  }
sigmask=which(pval<=0.05 & !is.na(tmp),arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=4,lwd=2,cex=2,xlim=c(110,175),ylim=c(-45,-10))
map(xlim=c(110,175),ylim=c(-45,-10),add=T,lwd=2)

ColorBar(bb1,cm1)
ColorBar(bb2,cm2)
dev.off()

### Entered/Formed v season

Wnum=c(1,2,NA,NA,3,4)

cvthresh=1
count2=array(NaN,c(20,12,6,4))
for(i in 1:20)
  for(j in 1:12)
    for(k in 1:6)
    {
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j & events[[k]]$CV2>=cvthresh & events[[k]]$EnteredFormed==1)
      count2[i,j,k,1]=length(I)
      
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j & events[[k]]$CV2>=cvthresh & events[[k]]$EnteredFormed>1)
      count2[i,j,k,2]=length(I)      
      
      if(!is.na(Wnum[k]))
      {
      I=which(eventsW[[Wnum[k]]]$Year==years[i] & eventsW[[Wnum[k]]]$Month==j & eventsW[[Wnum[k]]]$CV2>=cvthresh & eventsW[[Wnum[k]]]$EnteredFormed==1)
      count2[i,j,k,3]=length(I)      
      I=which(eventsW[[Wnum[k]]]$Year==years[i] & eventsW[[Wnum[k]]]$Month==j & eventsW[[Wnum[k]]]$CV2>=cvthresh & eventsW[[Wnum[k]]]$EnteredFormed>1)
      count2[i,j,k,4]=length(I)            
      }
    }

count2=abind(count2[,,,3],count2[,,,1]-count2[,,,3],count2[,,,2],apply(count2[,,,1:2],c(1,2,3),sum),along=4)

library(abind)
count3=abind(apply(count2,c(1,3,4),sum),apply(count2[,3:5,,],c(1,3,4),sum),apply(count2[,6:8,,],c(1,3,4),sum),
             apply(count2[,9:11,,],c(1,3,4),sum),apply(count2[,c(12,1:2),,],c(1,3,4),sum),
             apply(count2[,5:10,,],c(1,3,4),sum),apply(count2[,5:10,,],c(1,3,4),sum),along=4)
count3[1:19,,,5]=apply(abind(count2[1:19,12,,],count2[2:20,1:2,,],along=2),c(1,3,4),sum)
count3[1:19,,,7]=apply(abind(count2[1:19,11:12,,],count2[2:20,1:4,,],along=2),c(1,3,4),sum)
count3[20,,,c(5,7)]=NaN
dimnames(count3)[[4]]=c("Ann","MAM","JJA","SON","DJF","MJJASO","NDJFMA")

## Significant difference - seasonal

sigM=array(0,c(5,7,4))
dimnames(sigM)[[3]]=c("Formed 0.5","Entered 0.5","Entered 1","All")
dimnames(sigM)[[2]]=c("Ann","MAM","JJA","SON","DJF","MJJASO","NDJFMA")
dimnames(sigM)[[1]]=c("Count","Change NN","p NN","Change N","p N")
for(i in 1:7)
  for(e in 1:4)
  {
  sigM[1,i,e]=mean(count3[,1,e,i],na.rm=T)
  sigM[2,i,e]=100*((mean(count3[,5,e,i],na.rm=T)/mean(count3[,1,e,i],na.rm=T))-1)
  a=t.test(count3[,1,e,i],count3[,5,e,i])
  sigM[3,i,e]=a$p.value
  sigM[4,i,e]=100*((mean(count3[,6,e,i],na.rm=T)/mean(count3[,2,e,i],na.rm=T))-1)
  a=t.test(count3[,6,e,i],count3[,2,e,i])
  sigM[5,i,e]=a$p.value  
}


#### FIXED UP TO HERE
#### Another quick thing - Dowdy stuff comparison

wrfdir=c("/srv/ccrc/data34/z3478332/WRF/ERA-nonudge/","/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nonudging_notopo/out/impact/")
GV=matrix(NaN,20*365.25*4,2)

for(n in 1:2)
{
tmp=read.csv(paste(wrfdir[n],"GV_6hrly_timeseries.txt",sep=""),header=F) ## Same for d01 or d02, because of region
len1=dim(GV)[1]
len2=dim(tmp)[1]

for(i in 5:(len1-4)) {
  if(len1==len2) j=i else j=len2-len1+i
  GV[i,n]=mean(tmp[(j-4):(j+4),1])
}
}

apply(GV,2,quantile,0.9,na.rm=T)


makePDF(GV[,1],GV[,2],xlabel="500hPa GV",leg=c("NoNudge","NoNudge_NoTopo"))

###Change in various types?

tlistSB=c("ET","IT","CL","SSL")
tlistH=c("EC","SC","TC","Mixed")

countT=array(0,c(20,12,6,9))
dimnames(countT)[[4]]=c("All",tlistH,tlistSB)

cvthresh=1
count2=array(NaN,c(20,12,6))
for(i in 1:20)
  for(j in 1:12)
    for(k in 1:6)
    {
      I=which(events[[k]]$Year==years[i] & events[[k]]$Month==j & events[[k]]$CV2>=cvthresh)
      countT[i,j,k,1]=length(I)
      
      for(t in 1:4)
      {
        countT[i,j,k,t+1]=length(which(events[[k]]$HartType[I]==tlistH[t]))
        countT[i,j,k,t+5]=length(which(events[[k]]$TypeSB[I]==tlistSB[t]))
      }
    }

monthC=apply(countT,c(2,3,4),sum)
yearC=apply(countT,c(1,3,4),sum)

sigM=array(0,c(5,9))
dimnames(sigM)[[2]]=dimnames(countT)[[4]]
dimnames(sigM)[[1]]=c("Count","Change d01","p d01","Change d02","p d02")
for(i in 1:9)
  {
    sigM[1,i]=mean(yearC[,2,i],na.rm=T)
    sigM[2,i]=100*((mean(yearC[,3,i],na.rm=T)/mean(yearC[,2,i],na.rm=T))-1)
    a=t.test(yearC[,3,i],yearC[,2,i])
    sigM[3,i]=a$p.value
    sigM[4,i]=100*((mean(yearC[,6,i],na.rm=T)/mean(yearC[,5,i],na.rm=T))-1)
    a=t.test(yearC[,6,i],yearC[,5,i])
    sigM[5,i]=a$p.value  
}


########### Average rainfall change
library(RNetCDF)

a=open.nc("/srv/ccrc/data45/z3478332/WRF/output/extracted_data/WRF_d01_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
a=open.nc("/srv/ccrc/data34/z3478332/WRF_d01_ESB_mask.nc")
maskE=var.get.nc(a,"ESB")
maskE[maskE==0]=NaN

# CONTDir="/srv/ccrc/data34/z3478332/WRF/ERA-nonudge/"
#  DATADir="/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nonudging_notopo/out/impact/"

dirs=c("/srv/ccrc/data34/z3478332/WRF/ERA-nonudge/","/srv/ccrc/data34/z3478332/WRF/ERA-nudge/",
       "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nonudging_notopo/out/impact/",
       "/srv/ccrc/data34/z3478332/WRF/output/ERAI_R2_nudging_notopo_19902009/out/impact/")

precip=matrix(0,20,4)

for(i in 1:4)
{
  f1=open.nc(paste(dirs[i],"WRF_d01_monthly_uvpr.nc",sep=""))
  r1=var.get.nc(f1,"PRCP_d01")
  
  for(y in 1:20)
  {
    x=seq(((y-1)*12+1),y*12)
    tmp=apply(r1[x,,],c(2,3),sum)
    precip[y,i]=mean(tmp*maskE,na.rm=T)
  }
}

  
r1=f1->PRCP_d01
f2=addfile(DATADir+"WRF_d01_monthly_uvpr.nc","r")
r2=f2->PRCP_d01


