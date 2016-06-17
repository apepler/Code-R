######### Intial stuff
rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing')
library(abind)
library("R.matlab")
library(fields)
library(maps)
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")
c=3 ## Which version do I want? 3 is the default

##########
## Step 1: load all the data for my category of choice

events<-events_notopo<-fixes<-fixes_notopo<-list()

n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    fixes[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    fixes[[n]]$Location2<-0
    I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    fixes[[n]][fixes[[n]]==-Inf]=NaN    
    fixes[[n]][fixes[[n]]==Inf]=NaN  
    
    events[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    events[[n]][events[[n]]==-Inf]=NaN    
    events[[n]]$Location2=0
    for(i in 1:length(events[[n]]$ID))
    {
      I=which(fixes[[n]]$ID==events[[n]]$ID[i])
      events[[n]]$Location2[i]=sum(fixes[[n]]$Location2[I])
    }
    
    
    fixes_notopo[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_notopo_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    fixes_notopo[[n]]$Year=floor(fixes_notopo[[n]]$Date/10000)
    fixes_notopo[[n]]$Month=floor(fixes_notopo[[n]]$Date/100)%%100
    fixes_notopo[[n]]$Location2<-0
    I<-which(fixes_notopo[[n]][,7]>=149 & fixes_notopo[[n]][,7]<=154 & fixes_notopo[[n]][,8]<(-37) & fixes_notopo[[n]][,8]>=-41)
    fixes_notopo[[n]]$Location2[I]<-1
    I<-which(fixes_notopo[[n]][,7]>=(149+(37+fixes_notopo[[n]][,8])/2) & fixes_notopo[[n]][,7]<=(154+(37+fixes_notopo[[n]][,8])/2) & fixes_notopo[[n]][,8]<(-31) & fixes_notopo[[n]][,8]>=-37)
    fixes_notopo[[n]]$Location2[I]<-1
    I<-which(fixes_notopo[[n]][,7]>=152 & fixes_notopo[[n]][,7]<=157 & fixes_notopo[[n]][,8]<=(-24) & fixes_notopo[[n]][,8]>=-31)
    fixes_notopo[[n]]$Location2[I]<-1
    fixes_notopo[[n]]$Date2=as.POSIXct(paste(as.character(fixes_notopo[[n]]$Date),substr(fixes_notopo[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    fixes_notopo[[n]][fixes_notopo[[n]]==-Inf]=NaN    
    fixes_notopo[[n]][fixes_notopo[[n]]==Inf]=NaN  
    
    events_notopo[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_notopo_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    events_notopo[[n]]$Year=floor(events_notopo[[n]]$Date1/10000)
    events_notopo[[n]]$Month=floor(events_notopo[[n]]$Date1/100)%%100
    events_notopo[[n]][events_notopo[[n]]==-Inf]=NaN
    
    events_notopo[[n]]$Location2=0
    for(i in 1:length(events_notopo[[n]]$ID))
    {
      I=which(fixes_notopo[[n]]$ID==events_notopo[[n]]$ID[i])
      events_notopo[[n]]$Location2[i]=sum(fixes_notopo[[n]]$Location2[I])
    }
    
    n=n+1     
  }

####### Some testing - match between d01 & d02
matchtest<-array(0,c(3,2,2))
for(i in 1:3)
{
  a=eventmatch(events[[i]],fixes[[i]],events[[i+3]],fixes[[i+3]])
  matchtest[i,1,1]=length(which(a[,4]>0))/length(a[,4])
  
  a=eventmatch(events_notopo[[i]],fixes_notopo[[i]],events_notopo[[i+3]],fixes_notopo[[i+3]])
  matchtest[i,2,1]=length(which(a[,4]>0))/length(a[,4])
  
  a=eventmatch(events[[i+3]],fixes[[i+3]],events[[i]],fixes[[i]])
  matchtest[i,1,2]=length(which(a[,4]>0))/length(a[,4])
  
  a=eventmatch(events_notopo[[i+3]],fixes_notopo[[i+3]],events_notopo[[i]],fixes_notopo[[i]])
  matchtest[i,2,2]=length(which(a[,4]>0))/length(a[,4])
}

matchtest2<-array(0,c(6,2))
for(i in 1:6)
{
  a=eventmatch(events_notopo[[i]],fixes_notopo[[i]],events_notopo[[5]],fixes_notopo[[5]])
  matchtest2[i,1]=length(which(a[,4]>0))/length(a[,4])
  a=eventmatch(events_notopo[[5]],fixes_notopo[[5]],events_notopo[[i]],fixes_notopo[[i]])
  matchtest2[i,2]=length(which(a[,4]>0))/length(a[,4])
}


count=matrix(NaN,6,2)
colnames(count)=c("Control","NoTopo")
rownames(count)=c("R1 50km","R2 50km","R2 50km","R1 10km","R2 10km","R3 10km")
cvthresh=NaN ## Change this if I want a minimum CV threshold
countthresh=NaN ## Change this if I want to use a CV threshold that gives X events p.a. for the BRAN (control) case

for(i in 1:6)
{
  if(!is.na(cvthresh))
  {
    count[i,1]=length(which(events[[i]]$CV2>=cvthresh))
    count[i,2]=length(which(events_notopo[[i]]$CV2>=cvthresh))
  } else if(!is.na(countthresh)) {
    a=order(events[[i]]$CV2,decreasing=T)
    if(length(a)>=countthresh) b=events[[i]]$CV2[a[countthresh]] else b=min(events[[i]]$CV2)
    
    count[i,1]=length(which(events[[i]]$CV2>=b))
    count[i,2]=length(which(events_notopo[[i]]$CV2>=b))
    
  } else {
    count[i,1]=length(events[[i]]$ID)
    count[i,2]=length(events_notopo[[i]]$ID)
  }
}
median(count[,2]/count[,1])

##### What about count of Location2 events

count=matrix(NaN,6,2)
colnames(count)=c("Control","NoTopo")
rownames(count)=c("R1 50km","R1 20km","R1 30km","R1 10km","R2 10km","R3 10km")
for(i in 1:6)
{
  count[i,1]=length(which(events[[i]]$Location2>1))
  count[i,2]=length(which(events_notopo[[i]]$Location2>1))
}
median(count[,2]/count[,1])

##### No change in overall freq

count2=matrix(NaN,6,4)
colnames(count2)=c("Control","NoTopo","Control L2","NoTopo L2")
rownames(count2)=c("R1 50km","R2 50km","R3 50km","R1 10km","R2 10km","R3 10km")
cvthresh=2 ## Change this if I want a minimum CV threshold
countthresh=NaN ## Change this if I want to use a CV threshold that gives X events p.a. for the BRAN (control) case

for(i in 1:6)
{
  count2[i,1]=length(which(fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh))
  count2[i,2]=length(which(fixes_notopo[[i]]$Location==1 & fixes_notopo[[i]]$CV>=cvthresh))
  count2[i,3]=length(which(fixes[[i]]$Location2==1 & fixes[[i]]$CV>=cvthresh))
  count2[i,4]=length(which(fixes_notopo[[i]]$Location2==1 & fixes_notopo[[i]]$CV>=cvthresh))
}
median(count[,2]/count[,1])

### What about extreme events?
count=array(NaN,c(6,2,12))
dimnames(count)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(count)[[2]]=c("Control","NoTopo")
dimnames(count)[[3]]=c("All","CV>=2.5","Bombs","Formed in region","EC","SC","Mixed","Mean rain>6 mm/6hr","Mean rain>12 mm/6hr",
                       "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")

for(j in 1:2)
  for(i in 1:6)
  {
    if(j==1) ev=events[[i]] else
      if(j==2) ev=events_notopo[[i]] 
      
      count[i,j,1]=length(ev$CV2)
      count[i,j,2]=length(which(ev$CV2>=2.5))
      count[i,j,3]=sum(ev$Bomb)
      count[i,j,4]=length(which(ev$EnteredFormed==1))
      count[i,j,5]=length(which(ev$HartType=="EC"))
      count[i,j,6]=length(which(ev$HartType=="SC"))
      count[i,j,7]=length(which(ev$HartType=="Mixed"))
      count[i,j,8]=length(which(ev$MaxMeanRain250>=6))
      count[i,j,9]=length(which(ev$MaxMeanRain250>=12))
      count[i,j,10]=length(which(ev$MaxPointRain250>=50))
      count[i,j,11]=length(which(ev$MaxMeanWind250>=13.9))
      count[i,j,12]=length(which(ev$MaxPointWind250>=22.2))
  }
count[,2,]/count[,1,]

apply(count[,2,]/count[,1,],2,median)


#### What about Loc?

count2=array(NaN,c(6,2,12))
dimnames(count2)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(count2)[[2]]=c("Control","NoTopo")
dimnames(count2)[[3]]=c("All","CV>=2.5","Bombs","Formed in region","EC","SC","Mixed","Mean rain>6 mm/6hr","Mean rain>12 mm/6hr",
                        "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")

for(j in 1:2)
  for(i in 1:6)
  {
    if(j==1) ev=events[[i]][events[[i]]$Location2>0,] else if(j==2) ev=events_notopo[[i]][events_notopo[[i]]$Location2>0,] 
    
    if(j==1) data=fixes[[i]][fixes[[i]]$Location2>0,] else if(j==2) data=fixes_notopo[[i]][fixes_notopo[[i]]$Location2>0,] 
    
    for(k in 1:length(ev$ID))
    {
      I=which(data$ID==ev$ID[k])
      ev$TotalRain500[k]=sum(data$MeanRain500[I],na.rm=T)
      ev$MaxMeanRain500[k]=max(data$MeanRain500[I],na.rm=T)
      ev$MaxPointRain500[k]=max(data$MaxRain500[I],na.rm=T)
      ev$MaxMeanWind500[k]=max(data$MeanWind500[I],na.rm=T)
      ev$MaxPointWind500[k]=max(data$MaxWind500[I],na.rm=T)
      ev$TotalRain250[k]=sum(data$MeanRain250[I],na.rm=T)
      ev$MaxMeanRain250[k]=max(data$MeanRain250[I],na.rm=T)
      ev$MaxPointRain250[k]=max(data$MaxRain250[I],na.rm=T)
      ev$MaxMeanWind250[k]=max(data$MeanWind250[I],na.rm=T)
      ev$MaxPointWind250[k]=max(data$MaxWind250[I],na.rm=T)
    
    }
    count2[i,j,1]=length(ev$CV_loc2)
    count2[i,j,2]=length(which(ev$CV_loc2>=2))
    count2[i,j,3]=sum(ev$Bomb)
    count2[i,j,4]=length(which(ev$EnteredFormed==1))
    count2[i,j,5]=length(which(ev$HartType=="EC"))
    count2[i,j,6]=length(which(ev$HartType=="SC"))
    count2[i,j,7]=length(which(ev$HartType=="Mixed"))
    count2[i,j,8]=length(which(ev$MaxMeanRain250>=6))
    count2[i,j,9]=length(which(ev$MaxMeanRain250>=12))
    count2[i,j,10]=length(which(ev$MaxPointRain250>=50))
    count2[i,j,11]=length(which(ev$MaxMeanWind250>=13.9))
    count2[i,j,12]=length(which(ev$MaxPointWind250>=22.2))
    
  }
apply(count2[,2,]/count2[,1,],2,median)



for(i in 1:6) print(makePDF(fixes[[i]]$CV[fixes[[i]]$Location2==1],fixes_notopo[[i]]$CV[fixes_notopo[[i]]$Location2==1]))


######## Matching

match<-list()
n=1
for(n in 1:6)
  match[[n]]=eventmatch(events[[n]],fixes[[n]],events_notopo[[n]],fixes_notopo[[n]],T)

### Compare - what happens to cv, slp etc when matched/unmatched?

comp=array(NaN,c(6,5))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2","GV2")

for(n in 1:6)
{
  comp[n,1]=length(which(match[[n]][,5]>0))/length(match[[n]][,5])
  comp[n,2]=mean(match[[n]][,6],na.rm=T)
  comp[n,3]=mean(match[[n]][,7]-match[[n]][,1],na.rm=T)
  comp[n,4]=mean(match[[n]][,8]-match[[n]][,2],na.rm=T)
  comp[n,5]=mean(match[[n]][,10]-match[[n]][,4],na.rm=T)
}
apply(comp,2,mean)


for(n in 1:6) print(t.test(match[[n]][,10],match[[n]][,4]))

cvthresh=c(seq(1,3,0.5),Inf)
CVchange<-array(NaN,c(5,6,3))

for(i in 1:5)
  for(j in 1:6)
  {
    I=which(match[[j]][,1]>=cvthresh[i] & match[[j]][,1]<cvthresh[i+1])
    CVchange[i,j,1]=length(I)
    CVchange[i,j,2]=length(which(!is.na(match[[j]][I,7])))/length(I)
    CVchange[i,j,3]=mean(match[[j]][I,7]-match[[j]][I,1],na.rm=T)
  }

cvthresh=seq(1,2.5,0.25)
CVchange2<-array(NaN,c(7,6,3))

##This time, cumulative (or anti-cumulative?)
for(i in 1:7)
  for(j in 1:6)
  {
    I=which(match[[j]][,1]>=cvthresh[i])
    CVchange2[i,j,1]=length(I)
    CVchange2[i,j,2]=length(which(!is.na(match[[j]][I,7])))/length(I)
    CVchange2[i,j,3]=mean(match[[j]][I,7]-match[[j]][I,1],na.rm=T)
  }

pnt=c(1,2,4,1,2,4)
cm=c("black","black","black","darkgrey","darkgrey","darkgrey")
plot(cvthresh,CVchange2[,1,3],"b",ylim=c(-1.5,0.5),lwd=2,cex=1.5,pch=pnt[1],col=cm[1],
     xlab="Minimum intensity in control",ylab="Average CV change")
abline(h=0,col="red")
for(i in 2:6) lines(cvthresh,CVchange2[,i,3],"b",lwd=2,cex=1.5,pch=pnt[i],col=cm[i])
legend("bottomleft",c("R1","R2","R3","50km","10km"),pch=c(1,2,4,1,1),pt.lwd=2,col=c(1,1,1,1,"darkgrey"),pt.cex=1.5,bty="n",ncol=2)

### Repeat for GV change?

gvthresh=c(seq(0,50,10),Inf)
GVchange<-array(NaN,c(6,6,3))

for(i in 1:6)
  for(j in 1:6)
  {
    I=which(match[[j]][,4]>=gvthresh[i] & match[[j]][,4]<gvthresh[i+1])
    GVchange[i,j,1]=length(I)
    GVchange[i,j,2]=length(which(!is.na(match[[j]][I,10])))/length(I)
    GVchange[i,j,3]=mean(match[[j]][I,10]-match[[j]][I,4],na.rm=T)
  }

cvthresh=seq(1,2.5,0.25)
CVchange2<-array(NaN,c(7,6,3))

##This time, cumulative (or anti-cumulative?)
for(i in 1:7)
  for(j in 1:6)
  {
    I=which(match[[j]][,4]>=gvthresh[i])
    GVchange2[i,j,1]=length(I)
    GVchange2[i,j,2]=length(which(!is.na(match[[j]][I,10])))/length(I)
    GVchange2[i,j,3]=mean(match[[j]][I,10]-match[[j]][I,4],na.rm=T)
  }

####### Changes in mean rain rate/rain distribution

extremeF=array(NaN,c(6,8,3,6))
dimnames(extremeF)[[1]]=c("R1 50km","R2 50km","R2 50km","R1 10km","R2 10km","R3 10km")
dimnames(extremeF)[[2]]=colnames(fixes[[1]])[24:31]
dimnames(extremeF)[[3]]=c("All","Loc1","Loc2")
dimnames(extremeF)[[4]]=c("Count Control","Count NT","Mean Control","Mean NT","T.test","KS.test")

for(i in 1:6)
  for(j in 1:8)
  {
    extremeF[i,j,1,1]=length(fixes[[i]]$Location)
    extremeF[i,j,1,2]=length(fixes_notopo[[i]]$Location)
    a=t.test(fixes[[i]][,23+j],fixes_notopo[[i]][,23+j])
    extremeF[i,j,1,3:4]=a$estimate
    extremeF[i,j,1,5]=a$p.value
    a=ks.test(fixes[[i]][,23+j],fixes_notopo[[i]][,23+j])
    extremeF[i,j,1,6]=a$p.value
    
    I=which(fixes[[i]]$Location==1)
    J=which(fixes_notopo[[i]]$Location==1)
    extremeF[i,j,2,1]=length(I)
    extremeF[i,j,2,2]=length(J)
    a=t.test(fixes[[i]][I,23+j],fixes_notopo[[i]][J,23+j])
    extremeF[i,j,2,3:4]=a$estimate
    extremeF[i,j,2,5]=a$p.value
    a=ks.test(fixes[[i]][I,23+j],fixes_notopo[[i]][J,23+j])
    extremeF[i,j,2,6]=a$p.value
    
    I=which(fixes[[i]]$Location2==1)
    J=which(fixes_notopo[[i]]$Location2==1)
    extremeF[i,j,3,1]=length(I)
    extremeF[i,j,3,2]=length(J)
    a=t.test(fixes[[i]][I,23+j],fixes_notopo[[i]][J,23+j])
    extremeF[i,j,3,3:4]=a$estimate
    extremeF[i,j,3,5]=a$p.value
    a=ks.test(fixes[[i]][I,23+j],fixes_notopo[[i]][J,23+j])
    extremeF[i,j,3,6]=a$p.value
  }

apply(extremeF[,,,4]-extremeF[,,,3],c(2,3),mean,na.rm=T)
apply(extremeF[,,,5],c(2,3),median,na.rm=T)


######### Another test- look at ECL days

dates=seq(as.Date("20070101",format="%Y%m%d"),
          as.Date("20081231",format="%Y%m%d"),
          by="day",format="%Y%m%d")
dates=as.numeric(format.Date(dates,"%Y%m%d"))
dates2=array(0,c(length(dates),6,4))
dimnames(dates2)[[3]]=c("Control","NoTopo","Control Loc2","NoTopo Loc2")

for(i in 1:6)
{
  tmp=sort(unique(fixes[[i]]$Date[fixes[[i]]$Location==1]))
  for(j in 1:length(tmp))
  {
    I=which(dates==tmp[j])
    J=which(fixes[[i]]$Date==tmp[j] & fixes[[i]]$Location==1)
    dates2[I,i,1]=max(fixes[[i]]$CV[J])
  }
  tmp=sort(unique(fixes_notopo[[i]]$Date[fixes_notopo[[i]]$Location==1]))
  for(j in 1:length(tmp))
  {
    I=which(dates==tmp[j])
    J=which(fixes_notopo[[i]]$Date==tmp[j] & fixes_notopo[[i]]$Location==1)
    dates2[I,i,2]=max(fixes_notopo[[i]]$CV[J])
  }
  
  tmp=sort(unique(fixes[[i]]$Date[fixes[[i]]$Location2==1]))
  for(j in 1:length(tmp))
  {
    I=which(dates==tmp[j])
    J=which(fixes[[i]]$Date==tmp[j] & fixes[[i]]$Location2==1)
    dates2[I,i,3]=max(fixes[[i]]$CV[J])
  }
  tmp=sort(unique(fixes_notopo[[i]]$Date[fixes_notopo[[i]]$Location2==1]))
  for(j in 1:length(tmp))
  {
    I=which(dates==tmp[j])
    J=which(fixes_notopo[[i]]$Date==tmp[j] & fixes_notopo[[i]]$Location2==1)
    dates2[I,i,4]=max(fixes_notopo[[i]]$CV[J])
  }
  
}

dates3=dates2>0

CSI=array(0,c(6,6,4))

for(i in 1:6)
  for(j in 1:6)
    for(k in 1:4)
      CSI[i,j,k]=length(which(dates3[,i,k]==1 & dates3[,j,k]==1))/length(which(dates3[,i,k]==1 | dates3[,j,k]==1)) ## CSI


######## Location location

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

loc<-array(NaN,c(11,17,6,2,2))
for(i in 1:11)
  for(j in 1:17)
    for(k in 1:6)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1)
      loc[i,j,k,1,1]=length(I)
      loc[i,j,k,1,2]=median(fixes[[k]]$MeanRain250[I],na.rm=T)
      I=which(fixes_notopo[[k]]$Lat>=lat[i]-2.5 & fixes_notopo[[k]]$Lat<lat[i]+2.5 & fixes_notopo[[k]]$Lon>=lon[j]-2.5 & fixes_notopo[[k]]$Lon<lon[j]+2.5 & fixes_notopo[[k]]$Location==1)
      loc[i,j,k,2,1]=length(I)
      loc[i,j,k,2,2]=median(fixes_notopo[[k]]$MeanRain250[I],na.rm=T)
    }

loc2=apply(loc,c(1,2,4,5),mean,na.rm=T)


#cols=gray(seq(1,0.1,-0.15))
#bb=c(-0.5,0,1,2,4,8,12,100)

bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.3))
image(lon,lat,100*t(apply((loc[,,,2,2]/loc[,,,1,2])-1,c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),cex.axis=1.5,cex.main=1.5)
tmp=apply(loc[,,,2,2]/loc[,,,1,2]>1,c(1,2),mean) ## Proportion with a positive change
sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=2,lwd=2) ## Add some points
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)