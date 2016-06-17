######### Intial stuff
rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing')
library(abind)
library("R.matlab")
library(fields)
library(maps)
figdir="~/Documents/ECLs/WRFruns/0708/EACpaper/"
source("~/Documents/R/ECL_functions.R")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")
c=3 ## Which version do I want? 3 is the default

##########
## Step 1: load all the data for my category of choice

eventsBRAN<-events_noeac<-events_2eac<-fixesBRAN<-fixes_noeac<-fixes_2eac<-events_2eac2<-fixes_2eac2<-eventsLR<-fixesLR<-list()

n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    eventsBRAN[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    eventsBRAN[[n]]$Year=floor(eventsBRAN[[n]]$Date1/10000)
    eventsBRAN[[n]]$Month=floor(eventsBRAN[[n]]$Date1/100)%%100
    eventsBRAN[[n]][eventsBRAN[[n]]==-Inf]=NaN
    
    eventsLR[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    eventsLR[[n]]$Year=floor(eventsLR[[n]]$Date1/10000)
    eventsLR[[n]]$Month=floor(eventsLR[[n]]$Date1/100)%%100
    eventsLR[[n]][eventsLR[[n]]==-Inf]=NaN    
    fixesLR[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixesLR[[n]]$Year=floor(fixesLR[[n]]$Date/10000)
    fixesLR[[n]]$Month=floor(fixesLR[[n]]$Date/100)%%100
    fixesLR[[n]]$Date2=as.POSIXct(paste(as.character(fixesLR[[n]]$Date),substr(fixesLR[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    
    events_noeac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    events_noeac[[n]]$Year=floor(events_noeac[[n]]$Date1/10000)
    events_noeac[[n]]$Month=floor(events_noeac[[n]]$Date1/100)%%100
    events_noeac[[n]][events_noeac[[n]]==-Inf]=NaN
    
    fixesBRAN[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixesBRAN[[n]]$Year=floor(fixesBRAN[[n]]$Date/10000)
    fixesBRAN[[n]]$Month=floor(fixesBRAN[[n]]$Date/100)%%100
    fixesBRAN[[n]]$Location2<-0
    I<-which(fixesBRAN[[n]][,7]>=149 & fixesBRAN[[n]][,7]<=154 & fixesBRAN[[n]][,8]<(-37) & fixesBRAN[[n]][,8]>=-41)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=(149+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,7]<=(154+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,8]<(-31) & fixesBRAN[[n]][,8]>=-37)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=152 & fixesBRAN[[n]][,7]<=157 & fixesBRAN[[n]][,8]<=(-24) & fixesBRAN[[n]][,8]>=-31)
    fixesBRAN[[n]]$Location2[I]<-1
    fixesBRAN[[n]]$Date2=as.POSIXct(paste(as.character(fixesBRAN[[n]]$Date),substr(fixesBRAN[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    

    fixes_noeac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixes_noeac[[n]]$Year=floor(fixes_noeac[[n]]$Date/10000)
    fixes_noeac[[n]]$Month=floor(fixes_noeac[[n]]$Date/100)%%100
    fixes_noeac[[n]]$Location2<-0
    I<-which(fixes_noeac[[n]][,7]>=149 & fixes_noeac[[n]][,7]<=154 & fixes_noeac[[n]][,8]<(-37) & fixes_noeac[[n]][,8]>=-41)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=(149+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,7]<=(154+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,8]<(-31) & fixes_noeac[[n]][,8]>=-37)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=152 & fixes_noeac[[n]][,7]<=157 & fixes_noeac[[n]][,8]<=(-24) & fixes_noeac[[n]][,8]>=-31)
    fixes_noeac[[n]]$Location2[I]<-1
    fixes_noeac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_noeac[[n]]$Date),substr(fixes_noeac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
# 
#     if(dom=="d01")
#     {
#       fixes_2eac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
#       fixes_2eac[[n]]$Year=floor(fixes_2eac[[n]]$Date/10000)
#       fixes_2eac[[n]]$Month=floor(fixes_2eac[[n]]$Date/100)%%100
#       fixes_2eac[[n]]$Location2<-0
#       I<-which(fixes_2eac[[n]][,7]>=149 & fixes_2eac[[n]][,7]<=154 & fixes_2eac[[n]][,8]<(-37) & fixes_2eac[[n]][,8]>=-41)
#       fixes_2eac[[n]]$Location2[I]<-1
#       I<-which(fixes_2eac[[n]][,7]>=(149+(37+fixes_2eac[[n]][,8])/2) & fixes_2eac[[n]][,7]<=(154+(37+fixes_2eac[[n]][,8])/2) & fixes_2eac[[n]][,8]<(-31) & fixes_2eac[[n]][,8]>=-37)
#       fixes_2eac[[n]]$Location2[I]<-1
#       I<-which(fixes_2eac[[n]][,7]>=152 & fixes_2eac[[n]][,7]<=157 & fixes_2eac[[n]][,8]<=(-24) & fixes_2eac[[n]][,8]>=-31)
#       fixes_2eac[[n]]$Location2[I]<-1
#       fixes_2eac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_2eac[[n]]$Date),substr(fixes_2eac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
#       
# 
#       events_2eac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
#       events_2eac[[n]]$Year=floor(events_2eac[[n]]$Date1/10000)
#       events_2eac[[n]]$Month=floor(events_2eac[[n]]$Date1/100)%%100  
#     }
    
    fixes_2eac2[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixes_2eac2[[n]]$Year=floor(fixes_2eac2[[n]]$Date/10000)
    fixes_2eac2[[n]]$Month=floor(fixes_2eac2[[n]]$Date/100)%%100
    fixes_2eac2[[n]]$Location2<-0
    I<-which(fixes_2eac2[[n]][,7]>=149 & fixes_2eac2[[n]][,7]<=154 & fixes_2eac2[[n]][,8]<(-37) & fixes_2eac2[[n]][,8]>=-41)
    fixes_2eac2[[n]]$Location2[I]<-1
    I<-which(fixes_2eac2[[n]][,7]>=(149+(37+fixes_2eac2[[n]][,8])/2) & fixes_2eac2[[n]][,7]<=(154+(37+fixes_2eac2[[n]][,8])/2) & fixes_2eac2[[n]][,8]<(-31) & fixes_2eac2[[n]][,8]>=-37)
    fixes_2eac2[[n]]$Location2[I]<-1
    I<-which(fixes_2eac2[[n]][,7]>=152 & fixes_2eac2[[n]][,7]<=157 & fixes_2eac2[[n]][,8]<=(-24) & fixes_2eac2[[n]][,8]>=-31)
    fixes_2eac2[[n]]$Location2[I]<-1
    fixes_2eac2[[n]]$Date2=as.POSIXct(paste(as.character(fixes_2eac2[[n]]$Date),substr(fixes_2eac2[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    
    events_2eac2[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    events_2eac2[[n]]$Year=floor(events_2eac2[[n]]$Date1/10000)
    events_2eac2[[n]]$Month=floor(events_2eac2[[n]]$Date1/100)%%100  
    events_2eac2[[n]][events_2eac2[[n]]==-Inf]=NaN
    n=n+1     
  }

##########
## Step 2: Number of total events for a certain # events threshold.

count=matrix(NaN,6,5)
countthresh=NaN ##Number of total events
cvthresh=1.5 ##Minimum CV threshold

for(i in 1:6)
{
  if(!is.na(cvthresh))
  {
    count[i,1]=length(which(eventsBRAN[[i]]$CV2>=cvthresh))
    count[i,2]=length(which(events_noeac[[i]]$CV2>=cvthresh))
    #if(i<4) count[i,3]=length(which(events_2eac[[i]]$CV2>=cvthresh))
    count[i,4]=length(which(events_2eac2[[i]]$CV2>=cvthresh))
  } else if(!is.na(countthresh)) {
    a=order(eventsBRAN[[i]]$CV2,decreasing=T)
    if(length(a)>=countthresh) b=eventsBRAN[[i]]$CV2[a[countthresh]] else b=min(eventsBRAN[[i]]$CV2)
    
    count[i,1]=length(which(eventsBRAN[[i]]$CV2>=b))
    count[i,2]=length(which(events_noeac[[i]]$CV2>=b))
    #if(i<4) count[i,3]=length(which(events_2eac[[i]]$CV2>=b))
    count[i,4]=length(which(events_2eac2[[i]]$CV2>=b))
  } else {
    count[i,1]=length(eventsBRAN[[i]]$ID)
    count[i,2]=length(events_noeac[[i]]$ID)
    #if(i<4) count[i,3]=length(events_2eac[[i]]$ID)
    count[i,4]=length(events_2eac2[[i]]$ID)
    count[i,5]=length(eventsLR[[i]]$ID)
  }
}

count2=cbind(count[,2],count[,1],count[,4])
pnt=c(1,2,4,1,2,4)
cm=c("black","black","black","darkgrey","darkgrey","darkgrey")

pdf(paste(figdir,"ECLcount_change.pdf",sep=""),width=4,height=4)
par(mar=c(3,4,1,1))
plot(1:3,count2[1,],"b",axes=F,ylim=c(20,80),lwd=2,cex=1.5,pch=pnt[1],col=cm[1],
     xlab="",ylab="Number of ECLs")
axis(1,at=1:3,labels=c("No EAC","Control","2EAC"))
axis(2)
for(i in 2:6) lines(1:3,count2[i,],"b",lwd=2,cex=1.5,pch=pnt[i],col=cm[i])
points(rep(2,6),count[,5],pch=pnt,col=cm,cex=1.5,lwd=2)
legend("topleft",c("R1","R2","R3","50km","10km"),pch=c(1,2,4,1,1),pt.lwd=2,col=c(1,1,1,1,"darkgrey"),pt.cex=1.5,bty="n",ncol=2)
dev.off()

tmp=cbind(count[,1],count[,5],((count[,2]/count[,1])-1),((count[,4]/count[,1])-1))

median(count[,4]/count[,1],na.rm=T)

count2=array(NaN,c(6,4,2))
dimnames(count2)[[2]]=c("BRAN","NoEAC","2EAC","2EAC2")
dimnames(count2)[[3]]=c("Cool","Warm")
countthresh=NaN ##Number of total events
cvthresh=NaN ##Minimum CV threshold

for(i in 1:6)
{
  if(i<4) tmp=list(eventsBRAN[[i]],events_noeac[[i]],events_2eac[[i]],events_2eac2[[i]]) else tmp=list(eventsBRAN[[i]],events_noeac[[i]],events_2eac[[i-3]],events_2eac2[[i]])
  if(!is.na(cvthresh))
    for(j in 1:4)
    {
      tmp2=tmp[[j]][tmp[[j]]$CV2>=cvthresh,]
      mm=floor(tmp2$Date1/100)%%100
      count2[i,j,1]=length(which(mm>=5 & mm<=10))
      count2[i,j,2]=length(which(mm>=11 | mm<=4))
    }
  else if(!is.na(countthresh)) 
    for(j in 1:4)
    {
      a=order(tmp[[1]]$CV2,decreasing=T)
      if(length(a)>=countthresh) cvthresh=tmp[[1]]$CV2[a[countthresh]] else cvthresh=min(tmp[[1]]$CV2)
      
      tmp2=tmp[[j]][tmp[[j]]$CV2>=cvthresh,]
      mm=floor(tmp2$Date1/100)%%100
      count2[i,j,1]=length(which(mm>=5 & mm<=10))
      count2[i,j,2]=length(which(mm>=11 | mm<=4))
      
    } 
  else
      for(j in 1:4)
      {
        tmp2=tmp[[j]]
        mm=floor(tmp2$Date1/100)%%100
        count2[i,j,1]=length(which(mm>=5 & mm<=10))
        count2[i,j,2]=length(which(mm>=11 | mm<=4))
      }
}

apply(count2[,2,]/count2[,1,],2,mean)
apply(count2[,4,]/count2[,1,],2,mean)

#### Count of hours with location 1/2

count3=array(NaN,c(6,4,3))
dimnames(count3)[[2]]=c("BRAN","NoEAC","2EAC","2EAC2")
dimnames(count3)[[3]]=c("All","Location2","Outside Location 2")

for(i in 1:6)
{
  if(i<4) tmp=list(fixesBRAN[[i]],fixes_noeac[[i]],fixes_2eac[[i]],fixes_2eac2[[i]]) else tmp=list(fixesBRAN[[i]],fixes_noeac[[i]],fixes_2eac[[i-3]],fixes_2eac2[[i]])

    for(j in 1:4)
    {
      tmp2=tmp[[j]]
      count3[i,j,1]=length(unique(tmp2$Date2[tmp2$Location==1]))
      count3[i,j,2]=length(unique(tmp2$Date2[tmp2$Location2==1]))
      count3[i,j,3]=length(unique(tmp2$Date2[tmp2$Location==1 & tmp2$Location2==0]))
    }
}

apply(count3[,2,]/count3[,1,],2,mean)
apply(count3[,4,]/count3[,1,],2,mean)


####### Change in mean/distribution of different parameters
###### Based on events

change=array(NaN,c(6,3,4,3))
dimnames(change)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(change)[[2]]=c("NoEAC","2EAC","2EAC2")
dimnames(change)[[3]]=colnames(eventsBRAN[[1]])[8:11]
dimnames(change)[[4]]=c("Mean difference","t-test sig","ks-test sig")

for(j in 1:3)
  for(i in 1:6)
  {
    control=eventsBRAN[[i]] 
      if(j==1) events=events_noeac[[i]] else
        #if(i<4 & j==2) events=events_2eac[[i]] else 
          if(j==3) events=events_2eac2[[i]] else next
          
          
    for(k in 1:4)
    {
      change[i,j,k,1]=mean(events[,k+7])-mean(control[,k+7])
      a=t.test(control[,k+7],events[,k+7])
      change[i,j,k,2]=a$p.value
      a=ks.test(control[,k+7],events[,k+7])
      change[i,j,k,3]=a$p.value
    }
  }




###Histogram of change vs intensity threshold

cvthresh=c(seq(1,4,0.5),NaN)
cvcount=array(0,c(7,6,4))
for(i in 1:6)
{
  if(i<4) tmp=list(eventsBRAN[[i]],events_noeac[[i]],events_2eac[[i]],events_2eac2[[i]]) else tmp=list(eventsBRAN[[i]],events_noeac[[i]],events_2eac[[i-3]],events_2eac2[[i]])
  for(j in 1:4)
    for(x in 1:7)
      cvcount[x,i,j]=length(which(tmp[[j]]$CV2>=cvthresh[x] & tmp[[j]]$CV2<cvthresh[x+1]))
}
cvcount[,4:6,3]=NaN

a=apply(cvcount[1:2,,],c(2,3),sum)
b=apply(cvcount[3:7,,],c(2,3),sum)

apply(cvcount[,,2]-cvcount[,,1],1,mean,na.rm=T)
mean(apply(cvcount[1:2,,4],2,sum,na.rm=T)-apply(cvcount[1:2,,1],2,sum,na.rm=T))
mean(apply(cvcount[3:13,,4],2,sum,na.rm=T)-apply(cvcount[3:13,,1],2,sum,na.rm=T))

cvcount2=apply(cvcount,c(1,3),mean)

a=hist(eventsBRAN[[1]]$CV2,breaks=seq(1,4,0.5),plot=F)
tmp=a
tmp$counts=apply(cvcount[,,2]-cvcount[,,1],1,mean,na.rm=T)
plot(tmp,col=rgb(0,0,1,1/4),xlim=c(1,4),ylim=c(-8,8),
     xlab="Intensity (hPa/(deg.lat)^2)",ylab="Average number of ECLs",main="")
tmp$counts=apply(cvcount[,,4]-cvcount[,,1],1,mean,na.rm=T)
plot(tmp,col=rgb(1,0,0,1/4),add=T)
legend("topright",legend=c("No EAC","Doubled EAC"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4)),lwd=8,cex=1,bty="n")   


####################
#### Step 3: Location of events

## Did it start in the ECL domain?
Gprop=matrix(NaN,6,4)
for(i in 1:6)
{
  I=which(fixesBRAN[[i]]$Fix==1)
  J=which(fixesBRAN[[i]]$Location[I]==1)
  Gprop[i,1]=length(J)/length(I)
  
  I=which(fixes_noeac[[i]]$Fix==1)
  J=which(fixes_noeac[[i]]$Location[I]==1)
  Gprop[i,2]=length(J)/length(I)
  
  if(i<4)
  {
    I=which(fixes_2eac[[i]]$Fix==1)
    J=which(fixes_2eac[[i]]$Location[I]==1)
    Gprop[i,3]=length(J)/length(I)
  }
  
  I=which(fixes_2eac2[[i]]$Fix==1)
  J=which(fixes_2eac2[[i]]$Location[I]==1)
  Gprop[i,4]=length(J)/length(I)
  
}

##### Maps of event locations

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

loc<-array(NaN,c(11,17,6,4))
for(i in 1:11)
  for(j in 1:17)
    for(k in 1:6)
    {
      I=which(fixesBRAN[[k]]$Lat>=lat[i]-2.5 & fixesBRAN[[k]]$Lat<lat[i]+2.5 & fixesBRAN[[k]]$Lon>=lon[j]-2.5 & fixesBRAN[[k]]$Lon<lon[j]+2.5 & fixesBRAN[[k]]$Location==1)
      loc[i,j,k,1]=length(I)
      I=which(fixes_noeac[[k]]$Lat>=lat[i]-2.5 & fixes_noeac[[k]]$Lat<lat[i]+2.5 & fixes_noeac[[k]]$Lon>=lon[j]-2.5 & fixes_noeac[[k]]$Lon<lon[j]+2.5 & fixes_noeac[[k]]$Location==1)
      loc[i,j,k,2]=length(I)
      if(k<4)
      {
        I=which(fixes_2eac[[k]]$Lat>=lat[i]-2.5 & fixes_2eac[[k]]$Lat<lat[i]+2.5 & fixes_2eac[[k]]$Lon>=lon[j]-2.5 & fixes_2eac[[k]]$Lon<lon[j]+2.5 & fixes_2eac[[k]]$Location==1)
        loc[i,j,k,3]=length(I)
      }
      I=which(fixes_2eac2[[k]]$Lat>=lat[i]-2.5 & fixes_2eac2[[k]]$Lat<lat[i]+2.5 & fixes_2eac2[[k]]$Lon>=lon[j]-2.5 & fixes_2eac2[[k]]$Lon<lon[j]+2.5 & fixes_2eac2[[k]]$Location==1)
      loc[i,j,k,4]=length(I)
    }

loc2=apply(loc,c(1,2,4),mean,na.rm=T)

#### Genesis Locations

locG<-array(NaN,c(11,17,6,4))

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:6)
    {
      I=which(fixesBRAN[[k]]$Lat>=lat[i]-2.5 & fixesBRAN[[k]]$Lat<lat[i]+2.5 & fixesBRAN[[k]]$Lon>=lon[j]-2.5 & fixesBRAN[[k]]$Lon<lon[j]+2.5 & fixesBRAN[[k]]$Fix==1)
      locG[i,j,k,1]=length(I)
      I=which(fixes_noeac[[k]]$Lat>=lat[i]-2.5 & fixes_noeac[[k]]$Lat<lat[i]+2.5 & fixes_noeac[[k]]$Lon>=lon[j]-2.5 & fixes_noeac[[k]]$Lon<lon[j]+2.5 & fixes_noeac[[k]]$Fix==1)
      locG[i,j,k,2]=length(I)
      if(k<4)
      {
        I=which(fixes_2eac[[k]]$Lat>=lat[i]-2.5 & fixes_2eac[[k]]$Lat<lat[i]+2.5 & fixes_2eac[[k]]$Lon>=lon[j]-2.5 & fixes_2eac[[k]]$Lon<lon[j]+2.5 & fixes_2eac[[k]]$Fix==1)
        locG[i,j,k,3]=length(I)
      }
      I=which(fixes_2eac2[[k]]$Lat>=lat[i]-2.5 & fixes_2eac2[[k]]$Lat<lat[i]+2.5 & fixes_2eac2[[k]]$Lon>=lon[j]-2.5 & fixes_2eac2[[k]]$Lon<lon[j]+2.5 & fixes_2eac2[[k]]$Fix==1)
      locG[i,j,k,4]=length(I)
    }

locG2=apply(locG,c(1,2,4),mean,na.rm=T)

####### Some mapping code

#cols=gray(seq(1,0.1,-0.15))
#bb=c(-0.5,0,1,2,4,8,12,100)

bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
names=c("BRAN","NoEAC","2EAC","2EAC")

pdf(file=paste(figdir,"/ECL_locations_0708_",cat[c],"_SST_change_v2_median.pdf",sep=""),width=8,height=4,pointsize=12)
layout(cbind(1,2,3),c(1,1,0.3))
par(mar=c(3,3,3,2))
for(i in c(2,4))
{
  image(lon,lat,100*t(apply((loc[,,,i]/loc[,,,1])-1,c(1,2),median,na.rm=T)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=names[i],cex.axis=1.5,cex.main=1.5)
  map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
}
ColorBar(bb2,cm)
dev.off()

pdf(file=paste(figdir,"/ECL_locations_0708_",cat[c],"_SST_change_v2_points.pdf",sep=""),width=8,height=4,pointsize=12)
layout(cbind(1,2,3),c(1,1,0.3))
par(mar=c(3,3,3,2))
for(i in c(2,4))
{
  image(lon,lat,100*t(apply((loc[,,,i]/loc[,,,1])-1,c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=names[i],cex.axis=1.5,cex.main=1.5)
  
  tmp=apply(loc[,,,i]/loc[,,,1]>1,c(1,2),mean)
  sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T)
  points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=2,lwd=2)
  
  map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
}
ColorBar(bb2,cm)
dev.off()


###########################
### Chanegs in general statistics of ECLs

## New count matrix

count=array(NaN,c(6,4,11))
dimnames(count)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(count)[[2]]=c("Control","NoEAC","2EAC","2EAC2")
dimnames(count)[[3]]=c("All","CV>=2","Bombs","Formed in region","EC","SC","Mixed","Mean rain>6 mm/6hr",
                       "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")

for(j in 1:4)
  for(i in 1:6)
  {
    if(j==1) events=eventsBRAN[[i]] else
      if(j==2) events=events_noeac[[i]] else
        #if(i<4 & j==3) events=events_2eac[[i]] else 
          if(j==4) events=events_2eac2[[i]] else next
        
    count[i,j,1]=length(events$CV2)
    count[i,j,2]=length(which(events$CV2>=2))
    count[i,j,3]=sum(events$Bomb)
    count[i,j,4]=length(which(events$EnteredFormed==1))
    count[i,j,5]=length(which(events$HartType=="EC"))
    count[i,j,6]=length(which(events$HartType=="SC"))
    count[i,j,7]=length(which(events$HartType=="Mixed"))
    count[i,j,8]=length(which(events$MaxMeanRain250>=6))
    count[i,j,9]=length(which(events$MaxPointRain250>=50))
    count[i,j,10]=length(which(events$MaxMeanWind250>=13.9))
    count[i,j,11]=length(which(events$MaxPointWind250>=22.2))
  }
apply(count[,1,],2,mean)
apply(count[,4,]/count[,1,],2,mean,na.rm=T)
count[,2,]/count[,1,]

HartChange<-array(0,c(3,6,2))
dimnames(HartChange)[[1]]=c("EC","SC","Mixed")
dimnames(HartChange)[[3]]=c("NoEAC","2EAC")

for(i in 1:3)
  for(j in 1:6)
  {
    HartChange[i,j,1]=count[j,2,4+i]-count[j,1,4+i]
    HartChange[i,j,2]=count[j,4,4+i]-count[j,1,4+i]
  }


Form=count[,,4]
Out=count[,,1]-count[,,4]
median(Form[,2]/Form[,1])
median(Out[,2]/Out[,1],na.rm=T)


### Change in average extremeness
extreme=array(NaN,c(6,4,10))
dimnames(extreme)[[3]]=colnames(events)[18:27]
for(i in 1:6)
    for(j in 1:4)
    {
      if(j==1) events=eventsBRAN[[i]] else
        if(j==2) events=events_noeac[[i]] else
          #if(i<4 & j==3) events=events_2eac[[i]] else 
            if(j==4) events=events_2eac2[[i]] else next
          
            
      events[events==-Inf]=NaN
      extreme[i,j,1:10]=apply(events[,18:27],2,mean,na.rm=T)
    }
  
apply((extreme[,2,]-extreme[,1,]),2,mean,na.rm=T) #% change

testex=array(NaN,c(6,3,10))
dimnames(testex)[[3]]=colnames(eventsBRAN[[1]])[18:27]
for(i in 1:6)
  for(j in 1:10)
  {
    a=t.test(eventsBRAN[[i]][,j+17],events_noeac[[i]][,j+17])
    testex[i,1,j]=a$p.value
#     if(i<4)
#     {
#     a=t.test(eventsBRAN[[i]][,j+17],events_2eac[[i]][,j+17])
#     testex[i,2,j]=a$p.value
#     }
    a=t.test(eventsBRAN[[i]][,j+17],events_2eac2[[i]][,j+17])
    testex[i,3,j]=a$p.value
  }

##### Comparing all low centres w/in ECL region

extremeF=array(NaN,c(6,4,8))
dimnames(extremeF)[[3]]=colnames(fixesBRAN[[1]])[23:30]
for(i in 1:6)
  for(j in 1:4)
  {
    if(j==1) fixes=fixesBRAN[[i]] else
      if(j==2) fixes=fixes_noeac[[i]] else
        if(i<4 & j==3) fixes=fixes_2eac[[i]] else 
          if(j==4) fixes=fixes_2eac2[[i]] else next
          
          I=which(fixes$Location==1)
          extremeF[i,j,]=apply(fixes[I,23:30],2,mean,na.rm=T)
  }
apply((extremeF[,1,]),2,mean) #% change
apply((extremeF[,4,]-extremeF[,1,]),2,mean) #% change
apply((extremeF[,2,]-extremeF[,1,]),2,mean) #% change

testexF=array(NaN,c(6,3,8))
dimnames(testexF)[[3]]=colnames(fixesBRAN[[1]])[23:30]
for(i in 1:6)
  for(j in 1:8)
  {
    a=t.test(fixesBRAN[[i]][fixesBRAN[[i]]$Location==1,j+22],fixes_noeac[[i]][fixes_noeac[[i]]$Location==1,j+22])
    testex[i,1,j]=a$p.value
#     if(i<4)
#     {
#       a=t.test(fixesBRAN[[i]][fixesBRAN[[i]]$Location==1,j+22],fixes_2eac[[i]][fixes_2eac[[i]]$Location==1,j+22])
#       testex[i,2,j]=a$p.value
#     }
    a=t.test(fixesBRAN[[i]][fixesBRAN[[i]]$Location==1,j+22],fixes_2eac2[[i]][fixes_2eac2[[i]]$Location==1,j+22])
    testex[i,3,j]=a$p.value
  }

#############################################
####
### ECl matching

## 1 - Comparison of control across different WRF databases & ERAI

erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)
erai_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)

erai=erai[,-c(1,2)]
erai_E=erai_E[,-c(1,2)]
erai$Date2=as.POSIXct(paste(as.character(erai$Date),substr(erai$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

firstfixes=erai[erai$Fix==1,]
sum(firstfixes$Location)/length(firstfixes[,1])

match=array(0,c(7,6))
for(x in 1:6)
  for(y in 1:6)
  {
    matches=rep(NaN,length(eventsBRAN[[x]]$ID))
    for(i in 1:length(eventsBRAN[[x]]$ID))
    {
      tmp=fixesBRAN[[x]][(fixesBRAN[[x]]$ID==eventsBRAN[[x]]$ID[i] & fixesBRAN[[x]]$Location==1),]
      rn=range(tmp$Date2)
      I=which(fixesBRAN[[y]]$Date2<=rn[2]+(60*60*6) & fixesBRAN[[y]]$Date2>=rn[1]-(60*60*6) & fixesBRAN[[y]]$Location==1)
      if(length(I)>0)
      {
        J=unique(fixesBRAN[[y]]$ID[I])
        matches[i]=length(J) #All events that match
      } 
      
    }
    match[x,y]=length(which(!is.na(matches)))/length(matches)
  }

for(y in 1:6)
{
  matches=rep(NaN,length(erai_E$ID))
  for(i in 1:length(erai_E$ID))
  {
    tmp=erai[(erai$ID==erai_E$ID[i] & erai$Location==1),]
    rn=range(tmp$Date2)
    I=which(fixesBRAN[[y]]$Date2<=rn[2]+(60*60*6) & fixesBRAN[[y]]$Date2>=rn[1]-(60*60*6) & fixesBRAN[[y]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixesBRAN[[y]]$ID[I])
      matches[i]=length(J) #All events that match
    } 
    
  }
  match[7,y]=length(which(!is.na(matches)))/length(matches)
}

match[match==1]=NaN

mean(match[7,]) # With ERAI
mean(cbind(match[1:3,1:3],match[4:6,4:6]),na.rm=T)
mean(c(match[1,4],match[2,5],match[3,6]))


######################
##
## Matching 

match<-cv<-slp<-len<-list()
n=1
for(n in 1:6)
  {
    match[[n]]=array(NaN,c(length(eventsBRAN[[n]]$ID),5,4))
    dimnames(match[[n]])[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2","Length2")
    dimnames(match[[n]])[[3]]=c("NoEac","2EAC","2EAC2","LR")
    
    for(i in 1:length(eventsBRAN[[n]]$ID))
    {
      tmp=fixesBRAN[[n]][(fixesBRAN[[n]]$ID==eventsBRAN[[n]]$ID[i] & fixesBRAN[[n]]$Location==1),]
      rn=range(tmp$Date2)
      
      I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_noeac[[n]]$Location==1)
      if(length(I)>0)
      {
        J=unique(fixes_noeac[[n]]$ID[I])
        match[[n]][i,1,1]=length(J) #All events that match
        match[[n]][i,2,1]=length(which(fixes_noeac[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
        
        K=which(events_noeac[[n]]$ID %in% J)
        match[[n]][i,3,1]=max(events_noeac[[n]]$CV2[K])
        match[[n]][i,4,1]=min(events_noeac[[n]]$MSLP2[K])
        match[[n]][i,5,1]=min(events_noeac[[n]]$Length2[K])
      } else match[[n]][i,1,1]=0
      
      if(n<4)
      {
        I=which(fixes_2eac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_2eac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_2eac[[n]]$Location==1)
        if(length(I)>0)
        {
          J=unique(fixes_2eac[[n]]$ID[I])
          match[[n]][i,1,2]=length(J) #All events that match
          match[[n]][i,2,2]=length(which(fixes_2eac[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
          
          K=which(events_2eac[[n]]$ID %in% J)
          match[[n]][i,3,2]=max(events_2eac[[n]]$CV2[K])
          match[[n]][i,4,2]=min(events_2eac[[n]]$MSLP2[K])
          match[[n]][i,5,2]=min(events_2eac[[n]]$Length2[K])
        } else match[[n]][i,1,2]=0
      }
      
      I=which(fixes_2eac2[[n]]$Date2<=rn[2]+(60*60*6) & fixes_2eac2[[n]]$Date2>=rn[1]-(60*60*6) & fixes_2eac2[[n]]$Location==1)
      if(length(I)>0)
      {
        J=unique(fixes_2eac2[[n]]$ID[I])
        match[[n]][i,1,3]=length(J) #All events that match
        match[[n]][i,2,3]=length(which(fixes_2eac2[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
        
        K=which(events_2eac2[[n]]$ID %in% J)
        match[[n]][i,3,3]=max(events_2eac2[[n]]$CV2[K])
        match[[n]][i,4,3]=min(events_2eac2[[n]]$MSLP2[K])
        match[[n]][i,5,3]=min(events_2eac2[[n]]$Length2[K])
      }
      
      I=which(fixesLR[[n]]$Date2<=rn[2]+(60*60*6) & fixesLR[[n]]$Date2>=rn[1]-(60*60*6) & fixesLR[[n]]$Location==1)
      if(length(I)>0)
      {
        J=unique(fixesLR[[n]]$ID[I])
        match[[n]][i,1,4]=length(J) #All events that match
        match[[n]][i,2,4]=length(which(fixesLR[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
        
        K=which(eventsLR[[n]]$ID %in% J)
        match[[n]][i,3,4]=max(eventsLR[[n]]$CV2[K])
        match[[n]][i,4,4]=min(eventsLR[[n]]$MSLP2[K])
        match[[n]][i,5,4]=min(eventsLR[[n]]$Length2[K])
      }
      
    }
    
    cv[[n]]=cbind(eventsBRAN[[n]]$CV2,match[[n]][,3,])
    slp[[n]]=cbind(eventsBRAN[[n]]$MSLP2,match[[n]][,4,])
    len[[n]]=cbind(eventsBRAN[[n]]$Length2,match[[n]][,5,])
    
    n=n+1
}

### Compare - what happens to cv, slp etc when matched/unmatched?

comp=array(NaN,c(6,4,4))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("NoEAC","2EAC","2EAC2","LR")
dimnames(comp)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2")

for(n in 1:6)
  for(x in 1:4)
  {
    comp[n,x,1]=length(which(match[[n]][,1,x]>0))/length(match[[n]][,1,x])
    comp[n,x,2]=mean(match[[n]][,2,x],na.rm=T)
    comp[n,x,3]=mean(cv[[n]][,x+1]-cv[[n]][,1],na.rm=T)
    comp[n,x,4]=mean(slp[[n]][,x+1]-slp[[n]][,1],na.rm=T)
  }

comp[4:6,2,]=NaN

cvcomp=matrix(0,6,3)
dimnames(cvcomp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                        "d02 R1","d02 R2","d02 R3","All")
dimnames(cvcomp)[[2]]=c("NoEAC","2EAC","2EAC2")
for(n in 1:6)
  for(x in 1:3)
  {
    if(x==2)
    {
      if(n<=3)
      {
        b=t.test(cv[[n]][,x+1]-cv[[n]][,1])
        cvcomp[n,x]=b$p.value
      }
    } else
    {
      
      b=t.test(cv[[n]][,x+1]-cv[[n]][,1])
      cvcomp[n,x]=b$p.value
    }
  }

###### What about matched/unmatched by type - EC/SC/Mixed

Type=c("EC","SC","Mixed")

matchtypes=array(NaN,c(6,3,4))
dimnames(matchtypes)[[2]]=Type

for(n in 1:6)
  for(i in 1:3)
  {
    I=which(eventsBRAN[[n]]$HartType==Type[i])
    matchtypes[n,i,1]=length(I)
    for(x in 1:3) matchtypes[n,i,x+1]=length(which(match[[n]][I,1,x]>0))
  }

### Histogram of matched/unmatched events

cvthresh=seq(1,4,0.25)
dens2<-array(0,c(12,6,2))
for(r in 1:6)
  for(i in 1:12)
  {
    tmp=cv[[r]]
    I=which(is.na(tmp[,2]))
    J=which(tmp[I,1]>=cvthresh[i] & tmp[I,1]<cvthresh[i+1])
    dens2[i,r,1]=length(J)
    J=which(tmp[-I,1]>=cvthresh[i] & tmp[-I,1]<cvthresh[i+1])
    dens2[i,r,2]=length(J)
  }

a=apply(dens2[1:4,,],c(2,3),sum)
mean(a[,1]/(a[,1]+a[,2]))

###Stacked histogram
pdf(file=paste(figdir,"/ECL_match_SSTchange_CVdist.pdf",sep=""),width=6,height=3.5,pointsize=10)
par(mar=c(4,4,2,2))
counts <- apply(dens2,c(1,3),mean,na.rm=T)
colnames(counts)=c("Unmatched","Matched")
mp<-barplot(t(counts), col=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)),
        ylab="Number of ECLs",xlab="Intensity (hPa. (deg lat)^2))",legend = colnames(counts))
axis(1,at=c(mp-0.6,mp[12]+0.6),cvthresh)
dev.off()

##########matching extremes?

matchE<-list()
n=1
for(n in 1:6)
{
  matchE[[n]]=array(NaN,c(length(eventsBRAN[[n]]$CV2),10,3))
  matchE[[n]][,,1]=as.matrix(eventsBRAN[[n]][,18:27])
    
  for(i in 1:length(eventsBRAN[[n]]$ID))
  {
    tmp=fixesBRAN[[n]][(fixesBRAN[[n]]$ID==eventsBRAN[[n]]$ID[i] & fixesBRAN[[n]]$Location==1),]
    rn=range(tmp$Date2)
    
    I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_noeac[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_noeac[[n]]$ID[I])
      K=which(events_noeac[[n]]$ID %in% J)
      matchE[[n]][i,,2]=apply(events_noeac[[n]][K,18:27],2,max)
    } 
    
    I=which(fixes_2eac2[[n]]$Date2<=rn[2]+(60*60*6) & fixes_2eac2[[n]]$Date2>=rn[1]-(60*60*6) & fixes_2eac2[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_2eac2[[n]]$ID[I])
      K=which(events_2eac2[[n]]$ID %in% J)
      matchE[[n]][i,,3]=apply(events_2eac2[[n]][K,18:27],2,max)
    }
  }
  n=n+1
}

Rdiff<-array(0,c(6,10,2))
dimnames(Rdiff)[[2]]<-colnames(eventsBRAN[[1]])[18:27]
dimnames(Rdiff)[[3]]<-c("NoEac","2EAC")
for(n in 1:6)
{
  Rdiff[n,,1]=apply(matchE[[n]][,,2]-matchE[[n]][,,1],2,mean,na.rm=T)
  Rdiff[n,,2]=apply(matchE[[n]][,,3]-matchE[[n]][,,1],2,mean,na.rm=T)
}
apply(Rdiff,c(2,3),mean)

Rtest<-Rdiff
for(n in 1:6)
  for(i in 1:10)
{
  a=t.test(matchE[[n]][,i,2]-matchE[[n]][,i,1])
  Rtest[n,i,1]=a$p.value
  a=t.test(matchE[[n]][,i,3]-matchE[[n]][,i,1])
  Rtest[n,i,2]=a$p.value
}


############ Comparing to the cv0.5 case

match2<-cv2<-slp2<-len2<-list()
n=1
for(dom in c("d01","d02"))
  for(r in 1:3) 
{
  match2[[n]]=array(NaN,c(length(eventsBRAN[[n]]$ID),5,2))
  dimnames(match2[[n]])[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2","Length2")
  dimnames(match2[[n]])[[3]]=c("NoEac","2EAC")
  
  for(i in 1:length(eventsBRAN[[n]]$ID))
  {
    tmp=fixesBRAN[[n]][(fixesBRAN[[n]]$ID==eventsBRAN[[n]]$ID[i] & fixesBRAN[[n]]$Location==1),]
    rn=range(tmp$Date2)
    
    fixes=rbind(read.csv(paste("../ECLfixes_",dom,"_2007_R",r,"_BRAN_noeac_",cat[5],".csv",sep="")),
                     read.csv(paste("../ECLfixes_",dom,"_2008_R",r,"_BRAN_noeac_",cat[5],".csv",sep="")))
    fixes$ID=floor(fixes$Date/10000)*1000+fixes$ID
    fixes$Date2=as.POSIXct(paste(as.character(fixes$Date),substr(fixes$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events=rbind(read.csv(paste("../ECLevents_",dom,"_2007_R",r,"_BRAN_noeac_",cat[5],".csv",sep="")),
                read.csv(paste("../ECLevents_",dom,"_2008_R",r,"_BRAN_noeac_",cat[5],".csv",sep="")))  
    events$ID=floor(events$Date1/10000)*1000+events$ID
    
    I=which(fixes$Date2<=rn[2]+(60*60*6) & fixes$Date2>=rn[1]-(60*60*6) & fixes$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes$ID[I])
      match2[[n]][i,1,1]=length(J) #All events that match
      match2[[n]][i,2,1]=length(which(fixes$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
      
      K=which(events$ID %in% J)
      match2[[n]][i,3,1]=max(events$CV2[K])
      match2[[n]][i,4,1]=min(events$MSLP2[K])
      match2[[n]][i,5,1]=min(events$Length2[K])
    } else match2[[n]][i,1,1]=0
    
  }
  
  cv2[[n]]=cbind(eventsBRAN[[n]]$CV2,match2[[n]][,3,])
  slp2[[n]]=cbind(eventsBRAN[[n]]$MSLP2,match2[[n]][,4,])
  len2[[n]]=cbind(eventsBRAN[[n]]$Length2,match2[[n]][,5,])
  
  n=n+1
}

comp2=array(NaN,c(6,2,4))
dimnames(comp2)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp2)[[2]]=c("All","Unmatched by CV1")
dimnames(comp2)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2")

for(n in 1:6)
  {
    comp2[n,1,1]=length(which(match2[[n]][,1,1]>0))/length(match2[[n]][,1,1])
    comp2[n,1,2]=mean(match2[[n]][,2,1],na.rm=T)
    comp2[n,1,3]=mean(cv2[[n]][,2]-cv2[[n]][,1],na.rm=T)
    comp2[n,1,4]=mean(slp2[[n]][,2]-slp2[[n]][,1],na.rm=T)
    
    I=which(match[[n]][,1,1]==0)
    comp2[n,2,1]=length(which(match2[[n]][I,1,1]>0))/length(match2[[n]][I,1,1])
    comp2[n,2,2]=mean(match2[[n]][I,2,1],na.rm=T)
    comp2[n,2,3]=mean(cv2[[n]][I,2]-cv2[[n]][I,1],na.rm=T)
    comp2[n,2,4]=mean(slp2[[n]][I,2]-slp2[[n]][I,1],na.rm=T)
}

apply(comp2,c(2,3),mean)

#############
## What about average cv change as a function of cv?

cvthresh=c(1,1.5,5)

cvchange=array(NaN,c(6,2,7))
dimnames(cvchange)[[2]]=cvthresh[1:3]
dimnames(cvchange)[[3]]=c("Count","NoEAC Match","NoEAC CV change","2EAC Match","2EAC CV change","2EAC2 Match","2EAC2 CV change")
for(n in 1:6)
  for(j in 1:2)
  {
    I=which(cv[[n]][,1]>=cvthresh[j] & cv[[n]][,1]<cvthresh[j+1])
    cvchange[n,j,1]=length(I)
    
    cvchange[n,j,2]=length(which(!is.na(cv[[n]][I,2])))
    cvchange[n,j,3]=mean(cv[[n]][I,2]-cv[[n]][I,1],na.rm=T)
    if(n<4)
    {
    cvchange[n,j,4]=length(which(!is.na(cv[[n]][I,3])))
    cvchange[n,j,5]=mean(cv[[n]][I,3]-cv[[n]][I,1],na.rm=T)
    }
    cvchange[n,j,6]=length(which(!is.na(cv[[n]][I,4])))
    cvchange[n,j,7]=mean(cv[[n]][I,4]-cv[[n]][I,1],na.rm=T)
  }


#################
## What about change relative to GV - only for d01 though as GV region is big

eventsBRAN<-events_noeac<-events_2eac<-fixesBRAN<-fixes_noeac<-fixes_2eac<-fixes_2eac2<-events_2eac2<-list()

n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    eventsBRAN[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
    eventsBRAN[[n]]$Year=floor(eventsBRAN[[n]]$Date1/10000)
    eventsBRAN[[n]]$Month=floor(eventsBRAN[[n]]$Date1/100)%%100
    
    events_noeac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
    events_noeac[[n]]$Year=floor(events_noeac[[n]]$Date1/10000)
    events_noeac[[n]]$Month=floor(events_noeac[[n]]$Date1/100)%%100
    
    fixesBRAN[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
    fixesBRAN[[n]]$Year=floor(fixesBRAN[[n]]$Date/10000)
    fixesBRAN[[n]]$Month=floor(fixesBRAN[[n]]$Date/100)%%100
    fixesBRAN[[n]]$Location2<-0
    I<-which(fixesBRAN[[n]][,7]>=149 & fixesBRAN[[n]][,7]<=154 & fixesBRAN[[n]][,8]<(-37) & fixesBRAN[[n]][,8]>=-41)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=(149+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,7]<=(154+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,8]<(-31) & fixesBRAN[[n]][,8]>=-37)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=152 & fixesBRAN[[n]][,7]<=157 & fixesBRAN[[n]][,8]<=(-24) & fixesBRAN[[n]][,8]>=-31)
    fixesBRAN[[n]]$Location2[I]<-1
    fixesBRAN[[n]]$Date2=as.POSIXct(paste(as.character(fixesBRAN[[n]]$Date),substr(fixesBRAN[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    
    fixes_noeac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
    fixes_noeac[[n]]$Year=floor(fixes_noeac[[n]]$Date/10000)
    fixes_noeac[[n]]$Month=floor(fixes_noeac[[n]]$Date/100)%%100
    fixes_noeac[[n]]$Location2<-0
    I<-which(fixes_noeac[[n]][,7]>=149 & fixes_noeac[[n]][,7]<=154 & fixes_noeac[[n]][,8]<(-37) & fixes_noeac[[n]][,8]>=-41)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=(149+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,7]<=(154+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,8]<(-31) & fixes_noeac[[n]][,8]>=-37)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=152 & fixes_noeac[[n]][,7]<=157 & fixes_noeac[[n]][,8]<=(-24) & fixes_noeac[[n]][,8]>=-31)
    fixes_noeac[[n]]$Location2[I]<-1
    fixes_noeac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_noeac[[n]]$Date),substr(fixes_noeac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    
    if(dom=="d01")
    {
      fixes_2eac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
      fixes_2eac[[n]]$Year=floor(fixes_2eac[[n]]$Date/10000)
      fixes_2eac[[n]]$Month=floor(fixes_2eac[[n]]$Date/100)%%100
      fixes_2eac[[n]]$Location2<-0
      I<-which(fixes_2eac[[n]][,7]>=149 & fixes_2eac[[n]][,7]<=154 & fixes_2eac[[n]][,8]<(-37) & fixes_2eac[[n]][,8]>=-41)
      fixes_2eac[[n]]$Location2[I]<-1
      I<-which(fixes_2eac[[n]][,7]>=(149+(37+fixes_2eac[[n]][,8])/2) & fixes_2eac[[n]][,7]<=(154+(37+fixes_2eac[[n]][,8])/2) & fixes_2eac[[n]][,8]<(-31) & fixes_2eac[[n]][,8]>=-37)
      fixes_2eac[[n]]$Location2[I]<-1
      I<-which(fixes_2eac[[n]][,7]>=152 & fixes_2eac[[n]][,7]<=157 & fixes_2eac[[n]][,8]<=(-24) & fixes_2eac[[n]][,8]>=-31)
      fixes_2eac[[n]]$Location2[I]<-1
      fixes_2eac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_2eac[[n]]$Date),substr(fixes_2eac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
      
      
      events_2eac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
      events_2eac[[n]]$Year=floor(events_2eac[[n]]$Date1/10000)
      events_2eac[[n]]$Month=floor(events_2eac[[n]]$Date1/100)%%100  
    }
    
    fixes_2eac2[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_GV.csv",sep=""),stringsAsFactors = F)
    fixes_2eac2[[n]]$Year=floor(fixes_2eac2[[n]]$Date/10000)
    fixes_2eac2[[n]]$Month=floor(fixes_2eac2[[n]]$Date/100)%%100
    fixes_2eac2[[n]]$Location2<-0
    I<-which(fixes_2eac2[[n]][,7]>=149 & fixes_2eac2[[n]][,7]<=154 & fixes_2eac2[[n]][,8]<(-37) & fixes_2eac2[[n]][,8]>=-41)
    fixes_2eac2[[n]]$Location2[I]<-1
    I<-which(fixes_2eac2[[n]][,7]>=(149+(37+fixes_2eac2[[n]][,8])/2) & fixes_2eac2[[n]][,7]<=(154+(37+fixes_2eac2[[n]][,8])/2) & fixes_2eac2[[n]][,8]<(-31) & fixes_2eac2[[n]][,8]>=-37)
    fixes_2eac2[[n]]$Location2[I]<-1
    I<-which(fixes_2eac2[[n]][,7]>=152 & fixes_2eac2[[n]][,7]<=157 & fixes_2eac2[[n]][,8]<=(-24) & fixes_2eac2[[n]][,8]>=-31)
    fixes_2eac2[[n]]$Location2[I]<-1
    fixes_2eac2[[n]]$Date2=as.POSIXct(paste(as.character(fixes_2eac2[[n]]$Date),substr(fixes_2eac2[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    
    events_2eac2[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_GV.csv",sep=""),stringsAsFactors = F)
    events_2eac2[[n]]$Year=floor(events_2eac2[[n]]$Date1/10000)
    events_2eac2[[n]]$Month=floor(events_2eac2[[n]]$Date1/100)%%100  
    n=n+1     
  }

match<-cv<-slp<-len<-gv<-list()
n=1
for(n in 1:3)
{
  match[[n]]=array(NaN,c(length(eventsBRAN[[n]]$ID),6,3))
  dimnames(match[[n]])[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2","Length2","GV")
  dimnames(match[[n]])[[3]]=c("NoEac","2EAC","2EAC2")
  
  for(i in 1:length(eventsBRAN[[n]]$ID))
  {
    tmp=fixesBRAN[[n]][(fixesBRAN[[n]]$ID==eventsBRAN[[n]]$ID[i] & fixesBRAN[[n]]$Location==1),]
    rn=range(tmp$Date2)
    
    I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_noeac[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_noeac[[n]]$ID[I])
      match[[n]][i,1,1]=length(J) #All events that match
      match[[n]][i,2,1]=length(which(fixes_noeac[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
      
      K=which(events_noeac[[n]]$ID %in% J)
      match[[n]][i,3,1]=max(events_noeac[[n]]$CV2[K])
      match[[n]][i,4,1]=min(events_noeac[[n]]$MSLP2[K])
      match[[n]][i,5,1]=min(events_noeac[[n]]$Length2[K])
      match[[n]][i,6,1]=max(events_noeac[[n]]$GV[K])
    } else match[[n]][i,1,1]=0
    
    if(n<4)
    {
      I=which(fixes_2eac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_2eac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_2eac[[n]]$Location==1)
      if(length(I)>0)
      {
        J=unique(fixes_2eac[[n]]$ID[I])
        match[[n]][i,1,2]=length(J) #All events that match
        match[[n]][i,2,2]=length(which(fixes_2eac[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
        
        K=which(events_2eac[[n]]$ID %in% J)
        match[[n]][i,3,2]=max(events_2eac[[n]]$CV2[K])
        match[[n]][i,4,2]=min(events_2eac[[n]]$MSLP2[K])
        match[[n]][i,5,2]=min(events_2eac[[n]]$Length2[K])
        match[[n]][i,6,2]=max(events_2eac[[n]]$GV[K])
      } else match[[n]][i,1,2]=0
    }
    
    I=which(fixes_2eac2[[n]]$Date2<=rn[2]+(60*60*6) & fixes_2eac2[[n]]$Date2>=rn[1]-(60*60*6) & fixes_2eac2[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_2eac2[[n]]$ID[I])
      match[[n]][i,1,3]=length(J) #All events that match
      match[[n]][i,2,3]=length(which(fixes_2eac2[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
      
      K=which(events_2eac2[[n]]$ID %in% J)
      match[[n]][i,3,3]=max(events_2eac2[[n]]$CV2[K])
      match[[n]][i,4,3]=min(events_2eac2[[n]]$MSLP2[K])
      match[[n]][i,5,3]=min(events_2eac2[[n]]$Length2[K])
      match[[n]][i,6,3]=max(events_2eac2[[n]]$GV[K])
    } else match[[n]][i,1,3]=0
  }
  
  cv[[n]]=cbind(eventsBRAN[[n]]$CV2,match[[n]][,3,])
  slp[[n]]=cbind(eventsBRAN[[n]]$MSLP2,match[[n]][,4,])
  len[[n]]=cbind(eventsBRAN[[n]]$Length2,match[[n]][,5,])
  gv[[n]]=cbind(eventsBRAN[[n]]$GV,eventsBRAN[[n]]$GVthresh,match[[n]][,6,])
  
  n=n+1
}


comp=array(NaN,c(3,2,4))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3")
dimnames(comp)[[2]]=c("NoEAC","2EAC")
dimnames(comp)[[3]]=c("MatchEvents","Match NotGV","Match GV","GV diff")

for(n in 1:3)
  for(x in 1:2)
  {
    comp[n,x,1]=length(which(gv[[n]][,x+2]>0))/length(gv[[n]][,x+2])
    comp[n,x,2]=length(which(gv[[n]][,x+2]>0 & gv[[n]][,2]==0))/length(which(gv[[n]][,2]==0))
    comp[n,x,3]=length(which(gv[[n]][,x+2]>0 & gv[[n]][,2]==1))/length(which(gv[[n]][,2]==1))
    comp[n,x,4]=mean(gv[[n]][,x+2]-gv[[n]][,1],na.rm=T)
  }

### Histogram of matched/unmatched events

gvthresh=seq(0,60,5)
dens2<-array(0,c(12,6,2))
for(r in 1:6)
  for(i in 1:12)
  {
    tmp=gv[[r]]
    I=which(is.na(tmp[,2]))
    J=which(tmp[I,1]>=gvthresh[i] & tmp[I,1]<gvthresh[i+1])
    dens2[i,r,1]=length(J)
    J=which(tmp[-I,1]>=gvthresh[i] & tmp[-I,1]<gvthresh[i+1])
    dens2[i,r,2]=length(J)
  }

a=apply(dens2[1:4,,],c(2,3),sum)
mean(a[,1]/(a[,1]+a[,2]))

counts <- apply(dens2,c(1,3),mean,na.rm=T)
colnames(counts)=c("Unmatched","Matched")
mp<-barplot(t(counts), col=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)),
            ylab="Number of ECLs",xlab="Maximum geostrophic vorticity)",legend = colnames(counts))
axis(1,at=c(mp-0.6,mp[12]+0.6),gvthresh)

########
## Overall changes in GV distribution?

gvthresh=c(0,20,Inf)
count<-array(0,c(6,3,2))

for(i in 1:6)
  for(j in 1:2)
{
  count[i,1,j]=length(which(eventsBRAN[[i]]$GV>=gvthresh[j] & eventsBRAN[[i]]$GV<gvthresh[j+1]))
  count[i,2,j]=length(which(events_noeac[[i]]$GV>=gvthresh[j] & events_noeac[[i]]$GV<gvthresh[j+1]))
  count[i,3,j]=length(which(events_2eac2[[i]]$GV>=gvthresh[j] & events_2eac2[[i]]$GV<gvthresh[j+1]))
  }

a=apply(count[,,1:2],c(1,2),sum)
b=apply(count[,,3:6],c(1,2),sum)

apply(count[,2,]/count[,1,],2,mean)

cvthresh=c(1,2,Inf)
count<-array(0,c(3,3,2))

for(i in 1:3)
  for(j in 1:2)
  {
    count[i,1,j]=length(which(eventsBRAN[[i]]$CV2>=cvthresh[j] & eventsBRAN[[i]]$CV2<cvthresh[j+1]))
    count[i,2,j]=length(which(events_noeac[[i]]$CV2>=cvthresh[j] & events_noeac[[i]]$CV2<cvthresh[j+1]))
    count[i,3,j]=length(which(events_2eac[[i]]$CV2>=cvthresh[j] & events_2eac[[i]]$CV2<cvthresh[j+1]))
  }

apply(count[,3,]/count[,1,],2,mean)

###
## Combined thresh - GV & CV

cvthresh=c(1,1.5,Inf)
gvthresh=c(0,15,Inf)

match<-count<-array(0,c(3,2,2))
dimnames(match)[[2]]<-dimnames(count)[[2]]<-c("GV<15","GV>15")
dimnames(match)[[3]]<-dimnames(count)[[3]]<-c("CV<1.5","CV>1.5")

for(n in 1:3)
for(i in 1:2)
  for(j in 1:2)
  {
    I=which(gv[[n]][,1]>=gvthresh[i] & cv[[n]][,1]>=cvthresh[j] & gv[[n]][,1]<gvthresh[i+1] & cv[[n]][,1]<cvthresh[j+1])
    match[n,i,j]=length(which(!is.na(cv[[n]][I,2])))
    count[n,i,j]=length(I)
  }

a=(match[,1,2]+match[,2,1])/(count[,1,2]+count[,2,1])


####Compare GV percentile

######### Intial stuff
rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing')
library(abind)
library("R.matlab")
library(fields)
library(maps)
figdir="~/Documents/ECLs/WRFruns/0708/EACpaper"

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")
c=3 ## Which version do I want? 3 is the default



############# Doing Fei's GV
###
###

## First, make 9-point smoothed running average GV

wrfdirs=c("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN_2eac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_2eac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN_2eac/out/")

cat="rad2_p100"

GV_p90<-GV_count20<-matrix(0,15,4)
dom=c("d01","d02")
dates=seq.POSIXt(from=as.POSIXct("2007010100",format="%Y%m%d%H",tz="GMT"),
                 to=as.POSIXct("2008123118",format="%Y%m%d%H",tz="GMT"),
                 by="6 hours")

for(w in 1:length(wrfdirs))
  for(d in 1:2)
  {
    tmp=read.csv(paste(wrfdirs[w],"GV_6hrly_timeseries.txt",sep=""),header=F)
    GV=rep(NaN,length(tmp[,1]))
    for(i in 5:(length(GV)-4)) GV[i]=mean(tmp[(i-4):(i+4),1])
    GV_p90[w,d]=quantile(GV,0.9,na.rm=T)
    GV_count20[w,d]=length(which(GV>=20))
    
    dayGV=aggregate(GV,by=list(day = cut(dates, "days", right = TRUE)),max)
    GV_p90[w,d+2]=quantile(dayGV[,2],0.9,na.rm=T)
    GV_count20[w,d+2]=length(which(dayGV[,2]>=20))
  }

ks_pval<-matrix(0,6,2)
dom=c("","_d02")

n=1
for(w in 1:3)
  for(d in 1:2)
  {
    gv_cont=read.csv(paste(wrfdirs[w+6],"GV_6hrly_timeseries",dom[d],".txt",sep=""),header=F)
    gv_noeac=read.csv(paste(wrfdirs[w+9],"GV_6hrly_timeseries",dom[d],".txt",sep=""),header=F)
    gv_2eac=read.csv(paste(wrfdirs[w+12],"GV_6hrly_timeseries",dom[d],".txt",sep=""),header=F)
    
    a=ks.test(gv_cont[,1],gv_noeac[,1])
    ks_pval[n,1]=a$p.value
    a=ks.test(gv_cont[,1],gv_2eac[,1])
    ks_pval[n,2]=a$p.value
    n=n+1
  }


####### Playing with statistical significance
se <- function(x) sqrt(var(x)/length(x))
tmp=count[1:3,4]-count[1:3,1]

count2<-array(0,c(6,3,2))
year=c(2007,2008)
n=1
for(y in 1:2)
  for(i in 1:3)
  {
    count2[n,1,1]=length(eventsBRAN[[i]]$ID[eventsBRAN[[i]]$Year==year[y]])
    count2[n,2,1]=length(events_noeac[[i]]$ID[events_noeac[[i]]$Year==year[y]])
    count2[n,3,1]=length(events_2eac2[[i]]$ID[events_2eac2[[i]]$Year==year[y]])
    count2[n,1,2]=length(eventsBRAN[[i+3]]$ID[eventsBRAN[[i+3]]$Year==year[y]])
    count2[n,2,2]=length(events_noeac[[i+3]]$ID[events_noeac[[i+3]]$Year==year[y]])
    count2[n,3,2]=length(events_2eac2[[i+3]]$ID[events_2eac2[[i+3]]$Year==year[y]])
    n=n+1
  }

tmp=count2[,2,]/count2[,1,]
mean(tmp)
se(tmp)
c(mean(tmp)-2*se(tmp),mean(tmp)+2*se(tmp))
t.test(tmp)
wilcox.test(tmp)

################
######## Change in ECL rainfall - all ECLs (unmatched)

###Takes in a data column
### If th

library(boot)
bthresh <- function(data, indices,thresh) {
  return(length(which(data[indices]>=thresh)))
} 
a=boot(eventsBRAN[[1]]$MaxMeanRain250,bthresh,1000,thresh=6)
boot.ci(a)
b=boot(events_2eac[[1]]$MaxMeanRain250,bthresh,1000,thresh=6)
boot.ci(b)

bmean <- function(data, indices,thresh) {
  return(length(which(data[indices]>=thresh)))
} 
a=boot(eventsBRAN[[1]]$MaxMeanRain250,bthresh,1000,thresh=6)
boot.ci(a)
b=boot(events_2eac[[1]]$MaxMeanRain250,bthresh,1000,thresh=6)
boot.ci(b)


######### Matching vs 2EAC

match<-list()
n=1
for(n in 1:6)
{
  match[[n]]=array(NaN,c(length(events_2eac2[[n]]$ID),5,2))
  dimnames(match[[n]])[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2","Length2")
  dimnames(match[[n]])[[3]]=c("NoEac","BRAN")
  
  for(i in 1:length(events_2eac2[[n]]$ID))
  {
    tmp=fixes_2eac2[[n]][(fixes_2eac2[[n]]$ID==events_2eac2[[n]]$ID[i] & fixes_2eac2[[n]]$Location==1),]
    rn=range(tmp$Date2)
    
    I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_noeac[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_noeac[[n]]$ID[I])
      match[[n]][i,1,1]=length(J) #All events that match
      match[[n]][i,2,1]=length(which(fixes_noeac[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
      
      K=which(events_noeac[[n]]$ID %in% J)
      match[[n]][i,3,1]=max(events_noeac[[n]]$CV2[K])
      match[[n]][i,4,1]=min(events_noeac[[n]]$MSLP2[K])
      match[[n]][i,5,1]=min(events_noeac[[n]]$Length2[K])
    } else match[[n]][i,1,1]=0
    
    I=which(fixesBRAN[[n]]$Date2<=rn[2]+(60*60*6) & fixesBRAN[[n]]$Date2>=rn[1]-(60*60*6) & fixesBRAN[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixesBRAN[[n]]$ID[I])
      match[[n]][i,1,2]=length(J) #All events that match
      match[[n]][i,2,2]=length(which(fixesBRAN[[n]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
      
      K=which(eventsBRAN[[n]]$ID %in% J)
      match[[n]][i,3,2]=max(eventsBRAN[[n]]$CV2[K])
      match[[n]][i,4,2]=min(eventsBRAN[[n]]$MSLP2[K])
      match[[n]][i,5,2]=min(eventsBRAN[[n]]$Length2[K])
    } else match[[n]][i,1,2]=0
  }
  
  n=n+1
}

### Compare - what happens to cv, slp etc when matched/unmatched?

comp=array(NaN,c(6,2))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("NoEAC","BRAN")

for(n in 1:6)
  for(x in 1:2)
    comp[n,x]=mean(match[[n]][,1,x]>0)

############ Comparing to the cv0.5 case

match2<-cv2<-slp2<-len2<-list()
n=1
for(dom in c("d01","d02"))
  for(r in 1:3) 
  {
    match2[[n]]=array(NaN,c(length(events_2eac2[[n]]$ID),5,2))
    dimnames(match2[[n]])[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2","Length2")
    dimnames(match2[[n]])[[3]]=c("BRAN","NoEAC")
    
    for(i in 1:length(events_2eac2[[n]]$ID))
    {
      tmp=fixes_2eac2[[n]][(fixes_2eac2[[n]]$ID==events_2eac2[[n]]$ID[i] & fixes_2eac2[[n]]$Location==1),]
      rn=range(tmp$Date2)
      
      fixes=rbind(read.csv(paste("../ECLfixes_",dom,"_2007_R",r,"_BRAN_",cat[5],".csv",sep="")),
                  read.csv(paste("../ECLfixes_",dom,"_2008_R",r,"_BRAN_",cat[5],".csv",sep="")))
      fixes$ID=floor(fixes$Date/10000)*1000+fixes$ID
      fixes$Date2=as.POSIXct(paste(as.character(fixes$Date),substr(fixes$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
      events=rbind(read.csv(paste("../ECLevents_",dom,"_2007_R",r,"_BRAN_",cat[5],".csv",sep="")),
                   read.csv(paste("../ECLevents_",dom,"_2008_R",r,"_BRAN_",cat[5],".csv",sep="")))  
      events$ID=floor(events$Date1/10000)*1000+events$ID
      
      I=which(fixes$Date2<=rn[2]+(60*60*6) & fixes$Date2>=rn[1]-(60*60*6) & fixes$Location==1)
      if(length(I)>0)
      {
        J=unique(fixes$ID[I])
        match2[[n]][i,1,1]=length(J) #All events that match
        match2[[n]][i,2,1]=length(which(fixes$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
        
        K=which(events$ID %in% J)
        match2[[n]][i,3,1]=max(events$CV2[K])
        match2[[n]][i,4,1]=min(events$MSLP2[K])
        match2[[n]][i,5,1]=min(events$Length2[K])
      } else match2[[n]][i,1,1]=0
      
      
    }
    
    cv2[[n]]=cbind(events_2eac2[[n]]$CV2,match2[[n]][,3,])
    slp2[[n]]=cbind(events_2eac2[[n]]$MSLP2,match2[[n]][,4,])
    len2[[n]]=cbind(events_2eac2[[n]]$Length2,match2[[n]][,5,])
    
    n=n+1
  }

comp2=array(NaN,c(6,2,4))
dimnames(comp2)[[1]]=c("d01 R1","d01 R2","d01 R3",
                       "d02 R1","d02 R2","d02 R3")
dimnames(comp2)[[2]]=c("All","Unmatched by CV1")
dimnames(comp2)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2")

for(n in 1:6)
{
  comp2[n,1,1]=length(which(match2[[n]][,1,1]>0))/length(match2[[n]][,1,1])
  comp2[n,1,2]=mean(match2[[n]][,2,1],na.rm=T)
  comp2[n,1,3]=mean(cv2[[n]][,2]-cv2[[n]][,1],na.rm=T)
  comp2[n,1,4]=mean(slp2[[n]][,2]-slp2[[n]][,1],na.rm=T)
  
  I=which(match[[n]][,1,1]==0)
  comp2[n,2,1]=length(which(match2[[n]][I,1,1]>0))/length(match2[[n]][I,1,1])
  comp2[n,2,2]=mean(match2[[n]][I,2,1],na.rm=T)
  comp2[n,2,3]=mean(cv2[[n]][I,2]-cv2[[n]][I,1],na.rm=T)
  comp2[n,2,4]=mean(slp2[[n]][I,2]-slp2[[n]][I,1],na.rm=T)
}

apply(comp2,c(2,3),mean)
