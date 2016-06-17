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

eventsBRAN<-events_noeac<-events_2eac<-fixesBRAN<-fixes_noeac<-fixes_2eac<-eventsLR<-fixesLR<-list()

n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    eventsBRAN[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    eventsBRAN[[n]]$Year=floor(eventsBRAN[[n]]$Date1/10000)
    eventsBRAN[[n]]$Month=floor(eventsBRAN[[n]]$Date1/100)%%100
    eventsBRAN[[n]][eventsBRAN[[n]]==-Inf]=NaN
    
    fixesBRAN[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
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
    
    eventsLR[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    eventsLR[[n]]$Year=floor(eventsLR[[n]]$Date1/10000)
    eventsLR[[n]]$Month=floor(eventsLR[[n]]$Date1/100)%%100
    eventsLR[[n]][eventsLR[[n]]==-Inf]=NaN    
    
    fixesLR[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    fixesLR[[n]]$Year=floor(fixesLR[[n]]$Date/10000)
    fixesLR[[n]]$Month=floor(fixesLR[[n]]$Date/100)%%100
    fixesLR[[n]]$Date2=as.POSIXct(paste(as.character(fixesLR[[n]]$Date),substr(fixesLR[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    events_noeac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    events_noeac[[n]]$Year=floor(events_noeac[[n]]$Date1/10000)
    events_noeac[[n]]$Month=floor(events_noeac[[n]]$Date1/100)%%100
    events_noeac[[n]][events_noeac[[n]]==-Inf]=NaN
    
    fixes_noeac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
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
    
    events_2eac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
    events_2eac[[n]]$Year=floor(events_2eac[[n]]$Date1/10000)
    events_2eac[[n]]$Month=floor(events_2eac[[n]]$Date1/100)%%100  
    events_2eac[[n]][events_2eac[[n]]==-Inf]=NaN
    
    fixes_2eac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_impactsC2.csv",sep=""),stringsAsFactors = F)
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
    
    n=n+1     
  }

##########

######## Figure 2
####### Changes in total number of events

count=matrix(NaN,6,4)
colnames(count)=c("NoEAC","Control","2EAC","Low Resolution")
rownames(count)=c("R1 50km","R1 20km","R1 30km","R1 10km","R2 10km","R3 10km")
cvthresh=NaN ## Change this if I want a minimum CV threshold
countthresh=NaN ## Change this if I want to use a CV threshold that gives X events p.a. for the BRAN (control) case

for(i in 1:6)
{
  if(!is.na(cvthresh))
  {
    count[i,1]=length(which(events_noeac[[i]]$CV2>=cvthresh))
    count[i,2]=length(which(eventsBRAN[[i]]$CV2>=cvthresh))
    count[i,3]=length(which(events_2eac[[i]]$CV2>=cvthresh))
    count[i,4]=length(which(eventsLR[[i]]$CV2>=cvthresh))
  } else if(!is.na(countthresh)) {
    a=order(eventsBRAN[[i]]$CV2,decreasing=T)
    if(length(a)>=countthresh) b=eventsBRAN[[i]]$CV2[a[countthresh]] else b=min(eventsBRAN[[i]]$CV2)
    
    count[i,1]=length(which(events_noeac[[i]]$CV2>=b))
    count[i,2]=length(which(eventsBRAN[[i]]$CV2>=b))
    count[i,3]=length(which(events_2eac[[i]]$CV2>=b))
    count[i,4]=length(which(eventsLR[[i]]$CV2>=b))
  } else {
    count[i,1]=length(events_noeac[[i]]$ID)
    count[i,2]=length(eventsBRAN[[i]]$ID)
    count[i,3]=length(events_2eac[[i]]$ID)
    count[i,4]=length(eventsLR[[i]]$ID)
  }
}

##Changes in ECL frequency for thresholds above

change=array(0,c(3,2))
dimnames(change)[[1]]=c("NoEAC","2EAC","LR")
dimnames(change)[[2]]=c("Median % change","t.test")

tmp=c(1,3,4)
for(j in 1:3)
  {
    change[j,1]=100*(median(count[,tmp[j]]/count[,2])-1)
    a=t.test(count[,tmp[j]],count[,2])
    change[j,2]=a$p.value # Yes I know this is mostly meaningless for effective n=3
  }

### Make Figure 2

pnt=c(1,2,4,1,2,4)
cm=c("black","black","black","darkgrey","darkgrey","darkgrey")

pdf(paste(figdir,"ECLcount_change.pdf",sep=""),width=4,height=4)
par(mar=c(3,4,1,1))
plot(1:3,count[1,1:3],"b",axes=F,ylim=c(20,80),lwd=2,cex=1.5,pch=pnt[1],col=cm[1],
     xlab="",ylab="Number of ECLs")
axis(1,at=1:3,labels=c("No EAC","Control","2EAC"))
axis(2)
for(i in 2:6) lines(1:3,count[i,1:3],"b",lwd=2,cex=1.5,pch=pnt[i],col=cm[i])
points(rep(2,6),count[,4],pch=pnt,col=cm,cex=1.5,lwd=2)
legend("topleft",c("R1","R2","R3","50km","10km"),pch=c(1,2,4,1,1),pt.lwd=2,col=c(1,1,1,1,"darkgrey"),pt.cex=1.5,bty="n",ncol=2)
dev.off()

############## Figure 3
########## Changes in events by location

lat=seq(-60,-10,5)
lon=seq(100,180,5)

#### Location of ECLs when in the region
#### Edit slightly - needs to be unique time (in case two low centre same cell, unlikely)

loc<-array(NaN,c(11,17,6,3))
for(i in 1:11)
  for(j in 1:17)
    for(k in 1:6)
    {
      I=which(fixes_noeac[[k]]$Lat>=lat[i]-2.5 & fixes_noeac[[k]]$Lat<lat[i]+2.5 & fixes_noeac[[k]]$Lon>=lon[j]-2.5 & fixes_noeac[[k]]$Lon<lon[j]+2.5 & fixes_noeac[[k]]$Location==1)
      loc[i,j,k,1]=length(I)
      I=which(fixesBRAN[[k]]$Lat>=lat[i]-2.5 & fixesBRAN[[k]]$Lat<lat[i]+2.5 & fixesBRAN[[k]]$Lon>=lon[j]-2.5 & fixesBRAN[[k]]$Lon<lon[j]+2.5 & fixesBRAN[[k]]$Location==1)
      loc[i,j,k,2]=length(I)
      I=which(fixes_2eac[[k]]$Lat>=lat[i]-2.5 & fixes_2eac[[k]]$Lat<lat[i]+2.5 & fixes_2eac[[k]]$Lon>=lon[j]-2.5 & fixes_2eac[[k]]$Lon<lon[j]+2.5 & fixes_2eac[[k]]$Location==1)
      loc[i,j,k,3]=length(I)
    }

loc2=apply(loc,c(1,2,4),mean,na.rm=T)


#cols=gray(seq(1,0.1,-0.15))
#bb=c(-0.5,0,1,2,4,8,12,100)

bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
names=c("NoEAC","Control","2EAC")

pdf(file=paste(figdir,"/ECL_locations_0708_",cat[c],"_SST_change_v2_points.pdf",sep=""),width=8,height=4,pointsize=12)
layout(cbind(1,2,3),c(1,1,0.3))
par(mar=c(3,3,3,2))
for(i in c(1,3))
{
  image(lon,lat,100*t(apply((loc[,,,i]/loc[,,,2])-1,c(1,2),mean,na.rm=T)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=names[i],cex.axis=1.5,cex.main=1.5)
  
  tmp=apply(loc[,,,i]/loc[,,,2]>1,c(1,2),mean) ## Proportion with a positive change
  sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
  points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=2,lwd=2) ## Add some points
  map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
}
ColorBar(bb2,cm)
dev.off()

####### Change in mean/distribution of different parameters
###### Based on events

change=array(NaN,c(6,2,4,3))
dimnames(change)[[1]]=rownames(count)
dimnames(change)[[2]]=c("NoEAC","2EAC")
dimnames(change)[[3]]=colnames(eventsBRAN[[1]])[8:11]
dimnames(change)[[4]]=c("Mean difference","t-test sig","ks-test sig")

for(j in 1:2)
  for(i in 1:6)
  {
    control=eventsBRAN[[i]] 
      if(j==1) events=events_noeac[[i]] else
          if(j==2) events=events_2eac[[i]] else next
          
          
    for(k in 1:4)
    {
      change[i,j,k,1]=mean(events[,k+7])-mean(control[,k+7])
      a=t.test(control[,k+7],events[,k+7])
      change[i,j,k,2]=a$p.value
      a=ks.test(control[,k+7],events[,k+7])
      change[i,j,k,3]=a$p.value
    }
  }

###########################
### Chanegs in general statistics of ECLs

## GV distribution
## Overall changes in GV distribution?

gvthresh=c(0,20,Inf)
count<-array(0,c(6,3,2))

for(i in 1:6)
  for(j in 1:2)
  {
    count[i,1,j]=length(which(eventsBRAN[[i]]$GV>=gvthresh[j] & eventsBRAN[[i]]$GV<gvthresh[j+1]))
    count[i,2,j]=length(which(events_noeac[[i]]$GV>=gvthresh[j] & events_noeac[[i]]$GV<gvthresh[j+1]))
    count[i,3,j]=length(which(events_2eac[[i]]$GV>=gvthresh[j] & events_2eac[[i]]$GV<gvthresh[j+1]))
  }

a=apply(count[,,1:2],c(1,2),sum)
b=apply(count[,,3:6],c(1,2),sum)

apply(count[,2,]/count[,1,],2,mean)

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


## New count matrix

count=array(NaN,c(6,4,12))
dimnames(count)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(count)[[2]]=c("Control","NoEAC","2EAC")
dimnames(count)[[3]]=c("All","CV>=2","Bombs","Formed in region","EC","SC","Mixed","Mean rain>6 mm/6hr","Mean rain>12 mm/6hr",
                       "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")

for(j in 1:3)
  for(i in 1:6)
  {
    if(j==1) events=eventsBRAN[[i]] else
      if(j==2) events=events_noeac[[i]] else
          if(j==3) events=events_2eac[[i]] 
        
    count[i,j,1]=length(events$CV2)
    count[i,j,2]=length(which(events$CV2>=2))
    count[i,j,3]=sum(events$Bomb)
    count[i,j,4]=length(which(events$EnteredFormed==1))
    count[i,j,5]=length(which(events$HartType=="EC"))
    count[i,j,6]=length(which(events$HartType=="SC"))
    count[i,j,7]=length(which(events$HartType=="Mixed"))
    count[i,j,8]=length(which(events$MaxMeanRain250>=6))
    count[i,j,9]=length(which(events$MaxMeanRain250>=12))
    count[i,j,10]=length(which(events$MaxPointRain250>=50))
    count[i,j,11]=length(which(events$MaxMeanWind250>=13.9))
    count[i,j,12]=length(which(events$MaxPointWind250>=22.2))
  }

######## Table 3

table3=matrix(0,5,5)
colnames(table3)=c("Mean rain>6 mm/6hr","Mean rain>12 mm/6hr",
                   "Max rain>50 mm/6hr","Mean wind> 50 km/hr","Max wind> 80 km/hr")
rownames(table3)=c("Average (Control)","Change in NoEAC (%)","Change in 2EAC (%)",
                   "Change in NoEAC (% positive)","Change in 2EAC (% positive)")
table3[1,]=apply(count[,1,8:12],2,mean)
table3[2,]=100*(apply(count[,2,8:12]/count[,1,8:12],2,mean)-1)
table3[3,]=100*(apply(count[,3,8:12]/count[,1,8:12],2,mean)-1)
table3[4,]=100*apply(((count[,2,8:12]/count[,1,8:12])>1),2,mean)
table3[5,]=100*apply(((count[,3,8:12]/count[,1,8:12])>1),2,mean)


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
extreme=array(NaN,c(6,3,10))
dimnames(extreme)[[3]]=colnames(events)[18:27]
for(i in 1:6)
    for(j in 1:4)
    {
      if(j==1) events=eventsBRAN[[i]] else
        if(j==2) events=events_noeac[[i]] else
            if(j==3) events=events_2eac[[i]] else next
          
      events[events==-Inf]=NaN
      extreme[i,j,1:10]=apply(events[,18:27],2,mean,na.rm=T)
    }
  
apply((extreme[,2,]-extreme[,1,]),2,mean,na.rm=T) #% change

testex=array(NaN,c(6,2,10))
dimnames(testex)[[3]]=colnames(eventsBRAN[[1]])[18:27]
for(i in 1:6)
  for(j in 1:10)
  {
    a=t.test(eventsBRAN[[i]][,j+17],events_noeac[[i]][,j+17])
    testex[i,1,j]=a$p.value
    a=t.test(eventsBRAN[[i]][,j+17],events_2eac[[i]][,j+17])
    testex[i,2,j]=a$p.value
  }


extremeF=array(NaN,c(6,3,8))
dimnames(extremeF)[[3]]=colnames(fixesBRAN[[1]])[23:30]
for(i in 1:6)
  for(j in 1:4)
  {
    if(j==1) fixes=fixesBRAN[[i]] else
      if(j==2) fixes=fixes_noeac[[i]] else
          if(j==3) fixes=fixes_2eac[[i]] else next
          
          I=which(fixes$Location==1)
          extremeF[i,j,]=apply(fixes[I,23:30],2,mean,na.rm=T)
  }
apply((extremeF[,1,]),2,mean) #% change
apply((extremeF[,4,]-extremeF[,1,]),2,mean) #% change
apply((extremeF[,2,]-extremeF[,1,]),2,mean) #% change

testexF=array(NaN,c(6,2,8))
dimnames(testexF)[[3]]=colnames(fixesBRAN[[1]])[23:30]
for(i in 1:6)
  for(j in 1:8)
  {
    a=t.test(fixesBRAN[[i]][fixesBRAN[[i]]$Location==1,j+22],fixes_noeac[[i]][fixes_noeac[[i]]$Location==1,j+22])
    testex[i,1,j]=a$p.value
    a=t.test(fixesBRAN[[i]][fixesBRAN[[i]]$Location==1,j+22],fixes_2eac[[i]][fixes_2eac[[i]]$Location==1,j+22])
    testex[i,2,j]=a$p.value
  }


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

