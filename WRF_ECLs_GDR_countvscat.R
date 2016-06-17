###### New play - "Real" ECLs, w/ at least 2 fixes in location 2
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

count<-array(0,c(5,6,2,6))
dimnames(count)[[1]]<-cat
dimnames(count)[[2]]<-c("R1 50km","R2 50km","R2 50km","R1 10km","R2 10km","R3 10km")
dimnames(count)[[3]]<-c("Control","NoTopo")
dimnames(count)[[4]]<-c("All","Loc2","All 22","Loc2 22","All 5","Loc2 5")
countthresh=22*2
countthresh2=5*2

for(c in 1:5)
{
events<-events_notopo<-fixes<-fixes_notopo<-list()

n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    fixes[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
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
    
    events[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    events[[n]][events[[n]]==-Inf]=NaN    
    events[[n]]$Location2=0
    for(i in 1:length(events[[n]]$ID))
    {
      I=which(fixes[[n]]$ID==events[[n]]$ID[i])
      events[[n]]$Location2[i]=sum(fixes[[n]]$Location2[I])
    }
    
    fixes_notopo[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_notopo_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
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
    
    events_notopo[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_notopo_",cat[c],"_typing_GV.csv",sep=""),stringsAsFactors = F)
    events_notopo[[n]]$Year=floor(events_notopo[[n]]$Date1/10000)
    events_notopo[[n]]$Month=floor(events_notopo[[n]]$Date1/100)%%100
    events_notopo[[n]][events_notopo[[n]]==-Inf]=NaN
    
    events_notopo[[n]]$Location2=0
    for(i in 1:length(events_notopo[[n]]$ID))
    {
      I=which(fixes_notopo[[n]]$ID==events_notopo[[n]]$ID[i])
      events_notopo[[n]]$Location2[i]=sum(fixes_notopo[[n]]$Location2[I])
    }
    
    
    a=order(events[[n]]$CV2,decreasing=T)
    if(length(a)>=countthresh) b=events[[n]]$CV2[a[countthresh]] else b=min(events[[n]]$CV2)
    if(length(a)>=countthresh2) d=events[[n]]$CV2[a[countthresh2]] else d=min(events[[n]]$CV2)
    
    count[c,n,1,1]=length(events[[n]]$Location2)
    count[c,n,1,2]=length(which(events[[n]]$Location2>1))
    count[c,n,1,3]=length(which(events[[n]]$CV2>=b))
    count[c,n,1,4]=length(which(events[[n]]$Location2>1 & events[[n]]$CV2>=b))
    count[c,n,1,5]=length(which(events[[n]]$CV2>=d))
    count[c,n,1,6]=length(which(events[[n]]$Location2>1 & events[[n]]$CV2>=d))
    
    count[c,n,2,1]=length(events_notopo[[n]]$Location2)
    count[c,n,2,2]=length(which(events_notopo[[n]]$Location2>1))
    count[c,n,2,3]=length(which(events_notopo[[n]]$CV2>=b))
    count[c,n,2,4]=length(which(events_notopo[[n]]$Location2>1 & events_notopo[[n]]$CV2>=b))
    count[c,n,2,5]=length(which(events_notopo[[n]]$CV2>=d))
    count[c,n,2,6]=length(which(events_notopo[[n]]$Location2>1 & events_notopo[[n]]$CV2>=d))
    
    n=n+1     
  }
}

tmp=count[,,2,]/count[,,1,]
apply(tmp,3,mean)
apply(tmp>1,3,mean)
