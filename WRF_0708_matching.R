
rm(list=ls())
library(R.matlab)
library(RNetCDF)
library(abind)
library(sp)
library(geosphere)
library(abind)

####### Part 1 - Which events in ERAI are matched within +- 1 day for the WRF versions, & how do stats differ? 

setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing")

fixes<-events<-list()

for(r in 1:3) {
  fixes[[r]]=read.csv(paste("ECLfixes_d01_0708_R",r,"_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
  fixes[[r]]$Date2=as.POSIXct(paste(as.character(fixes[[r]]$Date),substr(fixes[[r]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  events[[r]]=read.csv(paste("ECLevents_d01_0708_R",r,"_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
}

### Comparing events, using ERAI as the comparator - just if on same day +- 24 hours, count times with overlap

erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)
erai_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)

erai=erai[,-c(1,2)]
erai_E=erai_E[,-c(1,2)]
erai$Date2=as.POSIXct(paste(as.character(erai$Date),substr(erai$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

erai_match=array(NaN,c(length(erai_E$ID),4,3))
dimnames(erai_match)[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2")
dimnames(erai_match)[[3]]=c("R1","R2","R3")

for(i in 1:length(erai_E$ID))
{
  tmp=erai[(erai$ID==erai_E$ID[i] & erai$Location==1),]
  rn=range(tmp$Date2)
  for(r in 1:3)
  {
    I=which(fixes[[r]]$Date2<=rn[2]+(60*60*24) & fixes[[r]]$Date2>=rn[1]-(60*60*24) & fixes[[r]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes[[r]]$ID[I])
      erai_match[i,1,r]=length(J) #All events that match
      erai_match[i,2,r]=length(which(fixes[[r]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
      
      K=which(events[[r]]$ID %in% J)
      erai_match[i,3,r]=max(events[[r]]$CV2[K])
      erai_match[i,4,r]=min(events[[r]]$MSLP2[K])
    } else erai_match[i,1,r]=0
    
  }
}

matchcount=matrix(0,3,3)
for(i in 1:3)
  for(j in 1:3)
    matchcount[i,j]=length(which(erai_match[,1,j]==i-1))

cv=cbind(erai_E$CV2,erai_match[,3,])
clist=c("red","blue","purple")
plot(NA,xlim=c(0,4),ylim=c(0,4),xlab="ERAI",ylab="WRF")
abline(0,1,col="grey")
for(i in 1:3) points(cv[,1],cv[,i+1],col=clist[i],lwd=2,pch=4)
legend("topleft",c("R1","R2","R3"),col=clist,pt.lwd=2,pch=4,bty='n')

slp=cbind(erai_E$MSLP2,erai_match[,4,])
clist=c("red","blue","purple")
plot(NA,xlim=c(980,1021),ylim=c(980,1021),xlab="ERAI",ylab="WRF")
abline(0,1,col="grey")
for(i in 1:3) points(slp[,1],slp[,i+1],col=clist[i],lwd=2,pch=4)
legend("topleft",c("R1","R2","R3"),col=clist,pt.lwd=2,pch=4,bty='n')


###Match 3 WRF versions with each other

setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing")

fixes<-events<-match<-cv<-slp<-list()
wnames=c("_notopo","_BRAN","_BRAN_noeac","_BRAN_2eac")
n=1
for(dom in c("d01","d02"))
  for(r in 1:3) {
    fixes[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    n=n+1
  }

match=array(0,c(6,6))
for(x in 1:6)
  for(y in 1:6)
  {
    matches=rep(NaN,length(events[[x]]$ID))
    for(i in 1:length(events[[x]]$ID))
    {
      tmp=fixes[[x]][(fixes[[x]]$ID==events[[x]]$ID[i] & fixes[[x]]$Location==1),]
      rn=range(tmp$Date2)
      I=which(fixes[[y]]$Date2<=rn[2]+(60*60*24) & fixes[[y]]$Date2>=rn[1]-(60*60*24) & fixes[[y]]$Location==1)
        if(length(I)>0)
        {
          J=unique(fixes[[y]]$ID[I])
          matches[i]=length(J) #All events that match
        } 
        
      }
    match[x,y]=length(which(!is.na(matches)))/length(matches)
  }

##Okay, so ERAI and WRF events match fairly well. 
##Now, we want to take each WRF version (1-3) & match it to ERAI & the four other versions
##That makes it clear how idnvidual events change



setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing")

fixes<-events<-match<-cv<-slp<-list()
wnames=c("_notopo","_BRAN","_BRAN_noeac","_BRAN_2eac")
n=1
for(dom in c("d01","d02"))
  for(r in 1:3) {
    fixes[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    
    match[[n]]=array(NaN,c(length(events[[n]]$ID),4,5))
    dimnames(match[[n]])[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2")
    dimnames(match[[n]])[[3]]=c("ERAI","NoTopo","BRAN","NoEac","2EAC")
    
    fixes2<-events2<-list()
    fixes2[[1]]=erai
    events2[[1]]=erai_E
    
    if(dom=="d01") x2=1:4 else x2=1:3
    
    for(x in x2)
    {
      fixes2[[x+1]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,wnames[x],"_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
      fixes2[[x+1]]$Date2=as.POSIXct(paste(as.character(fixes2[[x+1]]$Date),substr(fixes2[[x+1]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
      events2[[x+1]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,wnames[x],"_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    }
    
    if(dom=="d01") x3=1:5 else x3=1:4
    
    for(i in 1:length(events[[n]]$ID))
    {
      tmp=fixes[[n]][(fixes[[n]]$ID==events[[n]]$ID[i] & fixes[[n]]$Location==1),]
      rn=range(tmp$Date2)
      for(x in x3)
      {
        I=which(fixes2[[x]]$Date2<=rn[2]+(60*60*24) & fixes2[[x]]$Date2>=rn[1]-(60*60*24) & fixes2[[x]]$Location==1)
        if(length(I)>0)
        {
          J=unique(fixes2[[x]]$ID[I])
          match[[n]][i,1,x]=length(J) #All events that match
          match[[n]][i,2,x]=length(which(fixes2[[x]]$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
          
          K=which(events2[[x]]$ID %in% J)
          match[[n]][i,3,x]=max(events2[[x]]$CV2[K])
          match[[n]][i,4,x]=min(events2[[x]]$MSLP2[K])
        } else match[[n]][i,1,x]=0
        
      }
    }
    
    cv[[n]]=cbind(events[[n]]$CV2,match[[n]][,3,])
    slp[[n]]=cbind(events[[n]]$MSLP2,match[[n]][,4,])
    
    n=n+1
  }

comp=array(NaN,c(6,5,4))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("ERAI","NoTopo","BRAN","NoEAC","2EAC")
dimnames(comp)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2")

for(n in 1:6)
  for(x in 1:5)
  {
    comp[n,x,1]=length(which(match[[n]][,1,x]>0))/length(match[[n]][,1,x])
    comp[n,x,2]=mean(match[[n]][,2,x],na.rm=T)
    comp[n,x,3]=mean(cv[[n]][,x+1]-cv[[n]][,1],na.rm=T)
    comp[n,x,4]=mean(slp[[n]][,x+1]-slp[[n]][,1],na.rm=T)
  }

comp[4:6,5,]=NaN

cvcomp=matrix(0,7,5)
dimnames(cvcomp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3","All")
dimnames(cvcomp)[[2]]=c("ERAI","NoTopo","BRAN","NoEAC","2EAC")
for(n in 1:6)
  for(x in 1:5)
  {
    if(x==5)
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

for(x in 1:5)
{
  tmp=cv[[1]][,x+1]-cv[[1]][,1]
  if(x==5) n2=2:3 else n2=2:6
  for(n in n2) tmp=c(tmp,cv[[n]][,x+1]-cv[[n]][,1])
  
  b=t.test(tmp)
  cvcomp[7,x]=b$p.value
}

##### Not significantly different to 0

### Now, let's try to look at the dates of the unmatched events

dates=seq.POSIXt(as.POSIXct("2007010100",format="%Y%m%d%H",tz="GMT"),as.POSIXct("2008123118",format="%Y%m%d%H",tz="GMT"),by="day")
dates2=as.numeric(as.character(dates,format="%Y%m%d"))
unmatchdates=array(0,c(length(dates2),6,5))
dimnames(unmatchdates)[[1]]=dates2
dimnames(unmatchdates)[[2]]=c("d01 R1","d01 R2","d01 R3",
                              "d02 R1","d02 R2","d02 R3")
dimnames(unmatchdates)[[3]]=c("ERAI","NoTopo","BRAN","NoEAC","2EAC")

for(i in 1:6)
  for(j in 1:5)
  {
    I=which(match[[i]][,1,j]==0)
    if(length(I)>1)
    {
      J=which(fixes[[i]]$ID %in% events[[i]]$ID[I] & fixes[[i]]$Location==1)
      K=which(dates2 %in% fixes[[i]]$Date[J])
      unmatchdates[K,i,j]=1
    }
  }

tmp=apply(unmatchdates,1,sum)
I=which(tmp>1)
unmatch2=unmatchdates[I,,]

###################
#######
###### Repeat matching, but using NoEAC/2EAC vs BRAN

setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing")

fixes<-fixes_noeac<-fixes_2eac<-events<-events_noeac<-events_2eac<-match<-cv<-slp<-len<-list()
wnames=c("_notopo","_BRAN","_BRAN_noeac","_BRAN_2eac")
n=1
for(dom in c("d01","d02"))
  for(r in 1:3) {
    fixes[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    
    match[[n]]=array(NaN,c(length(events[[n]]$ID),5,2))
    dimnames(match[[n]])[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2","Length2")
    dimnames(match[[n]])[[3]]=c("NoEac","2EAC")
    
    fixes_noeac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_noeac_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    fixes_noeac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_noeac[[n]]$Date),substr(fixes_noeac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events_noeac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_noeac_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    
    if(dom=="d01")
    {
      fixes_2eac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_2eac_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
      fixes_2eac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_2eac[[n]]$Date),substr(fixes_2eac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
      events_2eac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_2eac_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    }

    
    for(i in 1:length(events[[n]]$ID))
    {
      tmp=fixes[[n]][(fixes[[n]]$ID==events[[n]]$ID[i] & fixes[[n]]$Location==1),]
      rn=range(tmp$Date2)
      
      I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*24) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*24) & fixes_noeac[[n]]$Location==1)
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
        
      if(dom=="d01")
      {
        I=which(fixes_2eac[[n]]$Date2<=rn[2]+(60*60*24) & fixes_2eac[[n]]$Date2>=rn[1]-(60*60*24) & fixes_2eac[[n]]$Location==1)
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
    }
    
    cv[[n]]=cbind(events[[n]]$CV2,match[[n]][,3,])
    slp[[n]]=cbind(events[[n]]$MSLP2,match[[n]][,4,])
    len[[n]]=cbind(events[[n]]$Length2,match[[n]][,5,])
    
    n=n+1
  }

comp=array(NaN,c(6,2,4))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("NoEAC","2EAC")
dimnames(comp)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2")

for(n in 1:6)
  for(x in 1:2)
  {
    comp[n,x,1]=length(which(match[[n]][,1,x]>0))/length(match[[n]][,1,x])
    comp[n,x,2]=mean(match[[n]][,2,x],na.rm=T)
    comp[n,x,3]=mean(cv[[n]][,x+1]-cv[[n]][,1],na.rm=T)
    comp[n,x,4]=mean(slp[[n]][,x+1]-slp[[n]][,1],na.rm=T)
  }

comp[4:6,2,]=NaN

cvcomp=matrix(0,7,2)
dimnames(cvcomp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                        "d02 R1","d02 R2","d02 R3","All")
dimnames(cvcomp)[[2]]=c("NoEAC","2EAC")
for(n in 1:6)
  for(x in 1:2)
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

comp2=array(0,c(6,3,2))
dimnames(comp2)[[2]]=c("Mean difference","CV of matched","CV of unmatched")
for(i in 1:6)
  for(j in 1:2)
  {
    tmp=cv[[i]]
    comp2[i,1,j]=mean(tmp[,j+1]-tmp[,1],na.rm=T)
    I=which(is.na(tmp[,j+1]))
    comp2[i,2,j]=mean(tmp[-I,1])
    comp2[i,3,j]=mean(tmp[I,1])
  }

cv2=rbind(cv[[1]],cv[[2]],cv[[3]])
I=which(is.na(cv2[,2]))
t.test(cv2[I,1],cv2[-I,1])

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

a=hist(events[[1]]$CV2,breaks=cvthresh,plot=F)
tmp=a
tmp$counts=apply(dens2[,1:3,2],1,mean)
plot(tmp,col=rgb(0,0,1,1/4),xlim=c(1,4),ylim=c(0,10),
     xlab="Minimum pressure (hPa)",ylab="Number of events",main="")
tmp$counts=apply(dens2[,1:3,1],1,mean)
plot(tmp,col=rgb(1,0,0,1/4),add=T)
legend("topright",legend=c("Matched","Unmatched"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=8,cex=1,bty="n")   

####### Repeat, w/ slp

comp2=array(0,c(6,3,2))
dimnames(comp2)[[2]]=c("Mean difference","SLP of matched","SLP of unmatched")
for(i in 1:6)
  for(j in 1:2)
  {
    tmp=slp[[i]]
    comp2[i,1,j]=mean(tmp[,j+1]-tmp[,1],na.rm=T)
    I=which(is.na(tmp[,j+1]))
    comp2[i,2,j]=mean(tmp[-I,1])
    comp2[i,3,j]=mean(tmp[I,1])
  }

slp2=rbind(slp[[1]],slp[[2]],slp[[3]])
I=which(is.na(slp2[,2]))
t.test(slp2[I,1],slp2[-I,1])

slpthresh=seq(975,1035,5)
dens2<-array(0,c(12,6,2))
for(r in 1:6)
  for(i in 1:12)
  {
    tmp=slp[[r]]
    I=which(is.na(tmp[,2]))
    J=which(tmp[I,1]>=slpthresh[i] & tmp[I,1]<slpthresh[i+1])
    dens2[i,r,1]=length(J)
    J=which(tmp[-I,1]>=slpthresh[i] & tmp[-I,1]<slpthresh[i+1])
    dens2[i,r,2]=length(J)
    
  }

a=hist(events[[1]]$MSLP2,breaks=slpthresh,plot=F)
tmp=a
tmp$counts=apply(dens2[,,2],1,mean)
plot(tmp,col=rgb(0,0,1,1/4),xlim=c(975,1035),ylim=c(0,10),
     xlab="Minimum pressure (hPa)",ylab="Number of events",main="")
tmp$counts=apply(dens2[,,1],1,mean)
plot(tmp,col=rgb(1,0,0,1/4),add=T)
legend("topright",legend=c("Matched","Unmatched"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=8,cex=1,bty="n")   


##Double matched
##Check for each event if it's matched by an event in R2 & R3. &, if so, if those events are matched by their BRAN

setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing")

fixes<-fixes_noeac<-fixes_noeac2<-events<-events_noeac<-events_noeac2<-match<-cv<-slp<-list()
wnames=c("_notopo","_BRAN","_BRAN_noeac","_BRAN_2eac")
n=1
for(dom in c("d01","d02"))
  for(r in 1:3) {
    fixes[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    events[[n]]$MatchNoEAC==0
    
    match[[n]]=array(NaN,c(length(events[[n]]$ID),4,2))
    dimnames(match[[n]])[[2]]=c("MatchEvents","MatchHours","MSLP2","CV2")
    dimnames(match[[n]])[[3]]=c("NoEac","2EAC")
    
    fixes_noeac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_noeac_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    fixes_noeac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_noeac[[n]]$Date),substr(fixes_noeac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events_noeac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_noeac_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
    
    fixes_noeac2[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_noeac_rad2_p100_cv0.5_typing.csv",sep=""),stringsAsFactors=F)
    fixes_noeac2[[n]]$Date2=as.POSIXct(paste(as.character(fixes_noeac2[[n]]$Date),substr(fixes_noeac2[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events_noeac2[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_noeac_rad2_p100_cv0.5_typing.csv",sep=""),stringsAsFactors=F)
    
    
    ### Find if matched in either normal or 0.5 cv files
    for(i in 1:length(events[[n]]$ID))
    {
      tmp=fixes[[n]][(fixes[[n]]$ID==events[[n]]$ID[i] & fixes[[n]]$Location==1),]
      rn=range(tmp$Date2)
      
      I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*24) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*24) & fixes_noeac[[n]]$Location==1)
      if(length(I)>0) events[[n]]$MatchNoEAC[i]=1 else {
        I=which(fixes_noeac2[[n]]$Date2<=rn[2]+(60*60*24) & fixes_noeac2[[n]]$Date2>=rn[1]-(60*60*24) & fixes_noeac2[[n]]$Location==1)
        if(length(I)>0) events[[n]]$MatchNoEAC[i]=0.5
      }
    }
    
    n=n+1
  }

###Now, add the additional columns

compcol=cbind(c(2,3),c(1,3),c(1,2),c(5,6),c(4,6),c(4,5))
events[[1]]$MatchR3<-events[[1]]$MatchR2<-NaN
events[[2]]$MatchR3<-events[[2]]$MatchR1<-NaN
events[[3]]$MatchR2<-events[[3]]$MatchR1<-NaN
events[[4]]$MatchR3<-events[[4]]$MatchR2<-NaN
events[[5]]$MatchR3<-events[[5]]$MatchR1<-NaN
events[[6]]$MatchR2<-events[[6]]$MatchR1<-NaN

for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    for(i in 1:length(events[[n]]$ID))
    {
      tmp=fixes[[n]][(fixes[[n]]$ID==events[[n]]$ID[i] & fixes[[n]]$Location==1),]
      rn=range(tmp$Date2)
      
      ### Look at the two other BRAN files, & find if there's an ECL there and if that one's matched by its buddy.
      for(j in 1:2)
      {
        I=which(fixes[[compcol[j,n]]]$Date2<=rn[2]+(60*60*24) & fixes[[compcol[j,n]]]$Date2>=rn[1]-(60*60*24) & fixes[[compcol[j,n]]]$Location==1)
        if(length(I)>0) 
          {
          J=unique(fixes[[compcol[j,n]]]$ID[I])
          K=which(events[[compcol[j,n]]]$ID %in% J)
          events[[n]][i,16+j]=max(events[[compcol[j,n]]]$MatchNoEAC[K])
        }
      }
      
      I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*24) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*24) & fixes_noeac[[n]]$Location==1)
      if(length(I)>0) events[[n]]$MatchNoEAC[i]=1 else {
        I=which(fixes_noeac2[[n]]$Date2<=rn[2]+(60*60*24) & fixes_noeac2[[n]]$Date2>=rn[1]-(60*60*24) & fixes_noea2c[[n]]$Location==1)
        events[[n]]$MatchNoEAC[i]=0.5
      }
      
    }
    n=n+1
  }
