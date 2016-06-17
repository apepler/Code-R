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
    eventsBRAN[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_impactsC_GV.csv",sep=""),stringsAsFactors = F)
    eventsBRAN[[n]]$Year=floor(eventsBRAN[[n]]$Date1/10000)
    eventsBRAN[[n]]$Month=floor(eventsBRAN[[n]]$Date1/100)%%100
    eventsBRAN[[n]][eventsBRAN[[n]]==-Inf]=NaN
    
    fixesBRAN[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_impactsC_GV.csv",sep=""),stringsAsFactors = F)
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
    
    eventsLR[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC_GV.csv",sep=""),stringsAsFactors = F)
    eventsLR[[n]]$Year=floor(eventsLR[[n]]$Date1/10000)
    eventsLR[[n]]$Month=floor(eventsLR[[n]]$Date1/100)%%100
    eventsLR[[n]][eventsLR[[n]]==-Inf]=NaN    
    
    fixesLR[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC_GV.csv",sep=""),stringsAsFactors = F)
    fixesLR[[n]]$Year=floor(fixesLR[[n]]$Date/10000)
    fixesLR[[n]]$Month=floor(fixesLR[[n]]$Date/100)%%100
    fixesLR[[n]]$Date2=as.POSIXct(paste(as.character(fixesLR[[n]]$Date),substr(fixesLR[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    events_noeac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_impactsC_GV.csv",sep=""),stringsAsFactors = F)
    events_noeac[[n]]$Year=floor(events_noeac[[n]]$Date1/10000)
    events_noeac[[n]]$Month=floor(events_noeac[[n]]$Date1/100)%%100
    events_noeac[[n]][events_noeac[[n]]==-Inf]=NaN
    
    fixes_noeac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_impactsC_GV.csv",sep=""),stringsAsFactors = F)
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
    
    events_2eac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_impactsC_GV.csv",sep=""),stringsAsFactors = F)
    events_2eac[[n]]$Year=floor(events_2eac[[n]]$Date1/10000)
    events_2eac[[n]]$Month=floor(events_2eac[[n]]$Date1/100)%%100  
    events_2eac[[n]][events_2eac[[n]]==-Inf]=NaN
    
    fixes_2eac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_impactsC_GV.csv",sep=""),stringsAsFactors = F)
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

#############################################
####
### ECL matching

## 1 - Comparison of control across different WRF databases & ERAI

erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)
erai_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)
erai=erai[,-c(1,2)]
erai_E=erai_E[,-c(1,2)]
erai$Date2=as.POSIXct(paste(as.character(erai$Date),substr(erai$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

match=array(0,c(7,6))
for(x in 1:6)
  for(y in 1:6)
  {
    tmp=eventmatch(eventsBRAN[[x]],fixesBRAN[[x]],eventsBRAN[[y]],fixesBRAN[[y]])
    match[x,y]=length(which(tmp[,4]>0))/length(tmp[,4])
  }
for(y in 1:6)
{
  tmp=eventmatch(erai_E,erai,eventsBRAN[[y]],fixesBRAN[[y]])
  match[7,y]=length(which(tmp[,4]>0))/length(tmp[,4])
}
match[match==1]=NaN

mean(match[7,]) # With ERAI
mean(cbind(match[1:3,1:3],match[4:6,4:6]),na.rm=T)
mean(c(match[1,4],match[2,5],match[3,6]))

## 2 - matching Control with NoEAC/2EAC

match<-cv<-slp<-len<-gv<-list()
n=1
for(n in 1:6)
{
  match[[n]]=array(NaN,c(length(eventsBRAN[[n]]$ID),10,3))
  dimnames(match[[n]])[[2]]=c("CV","MSLP","Length","GV","MatchEvents","MatchHours","CV2","MSLP2","Length2","GV2")
  dimnames(match[[n]])[[3]]=c("NoEac","2EAC","LR")
  
  match[[n]][,,1]=eventmatch(eventsBRAN[[n]],fixesBRAN[[n]],events_noeac[[n]],fixes_noeac[[n]],T)
  match[[n]][,,2]=eventmatch(eventsBRAN[[n]],fixesBRAN[[n]],events_2eac[[n]],fixes_2eac[[n]],T)
  match[[n]][,,3]=eventmatch(eventsBRAN[[n]],fixesBRAN[[n]],eventsLR[[n]],fixesLR[[n]],T)
  
  cv[[n]]=cbind(eventsBRAN[[n]]$CV2,match[[n]][,7,])
  slp[[n]]=cbind(eventsBRAN[[n]]$MSLP2,match[[n]][,8,])
  len[[n]]=cbind(eventsBRAN[[n]]$Length2,match[[n]][,9,])
  gv[[n]]=cbind(eventsBRAN[[n]]$GV,match[[n]][,10,])
  
  n=n+1
}

### Compare - what happens to cv, slp etc when matched/unmatched?

comp=array(NaN,c(6,3,4))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("NoEAC","2EAC","LR")
dimnames(comp)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2","GV2")

for(n in 1:6)
  for(x in 1:3)
  {
    comp[n,x,1]=length(which(match[[n]][,4,x]>0))/length(match[[n]][,4,x])
    comp[n,x,2]=mean(match[[n]][,5,x],na.rm=T)
    comp[n,x,3]=mean(cv[[n]][,x+1]-cv[[n]][,1],na.rm=T)
    comp[n,x,4]=mean(slp[[n]][,x+1]-slp[[n]][,1],na.rm=T)
    comp[n,x,5]=mean(gv[[n]][,x+1]-gv[[n]][,1],na.rm=T)
  }
apply(comp,c(2,3),mean)

### Is there a difference in the intensity
## a) Between matched & unmatched events?
## b) Between the event in control & its match

cvcomp=array(0,c(6,3,2))
dimnames(cvcomp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                        "d02 R1","d02 R2","d02 R3")
dimnames(cvcomp)[[2]]=c("NoEAC","2EAC","LR")
dimnames(cvcomp)[[3]]=c("Match v Unmatch","Match v Corresp")
gvcomp=cvcomp

for(n in 1:6)
  for(x in 1:3)
  {
    I=which(match[[n]][,4,x]==0)
    b=t.test(cv[[n]][I,1],cv[[n]][-I,1])
    cvcomp[n,x,1]=b$p.value
    b=t.test(gv[[n]][I,1],gv[[n]][-I,1])
    gvcomp[n,x,1]=b$p.value
    
    b=t.test(cv[[n]][,x+1]-cv[[n]][,1])
    cvcomp[n,x,2]=b$p.value
    b=t.test(gv[[n]][,x+1]-gv[[n]][,1])
    gvcomp[n,x,2]=b$p.value
  }

apply(cvcomp<0.05,c(2,3),mean) # What prop are sig?
apply(gvcomp<0.05,c(2,3),mean) # What prop are sig?

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

####### Figure 4
### Histogram of matched/unmatched events

cvthresh=seq(1,4,0.25)
dens2<-array(0,c(12,6,2))
for(r in 1:6)
  for(i in 1:12)
  {
    tmp=cv[[r]] ## The CV matrix for this day
    I=which(is.na(tmp[,2])) ## The NoEAC matching
    J=which(tmp[I,1]>=cvthresh[i] & tmp[I,1]<cvthresh[i+1])
    dens2[i,r,1]=length(J)
    J=which(tmp[-I,1]>=cvthresh[i] & tmp[-I,1]<cvthresh[i+1])
    dens2[i,r,2]=length(J)
  }

a=apply(dens2[1:4,,],c(2,3),sum)
mean(a[,1]/(a[,1]+a[,2]))

pdf(file=paste(figdir,"/ECL_match_SSTchange_CVdist.pdf",sep=""),width=6,height=3.5,pointsize=10)
par(mar=c(4,4,2,2))
counts <- apply(dens2,c(1,3),mean,na.rm=T)
colnames(counts)=c("Unmatched","Matched")
mp<-barplot(t(counts), col=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)),
            ylab="Number of ECLs",xlab="Intensity (hPa. (deg lat)^2))",legend = colnames(counts))
axis(1,at=c(mp-0.6,mp[12]+0.6),cvthresh)
dev.off()

############ Comparing to the cv0.5 case - is there a missing event there?

match2<-cv2<-slp2<-len2<-list()
n=1
for(dom in c("d01","d02"))
  for(r in 1:3) 
{
  match2[[n]]=array(NaN,c(length(eventsBRAN[[n]]$ID),8))
  dimnames(match2[[n]])[[2]]=c("CV","MSLP","Length","MatchEvents","MatchHours","CV2","MSLP2","Length2")

  fixes=rbind(read.csv(paste("../ECLfixes_",dom,"_2007_R",r,"_BRAN_noeac_",cat[5],".csv",sep="")),
              read.csv(paste("../ECLfixes_",dom,"_2008_R",r,"_BRAN_noeac_",cat[5],".csv",sep="")))
  fixes$ID=floor(fixes$Date/10000)*1000+fixes$ID
  fixes$Date2=as.POSIXct(paste(as.character(fixes$Date),substr(fixes$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  events=rbind(read.csv(paste("../ECLevents_",dom,"_2007_R",r,"_BRAN_noeac_",cat[5],".csv",sep="")),
               read.csv(paste("../ECLevents_",dom,"_2008_R",r,"_BRAN_noeac_",cat[5],".csv",sep="")))  
  events$ID=floor(events$Date1/10000)*1000+events$ID
  
  match2[[n]]=eventmatch(eventsBRAN[[n]],fixesBRAN[[n]],events[[n]],fixes[[n]])
  n=n+1
  }

comp2=array(NaN,c(6,2,4))
dimnames(comp2)[[1]]=c("d01 R1","d01 R2","d01 R3",
                       "d02 R1","d02 R2","d02 R3")
dimnames(comp2)[[2]]=c("All","Unmatched by CV1")
dimnames(comp2)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2")

for(n in 1:6)
{
  comp2[n,1,1]=length(which(match2[[n]][,4]>0))/length(match2[[n]][,4])
  comp2[n,1,2]=mean(match2[[n]][,5],na.rm=T)
  comp2[n,1,3]=mean(match2[[n]][,6]-match2[[n]][,1],na.rm=T)
  comp2[n,1,4]=mean(match2[[n]][,7]-match2[[n]][,2],na.rm=T)
  
  I=which(match[[n]][,4,1]==0)
  comp2[n,2,1]=length(which(match2[[n]][I,4]>0))/length(match2[[n]][I,4])
  comp2[n,2,2]=mean(match2[[n]][I,5],na.rm=T)
  comp2[n,2,3]=mean(match2[[n]][I,6]-match2[[n]][I,1],na.rm=T)
  comp2[n,2,4]=mean(match2[[n]][I,7]-match2[[n]][I,2],na.rm=T)
}

apply(comp2,c(2,3),mean)

############ This is repeated for the 2EAC case treated as control

match3<-list()
n=1
for(n in 1:6)
{
  match3[[n]]=array(NaN,c(length(events_2eac[[n]]$ID),8,2))
  dimnames(match3[[n]])[[2]]=c("CV","MSLP","Length","MatchEvents","MatchHours","CV2","MSLP2","Length2")
  dimnames(match3[[n]])[[3]]=c("Control","Control cv0.5")
  
  match3[[n]][,,1]=eventmatch(events_2eac[[n]],fixes_2eac[[n]],eventsBRAN[[n]],fixesBRAN[[n]])

  fixes=rbind(read.csv(paste("../ECLfixes_",dom,"_2007_R",r,"_BRAN_",cat[5],".csv",sep="")),
              read.csv(paste("../ECLfixes_",dom,"_2008_R",r,"_BRAN_",cat[5],".csv",sep="")))
  fixes$ID=floor(fixes$Date/10000)*1000+fixes$ID
  fixes$Date2=as.POSIXct(paste(as.character(fixes$Date),substr(fixes$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  events=rbind(read.csv(paste("../ECLevents_",dom,"_2007_R",r,"_BRAN_",cat[5],".csv",sep="")),
               read.csv(paste("../ECLevents_",dom,"_2008_R",r,"_BRAN_",cat[5],".csv",sep="")))  
  events$ID=floor(events$Date1/10000)*1000+events$ID
  match3[[n]][,,2]=eventmatch(events_2eac[[n]],fixes_2eac[[n]],events,fixes)
  n=n+1
}

comp3=array(NaN,c(6,3,4))
dimnames(comp3)[[1]]=c("d01 R1","d01 R2","d01 R3",
                       "d02 R1","d02 R2","d02 R3")
dimnames(comp3)[[2]]=c("CV1","CV0.5","Unmatched by CV1")
dimnames(comp3)[[3]]=c("MatchEvents","MatchHours","CV2","MSLP2")

for(n in 1:6)
{
  comp3[n,1,1]=length(which(match3[[n]][,4,1]>0))/length(match3[[n]][,4,1])
  comp3[n,1,2]=mean(match3[[n]][,5,1],na.rm=T)
  comp3[n,1,3]=mean(match3[[n]][,6,1]-match3[[n]][,1],na.rm=T)
  comp3[n,1,4]=mean(match3[[n]][,7,1]-match3[[n]][,2],na.rm=T)
  
  comp3[n,2,1]=length(which(match3[[n]][,4,2]>0))/length(match3[[n]][,4,2])
  comp3[n,2,2]=mean(match3[[n]][,5,2],na.rm=T)
  comp3[n,2,3]=mean(match3[[n]][,6,2]-match3[[n]][,1],na.rm=T)
  comp3[n,2,4]=mean(match3[[n]][,7,2]-match3[[n]][,2],na.rm=T)
  
  I=which(match3[[n]][,4,1]==0)
  comp3[n,3,1]=length(which(match3[[n]][I,4,2]>0))/length(match3[[n]][I,4,2])
  comp3[n,3,2]=mean(match3[[n]][I,5,2],na.rm=T)
  comp3[n,3,3]=mean(match3[[n]][I,6,2]-match3[[n]][I,1,2],na.rm=T)
  comp3[n,3,4]=mean(match3[[n]][I,7,2]-match3[[n]][I,2,2],na.rm=T)
}

apply(comp3,c(2,3),mean)

#############
## What about average cv change as a function of cv?

cvthresh=c(1,1.5,5)

cvchange=array(NaN,c(6,2,5))
dimnames(cvchange)[[2]]=cvthresh[1:2]
dimnames(cvchange)[[3]]=c("Count","NoEAC Match","NoEAC CV change","2EAC Match","2EAC CV change")
for(n in 1:6)
  for(j in 1:2)
  {
    I=which(cv[[n]][,1]>=cvthresh[j] & cv[[n]][,1]<cvthresh[j+1])
    cvchange[n,j,1]=length(I)
    
    cvchange[n,j,2]=length(which(!is.na(cv[[n]][I,2])))
    cvchange[n,j,3]=mean(cv[[n]][I,2]-cv[[n]][I,1],na.rm=T)
      cvchange[n,j,4]=length(which(!is.na(cv[[n]][I,3])))
      cvchange[n,j,5]=mean(cv[[n]][I,3]-cv[[n]][I,1],na.rm=T)
  }


