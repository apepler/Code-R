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

###

match<-cv<-slp<-len<-gv<-list()
n=1
for(n in 1:6)
{
  match[[n]]=array(NaN,c(length(events[[n]]$ID),10,3))
  dimnames(match[[n]])[[2]]=c("CV","MSLP","Length","GV","MatchEvents","MatchHours","CV2","MSLP2","Length2","GV2")
  dimnames(match[[n]])[[3]]=c("NoEac","2EAC","LR")
  
  match[[n]]=eventmatch(events[[n]],fixes[[n]],events_notopo[[n]],fixes_notopo[[n]],T)

  cv[[n]]=cbind(events[[n]]$CV2,match[[n]][,7])
  slp[[n]]=cbind(events[[n]]$MSLP2,match[[n]][,8])
  len[[n]]=cbind(events[[n]]$Length2,match[[n]][,9])
  gv[[n]]=cbind(events[[n]]$GV,match[[n]][,10])
  n=n+1
}

comp=array(NaN,c(6,5))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("MatchEvents","MatchHours","CV2","MSLP2","GV2")

for(n in 1:6)
  {
    comp[n,1]=length(which(match[[n]][,4]>0))/length(match[[n]][,4])
    comp[n,2]=mean(match[[n]][,5],na.rm=T)
    comp[n,3]=mean(cv[[n]][,2]-cv[[n]][,1],na.rm=T)
    comp[n,4]=mean(slp[[n]][,2]-slp[[n]][,1],na.rm=T)
    comp[n,5]=mean(gv[[n]][,2]-gv[[n]][,1],na.rm=T)
}

cvthresh=c(1,1.5,2,2.5,5)
cvchange=array(NaN,c(6,4,3))
dimnames(cvchange)[[1]]=c("d01 R1","d01 R2","d01 R3",
                          "d02 R1","d02 R2","d02 R3")
dimnames(cvchange)[[2]]=cvthresh[1:4]
dimnames(cvchange)[[3]]=c("Count","NoTopo Match","NoTopo CV change")
for(n in 1:6)
  for(j in 1:4)
  {
    I=which(cv[[n]][,1]>=cvthresh[j] & cv[[n]][,1]<cvthresh[j+1])
    cvchange[n,j,1]=length(I)
    cvchange[n,j,2]=length(which(!is.na(cv[[n]][I,2])))
    cvchange[n,j,3]=mean(cv[[n]][I,2]-cv[[n]][I,1],na.rm=T)
  }

for(n in 1:6) print(cor(cv[[n]][,2]-cv[[n]][,1],cv[[n]][,1],use="pairwise.complete.obs"))

boxplot(cvchange[,,3],xlab="Intensity",ylab="Intensity change",ylim=c(-1,0.5))
abline(h=0,col="red")


### Match within 500km

match<-cv<-slp<-len<-gv<-list()
n=1
for(n in 1:6)
{
  match[[n]]=eventmatch2(events[[n]],fixes[[n]],events_notopo[[n]],fixes_notopo[[n]],F)
  cv[[n]]=cbind(events[[n]]$CV2,match[[n]][,5])
  slp[[n]]=cbind(events[[n]]$MSLP2,match[[n]][,6])
  n=n+1
}

comp=array(NaN,c(6,5))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("MatchEvents","MatchHours","CV2","MSLP2")

for(n in 1:6)
{
  comp[n,1]=length(which(match[[n]][,4]>0))/length(match[[n]][,4])
  comp[n,2]=mean(match[[n]][,5],na.rm=T)
  comp[n,3]=mean(cv[[n]][,2]-cv[[n]][,1],na.rm=T)
  comp[n,4]=mean(slp[[n]][,2]-slp[[n]][,1],na.rm=T)
}

cvthresh=c(1,2,5)
cvchange=array(NaN,c(6,2,3))
dimnames(cvchange)[[2]]=cvthresh[1:2]
dimnames(cvchange)[[3]]=c("Count","NoTopo Match","NoTopo CV change")
for(n in 1:6)
  for(j in 1:2)
  {
    I=which(cv[[n]][,1]>=cvthresh[j] & cv[[n]][,1]<cvthresh[j+1])
    cvchange[n,j,1]=length(I)
    cvchange[n,j,2]=length(which(!is.na(cv[[n]][I,2])))
    cvchange[n,j,3]=mean(cv[[n]][I,2]-cv[[n]][I,1],na.rm=T)
  }

for(n in 1:6) print(cor(match[[n]][,5]-match[[n]][,1],match[[n]][,1],use="pairwise.complete.obs"))

#### FixMatch

match<-list()
n=1
for(n in 1:6) match[[n]]=as.data.frame(fixmatch(fixes[[n]],fixes_notopo[[n]],timediff=1,dist=500,rain=T))

comp=array(NaN,c(6,7))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("MatchHours","MatchEvents","CV2","MSLP2","GV2","MeanR","MaxR")

for(n in 1:6)
{
  match[[n]]=match[[n]][match[[n]]$Location2==1,]
  
  comp[n,1]=length(which(match[[n]][,7]>0))/length(match[[n]][,1])
  comp[n,2]=length(unique(match[[n]]$ID[match[[n]]$MatchHours>0]))/length(unique(match[[n]]$ID))
  comp[n,3]=mean(match[[n]]$CV2-match[[n]]$CV,na.rm=T)
  comp[n,4]=mean(match[[n]]$MSLP2-match[[n]]$MSLP,na.rm=T)
  comp[n,5]=mean(match[[n]]$GV2-match[[n]]$GV,na.rm=T)
  comp[n,6]=mean(match[[n]][,13]-match[[n]][,11],na.rm=T)
  comp[n,7]=mean(match[[n]][,14]-match[[n]][,12],na.rm=T)
}

cvthresh=c(1,2,5)
cvchange=array(NaN,c(6,2,5))
dimnames(cvchange)[[1]]=c("d01 R1","d01 R2","d01 R3",
                          "d02 R1","d02 R2","d02 R3")
dimnames(cvchange)[[2]]=cvthresh[1:2]
dimnames(cvchange)[[3]]=c("Count","NoTopo Match","NoTopo CV change","NoTopo MeanRain change","NoTopo MaxRain Change")
match<-list()
for(n in 1:6) match[[n]]=as.data.frame(fixmatch(fixes[[n]][fixes[[n]]$Location2==1,],fixes_notopo[[n]],timediff=1,dist=500,rain=T))
for(n in 1:6)
  for(j in 1:2)
  {
    I=which(match[[n]]$CV>=cvthresh[j] & match[[n]]$CV<cvthresh[j+1])
    cvchange[n,j,1]=length(I)
    cvchange[n,j,2]=length(which(!is.na(match[[n]]$CV[I])))
    cvchange[n,j,3]=mean(match[[n]]$CV2[I]-match[[n]]$CV[I],na.rm=T)
    cvchange[n,j,4]=mean(match[[n]][I,13]/match[[n]][I,11],na.rm=T)
    cvchange[n,j,5]=mean(match[[n]][I,14]/match[[n]][I,12],na.rm=T)
  }

boxplot(cvchange[,,4],xlab="Intensity",ylab="Intensity change")
abline(h=0,col="red")



match2=match[[1]]
for(i in 2:6) match2=rbind(match2,match[[i]])

cvthresh=c(1,1.5,2,2.5,5)
cvchange=array(NaN,c(4,4))
dimnames(cvchange)[[1]]=cvthresh[1:4]
dimnames(cvchange)[[2]]=c("Count","NoTopo Match","NoTopo CV change","t.test")
  for(j in 1:4)
  {
    I=which(match2[,4]>=cvthresh[j] & match2[,4]<cvthresh[j+1])
    cvchange[j,1]=length(I)
    cvchange[j,2]=length(which(!is.na(match2[I,8])))
    cvchange[j,3]=mean(match2[I,13]/match2[I,11],na.rm=T)
    a=t.test(match2[I,13],match2[I,11])
    cvchange[j,4]=a$p.value
  }

match2$CVd=match2$CV2-match2$CV
match2$MeanRd=match2[,13]-match2[,11]
match2$MaxRd=match2[,14]-match2[,12]

I=which(match2$CV<2)
makePDF(match2$CVd[I],match2$CVd[-I],xlabel="Intensity difference",leg=c("Intensity<2","Intensity>=2"))
makePDF(match2$MaxRd[I],match2$MaxRd[-I],xlabel="meanRain500 difference",leg=c("Intensity<2","Intensity>=2"))

a=density(match2$CVd,na.rm=T,from=-3,to=3)
plot(NA,xlim=c(-10,10),ylim=c(0,0.5),ylab="Frequency",xlab="Rain Difference",main="")
cvthresh=c(1,1.5,2,2.5,5)
for(j in 1:4)
{
  I=which(match2[,4]>=cvthresh[j] & match2[,4]<cvthresh[j+1] & !is.na(match2$MeanRd))
  a=density(match2$MeanRd[I])
  lines(a,lwd=4,lty=j,col=j)
}

legend("topright",c("1 <= Intensity < 1.5","1.5 <= Intensity < 2","2 <= Intensity < 2.5","Intensity >= 2.5"), col=1:4,lty=1:4,lwd=4,cex=1,bty="n")

######## Where are all the cyclones with changes > 0.5?

plot(NA,xlim=c(145,162),ylim=c(-42,-23),xlab="",ylab="")
map(xlim=c(145,162),ylim=c(-42,-23),lwd=2,add=T,col="gray")

for(j in 1:6)
{
  match3=as.data.frame(fixmatch2(fixes[[j]],fixes_notopo[[j]],timediff=1,dist=500))
  I=which(match3$CV2-match3$CV<=-0.5)
  J=unique(fixes[[j]]$ID[I])
  
  for(k in 1:length(J))
  {
    K=which(fixes[[j]]$ID==J[k] & match3$CV2-match3$CV<=-0.5)
    points(fixes[[j]]$Lon[K[1]],fixes[[j]]$Lat[K[1]],pch=4)
  }

}


