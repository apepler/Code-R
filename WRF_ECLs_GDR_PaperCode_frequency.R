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

#### By season?

count2=array(NaN,c(12,6,2))
  for(j in 1:12)
    for(k in 1:6)
    {
      I=which(events[[k]]$Month==j)
      count2[j,k,1]=length(I)
      I=which(events_notopo[[k]]$Month==j)
      count2[j,k,2]=length(I)
    }


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
cvthresh=1 ## Change this if I want a minimum CV threshold
countthresh=NaN ## Change this if I want to use a CV threshold that gives X events p.a. for the BRAN (control) case

for(i in 1:6)
{
  count2[i,1]=length(which(fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh))
  count2[i,2]=length(which(fixes_notopo[[i]]$Location==1 & fixes_notopo[[i]]$CV>=cvthresh))
  count2[i,3]=length(which(fixes[[i]]$Location2==1 & fixes[[i]]$CV>=cvthresh))
  count2[i,4]=length(which(fixes_notopo[[i]]$Location2==1 & fixes_notopo[[i]]$CV>=cvthresh))
}
median(count2[,4]/count2[,3])

########## Prep for box plot

count=array(NaN,c(6,2,5,4))
dimnames(count)[[1]]=c("R1 50km","R2 50km","R3 50km","R1 10km","R2 10km","R3 10km")
dimnames(count)[[2]]=c("Control","NoTopo")
cvthresh=c(1,1.5,2,2.5,3)
dimnames(count)[[3]]=cvthresh
dimnames(count)[[4]]=c("Events","Fixes","Days","CoastDays")

for(i in 1:6)
  for(j in 1:5)
{
    count[i,1,j,1]=length(which(events[[i]]$CV2>=cvthresh[j]))
    count[i,2,j,1]=length(which(events_notopo[[i]]$CV2>=cvthresh[j]))
    
    count[i,1,j,2]=length(which(fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh[j]))
    count[i,2,j,2]=length(which(fixes_notopo[[i]]$Location==1 & fixes_notopo[[i]]$CV>=cvthresh[j]))
    
    count[i,1,j,3]=length(unique(fixes[[i]]$Date[fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh[j]]))
    count[i,2,j,3]=length(unique(fixes_notopo[[i]]$Date[fixes_notopo[[i]]$Location==1 & fixes_notopo[[i]]$CV>=cvthresh[j]]))
    
    count[i,1,j,4]=length(unique(fixes[[i]]$Date[fixes[[i]]$Location2==1 & fixes[[i]]$CV>=cvthresh[j]]))
    count[i,2,j,4]=length(unique(fixes_notopo[[i]]$Date[fixes_notopo[[i]]$Location2==1 & fixes_notopo[[i]]$CV>=cvthresh[j]]))
    
}

change=100*((count[,2,,]/count[,1,,])-1)
names(dimnames(change))<-c("Source","Intensity","Stat")

library(ggplot2)
library(reshape2)
data=melt(change[,c(1,3),c(1,3:4)])

pdf(paste(figdir,"ECL_change_boxplot_NoTopo_vstrength_cv2.pdf",sep=""),height=5,width=5)
ggplot(data, aes(x = Stat, y = value, fill = Intensity)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  scale_y_continuous(breaks=seq(-75, 75, 25)) +
  theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0)
dev.off()


### What about extreme events?
count=array(NaN,c(6,2,13))
dimnames(count)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(count)[[2]]=c("Control","NoTopo")
dimnames(count)[[3]]=c("All","CV>=2.5","Bombs","Formed in region","Entered","EC","SC","Mixed","Mean rain>6 mm/6hr","Mean rain>12 mm/6hr",
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
      count[i,j,5]=length(which(ev$EnteredFormed!=1))      
      count[i,j,6]=length(which(ev$HartType=="EC"))
      count[i,j,7]=length(which(ev$HartType=="SC"))
      count[i,j,8]=length(which(ev$HartType=="Mixed"))
      count[i,j,9]=length(which(ev$MaxMeanRain250>=6))
      count[i,j,10]=length(which(ev$MaxMeanRain250>=12))
      count[i,j,11]=length(which(ev$MaxPointRain250>=50))
      count[i,j,12]=length(which(ev$MaxMeanWind250>=13.9))
      count[i,j,13]=length(which(ev$MaxPointWind250>=22.2))
  }
count[,2,]/count[,1,]

apply(count[,2,]/count[,1,],2,median)

##Change in rain distribution?
for(k in 1:6) print(t.test(fixes[[k]]$MeanRain500[fixes[[k]]$Location2==1],fixes_notopo[[k]]$MeanRain500[fixes_notopo[[k]]$Location2==1]))

for(k in 1:6)
{
  if(k==1) tmp1=fixes[[k]]$MeanRain500[fixes[[k]]$Location2==1] else tmp1=c(tmp1,fixes[[k]]$MeanRain500[fixes[[k]]$Location2==1])
  if(k==1) tmp2=fixes_notopo[[k]]$MeanRain500[fixes_notopo[[k]]$Location2==1] else tmp2=c(tmp2,fixes_notopo[[k]]$MeanRain500[fixes_notopo[[k]]$Location2==1])
}

t.test(tmp1,tmp2)

dens<-array(NaN,c(512,8,2))
dimnames(dens)[[2]]=c("R1 50km","R2 50km","R3 50km","R1 10km","R2 10km","R3 10km","All 50km","All 10km")
dimnames(dens)[[3]]=c("Control","NoTopo")

for(k in 1:6)
{
  a=density(fixes[[k]]$MeanRain500[fixes[[k]]$Location2==1],from=0,to=20,na.rm=T)
  dens[,k,1]=a$y
  a=density(fixes_notopo[[k]]$MeanRain500[fixes_notopo[[k]]$Location2==1],from=0,to=20,na.rm=T)
  dens[,k,2]=a$y
  
  if(k==1 | k==4) tmp1=fixes[[k]]$MeanRain500[fixes[[k]]$Location2==1] else tmp1=c(tmp1,fixes[[k]]$MeanRain500[fixes[[k]]$Location2==1])
  if(k==1 | k==4) tmp2=fixes_notopo[[k]]$MeanRain500[fixes_notopo[[k]]$Location2==1] else tmp2=c(tmp2,fixes_notopo[[k]]$MeanRain500[fixes_notopo[[k]]$Location2==1])
  
  a=density(tmp1,from=0,to=20,na.rm=T)
  b=density(tmp2,from=0,to=20,na.rm=T)
  
  if(k==3)
  {
    dens[,7,1]=a$y
    dens[,7,2]=b$y
  }
  if(k==6)
  {
    dens[,8,1]=a$y
    dens[,8,2]=b$y
  }
}

ty=c(1,2)
plot(NA,xlim=c(0,20),ylim=c(0,0.2))
for(i in 1:2) lines(a$x,dens[,i+6,1],col="darkgrey",lwd=3,lty=ty[i])
for(i in 1:2) lines(a$x,dens[,i+6,2],col="black",lwd=3,lty=ty[i])
lines(a$x,apply(dens[,4:6,1],1,mean),col="darkgrey",lwd=3)

lims=range(data1,data2,na.rm=T)
if((lims[2]-lims[1])<10)
{
  lims[1]=floor(lims[1])
  lims[2]=ceiling(lims[2])
} else {
  lims[1]=floor(lims[1]/5)*5
  lims[2]=ceiling(lims[2]/5)*5
}

plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab=xlabel,ylab="Frequency",main=tit)
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend(labloc,legend=leg,
       col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   


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

count2[,2,]/count2[,1,]
apply(count2[,2,]/count2[,1,],2,median)

for(i in 1:6) print(makePDF(fixes[[i]]$CV[fixes[[i]]$Location2==1],fixes_notopo[[i]]$CV[fixes_notopo[[i]]$Location2==1]))

count3=array(NaN,c(6,2,8))
dimnames(count3)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(count3)[[2]]=c("Control","NoTopo")
dimnames(count3)[[3]]=c("All","CV>=2.5","Bombs","Mean rain>6 mm/6hr","Mean rain>9 mm/6hr",
                        "Max rain>50 mm/6hr","Mean wind> 40 km/hr","Max wind> 80 km/hr")

for(j in 1:2)
  for(i in 1:6)
  {
    if(j==1) ev=fixes[[i]][fixes[[i]]$Location2>0,] else if(j==2) ev=fixes_notopo[[i]][fixes_notopo[[i]]$Location2>0,] 

    count3[i,j,1]=length(ev$CV)
    count3[i,j,2]=length(which(ev$CV>=2))
    count3[i,j,3]=length(which(ev$NDR>=1))
    count3[i,j,4]=length(which(ev$MeanRain500>=6))
    count3[i,j,5]=length(which(ev$MeanRain500>=9))
    count3[i,j,6]=length(which(ev$MaxRain500>=50))
    count3[i,j,7]=length(which(ev$MeanWind500>=11.1))
    count3[i,j,8]=length(which(ev$MaxWind500>=22.2))
  }


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

loc<-array(NaN,c(11,17,6,2,11))
dimnames(loc)[[5]]=c("Count","Mean CV","Mean MSLP","Mean GV","Mean Radius","Mean CVchange","Mean NDR",
                     "Mean MeanRain","Mean MaxRain","Mean MeanWind","Mean MaxWind")
cols=c(10,9,32,12,15,16,26,27,30,31)

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:6)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1)
      J=which(fixes_notopo[[k]]$Lat>=lat[i]-2.5 & fixes_notopo[[k]]$Lat<lat[i]+2.5 & fixes_notopo[[k]]$Lon>=lon[j]-2.5 & fixes_notopo[[k]]$Lon<lon[j]+2.5 & fixes_notopo[[k]]$Location==1)
      loc[i,j,k,1,1]=length(I)
      loc[i,j,k,2,1]=length(J)
      
      for(x in 1:length(cols))
      {
      loc[i,j,k,1,x+1]=mean(fixes[[k]][I,cols[x]],na.rm=T)
      loc[i,j,k,2,x+1]=mean(fixes_notopo[[k]][J,cols[x]],na.rm=T)
      }
    }


loc2=apply((loc[,,,2,]-loc[,,,1,]),c(1,2,4),mean,na.rm=T)
locPC=100*apply((loc[,,,2,]/loc[,,,1,])-1,c(1,2,4),mean,na.rm=T)

fnames=c("freq","CV","MSLP","GV","Rad","CVdeep","NDR","MeanRain250","MaxRain250","MeanWind250","MaxWind250")
ranges=c(8,0.2,4,4,0.4,0.2,0.4,4,20,2,4)

for(x in 2:11)
{
pdf(file=paste(figdir,"ECL_location_",fnames[x],"_change.pdf",sep=""),width=5,height=4)
      bb2=c(-10000,seq(-ranges[x],ranges[x],length.out=9),10000)
  cm=pal(10)
  layout(cbind(1,2),c(1,0.35))
  par(mar=c(2,2,2,0))
  image(lon,lat,t(loc2[,,x]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main=dimnames(loc)[[5]][x],cex.axis=1,cex.main=1)
  tmp=apply(loc[,,,2,x]>loc[,,,1,x],c(1,2),mean) ## Proportion with a positive change
  sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
  points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=2,lwd=2) ## Add some points
  map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}

x=1
pdf(file=paste(figdir,"ECL_location_",fnames[x],"_PCchange.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(-40,40,10),10000)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
par(mar=c(2,2,2,0))
image(lon,lat,t(locPC[,,x]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main=dimnames(loc)[[5]][x],cex.axis=1,cex.main=1)
tmp=apply(loc[,,,2,x]/loc[,,,1,x]>1,c(1,2),mean) ## Proportion with a positive change
sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=2,lwd=2) ## Add some points
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

### What is the average intensity anyway?
pdf(file=paste(figdir,"ECL_location_CV.pdf",sep=""),width=8.5,height=4)
bb2=c(seq(0.9,1.9,0.1))
pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))
cm=pal1(10)
layout(cbind(1,2,3),c(1,1,0.3))
par(mar=c(3,3,3,1))
image(lon,lat,t(apply(loc[,,,1,2],c(1,2),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="Control",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
image(lon,lat,t(apply(loc[,,,2,2],c(1,2),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

####### Frequency of all lows vs CV2 lows

loc<-array(NaN,c(11,17,6,2,4))
cvthresh=c(1,1.5,2,2.5)

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:6)
      for(x in 1:4)
    {
      I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & 
                fixes[[k]]$Location==1 & fixes[[k]]$CV>=cvthresh[x])
      J=which(fixes_notopo[[k]]$Lat>=lat[i]-2.5 & fixes_notopo[[k]]$Lat<lat[i]+2.5 & fixes_notopo[[k]]$Lon>=lon[j]-2.5 & fixes_notopo[[k]]$Lon<lon[j]+2.5 & 
                fixes_notopo[[k]]$Location==1 & fixes_notopo[[k]]$CV>=cvthresh[x])
      loc[i,j,k,1,x]=length(I)
      loc[i,j,k,2,x]=length(J)
      }

pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))
ran=c(50,25,10,5)

for(x in 1:4)
{
  pdf(file=paste(figdir,"ECL_location_frequency_CV",cvthresh[x],".pdf",sep=""),width=5,height=4)
  bb2=c(-10000,seq(0,ran[x],length.out=11),10000)
  cm=pal1(12)
  layout(cbind(1,2),c(1,0.35))
  par(mar=c(2,2,2,0))
  image(lon,lat,t(apply(loc[,,,1,x],c(1,2),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
        main=dimnames(loc)[[5]][x],cex.axis=1,cex.main=1)
  map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
  ColorBar(bb2,cm)
  dev.off()
}

pdf(file=paste(figdir,"ECL_location_frequency_CV2_comp.pdf",sep=""),width=8.5,height=4)
bb2=c(-10000,seq(0,10,length.out=11),10000)
cm=pal1(12)
layout(cbind(1,2,3),c(1,1,0.3))
par(mar=c(3,3,3,1))
image(lon,lat,t(apply(loc[,,,1,3],c(1,2),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="Control",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
image(lon,lat,t(apply(loc[,,,2,3],c(1,2),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="NoTopo",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()

###### Genesis location, & first-fix-above-2 location

for(k in 1:6)
{
  if(k==1) tmp=fixes[[k]][fixes[[k]]$Fix==1,] else tmp=rbind(tmp,fixes[[k]][fixes[[k]]$Fix==1,])
  I=events[[k]]$ID[events[[k]]$CV2>=2]
  for(i in 1:length(I))
  {
    J=which(fixes[[k]]$ID==I[i] & fixes[[k]]$CV>=2)
    if(i==1 & k==1) tmp2=fixes[[k]][J[1],] else tmp2=rbind(tmp2,fixes[[k]][J[1],])
  }
}

plot(tmp$Lon,tmp$Lat,xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),pch=4,col="blue")
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
points(tmp2$Lon,tmp2$Lat,xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),pch=4,col="blue",lwd=2)

loc<-array(NaN,c(11,17,2))
for(i in 1:11)
  for(j in 1:17)
      {
        I=which(tmp$Lat>=lat[i]-2.5 & tmp$Lat<lat[i]+2.5 & tmp$Lon>=lon[j]-2.5 & tmp$Lon<lon[j]+2.5)
        J=which(tmp2$Lat>=lat[i]-2.5 & tmp2$Lat<lat[i]+2.5 & tmp2$Lon>=lon[j]-2.5 & tmp2$Lon<lon[j]+2.5)
        loc[i,j,1]=length(I)/6
        loc[i,j,2]=length(J)/6
      }

bb2=c(-10000,seq(0,8,length.out=9),10000)
cm=pal1(10)
layout(cbind(1,2),c(1,0.3))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,1]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="Gen",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)
bb2=c(-10000,seq(0,4,length.out=9),10000)
cm=pal1(10)
layout(cbind(1,2),c(1,0.3))
par(mar=c(3,3,3,1))
image(lon,lat,t(loc[,,2]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main="Gen CV2",cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm)

