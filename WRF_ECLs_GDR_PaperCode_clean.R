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
c=3 ## Several different intensity measures were calculated. This is the default.

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

########## Chanegs in ECL frequency

count=array(NaN,c(6,2,5,5))
dimnames(count)[[1]]=c("R1 50km","R2 50km","R3 50km","R1 10km","R2 10km","R3 10km")
dimnames(count)[[2]]=c("Control","NoTopo")
cvthresh=c(1,1.5,2,2.5,3)
dimnames(count)[[3]]=cvthresh
dimnames(count)[[4]]=c("Events","Fixes","CoastFixes","Days","CoastDays")

for(i in 1:6)
  for(j in 1:5)
{
    count[i,1,j,1]=length(which(events[[i]]$CV2>=cvthresh[j])) 
    count[i,2,j,1]=length(which(events_notopo[[i]]$CV2>=cvthresh[j]))
    
    count[i,1,j,2]=length(which(fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh[j])) 
    count[i,2,j,2]=length(which(fixes_notopo[[i]]$Location==1 & fixes_notopo[[i]]$CV>=cvthresh[j]))
    
    count[i,1,j,3]=length(which(fixes[[i]]$Location2==1 & fixes[[i]]$CV>=cvthresh[j])) 
    count[i,2,j,3]=length(which(fixes_notopo[[i]]$Location2==1 & fixes_notopo[[i]]$CV>=cvthresh[j]))
    
    count[i,1,j,4]=length(unique(fixes[[i]]$Date[fixes[[i]]$Location==1 & fixes[[i]]$CV>=cvthresh[j]]))
    count[i,2,j,4]=length(unique(fixes_notopo[[i]]$Date[fixes_notopo[[i]]$Location==1 & fixes_notopo[[i]]$CV>=cvthresh[j]]))
    
    count[i,1,j,5]=length(unique(fixes[[i]]$Date[fixes[[i]]$Location2==1 & fixes[[i]]$CV>=cvthresh[j]]))
    count[i,2,j,5]=length(unique(fixes_notopo[[i]]$Date[fixes_notopo[[i]]$Location2==1 & fixes_notopo[[i]]$CV>=cvthresh[j]]))
    
}

change=100*((count[,2,,]/count[,1,,])-1)

### Events by month?

count=array(0,c(12,6,2))
for(i in 1:6)
  for(m in 1:12)
  {
    count[m,i,1]=length(which(events[[i]]$Month==m)) 
    count[m,i,2]=length(which(events_notopo[[i]]$Month==m))
  }

a=apply(count[9:11,,],c(2,3),sum)/apply(count,c(2,3),sum)

### Change in frequency of extreme events?
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

a=count[,,4]/count[,,1]

#### Change in ECL rainfall 

cvthresh1=1
cvthresh2=2
count3=array(NaN,c(6,2,9))
dimnames(count3)[[1]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")
dimnames(count3)[[2]]=c("Control","NoTopo")
dimnames(count3)[[3]]=c("Count","Bombs","Mean Rain","Max Rain","Mean rain>6 mm/6hr","Mean rain>9 mm/6hr",
                        "Max rain>50 mm/6hr","Mean wind> 40 km/hr","Max wind> 80 km/hr")

for(j in 1:2)
  for(i in 1:6)
  {
    if(j==1) ev=fixes[[i]][fixes[[i]]$Location>0 & fixes[[i]]$CV>=cvthresh1 & fixes[[i]]$CV<cvthresh2,] else if(j==2) ev=fixes_notopo[[i]][fixes_notopo[[i]]$Location>0 & fixes_notopo[[i]]$CV>=cvthresh1 & fixes_notopo[[i]]$CV<cvthresh2,] 
    
    count3[i,j,1]=length(ev$CV)
    count3[i,j,2]=length(which(ev$NDR>=1))
    count3[i,j,3]=mean(ev$MeanRain500,na.rm=T)
    count3[i,j,4]=mean(ev$MaxRain500,na.rm=T)
    count3[i,j,5]=length(which(ev$MeanRain500>=6))
    count3[i,j,6]=length(which(ev$MeanRain500>=9))
    count3[i,j,7]=length(which(ev$MaxRain500>=50))
    count3[i,j,8]=length(which(ev$MeanWind500>=11.1))
    count3[i,j,9]=length(which(ev$MaxWind500>=22.2))
  }

apply(count3[,1,],2,mean)
apply(count3[,2,]/count3[,1,],2,mean)
count3[,2,]/count3[,1,]
######## Change in ECL locations

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

#### Total numbers in location - to compare to NoNudge

pal1 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

pdf(file=paste(figdir,"ECL_location_d01_nudge_0708_allwrf.pdf",sep=""),width=5,height=4)
bb2=c(-10000,seq(0,20,length.out=11),10000)
cm=pal1(12)
layout(cbind(1,2),width=c(1,0.35))
par(mar=c(3,3,3,1))
image(lon,lat,t(apply(loc[,,,1,1],c(1,2),mean))/2,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(145,161),ylim=c(-42,-23),
      main="Nudged - 2007-2008",cex.axis=1.5,cex.main=1.5)
map(xlim=c(145,161),ylim=c(-42,-23),add=T,lwd=2)
ColorBar(bb2,cm)
dev.off()


### Figure 3

pdf(file=paste(figdir,"ECL_location_freq_CV_PCchange.pdf",sep=""),width=10,height=4)
bb1=c(-10000,seq(-40,40,10),10000)
bb2=c(-10000,seq(-0.2,0.2,0.05),10000)
cm=pal(10)
layout(cbind(1,3,2,4),c(1,0.35,1,0.35))
par(mar=c(3,3,3,0))
image(lon,lat,t(locPC[,,1]),xlab="",ylab="",breaks=bb1,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main=dimnames(loc)[[5]][1],cex.axis=1.5,cex.main=1.5)
tmp=apply(loc[,,,2,1]/loc[,,,1,1]>1,c(1,2),mean) ## Proportion with a positive change
sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=3,lwd=2) ## Add some points
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)

image(lon,lat,t(loc2[,,2]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main=dimnames(loc)[[5]][2],cex.axis=1.5,cex.main=1.5)
tmp=apply(loc[,,,2,2]/loc[,,,1,2]>1,c(1,2),mean) ## Proportion with a positive change
sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=3,lwd=2) ## Add some points
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)

ColorBar(bb1,cm)
ColorBar(bb2,cm,subsampleg = 2)
dev.off()

####### Matched events

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

### Figure 4 - v1 - Difference in intensity vs intensity for events

boxplot(cvchange[,,3],xlab="Intensity",ylab="Intensity change",ylim=c(-1,0.5))
abline(h=0,col="red")

#### Matching instances

match<-list()
n=1
for(n in 1:6) match[[n]]=as.data.frame(fixmatch3(fixes[[n]],fixes_notopo[[n]],timediff=1,dist=500,rain=T))

comp=array(NaN,c(6,7))
dimnames(comp)[[1]]=c("d01 R1","d01 R2","d01 R3",
                      "d02 R1","d02 R2","d02 R3")
dimnames(comp)[[2]]=c("Matched?","Lon diff","CV diff","MSLP diff","GV diff","MeanR diff","MaxR diff")

for(n in 1:6)
{
  comp[n,1]=length(which(match[[n]][,7]>0))/length(match[[n]][,1])
  comp[n,2]=mean(match[[n]]$Lon2-match[[n]]$Lon,na.rm=T)
  comp[n,3]=mean(match[[n]]$CV2-match[[n]]$CV,na.rm=T)
  comp[n,4]=mean(match[[n]]$MSLP2-match[[n]]$MSLP,na.rm=T)
  comp[n,5]=mean(match[[n]]$GV2-match[[n]]$GV,na.rm=T)
  comp[n,6]=mean(match[[n]][,17]-match[[n]][,15],na.rm=T)
  comp[n,7]=mean(match[[n]][,18]-match[[n]][,16],na.rm=T)
}

cvthresh=c(1,2,5)
cvchange=array(NaN,c(6,2,6))
dimnames(cvchange)[[1]]=c("d01 R1","d01 R2","d01 R3",
                          "d02 R1","d02 R2","d02 R3")
dimnames(cvchange)[[2]]=cvthresh[1:2]
dimnames(cvchange)[[3]]=c("Count","NoTopo Match","NoTopo Lon Change","NoTopo CV change","NoTopo MeanRain change","NoTopo MaxRain Change")
for(n in 1:6)
  for(j in 1:2)
  {
    I=which(match[[n]]$CV>=cvthresh[j] & match[[n]]$CV<cvthresh[j+1])
    cvchange[n,j,1]=length(I)
    cvchange[n,j,2]=length(which(!is.na(match[[n]]$CV[I])))
    cvchange[n,j,3]=mean(match[[n]]$Lon2[I]-match[[n]]$Lon[I],na.rm=T)
    cvchange[n,j,4]=mean(match[[n]]$CV2[I]-match[[n]]$CV[I],na.rm=T)
    cvchange[n,j,5]=mean(match[[n]][I,17]/match[[n]][I,15],na.rm=T)
    cvchange[n,j,6]=mean(match[[n]][I,18]/match[[n]][I,16],na.rm=T)
  }

## Figure 4 - v2 - Difference in intensity vs intensity for fixes

boxplot(cvchange[,,4],xlab="Intensity",ylab="Intensity change")
abline(h=0,col="red")

## Figure 4 - v3 - all lumped together

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
match2$MeanRd=match2[,17]/match2[,15]
match2$MaxRd=match2[,18]/match2[,16]
match2$Lond=match2$Lon2-match2$Lon
match2$Latd=match2$Lat2-match2$Lat
apply(match2[,19:23],2,mean,na.rm=T)

I=which(match2$CVd<(-0.5))
apply(match2[I,19:23],2,mean,na.rm=T)
apply(match2[-I,19:23],2,mean,na.rm=T)

I=which(match2$CV<2)
makePDF(match2$CVd[I],match2$CVd[-I],xlabel="Intensity difference",leg=c("Intensity<2","Intensity>=2"))

I=which(match2$Lat<=-33 & match2$Lon<=157)
mean(match2$CVd[I],na.rm=T)
lm(match2$CVd[I]~match2$CV[I])
plot(match2$CV,match2$CVd)


######## Change vs location for matched cyclones
lat=seq(-60,-10,5)
lon=seq(100,180,5)
loc<-array(NaN,c(11,17,6,7))
dimnames(loc)[[4]]=c("Count","NoTopo Match","NoTopo Lon Change","NoTopo Lat Change","NoTopo CV change","NoTopo MeanRain change","NoTopo MaxRain Change")

for(i in 1:11)
  for(j in 1:17)
    for(k in 1:6)
    {
      I=which(match[[k]]$Lat>=lat[i]-2.5 & match[[k]]$Lat<lat[i]+2.5 & match[[k]]$Lon>=lon[j]-2.5 & match[[k]]$Lon<lon[j]+2.5)
      
      loc[i,j,k,1]=length(I)
      loc[i,j,k,2]=100*(sum(match[[k]]$MatchHours[I])/length(I))
      loc[i,j,k,3]=mean(match[[k]]$Lon2[I]-match[[k]]$Lon[I],na.rm=T)
      loc[i,j,k,4]=mean(match[[k]]$Lat2[I]-match[[k]]$Lat[I],na.rm=T)
      loc[i,j,k,5]=mean(match[[k]]$CV2[I]-match[[k]]$CV[I],na.rm=T)
      loc[i,j,k,6]=100*((mean(match[[k]][I,17]/match[[k]][I,15],na.rm=T))-1)
      loc[i,j,k,7]=100*((mean(match[[k]][I,18]/match[[k]][I,16],na.rm=T))-1)
    }


## Prop matched - higher in S than in N
bb1=seq(0,100,10)
cm=pal(10)
layout(cbind(1,2),c(1,0.35))
image(lon,lat,t(apply(loc[,,,2],c(1,2),mean)),xlab="",ylab="",breaks=bb1,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main=dimnames(loc)[[4]][2],cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb1,cm)

pdf(file=paste(figdir,"ECL_location_match_Lonchange_CVchange.pdf",sep=""),width=10,height=4)
bb1=c(-10000,seq(-0.8,0.8,0.2),10000)
bb2=c(-10000,seq(-0.4,0.4,0.1),10000)
cm=pal(10)
layout(cbind(1,3,2,4),c(1,0.35,1,0.35))
par(mar=c(3,3,3,0))
image(lon,lat,t(apply(loc[,,,3],c(1,2),mean)),xlab="",ylab="",breaks=bb1,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main=dimnames(loc)[[4]][3],cex.axis=1.5,cex.main=1.5)
tmp=apply(loc[,,,3]>0,c(1,2),mean) ## Proportion with a positive change
sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=3,lwd=2) ## Add some points
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)

image(lon,lat,t(apply(loc[,,,5],c(1,2),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),
      main=dimnames(loc)[[4]][5],cex.axis=1.5,cex.main=1.5)
tmp=apply(loc[,,,5]>0,c(1,2),mean) ## Proportion with a positive change
sigmask=which(tmp>0.75 | tmp<0.25,arr.ind=T) ## Which have at least 3/4 (~ 5/6) in same direction 
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=3,lwd=2) ## Add some points
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)

ColorBar(bb1,cm,subsampleg = 2)
ColorBar(bb2,cm,subsampleg = 2)
dev.off()

##########



##### What about matching against weaker events?

fixes_notopoW<-list()

c=5
n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    fixes_notopoW[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_notopo_",cat[c],"_typing.csv",sep=""),stringsAsFactors = F)
    fixes_notopoW[[n]]$Year=floor(fixes_notopoW[[n]]$Date/10000)
    fixes_notopoW[[n]]$Month=floor(fixes_notopoW[[n]]$Date/100)%%100
    fixes_notopoW[[n]]$Date2=as.POSIXct(paste(as.character(fixes_notopoW[[n]]$Date),substr(fixes_notopoW[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    fixes_notopoW[[n]][fixes_notopoW[[n]]==-Inf]=NaN    
    fixes_notopoW[[n]][fixes_notopoW[[n]]==Inf]=NaN  
    
    n=n+1     
  }

n=1
match3=as.data.frame(fixmatch3(fixes[[n]],fixes_notopoW[[n]],timediff=1,dist=500,rain=F))
for(n in 2:6) match3=rbind(match3,fixmatch3(fixes[[n]],fixes_notopoW[[n]],timediff=1,dist=500,rain=F))



match3$CVd=match3$CV2-match3$CV

I=which(match3$Lat<=-33 & match3$Lon<=157)
mean(match3$CVd[I],na.rm=T)
lm(match3$CVd[I]~match3$CV[I])
plot(match3$CV,match3$CVd)


match<-list()
for(n in 1:6) match[[n]]=as.data.frame(fixmatch3(fixes[[n]],fixes_notopoW[[n]],timediff=1,dist=500,rain=F))

cvthresh=c(1,2,5)
cvchange=array(NaN,c(6,2,4))
dimnames(cvchange)[[1]]=c("d01 R1","d01 R2","d01 R3",
                          "d02 R1","d02 R2","d02 R3")
dimnames(cvchange)[[2]]=cvthresh[1:2]
dimnames(cvchange)[[3]]=c("Count","NoTopo Match","NoTopo Lon Change","NoTopo CV change")
for(n in 1:6)
  for(j in 1:2)
  {
    I=which(match[[n]]$CV>=cvthresh[j] & match[[n]]$CV<cvthresh[j+1])
    cvchange[n,j,1]=length(I)
    cvchange[n,j,2]=length(which(!is.na(match[[n]]$CV[I])))
    cvchange[n,j,3]=mean(match[[n]]$Lon2[I]-match[[n]]$Lon[I],na.rm=T)
    cvchange[n,j,4]=mean(match[[n]]$CV2[I]-match[[n]]$CV[I],na.rm=T)
  }

boxplot(cvchange[,,4],xlab="Intensity",ylab="Intensity change")
abline(h=0,col="red")


## Looking at events that fail to intensify - what's happening at the moment?

for(k in 1:6)
{
  match[[k]]$CVd=match[[k]]$CV2-match[[k]]$CV
  I=which(match[[k]]$CV2-match[[k]]$CV<=-0.5)
  J=unique(match[[k]]$ID[I])
}

match[[k]][match[[k]]$ID%in%J,c(1:14,19)]

k=6



###Change in various types?

tlistSB=c("ET","IT","CL","SSL")
tlistH=c("EC","SC","TC","Mixed")

countT=array(0,c(6,9,2))
dimnames(countT)[[2]]=c("All",tlistH,tlistSB)

cvthresh=1
    for(k in 1:6)
    {
      I=which(events[[k]]$CV2>=cvthresh)
      countT[k,1,1]=length(I)
      J=which(events_notopo[[k]]$CV2>=cvthresh)
      countT[k,1,2]=length(J)
      
      for(t in 1:4)
      {
        countT[k,t+1,1]=length(which(events[[k]]$HartType[I]==tlistH[t]))
        countT[k,t+5,1]=length(which(events[[k]]$TypeSB[I]==tlistSB[t]))
        countT[k,t+1,2]=length(which(events_notopo[[k]]$HartType[J]==tlistH[t]))
        countT[k,t+5,2]=length(which(events_notopo[[k]]$TypeSB[J]==tlistSB[t]))
      }
    }

