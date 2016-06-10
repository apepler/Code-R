rm(list=ls())
setwd('/home/nfs/z3478332/output/outputUM_wrf_2007/')

makePDF = function(data1,data2,xlabel,labloc="topleft") {
  a=density(data1,na.rm=T)
  b=density(data2,na.rm=T)
  
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
       xlab=xlabel,ylab="Frequency",cex.main=1.2,
       main="2007-2008 ECL statistics")
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=c("Control","NoTopo"),
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
}


#events=list(T60=rbind(read.csv("ECLevents_2007_p60.csv"),read.csv("ECLevents_2008_p60.csv")),
#            NT60=rbind(read.csv("ECLevents_2007_notopo_p60.csv"),read.csv("ECLevents_2008_notopo_p60.csv")),
#            T240=rbind(read.csv("ECLevents_2007_p240.csv"),read.csv("ECLevents_2008_p240.csv")),
#            NT240=rbind(read.csv("ECLevents_2007_notopo_p240.csv"),read.csv("ECLevents_2008_notopo_p240.csv")))
# mm=list(T60=floor(events[[1]]$Date1/100)%%100,NT60=floor(events[[2]]$Date1/100)%%100,
#         T240=floor(events[[3]]$Date1/100)%%100,NT240=floor(events[[4]]$Date1/100)%%100)


events=list(T60=rbind(read.csv("ECLevents_2007_cv1_p60.csv"),read.csv("ECLevents_2008_cv1_p60.csv")),
            NT60=rbind(read.csv("ECLevents_2007_notopo_cv1_p60.csv"),read.csv("ECLevents_2008_notopo_cv1_p60.csv")),
            T240=rbind(read.csv("ECLevents_2007_cv1_p240.csv"),read.csv("ECLevents_2008_cv1_p240.csv")),
            NT240=rbind(read.csv("ECLevents_2007_notopo_cv1_p240.csv"),read.csv("ECLevents_2008_notopo_cv1_p240.csv")),
            T60_d02=rbind(read.csv("ECLevents_d02_2007_cv1_p60.csv"),read.csv("ECLevents_d02_2008_cv1_p60.csv")),
            NT60_d02=rbind(read.csv("ECLevents_d02_2007_notopo_cv1_p60.csv"),read.csv("ECLevents_d02_2008_notopo_cv1_p60.csv")),
            T240_d02=rbind(read.csv("ECLevents_d02_2007_cv1_p240.csv"),read.csv("ECLevents_d02_2008_cv1_p240.csv")),
            NT240_d02=rbind(read.csv("ECLevents_d02_2007_notopo_cv1_p240.csv"),read.csv("ECLevents_d02_2008_notopo_cv1_p240.csv")))

fixes=list(T60=rbind(read.csv("ECLfixes_2007_cv1_p60.csv"),read.csv("ECLfixes_2008_cv1_p60.csv")),
           NT60=rbind(read.csv("ECLfixes_2007_notopo_cv1_p60.csv"),read.csv("ECLfixes_2008_notopo_cv1_p60.csv")),
           T240=rbind(read.csv("ECLfixes_2007_cv1_p240.csv"),read.csv("ECLfixes_2008_cv1_p240.csv")),
           NT240=rbind(read.csv("ECLfixes_2007_notopo_cv1_p240.csv"),read.csv("ECLfixes_2008_notopo_cv1_p240.csv")),
           T60_d02=rbind(read.csv("ECLfixes_d02_2007_cv1_p60.csv"),read.csv("ECLfixes_d02_2008_cv1_p60.csv")),
           NT60_d02=rbind(read.csv("ECLfixes_d02_2007_notopo_cv1_p60.csv"),read.csv("ECLfixes_d02_2008_notopo_cv1_p60.csv")),
           T240_d02=rbind(read.csv("ECLfixes_d02_2007_cv1_p240.csv"),read.csv("ECLfixes_d02_2008_cv1_p240.csv")),
           NT240_d02=rbind(read.csv("ECLfixes_d02_2007_notopo_cv1_p240.csv"),read.csv("ECLfixes_d02_2008_notopo_cv1_p240.csv")))

for(i in 1:8)
{
  I=which(events[[i]]$Date1>20080000)
  a=events[[i]]$ID[I[1]-1]
  events[[i]]$ID[I]=events[[i]]$ID[I]+a
  I=which(fixes[[i]]$Date>20080000)
  fixes[[i]]$ID[I]=fixes[[i]]$ID[I]+a
  
  data=fixes[[i]]
  fixes[[i]]$Date2=as.Date(as.character(data$Date),"%Y%m%d")+(as.numeric(data$Time)-1)/4
}

# mm=list(T60=floor(fixes[[1]]$Date/100)%%100,NT60=floor(fixes[[2]]$Date/100)%%100,
#         T240=floor(fixes[[3]]$Date/100)%%100,NT240=floor(fixes[[4]]$Date/100)%%100,
#         T60_d02=floor(fixes[[5]]$Date/100)%%100,NT60_d02=floor(fixes[[6]]$Date/100)%%100,
#         T240_d02=floor(fixes[[7]]$Date/100)%%100,NT240_d02=floor(fixes[[8]]$Date/100)%%100)

mm=list(T60=floor(events[[1]]$Date1/100)%%100,NT60=floor(events[[2]]$Date1/100)%%100,
        T240=floor(events[[3]]$Date1/100)%%100,NT240=floor(events[[4]]$Date1/100)%%100,
        T60_d02=floor(events[[5]]$Date1/100)%%100,NT60_d02=floor(events[[6]]$Date1/100)%%100,
        T240_d02=floor(events[[7]]$Date1/100)%%100,NT240_d02=floor(events[[8]]$Date1/100)%%100)



cvt=c(seq(0.5,3,0.5),Inf)
mt=c(seq(970,1020,10),Inf)
lt=c(seq(0,20,4),Inf)
rt=c(seq(0,6,0.5),Inf)

CVcount<-MSLPcount<-Lcount<-matrix(0,6,8)
Rcount=matrix(0,14,8)
Mcount=matrix(0,12,8)
colnames(CVcount)<-colnames(MSLPcount)<-colnames(Lcount)<-colnames(Rcount)<-colnames(Mcount)<-names(events)


for(j in 1:8)
{
  for(i in 1:6)
  {
    I=which(events[[j]]$CV2>=cvt[i] & events[[j]]$CV2<cvt[i+1])
    CVcount[i,j]=length(I)
    I=which(events[[j]]$MSLP2>=mt[i] & events[[j]]$MSLP2<mt[i+1])
    MSLPcount[i,j]=length(I)
    I=which(events[[j]]$Length2>=lt[i] & events[[j]]$Length2<lt[i+1])
    Lcount[i,j]=length(I)
  }
  for(i in 1:14)
  {
    I=which(events[[j]]$Rad2>=rt[i] & events[[j]]$Rad2<rt[i+1])
    Rcount[i,j]=length(I)
  }
  for(i in 1:12) Mcount[i,j]=length(which(mm[[j]]==i))
}
Mcount2=Mcount
for(j in 1:8) Mcount2[,j]=100*Mcount2[,j]/length(mm[[j]])

apply(Mcount,2,sum)
apply(Mcount[5:10,],2,sum)/apply(Mcount,2,sum)

plot(1,NA,xlim=c(1,12),ylim=c(0,25),xlab="Month",ylab="Count")
for(i in 3:4) {
  lines(1:12,Mcount[,i],col=i)
  lines(1:12,Mcount[,i+4],col=i,lty=2)
}

Mcount3=matrix(0,4,8)
colnames(Mcount3)=colnames(Mcount2)
rownames(Mcount3)=c("DJF","MAM","JJA","SON")
Mcount3[1,]=apply(Mcount[c(12,1,2),],2,sum)
for(i in 1:3) Mcount3[i+1,]=apply(Mcount[seq(i*3,i*3+2),],2,sum)

abline(h=0,col="red")

########Now, try w/ Fixes



cvt=c(seq(0.5,3,0.5),Inf)
mt=c(seq(970,1020,10),Inf)
lt=c(seq(0,20,4),Inf)
rt=c(seq(0,6,0.5),Inf)

CVcount<-MSLPcount<-matrix(0,6,8)
Rcount=matrix(0,14,8)
Mcount=matrix(0,12,8)
colnames(CVcount)<-colnames(MSLPcount)<-colnames(Lcount)<-colnames(Rcount)<-colnames(Mcount)<-names(fixes)


for(j in 1:8)
{
  for(i in 1:6)
  {
    I=which(fixes[[j]]$CV>=cvt[i] & fixes[[j]]$CV<cvt[i+1] & fixes[[j]]$Location==1)
    CVcount[i,j]=length(I)
    I=which(fixes[[j]]$MSLP>=mt[i] & fixes[[j]]$MSLP<mt[i+1] & fixes[[j]]$Location==1)
    MSLPcount[i,j]=length(I)
  }
  for(i in 1:14)
  {
    I=which(fixes[[j]]$Radius>=rt[i] & fixes[[j]]$Radius<rt[i+1] & fixes[[j]]$Location==1)
    Rcount[i,j]=length(I)
  }
  for(i in 1:12) Mcount[i,j]=length(which(mm[[j]]==i  & fixes[[j]]$Location==1))
}
Mcount2=Mcount
for(j in 1:8) Mcount2[,j]=100*Mcount2[,j]/length(mm[[j]])

Mcount3=matrix(0,4,8)
colnames(Mcount3)=colnames(Mcount2)
rownames(Mcount3)=c("DJF","MAM","JJA","SON")
Mcount3[1,]=apply(Mcount[c(12,1,2),],2,sum)
for(i in 1:3) Mcount3[i+1,]=apply(Mcount[seq(i*3,i*3+2),],2,sum)

apply(Mcount,2,sum)
apply(Mcount[5:10,],2,sum)/apply(Mcount,2,sum)

plot(1,NA,xlim=c(1,12),ylim=c(0,25),xlab="Month",ylab="Count")
for(i in 3:4) {
  lines(1:12,Mcount[,i],col=i)
  lines(1:12,Mcount[,i+4],col=i,lty=2)
}
abline(h=0,col="red")

############ Before matching, need to fix event numbers


library(sp)
makematches<- function(fixes1,fixes2){
  matches<-matrix(NaN,max(fixes1$ID),max(fixes2$ID))
  for(i in 1:max(fixes1$ID))
  {
    data=fixes1[which(fixes1$ID==i & fixes1$Location==1),] # Only the instances of event 1 in location
    for(j in 1:max(fixes2$ID))
    {
      data2=fixes2[which(fixes2$ID==j),] # All instances of event 2
      match=rep(NaN,length(data[,1]))
      for(k in 1:length(data[,1]))
      {
        I=which(data2$Date==data$Date[k] & data2$Time==data$Time[k]) # If event is at same date/time
        if(length(I)>0) match[k]=min(spDistsN1(as.matrix(data2[I,7:8]),as.numeric(data[k,7:8]),longlat=TRUE)) # Radial distance
      }
      if(sum(!is.na(match))>0) matches[i,j]=min(match,na.rm=T) # Closest (spatial) distance between events
    } 
  }
  return(matches)
}

makematches2<- function(fixes1,fixes2,thresh=NA){
  if(!is.na(thresh)){
    fixes1=fixes1[fixes1$CV>=thresh,]
    fixes2=fixes2[fixes2$CV>=thresh,]
  }
  
  fixes1$CV2<-fixes1$MSLP2<-fixes1$Lat2<-fixes1$Lon2<-fixes1$ID2<-fixes1$Dist<-NaN
  for(i in 1:length(fixes1$ID))
  {      
    I=which(fixes2$Date==fixes1$Date[i] & fixes2$Time==fixes1$Time[i]) # If event is at same date/time
    if(length(I)>0) 
    {
      match=spDistsN1(as.matrix(fixes2[I,7:8]),as.numeric(fixes1[i,7:8]),longlat=TRUE) # Radial distance
      if(sum(!is.na(match))>0)
      {
        J=which(match==min(match))
        fixes1$Dist[i]=min(match)
        fixes1$ID2[i]=fixes2$ID[I[J[1]]]
        fixes1[i,17:20]=fixes2[I[J[1]],7:10]
      }  
    }
  }
  return(fixes1)
}

for(i in c(1,3,5,7))
{
  fixes1<-makematches2(fixes[[i]],fixes[[i+1]])
  I=which(fixes1$Dist<200)
  print(c(mean(fixes1$Lon2[I]-fixes1$Lon[I],na.rm=T),
          mean(fixes1$Lat2[I]-fixes1$Lat[I],na.rm=T),
          mean(fixes1$CV2[I]-fixes1$CV[I],na.rm=T),
          mean(fixes1$MSLP2[I]-fixes1$MSLP[I],na.rm=T)))
  
  #         a=density(fixes1$Lon2-fixes1$Lon,from=-10,to=10,na.rm=T)
  #         b=density(fixes1$Lat2-fixes1$Lat,from=-10,to=10,na.rm=T)
  #         plot(a,col=rgb(0,0,1,1/4),xlim=c(-10,10),ylim=range(0,a$y,b$y),
  #              xlab="Difference (deg)",ylab="Frequency",cex.main=1.2,
  #              main="Difference in low centre location (for lows at same time)")
  #         lines(a,col="blue",lwd=2)
  #         lines(b,col="red",lwd=2)
  #         abline(v=0,col="grey",lwd=2)
  #         legend("topleft",legend=c("Longitude","Latitude"),
  #                col=c("blue","red"),lwd=2,cex=1,bty="n")   
}


eventdiff=data.frame(ID=unique(fixes1$ID),Count=rep(NaN,length(unique(fixes1$ID))),
                     MeanLon=rep(NaN,length(unique(fixes1$ID))),MeanLat=rep(NaN,length(unique(fixes1$ID))),
                     MeanMSLP=rep(NaN,length(unique(fixes1$ID))),MeanCV=rep(NaN,length(unique(fixes1$ID))))
for(i in 1:length(eventdiff[,1]))
{
  I=which(fixes1$ID==i & fixes1$Location==1 & !is.na(fixes1$CV2) & fixes1$Dist<500)
  eventdiff[i,2]=length(I)
  if(length(I)>0) eventdiff[i,3:6]=apply(fixes1[I,17:20]-fixes1[I,7:10],2,mean)
}
apply(eventdiff[eventdiff[,2]>0,],2,mean)


lon=seq(150,160,2.5)
lat=seq(-40,-25,2.5)
meanlat<-meanlon<-count<-matrix(NaN,7,5)
rownames(meanlat)<-rownames(meanlon)<-lat
colnames(meanlat)<-colnames(meanlon)<-lon

for(i in 1:length(lat))
  for(j in 1:length(lon))
  {
    I=which(fixes1$Lat>=lat[i]-1.25 &  fixes1$Lat<lat[i]+1.25 &  
              fixes1$Lon>=lon[j]-1.25 &  fixes1$Lon<lon[j]+1.25 &
              !is.na(fixes1$CV2) & fixes1$Dist<200)
    if(length(I>0))
    {
      count[i,j]<-length(I)
      meanlat[i,j]<-mean(fixes1$Lat2[I]-fixes1$Lat[I])
      meanlon[i,j]<-mean(fixes1$Lon2[I]-fixes1$Lon[I])
    }
  }

library(sp)

makematches3<- function(fixes1,fixes2,thresh=NA,dthresh=0){
  if(!is.na(thresh)){
    fixes1=fixes1[fixes1$CV>=thresh,]
    fixes2=fixes2[fixes2$CV>=thresh,]
  }
  
  fixes1$CV2<-fixes1$MSLP2<-fixes1$Lat2<-fixes1$Lon2<-fixes1$ID2<-fixes1$TDiff<-fixes1$Dist<-NaN
  for(i in 1:length(fixes1$ID))
  {      
    I=which(fixes2$Date==fixes1$Date[i] & fixes2$Time==fixes1$Time[i]) # If event is at same date/time
    if(length(I)>0) 
    {
      match=spDistsN1(as.matrix(fixes2[I,7:8]),as.numeric(fixes1[i,7:8]),longlat=TRUE) # Radial distance
      if(sum(!is.na(match))>0)
      {
        J=which(match==min(match))
        fixes1$Dist[i]=min(match)
        fixes1$TDiff[i]=0
        fixes1$ID2[i]=fixes2$ID[I[J[1]]]
        fixes1[i,19:22]=fixes2[I[J[1]],7:10]
      }  
    } else {
      if(dthresh>0)
        for(t in seq(0.25,dthresh,0.25))
        {
          I=which(fixes2$Date2<=fixes1$Date2[i]+t & fixes2$Date2>=fixes1$Date2[i]-t) # If event is at similar date/time
          if(length(I)>0) 
          {
            match=spDistsN1(as.matrix(fixes2[I,7:8]),as.numeric(fixes1[i,7:8]),longlat=TRUE) # Radial distance
            if(sum(!is.na(match))>0)
            {
              J=which(match==min(match))
              fixes1$Dist[i]=min(match)
              fixes1$TDiff[i]=t
              fixes1$ID2[i]=fixes2$ID[I[J[1]]]
              fixes1[i,19:22]=fixes2[I[J[1]],7:10]
              break
            }  
          }
        }
    }
  }
  return(fixes1)
}

fixtest1<-makematches2(fixes[[3]][,-15],fixes[[4]][,-15])
fixtest2<-makematches3(fixes[[3]],fixes[[4]],dthresh=1) 
## Cool, all the 0s are the same as the results from last version, but now can have many more "matched"
## Only few extra fixes are "matched" within 200 km tho (for 3/4 -> 230/241 200km, or 382/432 500km)
## Turns into 4 extra events (38 v 42, of 52) matched in ESB region
b=data.frame(ID=unique(fixtest1$ID),fix1match=rep(0,length(unique(fixtest1$ID))),fix1match=rep(0,length(unique(fixtest1$ID))))
for(i in 1:length(b[,1]))
{
  I=which(fixtest1$ID==b[i,1] & fixtest1$Location==1)
  b[i,2]=length(which(fixtest1$Dist[I]<500))
  b[i,3]=length(which(fixtest2$Dist[I]<500))
}

mldb=read.csv("~/Documents/ECLs/MLDB_200708.csv")
mldb$Date2=as.Date(as.character(mldb$Date),"%Y%m%d")
mldb$Lat=-mldb$Lat #Convert to SH

mldbmatch=matrix(NaN,42,8)
for(j in 1:8)
{
  data=fixes[[j]]
  for(i in 1:42)
  {
    I=which(mldb$ID==i)
    for(k in 1:length(I))
    {
      J=which(data$Date2<=mldb$Date2[I[k]]+1 & data$Date2>=mldb$Date2[I[k]]-1) # If event is at similar date/time
      match=rep(NaN,length(I))
      if(length(J)>0) match[k]=min(spDistsN1(as.matrix(data[J,7:8]),as.numeric(mldb[I[k],3:4]),longlat=TRUE)) # Radial distance
    }
    if(sum(!is.na(match))>0) mldbmatch[i,j]=min(match,na.rm=T) 
  }
}

apply(!is.na(mldbmatch),2,sum)


########################
## Warm/cold core stuff? Using d01/p60 for now
## So cases 1/2
library(RNetCDF)
library(sp)

dir1="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007/out/slp/"
dir2="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007_notopo/out/slp/"
a=open.nc(paste(dir1,'WRF_d01_gph_200701.nc',sep=""))
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
lev=var.get.nc(a,"lev")
close.nc(a)

f1=make.WarmCold(fixes[[1]],events[[1]],dir1,500)
f2=make.WarmCold(fixes[[2]],events[[2]],dir2,500)

### All are cold or mixed, most are symmetric (lots with the B opposite, which suggests I've done something horribly wrong)

make.WarmCold <- function(fixes1,events1,dir1,dthresh)
{
  ftime1=fixes1$Date*100+(as.numeric(fixes1$Time)-1)*6
  yy=floor(fixes1$Date/10000)
  mm=floor(fixes1$Date/100)%%100
  fixes1$B2<-fixes1$VT<-fixes1$VL<-fixes1$VU<-fixes1$B<-fixes1$dir<-0
  
  dthresh2=round(dthresh/100,1)+1
  
  for(k in 1:length(fixes1[,1]))
  {  
    I=which(lat>=fixes1$Lat[k]-dthresh2 & lat<=fixes1$Lat[k]+dthresh2 &
              lon>=fixes1$Lon[k]-dthresh2 & lon<=fixes1$Lon[k]+dthresh2,arr.ind=T)
    lat2=lat[sort(unique(I[,1])),sort(unique(I[,2]))]
    lon2=lon[sort(unique(I[,1])),sort(unique(I[,2]))]
    
    a=open.nc(paste(dir1,'WRF_d01_gph_',yy[k],sprintf("%2.2d",mm[k]),'.nc',sep=""))
    ll=cbind(range(unique(I[,1])),range(unique(I[,2])))
    time=var.get.nc(a,"Times")
    tt=which(time==ftime1[k])
    gph2=var.get.nc(a,"Z",start=c(ll[1,],1,tt),count=c((ll[2,]-ll[1,]+1),3,1)) ##Only the bits we need
    close.nc(a)
    
    dist=array(NaN,dim(lat2))
    for(i in 1:dim(lat2)[1])
      for(j in 1:dim(lat2)[2])
        dist[i,j]=spDistsN1(cbind(lon2[i,j],lat2[i,j]),cbind(fixes1$Lon[k],fixes1$Lat[k]),longlat=T)
    
    I=which(dist<=dthresh,arr.ind=T) ## The radius used for identifying cyclone
    ##Need a NEW GPH which has three columns
    gph3<-matrix(0,length(I[,1]),3)
    for(i in 1:length(I[,1]))
      gph3[i,]=gph2[I[i,1],I[i,2],]
    
    Zpert=apply(gph3,2,max,na.rm=T)-apply(gph3,2,min,na.rm=T)
    fixes1$VU[k]=(Zpert[3]-Zpert[2])/(log(300)-log(600))
    fixes1$VL[k]=(Zpert[2]-Zpert[1])/(log(600)-log(900))
    fixes1$VT[k]=(Zpert[3]-Zpert[1])/(log(300)-log(900))
    
    dir=(270-atan2(fixes1[k,14],fixes1[k,13])*180/pi)%%360
    fixes1$dir[k]=dir
    
    J=which(dist==min(dist),arr.ind=T)
    W=which(I[,2]<J[2]) ##  Western (first) half
    E=which(I[,2]>J[2]) ##  Eastern half
    S=which(I[,1]<J[1]) ##  Southern (first) half
    N=which(I[,1]>J[1]) ##  Northern half 
    
    if(dir>=45 & dir<135) ## Heading to west, L = SOUTH
      fixes1$B[k]=mean(gph3[S,2]-gph3[S,1],na.rm=T)-mean(gph3[N,2]-gph3[N,1],na.rm=T) else 
        if(dir>=135 & dir<225) ## Heading to north, L=WEST
          fixes1$B[k]=mean(gph3[W,2]-gph3[W,1],na.rm=T)-mean(gph3[E,2]-gph3[E,1],na.rm=T) else 
            if(dir>=225 & dir<315) ## Heading to east, L=NORTH
              fixes1$B[k]=mean(gph3[N,2]-gph3[N,1],na.rm=T)-mean(gph3[S,2]-gph3[S,1],na.rm=T) else ## Heading to south, left=EAST
                fixes1$B[k]=mean(gph3[E,2]-gph3[E,1],na.rm=T)-mean(gph3[W,2]-gph3[W,1],na.rm=T)
    
    if(fixes1$Fix[k]!=1)
    {
      u=fixes1$Lon[k]-fixes1$Lon[k-1] # Positive = has moved east (L is N)
      v=fixes1$Lat[k]-fixes1$Lat[k-1] # Positive = has move north (L is E)
      
      if(abs(u)>abs(v))
      {
        if(u>0) fixes1$B2[k]=mean(gph3[N,2]-gph3[N,1],na.rm=T)-mean(gph3[S,2]-gph3[S,1],na.rm=T) else 
          fixes1$B2[k]=mean(gph3[S,2]-gph3[S,1],na.rm=T)-mean(gph3[N,2]-gph3[N,1],na.rm=T)
      }
      else
      {
        if(v>0) fixes1$B2[k]=mean(gph3[E,2]-gph3[E,1],na.rm=T)-mean(gph3[W,2]-gph3[W,1],na.rm=T) else 
          fixes1$B2[k]=mean(gph3[W,2]-gph3[W,1],na.rm=T)-mean(gph3[E,2]-gph3[E,1],na.rm=T)    
      }
      
    }
    
    
  }
  
  fixes1$Type="Mixed"
  fixes1$Type[fixes1$VU<(-5) & fixes1$VL<(-5)]="Cold"
  fixes1$Type[fixes1$VU>5 & fixes1$VL>5]="Warm"
  fixes1$Symm="Symmetric"
  fixes1$Symm[fixes1$B2>=10]="Asymmetric"
  
  events1$B2<-events1$VT<-events1$VL<-events1$VU<-events1$B<-NaN
  events1$Var=F
  for(i in 1:length(events1[,1]))
  {
    I=which(fixes1$ID==events1$ID[i])
    events1[i,12:16]=apply(fixes1[I,17:21],2,mean,na.rm=T)
    
    a=unique(fixes1[I,22])
    b=unique(fixes1[I,23])
    if(length(a)>1 | length(b)>1) events1$Var[i]=T 
  }
  events1$Type="Mixed"
  events1$Type[events1$VU<(-5) & events1$VL<(-5)]="Cold"
  events1$Type[events1$VU>5 & events1$VL>5]="Warm"
  events1$Symm="Symmetric"
  events1$Symm[events1$B2>=10]="Asymmetric"
  
  return(list(fixes1,events1))
}




