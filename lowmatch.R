##Script to match lows.
##FIRST - match individual low at a given time

setwd('~/Documents/ECLs')
library(sp)

HR<-read.csv('CSV/ECLfixes_unsw_cfsr_50.csv')
LR<-read.csv('CSV/ECLfixes_unsw_cfsr_250.csv')
matchlows(HR,LR,c(10,50,100,250))

###Try again with, say, NCEP1 and NCEP1
ncep1<-read.csv('CSV/ECLfixes_unsw_ncep1_250.csv')
ncep1$Time=sprintf("%2.2i:00",ncep1$Time)
ncep2<-read.csv('CSV/ECLfixes_unsw_ncep2_250.csv')
erai<-read.csv('CSV/ECLfixes_unsw_erai_300.csv')
jra25<-read.csv('CSV/ECLfixes_unsw_jra_250.csv')
merra<-read.csv('CSV/ECLfixes_unsw_merra_250.csv')
cfsr<-read.csv('CSV/ECLfixes_unsw_cfsr_250.csv')
wrf<-read.csv('CSV/ECLfixes_unsw_wrf_250.csv')
datasets=list(ncep1,ncep2,erai,jra25,merra,cfsr,wrf)

all<-matchtime<-matrix(0,7,7)
all2<-matchtime2<-matrix(0,7,7)
for(i in 1:7)
  for(j in 1:7)
  {
    a=matchlows(datasets[[i]],datasets[[j]],c(250,100))
    all[i,j]=a[1,2]
    matchtime[i,j]=a[1,3]
    all2[i,j]=a[2,2]
    matchtime2[i,j]=a[2,3]
  }

library(fields)

tiff(file=paste("matches_100_all.tiff",sep=""), height=500, width=500)
par(mar=c(5.1,2.1,4.1,4.1))
image(1:6,1:6,t(all2[1:6,1:6]),axes=F,col='transparent',xlab="",ylab="")
axis(1,at=1:6,labels=c("NCEP1","NCEP2","ERAI","JRA25","MERRA","CFSR"),cex=1.5)
axis(2,at=1:6,labels=c("NCEP1","NCEP2","ERAI","JRA25","MERRA","CFSR"),cex=1.5)
image.plot(1:6,1:6,t(all2[1:6,1:6]),add=T,legend.mar=3.1,zlim=c(0,1),xlab="",ylab="")
dev.off()
tiff(file=paste("matches_100_time.tiff",sep=""), height=500, width=500)
par(mar=c(5.1,2.1,4.1,4.1))
image(1:6,1:6,t(matchtime2[1:6,1:6]),axes=F,col='transparent',xlab="",ylab="")
axis(1,at=1:6,labels=c("NCEP1","NCEP2","ERAI","JRA25","MERRA","CFSR"),cex=1.5)
axis(2,at=1:6,labels=c("NCEP1","NCEP2","ERAI","JRA25","MERRA","CFSR"),cex=1.5)
image.plot(1:6,1:6,t(matchtime2[1:6,1:6]),add=T,legend.mar=3.1,zlim=c(0,1),xlab="",ylab="")
dev.off()

###Try again with, say, NCEP1 and WRF
ncep1<-read.csv('CSV/ECLfixes_unsw_ncep1_250.csv')
ncep1=ncep1[,c(-5,-6)]
wrf<-read.csv('CSV/ECLfixes_unsw_wrf_250.csv')
wrf=wrf[,-1]
wrf$Time=(as.numeric(wrf$Time)-1)*6 #Convert to integer
matchlows(ncep1,wrf,c(10,50,100,250,500,1000))
matchevents(ncep1[ncep1$CV>=0.25,],wrf[wrf$CV>=0.25,],c(10,50,100,250,500,1000))

##Now, mine and Alejandro's ERAI
mine<-read.csv('Algorithm Comparison/Mine_cv4.csv')
ale<-read.csv('Algorithm Comparison/Ale.csv')
mine$Time=(as.numeric(mine$Time)-1)*6 #Convert to integer
ale=ale[ale$Closed==-1,]
mldb<-read.csv('CSV/ECLfixes_mldb.csv')
names(mldb)=c("ID","Date","Lat","Lon","MSLP","Type","Sig","Bomb")
mldb$Time=0

matchlows2(mine,ale,c(10,50,100,150,250,500,1000))
matchlows2(mldb,mine,c(10,50,100,150,250,500,1000))
matchlows2(mldb,ale,c(10,50,100,150,250,500,1000))
matchlows2a(mldb,mine,c(10,50,100,150,250,500,1000))
matchlows2a(mldb,ale,c(10,50,100,150,250,500,1000))
matchevents(mine,ale,c(10,50,100,150,250,500,1000))
a<-matchevents_matrix(mine,mine$CV,ale,ale$Grad,150)

setwd("~/Documents/ECLs/Algorithm Comparison/")
#read.csv('Mine_cv4.csv')->mine
read.csv('Mine_rad2.csv')->mine
read.csv('Fei.csv')->Fei
read.csv('Ale_v15CTL.csv')->Ale
read.csv('MLDB.csv')->MLDB
#read.csv('UM_cv35.csv')->UM
read.csv('UM_rad2_p100.csv')->UM
mine$Time=(as.numeric(mine$Time)-1)*6 #Convert to integer
UM$Time=(as.numeric(UM$Time)-1)*6 #Convert to integer
datasets=list(MLDB,mine,UM,Ale,Fei)

all<-time<-events<-matrix(0,4,4)
rownames(all)<-rownames(time)<-rownames(events)<-c("MLDB","LAPB","LAPM","PG")
colnames(all)<-colnames(time)<-colnames(events)<-c("MLDB","LAPB","LAPM","PG")
for(i in 1:4)
  for(j in 1:4)
  {
    a=matchlows2(datasets[[i]],datasets[[j]],150)
    all[i,j]=a[2]
    time[i,j]=a[3]
    a=matchevents2a(datasets[[i]],datasets[[j]],150)
    events[i,j]=a[2]
  }

all<-matrix(0,4,4)
for(i in 1:4)
  {
    b=datasets[[i]]
    a=matchlows2(b[b$Time==0,],datasets[[1]],150)
    all[i,1]=a[2]
    a=matchlows2(b[b$Time==0,],datasets[[1]],500)
    all[i,2]=a[2]
    a=matchevents2a(b[b$Time==0,],datasets[[1]],150)
    all[i,3]=a[2]
    a=matchevents2a(b[b$Time==0,],datasets[[1]],500)
    all[i,4]=a[2]
  }


matchlows<-function(d1,d2,res)
{
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  ##So same time period & curvature threshold; only d1 is location-restricted
  d1=d1[y1>=ymin & y1<=ymax & d1$CV>=0.25 & d1$Location==1,]
  d2=d2[y2>=ymin & y2<=ymax & d2$CV>=0.25,]
  
  match=matrix(NaN,length(d1[,1]))
  

    for(i in 1:length(d1[,1])) 
    {
      I=which(d2$Date==d1$Date[i] & d2$Time==d1$Time[i])
      if(length(I)>0) match[i]=min(spDistsN1(as.matrix(cbind(d2$Lon[I],d2$Lat[I])),as.numeric(c(d1$Lon[i],d1$Lat[i])),longlat=TRUE))
    } 

  res=cbind(res,matrix(0,length(res),2))
  for(i in 1:length(res[,1])){
    res[i,2]=length(which(match<=res[i,1]))/length(match)
    res[i,3]=length(which(match<=res[i,1]))/sum(1-is.na(match))
  } 
  
  return(res)
}

matchlows2<-function(d1,d2,res)
{
  ##Where different datasets in terms of curvature etc, but need to create location
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  d1$Location=0
  I<-which(d1$Lon>=149 & d1$Lon<=161 & d1$Lat<(-37) & d1$Lat>=-41)
  d1$Location[I]<-1
  I<-which(d1$Lon>=(149+(37+d1$Lat)/2) & d1$Lon<=161 & d1$Lat<(-31) & d1$Lat>=-37)
  d1$Location[I]<-1
  I<-which(d1$Lon>=152 & d1$Lon<=161 & d1$Lat<=(-24) & d1$Lat>=-31)
  d1$Location[I]<-1
  
  d2$Location=0
  I<-which(d2$Lon>=149 & d2$Lon<=161 & d2$Lat<(-37) & d2$Lat>=-41)
  d2$Location[I]<-1
  I<-which(d2$Lon>=(149+(37+d2$Lat)/2) & d2$Lon<=161 & d2$Lat<(-31) & d2$Lat>=-37)
  d2$Location[I]<-1
  I<-which(d2$Lon>=152 & d2$Lon<=161 & d2$Lat<=(-24) & d2$Lat>=-31)
  d2$Location[I]<-1
  
  
  ll=min(length(d1[,1]),length(d2[,1]))
  
  ##So same time period & the main one is in location
#  if(length(d1[d1$Location==1,1])<=length(d2[d2$Location==1,1]))
#  {
    d1=d1[y1>=ymin & y1<=ymax & d1$Location==1,]
    d2=d2[y2>=ymin & y2<=ymax,]
    ll=length(d1[,1])
    match=matrix(NaN,ll)
    
    for(i in 1:length(d1[,1])) 
    {
      I=which(d2$Date==d1$Date[i] & d2$Time==d1$Time[i])
      if(length(I)>0) match[i]=min(spDistsN1(as.matrix(cbind(d2$Lon[I],d2$Lat[I])),as.numeric(c(d1$Lon[i],d1$Lat[i])),longlat=TRUE))
    } 
#  }
#   else
#   {
#     d1=d1[y1>=ymin & y1<=ymax,]
#     d2=d2[y2>=ymin & y2<=ymax & d2$Location==1,]
#     ll=length(d2[,1])
#     match=matrix(NaN,ll)
#     
#     for(i in 1:length(d2[,1])) 
#     {
#       I=which(d1$Date==d2$Date[i] & d1$Time==d2$Time[i])
#       if(length(I)>0) match[i]=min(spDistsN1(as.matrix(cbind(d1$Lon[I],d1$Lat[I])),as.numeric(c(d2$Lon[i],d2$Lat[i])),longlat=TRUE))
#     } 
#   }

  res=cbind(res,matrix(0,length(res),2))
  for(i in 1:length(res[,1])){
    res[i,2]=length(which(match<=res[i,1]))/ll
    res[i,3]=length(which(match<=res[i,1]))/sum(1-is.na(match))
  } 
  
  return(res)
}

matchlows2a<-function(d1,d2,res)
{
  ##Where different datasets in terms of curvature etc, but need to create location
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  d1$Location=0
  I<-which(d1$Lon>=149 & d1$Lon<=161 & d1$Lat<(-37) & d1$Lat>=-41)
  d1$Location[I]<-1
  I<-which(d1$Lon>=(149+(37+d1$Lat)/2) & d1$Lon<=161 & d1$Lat<(-31) & d1$Lat>=-37)
  d1$Location[I]<-1
  I<-which(d1$Lon>=152 & d1$Lon<=161 & d1$Lat<=(-24) & d1$Lat>=-31)
  d1$Location[I]<-1
  
  d2$Location=0
  I<-which(d2$Lon>=149 & d2$Lon<=161 & d2$Lat<(-37) & d2$Lat>=-41)
  d2$Location[I]<-1
  I<-which(d2$Lon>=(149+(37+d2$Lat)/2) & d2$Lon<=161 & d2$Lat<(-31) & d2$Lat>=-37)
  d2$Location[I]<-1
  I<-which(d2$Lon>=152 & d2$Lon<=161 & d2$Lat<=(-24) & d2$Lat>=-31)
  d2$Location[I]<-1
    
  ##So same DAY & the main one is in location
  d1=d1[y1>=ymin & y1<=ymax & d1$Location==1,]
  d2=d2[y2>=ymin & y2<=ymax,]
  ll=length(d1[,1])
  match=matrix(NaN,ll)
    
  for(i in 1:length(d1[,1])) 
    {
      I=which(d2$Date==d1$Date[i])
      if(length(I)>0) match[i]=min(spDistsN1(as.matrix(cbind(d2$Lon[I],d2$Lat[I])),as.numeric(c(d1$Lon[i],d1$Lat[i])),longlat=TRUE))
    } 
  
  res=cbind(res,matrix(0,length(res),2))
  for(i in 1:length(res[,1])){
    res[i,2]=length(which(match<=res[i,1]))/ll
    res[i,3]=length(which(match<=res[i,1]))/sum(1-is.na(match))
  } 
  
  return(res)
}


matchevents<-function(d1,d2,res)
{
  ##Where different datasets in terms of curvature etc, but need to create location
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  d1$Location=0
  I<-which(d1$Lon>=149 & d1$Lon<=161 & d1$Lat<(-37) & d1$Lat>=-41)
  d1$Location[I]<-1
  I<-which(d1$Lon>=(149+(37+d1$Lat)/2) & d1$Lon<=161 & d1$Lat<(-31) & d1$Lat>=-37)
  d1$Location[I]<-1
  I<-which(d1$Lon>=152 & d1$Lon<=161 & d1$Lat<=(-24) & d1$Lat>=-31)
  d1$Location[I]<-1
  
  d2$Location=0
  I<-which(d2$Lon>=149 & d2$Lon<=161 & d2$Lat<(-37) & d2$Lat>=-41)
  d2$Location[I]<-1
  I<-which(d2$Lon>=(149+(37+d2$Lat)/2) & d2$Lon<=161 & d2$Lat<(-31) & d2$Lat>=-37)
  d2$Location[I]<-1
  I<-which(d2$Lon>=152 & d2$Lon<=161 & d2$Lat<=(-24) & d2$Lat>=-31)
  d2$Location[I]<-1
  
  x<-rle(d1$ID)
  ev1<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),MinMatch=rep(NaN,length(x$values)))
  for(i in 1:length(ev1[,1])) ev1[i,3]=max(d1$Location[d1$ID==ev1[i,1]])
  ev1=ev1[ev1$Location==1,]
  
  x<-rle(d2$ID)
  ev2<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),MinMatch=rep(NaN,length(x$values)))
  for(i in 1:length(ev2[,1])) ev2[i,3]=max(d2$Location[d2$ID==ev2[i,1]])
  ev2=ev2[ev2$Location==1,]

  match1=matrix(NaN,length(d1[,1]))
  for(i in 1:length(d1[,1])) 
    {
      I=which(d2$Date==d1$Date[i] & d2$Time==d1$Time[i])
      if(length(I)>0) match1[i]=min(spDistsN1(as.matrix(d2[I,7:8]),as.numeric(d1[i,7:8]),longlat=TRUE))
    } 
  for(i in 1:length(ev1[,1])) ev1[i,4]=min(match1[d1$ID==ev1[i,1] & d1$Location==1,],na.rm=T)
  
  match2=matrix(NaN,length(d2[,1]))
  for(i in 1:length(d2[,1])) 
  {
    I=which(d1$Date==d2$Date[i] & d1$Time==d2$Time[i])
    if(length(I)>0) match2[i]=min(spDistsN1(as.matrix(d1[I,7:8]),as.numeric(d2[i,7:8]),longlat=TRUE))
  } 
  for(i in 1:length(ev2[,1])) ev2[i,4]=min(match2[d2$ID==ev2[i,1] & d2$Location==1,],na.rm=T)

  res=cbind(res,matrix(0,length(res),2))
  for(i in 1:length(res[,1]))
    {
    res[i,2]=length(which(ev1[,4]<=res[i,1]))/length(ev1[,4]) 
    res[i,3]=length(which(ev2[,4]<=res[i,1]))/length(ev2[,4]) 
  }
  
  return(res)
}

matchevents2a<-function(d1,d2,res)
{
  ##Where different datasets in terms of curvature etc, but need to create location
  ##Designed for comparing against mldb - only asks for same day
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  d1$Location=0
  I<-which(d1$Lon>=149 & d1$Lon<=161 & d1$Lat<(-37) & d1$Lat>=-41)
  d1$Location[I]<-1
  I<-which(d1$Lon>=(149+(37+d1$Lat)/2) & d1$Lon<=161 & d1$Lat<(-31) & d1$Lat>=-37)
  d1$Location[I]<-1
  I<-which(d1$Lon>=152 & d1$Lon<=161 & d1$Lat<=(-24) & d1$Lat>=-31)
  d1$Location[I]<-1
  
  d2$Location=0
  I<-which(d2$Lon>=149 & d2$Lon<=161 & d2$Lat<(-37) & d2$Lat>=-41)
  d2$Location[I]<-1
  I<-which(d2$Lon>=(149+(37+d2$Lat)/2) & d2$Lon<=161 & d2$Lat<(-31) & d2$Lat>=-37)
  d2$Location[I]<-1
  I<-which(d2$Lon>=152 & d2$Lon<=161 & d2$Lat<=(-24) & d2$Lat>=-31)
  d2$Location[I]<-1
  
  x<-rle(d1$ID)
  ev1<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),MinMatch=rep(NaN,length(x$values)))
  for(i in 1:length(ev1[,1])) ev1[i,3]=max(d1$Location[d1$ID==ev1[i,1]])
  ev1=ev1[ev1$Location==1,]
  
  match1=matrix(NaN,length(d1[,1]))
  for(i in 1:length(d1[,1])) 
  {
    I=which(d2$Date==d1$Date[i])
    if(length(I)>0) match1[i]=min(spDistsN1(as.matrix(cbind(d2$Lon[I],d2$Lat[I])),as.numeric(c(d1$Lon[i],d1$Lat[i])),longlat=TRUE))
  } 
  for(i in 1:length(ev1[,1])) ev1[i,4]=min(match1[d1$ID==ev1[i,1] & d1$Location==1,],na.rm=T)
  
  res=cbind(res,rep(0,length(res)))
  for(i in 1:length(res[,1])) res[i,2]=length(which(ev1[,4]<=res[i,1]))/length(ev1[,4]) 
 
  return(res)
}

matchevents_matrix<-function(d1,d1I,d2,d2I,res)
{
  ##Where different datasets in terms of curvature etc, but need to create location
  ##Compares whether the event is ever matched within res
  ##d1I and d2I are the intensity measures
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  d1$Location=0
  I<-which(d1$Lon>=149 & d1$Lon<=161 & d1$Lat<(-37) & d1$Lat>=-41)
  d1$Location[I]<-1
  I<-which(d1$Lon>=(149+(37+d1$Lat)/2) & d1$Lon<=161 & d1$Lat<(-31) & d1$Lat>=-37)
  d1$Location[I]<-1
  I<-which(d1$Lon>=152 & d1$Lon<=161 & d1$Lat<=(-24) & d1$Lat>=-31)
  d1$Location[I]<-1
  d1$Intensity=d1I
  
  d2$Location=0
  I<-which(d2$Lon>=149 & d2$Lon<=161 & d2$Lat<(-37) & d2$Lat>=-41)
  d2$Location[I]<-1
  I<-which(d2$Lon>=(149+(37+d2$Lat)/2) & d2$Lon<=161 & d2$Lat<(-31) & d2$Lat>=-37)
  d2$Location[I]<-1
  I<-which(d2$Lon>=152 & d2$Lon<=161 & d2$Lat<=(-24) & d2$Lat>=-31)
  d2$Location[I]<-1
  d2$Intensity=d2I
  
  x<-rle(d1$ID)
  ev1<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),Intensity=rep(0,length(x$values)),MinMatch=rep(NaN,length(x$values)))
  for(i in 1:length(ev1[,1]))
    {
    ev1[i,3]=max(d1$Location[d1$ID==ev1[i,1]])
    ev1[i,4]=max(d1$Intensity[d1$ID==ev1[i,1]])
  } 
  ev1=ev1[ev1$Location==1,]
  
  x<-rle(d2$ID)
  ev2<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),Intensity=rep(0,length(x$values)),MinMatch=rep(NaN,length(x$values)))
  for(i in 1:length(ev2[,1]))
  {
    ev2[i,3]=max(d2$Location[d2$ID==ev2[i,1]])
    ev2[i,4]=max(d2$Intensity[d2$ID==ev2[i,1]])
  } 
  ev2=ev2[ev2$Location==1,]
  
  lenmin=c(0,5,9)
  imin1=c(0,quantile(ev1[,4],c(0.5,0.9)))
  imin2=c(0,quantile(ev2[,4],c(0.5,0.9)))
  
  match<-matchT<-count<-count2<-matrix(NaN,3,3)
  for(i in 1:3)
    for(j in 1:3)
    {
      ev1a=ev1[(ev1$Length>=lenmin[j] & ev1$Intensity>=imin1[i]),]
      ev2a=ev2[(ev2$Length>=lenmin[j] & ev2$Intensity>=imin2[i]),]
      
      include<-match(d1$ID,ev1a$ID)
      d1a=d1[which(is.na(include)==0),]
      include<-match(d2$ID,ev2a$ID)
      d2a=d2[which(is.na(include)==0),]
      count[i,j]=length(ev1a[,1])
      count2[i,j]=length(ev2a[,1])
      
      match1=rep(NaN,length(d1a[,1]))
      for(k in 1:length(d1a[,1])) 
      {
        I=which(d2a$Date==d1a$Date[k] & d2a$Time==d1a$Time[k])
        if(length(I)>0) match1[k]=min(spDistsN1(as.matrix(d2a[I,7:8]),as.numeric(d1a[k,7:8]),longlat=TRUE))
      } 
      for(k in 1:length(ev1a[,1]))
        {
        I=is.finite(match1[(d1a$ID==ev1a[k,1] & d1a$Location==1)])
        if(sum(I)>0) ev1a[k,5]=min(match1[d1a$ID==ev1a[k,1] & d1a$Location==1],na.rm=T)
        }
      match[i,j]=length(which(ev1a[,5]<=res))/length(ev1a[,4]) 
      matchT[i,j]=length(which(ev1a[,5]<=res))/sum(1-is.na(ev1a[,5]))
    }
  return(list(imin1,imin2,match,matchT))
}

matchevents_matrix2<-function(d1,d1I,d2,d2I,res)
{
  ##Where different datasets in terms of curvature etc, but need to create location
  ##Compares whether the event is ever matched within res
  ##d1I and d2I are the intensity measures
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  d1$Location=0
  I<-which(d1$Lon>=149 & d1$Lon<=161 & d1$Lat<(-37) & d1$Lat>=-41)
  d1$Location[I]<-1
  I<-which(d1$Lon>=(149+(37+d1$Lat)/2) & d1$Lon<=161 & d1$Lat<(-31) & d1$Lat>=-37)
  d1$Location[I]<-1
  I<-which(d1$Lon>=152 & d1$Lon<=161 & d1$Lat<=(-24) & d1$Lat>=-31)
  d1$Location[I]<-1
  d1$Intensity=d1I
  
  d2$Location=0
  I<-which(d2$Lon>=149 & d2$Lon<=161 & d2$Lat<(-37) & d2$Lat>=-41)
  d2$Location[I]<-1
  I<-which(d2$Lon>=(149+(37+d2$Lat)/2) & d2$Lon<=161 & d2$Lat<(-31) & d2$Lat>=-37)
  d2$Location[I]<-1
  I<-which(d2$Lon>=152 & d2$Lon<=161 & d2$Lat<=(-24) & d2$Lat>=-31)
  d2$Location[I]<-1
  d2$Intensity=d2I
  
  x<-rle(d1$ID)
  ev1<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),Intensity=rep(0,length(x$values)),MinMatch=rep(NaN,length(x$values)))
  for(i in 1:length(ev1[,1]))
  {
    ev1[i,3]=max(d1$Location[d1$ID==ev1[i,1]])
    ev1[i,4]=max(d1$Intensity[d1$ID==ev1[i,1]])
  } 
  ev1=ev1[ev1$Location==1,]
  
  x<-rle(d2$ID)
  ev2<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),Intensity=rep(0,length(x$values)),MinMatch=rep(NaN,length(x$values)))
  for(i in 1:length(ev2[,1]))
  {
    ev2[i,3]=max(d2$Location[d2$ID==ev2[i,1]])
    ev2[i,4]=max(d2$Intensity[d2$ID==ev2[i,1]])
  } 
  ev2=ev2[ev2$Location==1,]
  
  lenmin=c(0,5,9)
  imin1=c(0,quantile(ev1[,4],c(0.5,0.9)))
  imin2=c(0,quantile(ev2[,4],c(0.5,0.9)))
  
  match_d1<-match_d2<-count1<-count2<-matrix(NaN,3,3)
  for(i in 1:3)
    for(j in 1:3)
    {
      ev1a=ev1[(ev1$Length>=lenmin[j] & ev1$Intensity>=imin1[i]),]
      ev2a=ev2[(ev2$Length>=lenmin[j] & ev2$Intensity>=imin2[i]),]
      
      include<-match(d1$ID,ev1a$ID)
      d1a=d1[which(is.na(include)==0),]
      include<-match(d2$ID,ev2a$ID)
      d2a=d2[which(is.na(include)==0),]
      count1[i,j]=length(ev1a[,1])
      count2[i,j]=length(ev2a[,1])
      
      match1=rep(NaN,length(d1a[,1]))
      for(k in 1:length(d1a[,1])) 
      {
        I=which(d2$Date==d1a$Date[k] & d2$Time==d1a$Time[k])
        if(length(I)>0) match1[k]=min(spDistsN1(as.matrix(d2[I,7:8]),as.numeric(d1a[k,7:8]),longlat=TRUE))
      } 
      for(k in 1:length(ev1a[,1]))
      {
        I=is.finite(match1[(d1a$ID==ev1a[k,1] & d1a$Location==1)])
        if(sum(I)>0) ev1a[k,5]=min(match1[d1a$ID==ev1a[k,1] & d1a$Location==1],na.rm=T)
      }
      match_d1[i,j]=length(which(ev1a[,5]<=res))/length(ev1a[,4]) 
      
      match2=rep(NaN,length(d2a[,1]))
      for(k in 1:length(d2a[,1])) 
      {
        I=which(d1$Date==d2a$Date[k] & d1$Time==d2a$Time[k])
        if(length(I)>0) match2[k]=min(spDistsN1(as.matrix(d1[I,7:8]),as.numeric(d2a[k,7:8]),longlat=TRUE))
      } 
      for(k in 1:length(ev2a[,1]))
      {
        I=is.finite(match2[(d2a$ID==ev2a[k,1] & d2a$Location==1)])
        if(sum(I)>0) ev2a[k,5]=min(match1[d2a$ID==ev2a[k,1] & d2a$Location==1],na.rm=T)
      }
      match_d2[i,j]=length(which(ev2a[,5]<=res))/length(ev2a[,4])
    }
  return(list(imin1,imin2,match_d1,match_d2,count1,count2))
}


