library(R.matlab)
library(RNetCDF)
library(abind)
setwd("~/Documents/ECLs/StuartBrowning/")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
readMat("ECC_DatabaseMWR2013.mat")->Stuart
events<-Stuart$Xevents[,c(1,5:7,9,10)]
colnames(events)<-c("ID","Date1","Date2","Type","Duration","PG max")
types<-c("ET","SSL","IT","CL","XTC")
fixes<-Stuart$Xloc[,c(1,3:9)]
colnames(fixes)<-c("ID","Year","Month","Day","Time","Lag","Lat","Lon")
for(i in 1:5) print(paste(types[i],"-",length(which(events[,4]==i))))

fixes<-as.data.frame(fixes)
events<-as.data.frame(events)
fixes$Date=fixes$Year*10000+fixes$Month*100+fixes$Day

fixes$Location=0
I<-which(fixes$Lon>=149 & fixes$Lon<=161 & fixes$Lat<(-37) & fixes$Lat>=-41)
fixes$Location[I]<-1
I<-which(fixes$Lon>=(149+(37+fixes$Lat)/2) & fixes$Lon<=161 & fixes$Lat<(-31) & fixes$Lat>=-37)
fixes$Location[I]<-1
I<-which(fixes$Lon>=152 & fixes$Lon<=161 & fixes$Lat<=(-24) & fixes$Lat>=-31)
fixes$Location[I]<-1

write.csv(events,"~/Documents/ECLs/StuartBrowning/Events_SB.csv")
write.csv(fixes,"~/Documents/ECLs/StuartBrowning/Fixes_SB.csv")

#######
# Okay, now we have the data, do two things:
# 1. Run my version of the tracking on Stuart's events, and see how well it matches his results
# 2. Match his events with MLDB events & report the typing 9as opposed to typing from scratch)

IDs=unique(fixes$ID)

i=1
I=which(fixes$ID==IDs[i] & fixes$Lag==0)
firstfixes=fixes[I[1],]
for(i in 2:length(IDs))
{
  I=which(fixes$ID==IDs[i] &  fixes$Lag>=0)
  if(length(I)>0) firstfixes=rbind(firstfixes,fixes[I[1],])
}
firstfixes$Date2=as.Date(as.character(firstfixes$Date*100+firstfixes$Time),"%Y%m%d%H")

type1=types[events$Type[firstfixes[,1]]]

I=which(fixes$ID==1 & fixes$Lag<=0)
backtrack1=as.matrix(fixes[I,8:7])
dim(backtrack1)<-c(1,9,2)
for(i in 2:length(IDs))
{
  I=which(fixes$ID==i & fixes$Lag<=0)
  backtrack1=abind(backtrack1,fixes[I,8:7],along=1)
}
backtrack1=backtrack1[,9:1,]
type1a=run_typing(backtrack1) ## Running my typing on Stuart's backtrack
type1b=run_typing_v2(backtrack1) ## Running my typing with Stuart's orog on Stuart's backtrack
backtrack2=run_backtrack(firstfixes)
type2=run_typing(backtrack2)  ## Running my typing on my backtrack
type2a=run_typing_v2(backtrack2)  ## Running my typing on my backtrack

typecomp=matrix(0,4,4)
for(i in 1:4)
  for(j in 1:4)
    typecomp[i,j]=length(which(type1==types[i] & type1b==types[j]))

print(length(which(type1b==type1))/length(which(type1!="XTC")))

J=which(type1!=Type)
for(i in 1:20)
{
I=which(fixes$ID==J[i] & fixes$Lag<=0)
lines(fixes$Lon[I],fixes$Lat[I],type="l",col=events$Type[J[i]],lwd=3)
}
legend("topleft",types[1:5],col=1:5,lwd=2)

### Testing different land mask
## Lazy GDR boundary

lat2<-seq(-10,-40,-1.5)
GDR<-rep(0,length(lat2))
I=which(lat2>(-25) & lat2<=-10)
GDR[I]=151-(25+lat2[I])
I=which(lat2<=(-25) & lat2>=-31)
GDR[I]=151
I=which(lat2>=(-39) & lat2<(-31))
GDR[I]=148+(37+lat2[I])/2
lines(GDR,lat2,col="grey",lwd=3)

## Stuart's mask

a=open.nc("~/Documents/ECLs/StuartBrowning/land_mask.nc")
landmask=var.get.nc(a,"lsm",unpack=T)
lat=var.get.nc(a,"latitude")
lon=var.get.nc(a,"longitude")
a=open.nc("~/Documents/ECLs/StuartBrowning/orography.nc")
orog=var.get.nc(a,"z",unpack=T)

I=which(lon>=140 & lon<=155)
J=which(lat<=-12 & lat>=-38)
orog2=orog[I,J]
GDR2=rep(0,length(J))
for(i in 1:length(J))
{
  k=which(orog2[,i]==max(orog2[,i]))
  GDR2[i]=lon[I[k]]
}

landmask2=landmask
landmask2[lon>=124.5 & lon<=139.5,lat>=-36 & lat<=-33]=1 #Bight is land
landmask2[lon>=120 & lon<=132,lat>=-15 & lat<=-9]=1 #NW WA is land

run_backtrack<-function(firstfixes,DL=4,DS=8)
{
  years=floor(firstfixes$Date/10000)
  backtrack=array(0,c(length(firstfixes[,1]),DS+1,2))
  dimnames(backtrack)[[3]]=c("Lon","Lat")
  backtrack[,1,]=as.matrix(round(c(firstfixes$Lon,firstfixes$Lat)/1.5)*1.5)
  
  ### Okay, now set up to loop through ALL
  ylist=c(1979,1984,1989,1994,1999,2005,2011,2015)
  yold=years[1]
  I=which(yold>=ylist[1:7] & yold<=ylist[2:8])
  file=open.nc(paste("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_",ylist[I],"-01_",ylist[I+1]-1,"-12.nc",sep=""))
  lat=var.get.nc(file,"latitude")
  lon=var.get.nc(file,"longitude")
  dates=seq.Date(as.Date(paste(ylist[I],"010100",sep=""),"%Y%m%d%H"),as.Date(paste(ylist[I+1]-1,"123118",sep=""),"%Y%m%d%H"),0.25)
  
  ####Do stuart's lazy backtracking
  
  for(i in 1:length(firstfixes[,1]))
  {
    ynew=years[i]
    I=which(yold>=ylist[1:7] & yold<ylist[2:8])
    J=which(ynew>=ylist[1:7] & ynew<ylist[2:8])  
    if(I!=J) ## Test if need to open a new ERAI file
    {
      close.nc(file)
      file=open.nc(paste("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_",ylist[J],"-01_",ylist[J+1]-1,"-12.nc",sep=""))
      dates=seq.Date(as.Date(paste(ylist[J],"010100",sep=""),"%Y%m%d%H"),as.Date(paste(ylist[J+1]-1,"123118",sep=""),"%Y%m%d%H"),0.25)
    }
    yold=ynew
    
    I=which(dates==firstfixes$Date2[i])
    slp=var.get.nc(file,"msl",start=c(1,1,I-8),count=c(length(lon),length(lat),8),unpack=T)/100
    
    for(j in 1:DS)
    {
      J=which(lon==backtrack[i,j,1])
      K=which(lat==backtrack[i,j,2])
      
      n=1
      tmp=array(0,c((DL*2+1)^2,5))
      for(p in (J-DL):(J+DL))
        for(q in (K-DL):(K+DL))
        {    
          if(p>2 & q>2 & p<length(lon)-1 & q<length(lat)-1)
          {
            CELL = slp[p,q,9-j];
            
            # MEASURE PRESSURE GRADIENTS AROUND TEST CELL
            test = c(slp[p-2,q,9-j]-CELL, slp[p-2,q+2,9-j]-CELL, slp[p-2,q-2,9-j]-CELL, 
                     slp[p+2,q,9-j]-CELL, slp[p+2,q+2,9-j]-CELL, slp[p+2,q-2,9-j]-CELL, 
                     slp[p,q-2,9-j]-CELL, slp[p,q+2,9-j]-CELL)
            
            tmp[n,] = c(length(which(test>0)),max(test),lon[p],lat[q],CELL)
          }
          n=n+1
        }  
      a=order(-tmp[,1],tmp[,5])
      backtrack[i,j+1,]=tmp[a[1],3:4]
    }
  } 
  return(backtrack)
}

run_typing<-function(backtrack,DL=4,DS=8)
{
  motion=array(0,c(length(backtrack[,1,1]),8,3))
  for(i in 1:8) 
  {
    motion[,i,1:2]=backtrack[,i,]-backtrack[,i+1,]
    
    ##I'll use my lazy GDR boundary for ECLs from before!
    I<-which(backtrack[,i+1,1]<(151-(25+backtrack[,i+1,2])) & backtrack[,i+1,2]>(-25) & backtrack[,i+1,2]<=-10) ## Over more land somehow?
    if(length(I)>0) motion[I,i,3]<-1
    I<-which(backtrack[,i+1,1]<151 & backtrack[,i+1,2]<=(-25) & backtrack[,i+1,2]>=-31)
    if(length(I)>0) motion[I,i,3]<-1
    I<-which(backtrack[,i+1,2]<(-31) & backtrack[,i+1,2]>=-39 & backtrack[,i+1,1]<(148+(37+backtrack[,i+1,2])/2))
    if(length(I)>0) motion[I,i,3]<-1
  }
  
  typing=cbind(apply(motion[,,1:2]*(motion[,,1:2]>1),c(1,3),sum),apply(motion[,,1:2]*(motion[,,1:2]<(-1)),c(1,3),sum)*-1,apply(motion[,,3],1,sum)/DS,matrix(0,length(backtrack[,1,1]),2))
  colnames(typing)=c("E","N","W","S","Land prop","N points","S points")
  for(i in 1:length(typing[,1])){
    typing[i,6]=length(which(backtrack[i,2:9,2]>=-27))
    typing[i,7]=length(which(backtrack[i,2:9,2]<=-39))
  } 
  
  Type=rep("NA",length(typing[,1]))
  I=which(typing[,5]<=0.5 & (typing[,4]>typing[,2] | typing[,7]<2) )
  Type[I]="ET"
  I=which(typing[,5]<=0.5 & typing[,2]>=typing[,4] & typing[,7]>=2)
  Type[I]="SSL"
  I=which(typing[,5]>0.5 & typing[,6]>2)
  Type[I]="IT"
  I=which(typing[,5]>0.5 & typing[,6]<=2)
  Type[I]="CL"
  
  return(Type)
}

run_typing_v2<-function(backtrack,DL=4,DS=8)
{
  motion=array(0,c(length(backtrack[,1,1]),8,3))
  for(i in 1:8) 
  {
    motion[,i,1:2]=backtrack[,i,]-backtrack[,i+1,]
    
    # Using Stuart's landmask
    for(k in 1:length(backtrack[,1,1]))
      motion[k,i,3]=landmask[lon==backtrack[k,i+1,1],lat==backtrack[k,i+1,2]]

    ## East/west of divide
    
    lat2=lat[lat<=-12 & lat>=-38]
    lon2=lon[lon>=140 & lon<=155]
    
    I=which((backtrack[,i+1,2]>=-38 & backtrack[,i+1,2]<=-12) & motion[,i,3]==1)
    for(k in 1:length(I))
    {
      gdr=GDR2[which(lat2==backtrack[I[k],i+1,2])]
      if(backtrack[I[k],i+1,1]>=(gdr-1.5)) motion[I[k],i,3]==0 ## If east of GDR, change land-sea to ocean
    }
  }
  
  typing=cbind(apply(motion[,,1:2]*(motion[,,1:2]>1),c(1,3),sum),apply(motion[,,1:2]*(motion[,,1:2]<(-1)),c(1,3),sum)*-1,apply(motion[,,3],1,sum)/DS,matrix(0,length(backtrack[,1,1]),2))
  colnames(typing)=c("E","N","W","S","Land prop","N points","S points")
  for(i in 1:length(typing[,1])){
    typing[i,6]=length(which(backtrack[i,2:9,2]>=-25))
    typing[i,7]=length(which(backtrack[i,2:9,2]<=-39))
  } 
  
  Type=rep("NA",length(typing[,1]))
  I=which(typing[,5]<=0.5 & (typing[,4]>typing[,2] | typing[,7]<2) )
  Type[I]="ET"
  I=which(typing[,5]<=0.5 & typing[,2]>=typing[,4] & typing[,7]>=2)
  Type[I]="SSL"
  I=which(typing[,5]>0.5 & typing[,6]>2)
  Type[I]="IT"
  I=which(typing[,5]>0.5 & typing[,6]<=2)
  Type[I]="CL"
  
  I=which(Type=="ET" & backtrack[,1,2]<=-35)
  Type[I]="SSL"
  I=which(Type=="SSL" & backtrack[,1,2]>=-25)
  Type[I]="ET"
  
  return(Type)
}


######  Version 2 - Actual matching of lows

mldb=read.csv("~/Documents/ECLs/Algorithm Comparison/MLDB.csv",stringsAsFactors=F)
mldb$Date2=as.POSIXct(as.character(mldb$Date*100+mldb$Time),format="%Y%m%d%H")
mldb$MatchID=NaN
fixes$Date2=as.POSIXct(as.character(fixes$Date*100+fixes$Time),format="%Y%m%d%H")
library(sp)

IDs=data.frame(IDs=unique(mldb$ID),Year=rep(NaN,length(unique(mldb$ID))),Month=rep(NaN,length(unique(mldb$ID))),Type=rep("aaa",length(unique(mldb$ID))),MatchID=rep(NaN,length(unique(mldb$ID))),MatchType=rep("aaa",length(unique(mldb$ID))),stringsAsFactors=F)
dist=500
tdist=24

for(i in 1:length(IDs$IDs))
{
  I=which(mldb$ID==IDs$IDs[i])
  IDs$Year[i]=floor(mldb$Date[I[1]]/10000)
  IDs$Month[i]=floor(mldb$Date[I[1]]/100)%%100
  IDs$Type[i]=mldb$Type[I[1]]
  
  for(j in 1:length(I))
  {
    td=difftime(fixes$Date2,mldb$Date2[I[j]],units="hours")
    J=which(abs(td)<=tdist)
    jj=0
    if(length(J)>0) 
      {
      jj=jj+1
      tmp=spDistsN1(as.matrix(cbind(fixes$Lon[J],fixes$Lat[J])),as.numeric(c(mldb$Lon[I[j]],mldb$Lat[I[j]])),longlat=TRUE) 
      if(jj==1) matches=data.frame(ID=fixes$ID[J],Date=fixes$Date2[J],TimeDiff=td[J],Dist=tmp) else 
        matches=rbind(matches,data.frame(ID=fixes$ID[J],Date=fixes$Date2[J],TimeDiff=td[J],Dist=tmp))
    }
  }
   
  if(exists("matches"))
  {
  J=unique(matches[which(matches[,4]<=dist),1])
  if(length(J)==1) IDs$MatchID[i]=J else {
    if(length(J)>1) {
    t=table(matches[,1])
    TT=which(t==max(t))
    IDs$MatchID[i]=as.numeric(names(t)[TT[1]])}}
   rm(matches)
  }
}

I=which(!is.na(IDs$MatchID))
IDs$MatchType[I]=types[events$Type[IDs$MatchID[I]]]

typesA=c("XTC","ET","IT","SSL","CL")
typesB=c("TC","ET","IT","WF","LW","DF")
typecomp=matrix(0,5,6)
for(i in 1:5)
  for(j in 1:6)
    typecomp[i,j]=length(which(IDs$Type==typesB[j] & IDs$MatchType==typesA[i]))
rownames(typecomp)=typesA
colnames(typecomp)=typesB

for(i in 1:6)
{
  I=which(IDs$Type==typesB[i])
  J=which(!is.na(IDs$MatchID[I]))
  print(paste(typesB[i],100*length(J)/length(I)))
}

library(reshape)
library(ggplot2)
test=melt(typecomp)

ggplot(test, aes(x=X2, fill=X1, y=value)) + geom_bar(stat="identity") +
  theme_bw() + ylab("Number of events") + xlab("Speer classification") + labs(fill="Browning classification")
