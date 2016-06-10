rm(list=ls())

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful2$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
setwd("~/Documents/ECLs")

dates=seq.Date(as.Date("2007060100","%Y%m%d%H"),as.Date("2007063018","%Y%m%d%H"),0.25)
cat=c("d01_p100_rad2cv05","d01_p100_rad5cv015","d02_p100_rad2cv05","d02_p100_rad5cv015")
SSTtype=c("ERAI","BRAN","NoEAC","DoubleEAC")

##From stewart:
#Set days to backtrack = 2 (i.e. 8 timesteps)
DS=8
#Set region for tracking - ~ surrounding 6 degrees, so +- 12 0.5 deg cells
DL=12


c=1
r=1
day=27
hour="00"

dir=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep="")
dir2=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_ensemble/out/200705",day,hour,"/",sep="")

data=read.csv(paste(dir,"ECLfixes_200705",day,hour,".csv",sep=""))
i=1
I=which(data$ID==i & data$Location==1)
firstfixes=data[I[1],2:10]
for(i in 2:length(unique(data$ID)))
{
  I=which(data$ID==i & data$Location==1)
  firstfixes=rbind(firstfixes,firstfixes=data[I[1],2:10])
}
firstfixes$Date2=as.Date(as.character(firstfixes$Date*100+(as.numeric(firstfixes$Time)-1)*6),"%Y%m%d%H")

a=open.nc(paste(dir2,"WRF_d01_slp_regrid.nc",sep=""))
lat=var.get.nc(a,"lat0")
lon=var.get.nc(a,"lon0")
slp=var.get.nc(a,"slp0")
time=var.get.nc(a,"Times")
time2=as.Date(as.character(time),"%Y%m%d%H")

backtrack=array(0,c(length(firstfixes[,1]),DS+1,2))
dimnames(backtrack)[[3]]=c("Lon","Lat")
backtrack[,1,]=as.matrix(round(firstfixes[,6:7]*2)/2)


####Do stuart's lazy backtracking

for(i in 1:length(firstfixes[,1]))
  for(j in 1:DS)
  {
    I=which(dates==firstfixes$Date2[i]-j*0.25)
    J=which(lon==backtrack[i,j,1])
    K=which(lat==backtrack[i,j,2])
    
    n=1
    tmp=array(0,c((DL*2+1)^2,5))
    for(p in (J-DL):(J+DL))
      for(q in (K-DL):(K+DL))
      {    
        if(p>2 & q>2 & p<length(lon)-1 & q<length(lat)-1)
        {
          CELL = slp[p,q,I];
          
          # MEASURE PRESSURE GRADIENTS AROUND TEST CELL
          test = c(slp[p-2,q,I]-CELL, slp[p-2,q+2,I]-CELL, slp[p-2,q-2,I]-CELL, 
                   slp[p+2,q,I]-CELL, slp[p+2,q+2,I]-CELL, slp[p+2,q-2,I]-CELL, 
                   slp[p,q-2,I]-CELL, slp[p,q+2,I]-CELL)
          
          tmp[n,] = c(length(which(test>0)),max(test),lon[p],lat[q],CELL)
        }
        n=n+1
      }  
    a=order(-tmp[,1],tmp[,5])
    backtrack[i,j+1,]=tmp[a[1],3:4]
  }

contour(Useful$x,Useful$y,mask,xlim=c(100,180),ylim=c(-50,-10),drawlabels=F)

for(i in 1:9)
{
  points(backtrack[i,1,1],backtrack[i,1,2],pch=4,cex=2,lwd=2,col=i+1)
  lines(backtrack[i,,1],backtrack[i,,2],lwd=2,col=i+1)
}
legend("bottomleft",legend=1:9,col=2:10,lwd=2,bty="n",ncol=3)


##Test how defining land

contour(Useful$x,Useful$y,mask,drawlabels=F)

y=seq(0,-45,-1.5)
x=rep(NaN,length(y))
I=which(y>(-25) & y<=(-10))
x[I]=151-(25+y[I])
I=which(y<=(-25) & y>=(-31))
x[I]=151
I=which(y<(-31) & y>=(-39))
x[I]=148+(37+y[I])/2
lines(x,y,col="red",lwd=3)

#### Redo for an effective 1.5 degree resolution

rm(list=ls())

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful2$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
setwd("~/Documents/ECLs")

dates=seq.Date(as.Date("2007060100","%Y%m%d%H"),as.Date("2007063018","%Y%m%d%H"),0.25)
cat=c("d01_p100_rad2cv05","d01_p100_rad5cv015","d02_p100_rad2cv05","d02_p100_rad5cv015")
SSTtype=c("ERAI","BRAN","NoEAC","DoubleEAC")

##From stewart:
#Set days to backtrack = 2 (i.e. 8 timesteps)
DS=8
#Set region for tracking - ~ surrounding 6 degrees, so +- 4 1.5 deg cells
DL=4


c=1
r=1
day=27
hour="00"

dir=paste("/home/nfs/z3478332/output/outputUM_wrf_2007-06/ERAI_R",r,"_ensemble/",cat[c],"/",sep="")
dir2=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_ensemble/out/200705",day,hour,"/",sep="")

data=read.csv(paste(dir,"ECLfixes_200705",day,hour,".csv",sep=""))
i=1
I=which(data$ID==i & data$Location==1)
firstfixes=data[I[1],2:10]
for(i in 2:length(unique(data$ID)))
{
  I=which(data$ID==i & data$Location==1)
  firstfixes=rbind(firstfixes,firstfixes=data[I[1],2:10])
}
firstfixes$Date2=as.Date(as.character(firstfixes$Date*100+(as.numeric(firstfixes$Time)-1)*6),"%Y%m%d%H")

backtrack=array(0,c(length(firstfixes[,1]),DS+1,2))
dimnames(backtrack)[[3]]=c("Lon","Lat")
backtrack[,1,]=as.matrix(round(firstfixes[,6:7]/1.5)*1.5)

a=open.nc(paste(dir2,"WRF_d01_slp_regrid.nc",sep=""))
lat=var.get.nc(a,"lat0")
lon=var.get.nc(a,"lon0")
slp=var.get.nc(a,"slp0")
time=var.get.nc(a,"Times")
time2=as.Date(as.character(time),"%Y%m%d%H")

## Regrid as lazy average
lat2=rev(seq(-1.5,-48.5,-1.5))
lon2=rev(seq(178.5,106,-1.5))
slp2=array(NaN,c(length(lon2),length(lat2),length(time)))

for(i in 1:length(lon2))
  for(j in 1:length(lat2))
  {
    I=which(lon==lon2[i])
    J=which(lat==lat2[j])
    slp2[i,j,]=apply(slp[(I-1):(I+1),(J-1):(J+1),],3,mean,na.rm=T)
  }


####Do stuart's lazy backtracking

for(i in 1:length(firstfixes[,1]))
  for(j in 1:DS)
  {
    I=which(dates==firstfixes$Date2[i]-j*0.25)
    J=which(lon2==backtrack[i,j,1])
    K=which(lat2==backtrack[i,j,2])
    
    n=1
    tmp=array(0,c((DL*2+1)^2,5))
    for(p in seq((J-DL),(J+DL)))
      for(q in seq((K-DL),(K+DL)))
      {    
        if(p>1 & q>1 & p<length(lon2) & q<length(lat2))
        {
          CELL = slp2[p,q,I];
          
          # MEASURE PRESSURE GRADIENTS AROUND TEST CELL
          test = c(slp2[p-1,q,I]-CELL, slp2[p-1,q+1,I]-CELL, slp2[p-1,q-1,I]-CELL, 
                   slp2[p+1,q,I]-CELL, slp2[p+1,q+1,I]-CELL, slp2[p+1,q-1,I]-CELL, 
                   slp2[p,q-1,I]-CELL, slp2[p,q+1,I]-CELL)
          
          tmp[n,] = c(length(which(test>0)),max(test),lon2[p],lat2[q],CELL)
        }
        n=n+1
      }  
    a=order(-tmp[,1],tmp[,5])
    backtrack[i,j+1,]=tmp[a[1],3:4]
  }

contour(Useful$x,Useful$y,mask,xlim=c(100,180),ylim=c(-50,-10),drawlabels=F)

for(i in 1:9)
{
  points(backtrack[i,1,1],backtrack[i,1,2],pch=4,cex=2,lwd=2,col=i+1)
  lines(backtrack[i,,1],backtrack[i,,2],lwd=2,col=i+1)
}
legend("bottomleft",legend=1:9,col=2:10,lwd=2,bty="n",ncol=3)

IDs=firstfixes$ID
motion=array(0,c(length(IDs),8,3))
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

typing=cbind(apply(motion[,,1:2]*(motion[,,1:2]>1),c(1,3),sum),apply(motion[,,1:2]*(motion[,,1:2]<(-1)),c(1,3),sum)*-1,apply(motion[,,3],1,sum)/DS,matrix(0,length(IDs),2))
colnames(typing)=c("E","N","W","S","Land prop","N points","S points")
for(i in 1:length(IDs)){
  typing[i,6]=length(which(backtrack[i,2:9,2]>=-27))
  typing[i,7]=length(which(backtrack[i,2:9,2]<=-39))
} 

Type=rep("NA",length(IDs))
I=which(typing[,5]<=0.5 & (typing[,4]>typing[,2] | typing[,7]<2) )
Type[I]="ET"
I=which(typing[,5]<=0.5 & typing[,2]>=typing[,4] & typing[,7]>=2)
Type[I]="SSL"
I=which(typing[,5]>0.5 & typing[,6]>=2)
Type[I]="IT"
I=which(typing[,5]>0.5 & typing[,6]<2)
Type[I]="CL"

############# 2007-2008

rm(list=ls())

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful2$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
setwd("~/Documents/ECLs")
dates=seq.Date(as.Date("2007010100","%Y%m%d%H"),as.Date("2008123118","%Y%m%d%H"),0.25)

##From stewart:
#Set days to backtrack = 2 (i.e. 8 timesteps)
DS=8
#Set region for tracking - ~ surrounding 6 degrees, so +- 4 1.5 deg cells
DL=4


rad=2
proj=100
dom="d01"
nt=""

types=c("ET","IT","SSL","CL")

tlist<-tlistNT<-data.frame(R1=rep("aaa",200),R3=rep("aaa",200),R3=rep("aaa",200),stringsAsFactors=F)
bomb<-bombNT<-matrix(NaN,200,3)

clist=c("red","purple","blue","green")
count=matrix(0,4,3)

nt="_notopo"

for(r in 1:3)
{
  dir="~/output/outputUM_wrf_2007_all/"
  dir2=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",nt,"/out/slp/",sep="")
  
  data=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,nt,"_rad",rad,"_p",proj,".csv",sep="")),
             read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,nt,"_rad",rad,"_p",proj,".csv",sep="")))
  year=floor(data$Date/10000)
  data$ID=data$ID+year*1000
  
  IDs=unique(data$ID)
  I=which(data$ID==IDs[1] & data$Location==1)
  firstfixes=data[I[1],2:10]
  for(i in 2:length(IDs))
  {
    I=which(data$ID==IDs[i] & data$Location==1)
    firstfixes=rbind(firstfixes,firstfixes=data[I[1],2:10])
  }
  firstfixes$Date2=as.Date(as.character(firstfixes$Date*100+(as.numeric(firstfixes$Time)-1)*6),"%Y%m%d%H")
  
  backtrack=array(0,c(length(firstfixes[,1]),DS+1,2))
  dimnames(backtrack)[[3]]=c("Lon","Lat")
  backtrack[,1,]=as.matrix(round(firstfixes[,6:7]/1.5)*1.5)
  
  file=open.nc(paste(dir2,"WRF_d01_0708_regrid_coarse.nc",sep=""))
  lat=var.get.nc(file,"lat")
  lon=var.get.nc(file,"lon")
  time=var.get.nc(file,"Times")
  time2=as.Date(as.character(time),"%Y%m%d%H")
  
  ####Do stuart's lazy backtracking
  
  for(i in 1:length(firstfixes[,1]))
  {
    I=which(dates==firstfixes$Date2[i])
    slp=var.get.nc(file,"slp0",start=c(1,1,I-8),count=c(length(lon),length(lat),8))
    
    for(j in 1:DS)
    {
      J=which(lon==backtrack[i,j,1])
      K=which(lat==backtrack[i,j,2])
      
      n=1
      tmp=array(0,c((DL*2+1)^2,5))
      for(p in seq((J-DL),(J+DL)))
        for(q in seq((K-DL),(K+DL)))
        {    
          if(p>1 & q>1 & p<length(lon) & q<length(lat))
          {
            CELL = slp[p,q,9-j];
            
            # MEASURE PRESSURE GRADIENTS AROUND TEST CELL
            test = c(slp[p-1,q,9-j]-CELL, slp[p-1,q+1,9-j]-CELL, slp[p-1,q-1,9-j]-CELL, 
                     slp[p+1,q,9-j]-CELL, slp[p+1,q+1,9-j]-CELL, slp[p+1,q-1,9-j]-CELL, 
                     slp[p,q-1,9-j]-CELL, slp[p,q+1,9-j]-CELL)
            
            tmp[n,] = c(length(which(test>0)),max(test),lon[p],lat[q],CELL)
          }
          n=n+1
        }  
      a=order(-tmp[,1],tmp[,5])
      backtrack[i,j+1,]=tmp[a[1],3:4]
    }
  }
  
  motion=array(0,c(length(IDs),8,3))
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
  
  typing=cbind(apply(motion[,,1:2]*(motion[,,1:2]>1),c(1,3),sum),apply(motion[,,1:2]*(motion[,,1:2]<(-1)),c(1,3),sum)*-1,apply(motion[,,3],1,sum)/DS,matrix(0,length(IDs),2))
  colnames(typing)=c("E","N","W","S","Land prop","N points","S points")
  for(i in 1:length(IDs)){
    typing[i,6]=length(which(backtrack[i,2:9,2]>=-27))
    typing[i,7]=length(which(backtrack[i,2:9,2]<=-39))
  } 
  
  Type=rep("NA",length(IDs))
  I=which(typing[,5]<=0.5 & (typing[,4]>typing[,2] | typing[,7]<2) )
  Type[I]="ET"
  I=which(typing[,5]<=0.5 & typing[,2]>=typing[,4] & typing[,7]>=2)
  Type[I]="SSL"
  I=which(typing[,5]>0.5 & typing[,6]>=2)
  Type[I]="IT"
  I=which(typing[,5]>0.5 & typing[,6]<2)
  Type[I]="CL"
  
  tlistNT[1:length(Type),r]=Type
  
  ####### Okay, add deepening
  
  data$NDR=NaN
  I=which(data$ID[2:length(data[,1])]==data$ID[1:length(data[,1])-1])+1
  data$NDR[I]= (data$MSLP[I]-data$MSLP[I-1])*sin(60*pi/180)/(6*sin(data$Lat[I]*pi/180))
  
  events=rbind(read.csv(paste(dir,"ECLevents_",dom,"_2007_R",r,"_rad",rad,"_p",proj,".csv",sep="")),
               read.csv(paste(dir,"ECLevents_",dom,"_2008_R",r,"_rad",rad,"_p",proj,".csv",sep="")))
  year=floor(events$Date1/10000)
  events$ID=events$ID+year*1000
  
  events$NDR=NaN
  for(i in 1:length(events[,1]))
  {
    I=which(data$ID==events$ID[i] & data$Location==1 & !is.na(data$NDR))
    if(length(I)>0) events$NDR[i]=max(data$NDR[I],na.rm=T)
  }
  
  I=which(events$NDR>=1)
  bombNT[I,r]=1
}

save(bomb,bombNT,tlist,tlistNT,file="~/output/outputUM_wrf_2007_all/typing.RData")

######## Something for NARCLIM
##
## Only do for my 0.8 CV 22 events subset
## And the Rad2 Proj100 for now




