####### Okay, taking bits from existing code for ECL typing.
### First, we need to load in ERAI and restrict to just one year where a range of MLDB ECLs
library(R.matlab)
library(RNetCDF)
library(abind)

erai=read.csv("~/output/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")
I=which(erai$Date>=19980000 & erai$Date<=19990000)
erai=erai[I,]

erai_E=read.csv("~/output/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
I=which(erai_E$Date2>=19980000 & erai_E$Date1<=19990000)
erai_E=erai_E[I,]

readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
types<-c("ET","SSL","IT","CL","XTC")
stuart_events=read.csv("~/Documents/ECLs/StuartBrowning/Events_SB.csv")
stuart=read.csv("~/Documents/ECLs/StuartBrowning/Fixes_SB.csv")

mldb=read.csv("~/Documents/ECLs/Algorithm Comparison/MLDB.csv",stringsAsFactors=F)
mldb$Date2=as.POSIXct(as.character(mldb$Date*100+mldb$Time),format="%Y%m%d%H")
mldb_events=read.csv("~/Documents/ECLs/mldb_events.csv")

# Step 1 - Calculate NDR/bomb for each event

erai$NDR=NaN
I=which(erai$ID[2:length(erai[,1])]==erai$ID[1:length(erai[,1])-1])+1
erai$NDR[I]= (erai$MSLP[I]-erai$MSLP[I-1])*sin(60*pi/180)/(6*sin(erai$Lat[I]*pi/180))

erai_E$Bomb=0
for(i in 1:length(erai_E))
{
  I=which(erai$ID==erai_E$ID[i] & erai$Location==1 & !is.na(erai$NDR))
  if(length(I)>0) if(max(erai$NDR[I],na.rm=T)>=1) erai_E$Bomb[i]=1
}

# Step 2: Identify the event as entered or formed
#1: Formed in region 
#2: Entered & didn't intensify
#3: Entered and intensified - look at NDR & cv?
#4: Formed in region, but may need to tie to a different event

erai$Date2=as.POSIXct(as.character(erai$Date*100+(as.numeric(erai$Time)-1)*6),format="%Y%m%d%H")
erai_E$EnteredFormed=0

library(sp)

for(i in 1:31)
{
  tmp=erai[which(erai$ID==erai_E$ID[i]),]
  
  if(tmp$Location[1]==1)
  {
    I=which(erai$ID!=erai_E$ID[i] & (erai$Date2>=tmp$Date2[1]-43200 & erai$Date2<=tmp$Date2[1]+43200))
    if(length(I)>1) {
      dist=spDistsN1(as.matrix(erai[I,7:8]),as.numeric(tmp[1,7:8]),longlat=T)
      if(min(dist)<500) 
        {
        k=which(erai_E$ID==min(erai$ID[I]))
        erai_E$EnteredFormed[i]=erai_E$EnteredFormed[k] } else erai_E$EnteredFormed[i]=1
    }  else erai_E$EnteredFormed[i]=1
  } else {
    J=which(tmp$Location==1)
    c2=max(tmp$CV[J])
    n2=max(tmp$NDR[J])
    
    K=max(1,J[1]-3):(J[1]-1)
    c1=max(tmp$CV[K])
    
    if(n2>=1 | c2>=c1+0.1)  erai_E$EnteredFormed[i]=3 else erai_E$EnteredFormed[i]=2 
  }
}



rm(list=ls())
library(RNetCDF)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

#Set days to backtrack = 2 (i.e. 8 timesteps)
DS=8
#Set region for tracking - ~ surrounding 6 degrees, so +- 4 1.5 deg cells
DL=4

file="~/output/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv"
#year=2007

ylist=c(1979,1984,1989,1994,1999,2005,2011,2014)

data=read.csv(file)
# I=which(data$Date>=year*10000 & data$Date<=(year+1)*10000)
# data=data[I,]
IDs=unique(data$ID)

i=1
I=which(data$ID==IDs[i] & data$Location==1)
firstfixes=data[I[1],2:10]
for(i in 2:length(IDs))
{
  I=which(data$ID==IDs[i] & data$Location==1)
  firstfixes=rbind(firstfixes,firstfixes=data[I[1],2:10])
}
firstfixes$Date2=as.Date(as.character(firstfixes$Date*100+(as.numeric(firstfixes$Time)-1)*6),"%Y%m%d%H")
years=floor(firstfixes$Date/10000)
backtrack=array(0,c(length(firstfixes[,1]),DS+1,2))
dimnames(backtrack)[[3]]=c("Lon","Lat")
backtrack[,1,]=as.matrix(round(firstfixes[,6:7]/1.5)*1.5)

### Okay, now set up to loop through ALL

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

# contour(Useful$x,Useful$y,mask,xlim=c(100,180),ylim=c(-50,-10),drawlabels=F)
# 
# for(i in 1:length(IDs))
# {
#   points(backtrack[i,1,1],backtrack[i,1,2],pch=4,cex=2,lwd=2,col=i+1)
#   lines(backtrack[i,,1],backtrack[i,,2],lwd=2,col=i+1)
# }

##### Resuls are still fairly messy, but let's try to apply Stuart's cats anyway
## Need to know:
## Whether primary movement is from N, S or W
## If (mostly) east or west of GDR

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

types=c("ET","IT","SSL","CL")
clist=c("red","purple","blue","green")

contour(Useful$x,Useful$y,mask,xlim=c(100,180),ylim=c(-50,-10),drawlabels=F)

for(t in 1:4)
{
  I=which(Type==types[t])
  print(paste(types[t],"=",length(I)))  
  points(mean(backtrack[I,1,1]),mean(backtrack[I,1,2]),pch=4,cex=2,lwd=2,col=clist[t])
  lines(apply(backtrack[I,,1],2,mean),apply(backtrack[I,,2],2,mean),lwd=2,col=clist[t])
}
legend("topright",types,col=clist,lwd=2)

meanfreq=matrix(0,12,4)
for(t in 1:4)
{
  I=which(Type==types[t])
  month=floor(firstfixes$Date[I]/100)%%100
  for(m in 1:12) meanfreq[m,t] = length(which(month==m))/length(unique(years))
}
colnames(meanfreq)=types

plot(1:12,rep(NaN,12),xlab="Month",ylab="Average frequency",ylim=range(0,meanfreq))
for(i in 1:4) lines(1:12,meanfreq[,i],lwd=4,col=clist[i])
legend("topright",types,col=clist,lwd=4)

yy=unique(years)
ycount=matrix(0,length(yy),4)
for(i in 1:length(yy))
  for(t in 1:4)
    ycount[i,t]=length(which(years==yy[i] & Type==types[t]))

events=read.csv(file="~/output/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")

mslp=cbind(seq(970,1030,length.out=512),matrix(0,512,4))
cv=cbind(seq(1,5,length.out=512),matrix(0,512,4))
len=cbind(seq(1,16),matrix(0,16,4))

for(t in 1:4)
{
  I=which(Type==types[t])
  a=density(events$CV2[I],from=1,to=5)
  cv[,t+1]=a$y
  a=density(events$MSLP2[I],from=970,to=1030)
  mslp[,t+1]=a$y
  for(l in 1:16)
  {
    J=which(events$Length2[I]==l)
    len[l,t+1]=length(J)/length(I)
  }
}

plot(1,NaN,xlab="Intensity",ylab="Density",ylim=range(0,cv[,2:5]),xlim=range(cv[,1]))
for(i in 1:4) lines(cv[,1],cv[,i+1],lwd=4,col=clist[i])
legend("topright",types,col=clist,lwd=4)
plot(1,NaN,xlab="Central pressure",ylab="Density",ylim=range(0,mslp[,2:5]),xlim=range(mslp[,1]))
for(i in 1:4) lines(mslp[,1],mslp[,i+1],lwd=4,col=clist[i])
legend("topright",types,col=clist,lwd=4)
plot(1,NaN,xlab="Days in ECL region",ylab="Density",ylim=range(0,len[,2:5]),xlim=range(len[,1]/4))
for(i in 1:4) lines(len[,1]/4,len[,i+1],lwd=4,col=clist[i])
legend("topright",types,col=clist,lwd=4)

####### Okay, add deepening

data$NDR=NaN
I=which(data$ID[2:length(data$ID)]==data$ID[1:length(data$ID)-1])+1
data$NDR[I]= (data$MSLP[I]-data$MSLP[I-1])*sin(60*pi/180)/(6*sin(data$Lat[I]*pi/180))

events=read.csv(file="~/output/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
events$NDR=NaN
for(i in 1:length(events[,1]))
{
  I=which(data$ID==events$ID[i] & data$Location==1 & !is.na(data$NDR))
  if(length(I)>0) events$NDR[i]=max(data$NDR[I],na.rm=T)
}

length(which(events$NDR>=1))/length(which(!is.na(events$NDR)))

thresh=5
b=order(events$CV2,decreasing=T)
thresh2=events$CV2[b[30*thresh]]
J=which(events$CV2>=thresh2)
length(which(events$NDR[J]>=1))/length(which(!is.na(events$NDR[J])))

I=which(events$NDR>=1)
types=c("ET","IT","SSL","CL")

for(t in 1:4)
{
  I=which(Type==types[t])
  J=which(events$NDR[I]>=1)
  print(paste(types[t],"=",round(100*(length(J)/length(I))),"%"))
}

############################################################
####   What about using ERAI/Browning to type MLDB?   ######
############################################################

rm(list=ls())
library(RNetCDF)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
setwd("~/Documents/ECLs/StuartBrowning/")

#Set days to backtrack = 2 (i.e. 8 timesteps)
DS=8
#Set region for tracking - ~ surrounding 6 degrees, so +- 4 1.5 deg cells
DL=4
ylist=c(1979,1984,1989,1994,1999,2005,2011,2014)

file="~/Documents/ECLs/Algorithm Comparison/MLDB.csv"
data=read.csv(file,stringsAsFactors=F)
IDs=unique(data$ID)

i=1
I=which(data$ID==IDs[i])
firstfixes=data[I[1],c(1:4,8,6,7)]
for(i in 2:length(IDs))
{
  I=which(data$ID==IDs[i])
  firstfixes=rbind(firstfixes,firstfixes=data[I[1],c(1:4,8,6,7)])
}
firstfixes$Date2=as.Date(as.character(firstfixes$Date*100+firstfixes$Time),"%Y%m%d%H")

years=floor(firstfixes$Date/10000)
backtrack=array(0,c(length(firstfixes[,1]),DS+1,2))
dimnames(backtrack)[[3]]=c("Lon","Lat")
backtrack[,1,]=as.matrix(round(firstfixes[,6:7]/1.5)*1.5)

### Okay, now set up to loop through ALL

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

# contour(Useful$x,Useful$y,mask,xlim=c(100,180),ylim=c(-50,-10),drawlabels=F)
# 
# for(i in 1:length(IDs))
# {
#   points(backtrack[i,1,1],backtrack[i,1,2],pch=4,cex=2,lwd=2,col=i+1)
#   lines(backtrack[i,,1],backtrack[i,,2],lwd=2,col=i+1)
# }

##### Resuls are still fairly messy, but let's try to apply Stuart's cats anyway
## Need to know:
## Whether primary movement is from N, S or W
## If (mostly) east or west of GDR

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

Type2=rep("NA",length(IDs))
for(i in 1:length(IDs))
{
  I=which(data$ID==IDs[i])
  Type2[i]=data$Type[I[1]]
}

types1=unique(Type)
types2=unique(Type2)

comp=matrix(0,length(types1),length(types2))
rownames(comp)=types1
colnames(comp)=types2
for(i in 1:length(types1))
  for(j in 1:length(types2))
  {
    I=which(Type==types1[i] & Type2==types2[j])
    if(length(I)>0) comp[i,j]=length(I)
  }

comp2=cbind(comp[,1]+comp[,4],comp[,3],comp[,2]+comp[,5]+comp[,6])
colnames(comp2)=c("ET/TC","IT","W")
rownames(comp2)=c("ET/TC","SSL","CL","IT")

comp4<-comp3<-comp2
for(i in 1:4) comp3[i,]=comp3[i,]/apply(comp2,2,sum)
for(i in 1:3) comp4[,i]=comp4[,i]/apply(comp2,1,sum)


library(reshape)
library(ggplot2)
test=melt(comp)

ggplot(test, aes(x=X2, fill=X1, y=value)) + geom_bar(stat="identity") +
  theme_bw() + ylab("Number of events") + xlab("Speer classification") + labs(fill="Browning classification")


means.barplot <- qplot(x=group, y=mean, fill=variable,
                       data=means, geom="bar", stat="identity",
                       position="dodge")



