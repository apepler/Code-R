####### Okay, taking bits from existing code for ECL typing.
### First, we need to load in ERAI and restrict to just one year where a range of MLDB ECLs
library(R.matlab)
library(RNetCDF)
library(abind)
library(sp)

setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/")

erai=read.csv("ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")
I=which(erai$Date>=19980000 & erai$Date<=19990000)
erai=erai[I,]

erai_E=read.csv("ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
I=which(erai_E$Date2>=19980000 & erai_E$Date1<=19990000)
erai_E=erai_E[I,]

readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
types<-c("ET","SSL","IT","CL","XTC")
stuart_events=read.csv("~/Documents/ECLs/StuartBrowning/Events_SB.csv",stringsAsFactors=F)
stuart=read.csv("~/Documents/ECLs/StuartBrowning/Fixes_SB.csv")
stuart$Date2=as.POSIXct(as.character(stuart$Date*100+stuart$Time),format="%Y%m%d%H",tz="GMT")

mldb=read.csv("~/Documents/ECLs/Algorithm Comparison/MLDB.csv",stringsAsFactors=F)
mldb$Date2=as.POSIXct(as.character(mldb$Date*100),format="%Y%m%d%H",tz="GMT")
mldb_events=read.csv("~/Documents/ECLs/mldb_events.csv",stringsAsFactors=F)

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

erai$Date2=as.POSIXct(as.character(erai$Date*100+(as.numeric(erai$Time)-1)*6),format="%Y%m%d%H",tz="GMT")
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

### Nest set: Browning classification

#Set days to backtrack = 2 (i.e. 8 timesteps)
DS=8
#Set region for tracking - ~ surrounding 6 degrees, so +- 4 1.5 deg cells
DL=4

IDs=unique(erai$ID)

i=1
I=which(erai$ID==IDs[i] & erai$Location==1)
firstfixes=erai[I[1],]
for(i in 2:length(IDs))
{
  I=which(erai$ID==IDs[i] & erai$Location==1)
  firstfixes=rbind(firstfixes,erai[I[1],])
}
backtrack=array(0,c(length(firstfixes[,1]),DS+1,2))
dimnames(backtrack)[[3]]=c("Lon","Lat")
backtrack[,1,]=as.matrix(round(firstfixes[,7:8]/1.5)*1.5)

### Okay, now set up to loop through ALL

file=open.nc(paste("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_1994-01_1998-12.nc",sep=""))
lat=var.get.nc(file,"latitude")
lon=var.get.nc(file,"longitude")
dates=seq.POSIXt(as.POSIXct("1994010100",format="%Y%m%d%H",tz="GMT"),as.POSIXct("1998123118",format="%Y%m%d%H",tz="GMT"),21600)

####Do stuart's lazy backtracking

for(i in 1:length(firstfixes[,1]))
{
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

erai_E$TypeSB=Type

contour(Useful$x,Useful$y,mask,xlim=c(100,180),ylim=c(-50,-10),drawlabels=F)

for(t in 1:4)
{
  I=which(Type==types[t])
  for(j in 1:length(I))
  {
  points(backtrack[I[j],1,1],backtrack[I[j],1,2],pch=4,cex=2,lwd=2,col=clist[t])
  lines(backtrack[I[j],,1],backtrack[I[j],,2],lwd=2,col=clist[t])
  }
}
legend("topright",types,col=clist,lwd=2)

### Next: Hart phase space

erai$VL<-erai$VU<-erai$B<-erai$Bearing<-0

a<-open.nc('/srv/ccrc/data34/z3478332/ERAI/WarmCold/ERAI_GPH_1998_50hPa.nc')
lat<-var.get.nc(a,'latitude')
lon<-var.get.nc(a,'longitude')
lev<-var.get.nc(a,'level')
time<-var.get.nc(a,'time')
time1=as.POSIXct(time*60*60,origin="1900-01-01-00",tz="GMT")
logp=log(lev)

library(geosphere)
erai$Lon[erai$Lon>180]=erai$Lon[erai$Lon>180]-360
erai$Bearing=NaN
library(sp)
lon[lon>180]=lon[lon>180]-360

for(k in 1:length(erai[,1]))
{
  dist=matrix(NaN,length(lon),length(lat))
  for(i in 1:length(lon))
    for(j in 1:length(lat))
      dist[i,j]=spDistsN1(cbind(lon[i],lat[j]),cbind(erai$Lon[k],erai$Lat[k]))
    
    loc=which(dist==min(dist),arr.ind=T) ##I.e. the lon/lat location
    I=which(time1==erai$Date2[k])             ##I.e. the time
    
    if(loc[1]>7 & loc[2]>7 & loc[1]<=(length(lon)-7) & loc[2]<=(length(lat)-7) & min(dist)<=100)
    {
    lon2=lon[(loc[1]-7):(loc[1]+7)]
    lat2=lat[(loc[2]-7):(loc[2]+7)]
    hgt<-var.get.nc(a,'z',start=c(loc[1]-7,loc[2]-7,1,I),count=c(15,15,13,1),unpack=T)/9.80665
    
    delZ=apply(hgt,3,max,na.rm=T)-apply(hgt,3,min,na.rm=T)
    
    a=lm(delZ[1:7]~logp[1:7])
    erai_E$VU[k]=a$coefficients[2]
    a=lm(delZ[7:13]~logp[7:13])
    erai_E$VL[k]=a$coefficients[2]
    
    ##Bearing of low - in past 6 hours if possible, otherwise next six hours
    if(erai$Fix[k]==1 | k==1)
      erai$Bearing[k]=bearingRhumb(erai[k,7:8],erai[k+1,7:8]) else
        erai$Bearing[k]=bearingRhumb(erai[k-1,7:8],erai[k,7:8])
    
    bearing2=matrix(0,15,15)
    for(i in 1:15)
      for(j in 1:15)
        bearing2[i,j]=bearingRhumb(erai[k,7:8],c(lon2[i],lat2[j]))
    
    L=which((bearing2<erai$Bearing[k] & bearing2>=erai$Bearing[k]-180) | (bearing2<erai$Bearing[k]+360 & bearing2>=erai$Bearing[k]+180))
    R=which((bearing2>=erai$Bearing[k]-360 & bearing2<erai$Bearing[k]-180) | (bearing2<erai$Bearing[k]+180 & bearing2>=erai$Bearing[k]))
     
    Thick=hgt[,,which(lev==600)]-hgt[,,which(lev==900)]
    erai$B[k]=mean(Thick[L],na.rm=T)-mean(Thick[R],na.rm=T)
    }
}

erai$Type="Mixed"
erai$Type[erai$VL<(-50)]="EC"
erai$Type[erai$VU>-10 & erai$VL>0]="TC"
erai$Type[erai$VU<(-10) & erai$VL>0]="SC"

########## Ran on storm, now we can play

erai$Date2=as.POSIXct(paste(as.character(erai$Date),substr(erai$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

I=which(erai$Location==1)
plot(erai$VL[I],erai$VU[I])

erai_E$HartType="Unsure"
for(k in 1:length(erai_E[,1]))
{
  I=which(erai$ID==erai_E$ID[k] & erai$Location==1)
  b=unique(erai$Type[I])
  
  if(length(b)==1) erai_E$HartType[k]=b else {
    if(length(b)==2 & "Mixed" %in% b) erai_E$HartType[k]=b[which(b!="Mixed")] else
    {
      count=rep(0,length(b))
      for(j in 1:length(b)) count[j]=length(which(erai$Type[I]==b[j]))
      J=which(count==max(count))
      if(length(J)==1) erai_E$HartType[k]=b[J] else erai_E$HartType[k]="Mixed"
    }
  }
}

erai_E$B<-0
for(k in 1:length(erai_E[,1]))
{
  I=which(erai$ID==erai_E$ID[k] & erai$Location==1)
  erai_E$B[k]=mean(erai$B[I])
}
### Now, we need to see if there's an MLDB or SB event that matches the event for day/location (w/in 1000 km) and 
### Check what type it was



erai_E$BGtype<-erai_E$Mtype<-"None"
for(i in 1:length(erai_E[,1]))
{
  I=which(erai$ID==erai_E$ID[i] & erai$Location==1)
  tmp=erai[I,]
  
  for(j in 1:length(tmp[,1]))
  {
    I=which(mldb$Date2>=(tmp$Date2[j]-86400) & mldb$Date2<=(tmp$Date2[j]+86400))
    if(length(I)>0)
    {
      dist=spDistsN1(as.matrix(mldb[I,6:7]),as.numeric(tmp[j,7:8]),longlat=T)
      if(min(dist)<500) 
      {
        K=which(dist==min(dist))
        k=which(mldb_events$evntnum==mldb$ID[I[K]])
        erai_E$Mtype[i]=mldb_events$category[k] 
        } 
    }
    }
}

for(i in 1:length(erai_E[,1]))
{
  I=which(erai$ID==erai_E$ID[i] & erai$Location==1)
  tmp=erai[I,]
  
  for(j in 1:length(tmp[,1]))
  {
    I=which(stuart$Date2>=(tmp$Date2[j]-43200) & stuart$Date2<=(tmp$Date2[j]+43200))
    if(length(I)>0)
    {
      dist=spDistsN1(as.matrix(stuart[I,8:7]),as.numeric(tmp[j,7:8]),longlat=T)
      if(min(dist)<500) 
      {
        K=which(dist==min(dist))
        k=which(stuart_events$ID==stuart$ID[I[K[1]]])
        erai_E$BGtype[i]=types[stuart_events$Type[k]] 
      } 
    }
  }
}
  
write.csv(erai,file="ECLfixes_umelb_erai_150_topo_rad2_proj100_typing.csv",stringsAsFactors=F)
write.csv(erai_E,file="ECLevents_umelb_erai_150_topo_rad2_proj100_typing.csv")
