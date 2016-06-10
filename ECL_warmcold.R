rm(list=ls())
setwd('~/Documents/ECLs')
erai<-read.csv('CSV/ECLfixes_unsw_erai_150.csv')

##Testing - only using lows in 1981 & when in the location
erai=erai[erai$Date>=19810000 & erai$Date<=19820000 & erai$Location==1,]
date1=erai$Date+erai$Time/24
erai$Type<-erai$VL<-erai$VU<-erai$VT<-erai$B<-erai$Bearing<-0

library(RNetCDF)
a<-open.nc('/srv/ccrc/data23/z3478332/ERAI/WarmCold/erai_upper_1981_150.nc')
lat<-var.get.nc(a,'latitude')
lon<-var.get.nc(a,'longitude')
lev<-var.get.nc(a,'level')
time<-var.get.nc(a,'time')
time1=as.Date((time/24),origin="1900-01-01")
date2=as.double(format(time1,"%Y%m%d"))+(time %% 24)/24
logp=log(lev)

library(geosphere)
dirs=data.frame(angle=c(0,45,135,225,315,360),dir=c("N","E","S","W","N",NA))
erai$Lon[erai$Lon>180]=erai$Lon[erai$Lon>180]-360
erai$Bearing="none"

library(sp)
for(k in 1:length(erai[,1]))
{
  dist=matrix(NaN,length(lon),length(lat))
  for(i in 1:length(lon))
    for(j in 1:length(lat))
      dist[i,j]=spDistsN1(cbind(lon[i],lat[j]),cbind(erai$Lon[k],erai$Lat[k]))
  
  loc=which(dist==min(dist),arr.ind=T) ##I.e. the lon/lat location
  I=which(date2==date1[k])             ##I.e. the time
  hgt<-var.get.nc(a,'z',start=c(loc[1]-7,loc[2]-7,1,I),count=c(15,15,14,1),unpack=T)/9.80665
  
  delZ=apply(hgt,3,max,na.rm=T)-apply(hgt,3,min,na.rm=T)
  
  dd=rep(0,12)
  for(i in 1:12) dd[1]=(delZ[i+1]-delZ[i])/(logp[i+1]-log[[i]]) 
  erai$VU[k]=sum(dd[1:6])
  erai$VL[k]=sum(dd[7:12])
  
  ##Bearing of low - in past 6 hours if possible, otherwise next six hours
  if(erai$Fix[k]==1 | k==1)
    erai$Bearing[k]=bearing(erai[k,6:7],erai[k+1,6:7])
  else
    erai$Bearing[k]=bearing(erai[k-1,6:7],erai[k,6:7])
}

erai$VU=(erai$d300-erai$d600)/(log(300)-log(600))
erai$VL=(erai$d600-erai$d900)/(log(600)-log(900))
erai$VT=(erai$d300-erai$d900)/(log(300)-log(900))

erai$Type="Mixed"
erai$Type[erai$VU<(-5) & erai$VL<(-5)]="Cold"
erai$Type[erai$VU>5 & erai$VL>5]="Warm"
length(which(erai$Type=="Cold"))

###########################################################
##Trying to do L/R asymmetry - so need to find the point that it's between
##Also need to find the direction of movement - simplified to be N,S,W,E
rm(list=ls())
erai<-read.csv('CSV/ECLfixes_unsw_erai_150.csv')
erai=erai[erai$Date>=19810000 & erai$Date<=19820000,]
date1=erai$Date+erai$Time/24

library(RNetCDF)
a<-open.nc('/srv/ccrc/data23/z3478332/ERAI/WarmCold/erai_upper_1981_150.nc')
lat<-var.get.nc(a,'latitude')
lon<-var.get.nc(a,'longitude')
lev<-var.get.nc(a,'level')
time<-var.get.nc(a,'time')
time1=as.Date((time/24),origin="1900-01-01")
date2=as.double(format(time1,"%Y%m%d"))+(time %% 24)/24
hgt<-var.get.nc(a,'z',unpack=T)/9.80665

##This version gives a bearing in E/N/S/W between the point at k and its subsequent time
library(geosphere)
dirs=data.frame(angle=c(0,45,135,225,315,360),dir=c("N","E","S","W","N",NA))
erai$Lon[erai$Lon>180]=erai$Lon[erai$Lon>180]-360
erai$Bearing="none"
for(k in 1:(length(erai[,1])-1))
{
  if(erai$ID[k+1]==erai$ID[k])
     erai$Bearing[k]=as.character(dirs[findInterval(bearing(erai[k,6:7],erai[k+1,6:7]),dirs[,1]),2])
  else
    erai$Bearing[k]=erai$Bearing[k-1]
}
erai$Bearing[k+1]=erai$Bearing[k]

# ##Actually, I'd rather do a bearing that's between point k-1 & k+1 if possible
# ##Although, Hart does seem to use the direction up to current time.
# erai$Bearing2="none"
# erai$Bearing2[1]=as.character(dirs[findInterval(bearing(erai[1,6:7],erai[2,6:7]),dirs[,1]),2])
# for(k in 2:(length(erai[,1])-1))
# {
#   if(erai$ID[k+1]==erai$ID[k-1]) {
#     erai$Bearing2[k]=as.character(dirs[findInterval(bearing(erai[k-1,6:7],erai[k+1,6:7]),dirs[,1]),2]) 
#   } else if(erai$ID[k-1]==erai$ID[k]) {
#     erai$Bearing2[k]=as.character(dirs[findInterval(bearing(erai[k-1,6:7],erai[k,6:7]),dirs[,1]),2])
#   } else {
#     erai$Bearing2[k]=as.character(dirs[findInterval(bearing(erai[k,6:7],erai[k+1,6:7]),dirs[,1]),2])
#   }
# }
# erai$Bearing2[length(erai[,1])]=erai$Bearing[k]=as.character(dirs[findInterval(bearing(erai[length(erai[,1])-1,6:7],erai[length(erai[,1]),6:7]),dirs[,1]),2])

erai=erai[erai$Location==1,]
erai$B<-erai$Type<-erai$VL<-erai$VU<-erai$d900<-erai$d600<-erai$d300<-0

for(k in 1:(length(erai[,1])-1))    
{
  I=which(date2==date1[k])     
  lonL=findInterval(erai$Lon[k],lon)
  latL=29-findInterval(erai$Lat[k],lat[29:1])
  erai[k,12:14]=apply(hgt[(lonL-1):(lonL+2),(latL-1):(latL+2),,I],3,max)-apply(hgt[(lonL-1):(lonL+2),(latL-1):(latL+2),,I],3,min)  

  ##Now, need to do B - depends on location
  ##If direction is north, want west-east
  ##If direction is east, want north-south
  
  switch(erai$Bearing[k],
         N={
           left=mean(hgt[(lonL-1):(lonL),(latL-1):(latL+2),2,I]-hgt[(lonL-1):(lonL),(latL-1):(latL+2),3,I])
           right=mean(hgt[(lonL+1):(lonL+2),(latL-1):(latL+2),2,I]-hgt[(lonL+1):(lonL+2),(latL-1):(latL+2),3,I])
           erai[k,18]=left-right
         },
         E={
           left=mean(hgt[(lonL-1):(lonL+2),(latL-1):(latL),2,I]-hgt[(lonL-1):(lonL+2),(latL-1):(latL),3,I])
           right=mean(hgt[(lonL-1):(lonL+2),(latL+1):(latL+2),2,I]-hgt[(lonL-1):(lonL+2),(latL+1):(latL+2),3,I])
           erai[k,18]=left-right
         },
         S={
           left=mean(hgt[(lonL+1):(lonL+2),(latL-1):(latL+2),2,I]-hgt[(lonL+1):(lonL+2),(latL-1):(latL+2),3,I])
           right=mean(hgt[(lonL-1):(lonL),(latL-1):(latL+2),2,I]-hgt[(lonL-1):(lonL),(latL-1):(latL+2),3,I])
           erai[k,18]=left-right
         },
         W={
            left=mean(hgt[(lonL-1):(lonL+2),(latL+1):(latL+2),2,I]-hgt[(lonL-1):(lonL+2),(latL+1):(latL+2),3,I])
            right=mean(hgt[(lonL-1):(lonL+2),(latL-1):(latL),2,I]-hgt[(lonL-1):(lonL+2),(latL-1):(latL),3,I])
            erai[k,18]=left-right
        },
         stop("Enter something that switches me!"))
}

erai$VU=(erai$d300-erai$d600)/(log(300)-log(600))
erai$VL=(erai$d600-erai$d900)/(log(600)-log(900))
erai$Type="Mixed"
erai$Type[erai$VU<(-5) & erai$VL<(-5)]="Cold"
erai$Type[erai$VU>5 & erai$VL>5]="Warm"







