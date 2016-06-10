rm(list=ls())
library(RNetCDF)
library(sp)
setwd("~/Documents/ECLs/Algorithm Comparison/")

fin=c('MLDB.csv','Mine_rad2.csv','UM_rad2_p100.csv','Ale_v15CTL.csv')
fout=c('MLDB_WC.csv','Mine_rad2_WC.csv','UM_rad2_p100_WC.csv','Ale_v15CTL_WC.csv')
i=1

for(i in 1:5)
{
  erai=read.csv(fin[i])
  if(class(erai$Time)=="factor") erai$Time = (as.numeric(erai$Time)-1)*6
  DD=dim(erai)

date1=erai$Date+erai$Time/24
erai$VL<-erai$VU<-erai$d900<-erai$d600<-erai$d300<-0
yy=floor(erai$Date/10000)
years=unique(year)

for(j in 1:length(year))
{
a<-open.nc(paste("/srv/ccrc/data34/z3478332/ERAI/WarmCold/erai_upper_",years[j],".nc",sep=""))
lat<-var.get.nc(a,'latitude')
lon<-var.get.nc(a,'longitude')
lev<-var.get.nc(a,'level')
time<-var.get.nc(a,'time')
time1=as.Date((time/24),origin="1900-01-01")
date2=as.double(format(time1,"%Y%m%d"))+(time %% 24)/24
hgt<-var.get.nc(a,'z',unpack=T)/9.80665

J=which(yy==years[j])
for(k in 1:length(J))
{
  lat2=lat[lat>=erai$Lat[J[k]]-1.5 & lat<=erai$Lat[J[k]]+1.5]
  lon2=lon[lon>=erai$Lon[J[k]]-1.5 & lon<=erai$Lon[J[k]]+1.5]
  
  dist=matrix(NaN,length(lon2),length(lat2))
  for(i in 1:length(lon2))
    for(j in 1:length(lat2))
      dist[i,j]=spDistsN1(cbind(lon2[i],lat2[j]),cbind(erai$Lon[J[k]],erai$Lat[J[k]]))
  loc=which(dist==min(dist),arr.ind=T) ##I.e. the lon/lat location
  loc=c(which(lon==lon2[loc[1]]),which(lat==lat2[loc[2]]))
  loc=cbind((loc[1]-1):(loc[1]+1),(loc[2]-1):(loc[2]+1))
  
  loc[loc[,1]>length(lon),1]=loc[loc[,1]>length(lon),1]-length(lon)
  loc[loc[,2]>length(lat),2]=NaN
  loc[loc<1]=NaN
  
  
  I=which(date2==date1[J[k]])             ##I.e. the time
  erai[J[k],(DD[2]+1):(DD[2]+3)]=apply(hgt[loc[,1],loc[,2],,I],3,max)-apply(hgt[loc[,1],loc[,2],,I],3,min)  
}
}

erai$VU=(erai$d300-erai$d600)/(log(300)-log(600))
erai$VL=(erai$d600-erai$d900)/(log(600)-log(900))
erai$VT=(erai$d300-erai$d900)/(log(300)-log(900))

erai$Type="Mixed"
erai$Type[erai$VU<(-5) & erai$VL<(-5)]="Cold"
erai$Type[erai$VU>5 & erai$VL>5]="Warm"
length(which(erai$Type=="Cold"))
write.csv(erai,file=fout)
}

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







