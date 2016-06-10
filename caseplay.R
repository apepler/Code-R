##This is to play with looking at the GPH & temperature profiles for four case study lows
library("RNetCDF")
setwd('~/Documents/ECLs')
read.csv('~/Documents/ECLs/cases.csv')->cases

f1=open.nc('~/Documents/Data/ncep.hgt.1998.nc')
f2=open.nc('~/Documents/Data/ncep.air.1998.nc')
lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
lev=var.get.nc(f1,'level')
time=var.get.nc(f1,'time')
hh=(time %% 24)
time=as.Date((time/24-2),origin="0001-01-01")
time=as.numeric(format(time,"%Y"))*10000+as.numeric(format(time,"%m"))*100+as.numeric(format(time,"%d"))

air=var.get.nc(f2,'air',unpack=T)
hgt=var.get.nc(f1,'hgt',unpack=T)
a<-h<-matrix(0,17,5)
a[,1]<-h[,1]<-lev
for(i in 1:4)
{
  I=which(lon==cases[i,5])
  J=which(lat==cases[i,6])
  K=which(time==cases[i,1] & hh==cases[i,2])
  a[,i+1]=air[I,J,,K]
  h[,i+1]=hgt[I,J,,K]
}

##Nope, that's overwhelmed by lapserates
for(i in 1:4)
{
  I=which(lon==cases[i,5])
  J=which(lat==cases[i,6])
  K=which(time==cases[i,1] & hh==cases[i,2])
  for(j in 1:17) h[j,i+1]=max(hgt[(I-2):(I+2),(J-2):(J+2),j,K])-min(hgt[(I-2):(I+2),(J-2):(J+2),j,K])
}


##Version 2 - full progression, but do just 925-300
read.csv('~/Documents/ECLs/cases2.csv')->cases
cases$Zdiff=0
for(i in 1:length(cases[,1]))
{
  I=which(lon==cases[i,5])
  J=which(lat==cases[i,6])
  K=which(time==cases[i,1] & hh==cases[i,2])
  for(j in 2:8) aa[j-1]=max(hgt[(I-2):(I+2),(J-2):(J+2),j,K])-min(hgt[(I-2):(I+2),(J-2):(J+2),j,K])
  aa2=lm(aa[2:4]~log(lev[3:5]))
  cases[i,9]=aa2$coefficients[2] ##Lower level, 600-850
  aa2=lm(aa[4:7]~log(lev[5:8]))
  cases[i,10]=aa2$coefficients[2] ##Upper level, 300-600
}

