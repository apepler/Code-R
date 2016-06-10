library(RNetCDF)
f1=open.nc('/srv/ccrc/data34/z3478332/NCEP/Monthly/slp.mon.mean.nc')
time=var.get.nc(f1,'time')
time=as.Date((time/24-2),origin="0001-01-01")
time=data.frame(yy=as.numeric(format(time,"%Y")), mm=as.numeric(format(time,"%m")))

lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
slp=var.get.nc(f1,'slp')

I=which(lat>=-45 & lat<=-10)
lat=lat[I]
J=which(lon>=145 & lon<=150)
slp2=apply(slp[J,I,],c(2,3),mean) ##Profile of mean 145-150 pressure

str=data.frame(Year=time$yy,Month=time$mm,STRI=rep(0,length(time[,1])),STRP=rep(0,length(time[,1])))
for(i in 1:length(time[,1]))
{
  str[i,3]=max(slp2[,i])
  I=which(slp2[,i]==str[i,3])
  str[i,4]=max(lat[I])
}

STR_bt=read.csv('~/Documents/Timeseries/str.csv')

str3=str[(str[,1]>=1979 & str[,1]<=2009),]
time3=time[(str[,1]>=1979 & str[,1]<=2009),]
cor(STR_bt[,3],str3[,3])
cor(STR_bt[,4],str3[,4])
mean(STR_bt[,3]-str3[,3])
aves=data.frame(mm=seq(1,12),I1=rep(0,12),P1=rep(0,12),I2=rep(0,12),P2=rep(0,12))
for(i in 1:12)
{
  I=which(time3[,2]==i)
  aves[i,2:3]=mean(str3[I,3:4])
  I=which(STR_bt[,2]==i)
  aves[i,4:5]=mean(STR_bt[I,3:4])
}
aves[,3]=-aves[,3]
aves=rbind(aves,aves[1,])


###What about a uwnd profile?
library(RNetCDF)
f1=open.nc('/srv/ccrc/data34/z3478332/NCEP/Monthly/uwnd.mon.mean.nc')
time=var.get.nc(f1,'time')
time=as.Date((time/24-2),origin="0001-01-01")
time=data.frame(yy=as.numeric(format(time,"%Y")), mm=as.numeric(format(time,"%m")))

lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(length(lon),length(lat),1,length(time[,1])),unpack=T)

I=which(lat>=-60 & lat<=-10)
lat=lat[I]
J=which(lon==150)
u1=uwnd[J,I,]
J=which(lon>=145 & lon<=150)
u2=apply(uwnd[J,I,],c(2,3),mean) ##Profile of mean 145-150 uwnd
J=which(lon>=145 & lon<=155)
u3=apply(uwnd[J,I,],c(2,3),mean) ##Profile of mean 145-155 uwnd

u1m<-u2m<-u3m<-matrix(0,length(lat),12)
for(i in 1:12)
{
  I=which(time[,2]==i)
  u1m[,i]=apply(u1[,I],1,mean)
  u2m[,i]=apply(u2[,I],1,mean)
  u3m[,i]=apply(u3[,I],1,mean)
}

Ucool<-matrix(0,length(lat),3)
I=which(time[,2]>=5 & time[,2]<=10)
Ucool[,1]=apply(u1[,I],1,mean)
Ucool[,2]=apply(u2[,I],1,mean)
Ucool[,3]=apply(u3[,I],1,mean)

plot(lat,Ucool[,1],col="black",type="l")
lines(lat,Ucool[,2],col="red")
lines(lat,Ucool[,3],col="blue")

dU<-matrix(0,length(lat)-2,3)
for(i in 1:length(dU[,1]))
  for(j in 1:3)
    dU[i,j]=(Ucool[i+2,j]-Ucool[i,j])/5

plot(lat[2:(length(lat)-1)],dU[,1],col="black",type="l")
lines(lat[2:(length(lat)-1)],dU[,2],col="red")
lines(lat[2:(length(lat)-1)],dU[,3],col="blue")

library(RNetCDF)
f1=open.nc('/srv/ccrc/data34/z3478332/NCEP/Monthly/slp.mon.mean.nc')
time2=var.get.nc(f1,'time')
time2=as.Date((time2/24-2),origin="0001-01-01")
time2=data.frame(yy=as.numeric(format(time2,"%Y")), mm=as.numeric(format(time2,"%m")))

lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
slp=var.get.nc(f1,'slp')

I=which(lat>=-60 & lat<=-10)
lat=lat[I]
J=which(lon>=145 & lon<=150)
slp2=apply(slp[J,I,],c(2,3),mean) ##Profile of mean 145-150 pressure
I=which(time2[,2]>=5 & time2[,2]<=10)
Pcool=apply(slp2[,I],1,mean)

dUa<-dU2a<-matrix(0,length(lat)-1,3)
for(i in 1:length(dUa[,1]))
  for(j in 1:3)
    dUa[i,j]=(Ucool[i+1,j]-Ucool[i,j])/2.5

lat1=seq(-11.25,-58.75,-2.5)
for(i in 1:3)
  dU2a[,i]=(dUa[,i]-mean(dUa[,i]))/sd(dUa[,i])


dU2<-matrix(0,length(lat)-2,3)
for(i in 1:3)
  dU2[,i]=(dU[,i]-mean(dU[,i]))/sd(dU[,i])
P2=(Pcool-mean(Pcool))/sd(Pcool)

plot(lat[2:(length(lat)-1)],dU2[,1],col="black",type="l")
lines(lat[2:(length(lat)-1)],dU2[,2],col="red")
lines(lat[2:(length(lat)-1)],dU2[,3],col="blue")
lines(lat1,dU2a[,1],col="black",lty=2)
lines(lat1,dU2a[,2],col="red",lty=2)
lines(lat1,dU2a[,3],col="blue",lty=2)
lines(lat,P2,col="purple")

lat[Pcool==max(Pcool)]

I=which(lat>=-45 & lat<=-20)
for(i in 1:3) {
  J=which(dU[I,i]==min(dU[I,i]))
  print(lat[I[J]])
}
  

##Okay, over time, let's get STR position, STR intensity, and position of minimum Uwind gradient
lat1=seq(-11.25,-58.75,-2.5)
du<-matrix(0,length(lat)-1,length(time[,1]))
for(i in 1:length(du[,1])) du[i,]=(u3[i+1,]-u3[i,])/2.5

I=which(lat>=-45 & lat<=-20)
I1=which(lat1>=-45 & lat1<=-20)
str=data.frame(Year=time$yy,Month=time$mm,STRI=rep(0,length(time[,1])),STRP=rep(0,length(time[,1])),duP=rep(0,length(time[,1])))
for(i in 1:length(time[,1]))
{
  str[i,3]=max(slp2[I,i])
  J=which(slp2[I,i]==str[i,3])
  str[i,4]=max(lat[I[J]])
  J=which(du[I1,i]==min(du[I1,i]))
  str[i,5]=max(lat1[I1[J]])
}


##Comp with other uwnd stuff

uwnd=read.csv('~/Documents/Data/CMIP5/uwnd.csv')
uwnd$STRPncep<-uwnd$STRIncep<-rep(0,60)
for(i in 1:length(uwnd[,1]))
{
  I=which(str[,1]==uwnd[i,1] & str[,2]>=5 & str[,2]<=10)
  uwnd[i,6:7]=apply(str[I,3:4],2,mean)
}

lmI=lm(uwnd[,5]~uwnd[,3])
lmP=lm(uwnd[,5]~uwnd[,4])
lmIP=lm(uwnd[,5]~uwnd[,3]+uwnd[,4])

plot(uwnd[,3],uwnd[,4])

summary(lm(uwnd[,5]~uwnd[,3]+uwnd[,4]))



