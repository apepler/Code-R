rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

wrfv=c("R1","R2","R3")
monthrain<-monthrainNT<-matrix(0,24,3)
colnames(monthrain)<-colnames(monthrainNT)<-wrfv

# for(i in 1:3)
# {
#   n=1
#   for(year in c(2007,2008))
#     for(month in 1:12)
#     {
#       dir1=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_",wrfv[i],"_nudging_default_2007/out/",sep="")
#       dir2=paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_",wrfv[i],"_nudging_default_2007_notopo/out/",sep="")
#       fname=paste("wrfhrly_d02_",year,"-",sprintf("%2.2d",month),"-01_00:00:00",sep="")
#       
#       a=open.nc(paste(dir1,fname,sep=""))
#       latW=var.get.nc(a,"XLAT",c(1,1,1),c(325,200,1))
#       lonW=var.get.nc(a,"XLONG",c(1,1,1),c(325,200,1))
#       rain1=apply(var.get.nc(a,"PREC_ACC_C"),c(1,2),sum)+apply(var.get.nc(a,"PREC_ACC_NC"),c(1,2),sum)
#       a=open.nc(paste(dir2,fname,sep=""))
#       rain2=apply(var.get.nc(a,"PREC_ACC_C"),c(1,2),sum)+apply(var.get.nc(a,"PREC_ACC_NC"),c(1,2),sum)
#          
#       lont=as.vector(lonW)
#       latt=as.vector(latW)
#       rr=as.vector(rain1)
#       rr[which(is.na(rr))]=0
#       rr2=interp(lont,latt,rr,Useful$x,Useful$y)
#       monthrain[n,i]=mean(rr2$z*t(Useful$mask),na.rm=T)
#       rr=as.vector(rain2)
#       rr[which(is.na(rr))]=0
#       rr2=interp(lont,latt,rr,Useful$x,Useful$y)
#       monthrainNT[n,i]=mean(rr2$z*t(Useful$mask),na.rm=T)
#       n=n+1
#     }
# }
# save(monthrain,monthrainNT,file="monthrain_d02_0708.RData")

load("monthrain_0708.RData")
diff=monthrainNT-monthrain
a=density(diff[,1],na.rm=T,from=-300,to=100)
dens=data.frame(X=a$x,R1=a$y,R2=rep(0,length(a$y)),R3=rep(0,length(a$y)))
for(i in 2:3){
  a=density(diff[,i],na.rm=T,from=-300,to=100)
  dens[,i+1]=a$y
} 

plot(NA,xlim=range(dens[,1]),ylim=range(dens[,2:4]),xlab="Monthly rainfall difference",
     ylab="Frequency",main="ESB rainfall differences in 2007-2008",cex.main=1.2)
for(i in 1:3) lines(dens[,1],dens[,i+1],col=i,lwd=2)

load("monthrain_d02_0708.RData")
diff=monthrainNT-monthrain
a=density(diff[,1],na.rm=T,from=-300,to=100)
dens=data.frame(X=a$x,R1=a$y,R2=rep(0,length(a$y)),R3=rep(0,length(a$y)))
for(i in 2:3){
  a=density(diff[,i],na.rm=T,from=-300,to=100)
  dens[,i+1]=a$y
} 
for(i in 1:3) lines(dens[,1],dens[,i+1],col=i,lwd=2,lty=2)
legend("topleft",legend=c(wrfv,"d01","d02"),col=c(1:3,1,1),lwd=2,lty=c(1,1,1,2,2),ncol=2)

rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

wrfv=c("R1","R2","R3")
dates=seq.Date(as.Date("2007-01-01"),as.Date("2008-12-31"),by=1)
load("dayrain_0708.RData")

diff=dayrainNT[,2:4]-dayrain[,2:4]
a=density(diff[,1],na.rm=T,from=-300,to=100)
dens=data.frame(X=a$x,R1=a$y,R2=rep(0,length(a$y)),R3=rep(0,length(a$y)))
for(i in 2:3){
  a=density(diff[,i],na.rm=T,from=-300,to=100)
  dens[,i+1]=a$y
} 

plot(NA,xlim=range(dens[,1]),ylim=range(dens[,2:4]),xlab="Monthly rainfall difference",
     ylab="Frequency",main="ESB rainfall differences in 2007-2008",cex.main=1.2)
for(i in 1:3) lines(dens[,1],dens[,i+1],col=i,lwd=2)

a=open.nc("/srv/ccrc/data13/z3393020/Analyses/share/geo_em_files/geo_em.d01.narclim.nc")
latW=var.get.nc(a,"XLAT_M",c(1,1,1),c(215,144,1))
lonW=var.get.nc(a,"XLONG_M",c(1,1,1),c(215,144,1))
lont=as.vector(lonW)
latt=as.vector(latW)
library(abind)
mask2=abind(t(Useful$mask),t(Useful$mask),t(Useful$mask),along=3)
range(anndiff2*mask2,na.rm=T)

anngrid<-anngridNT<-array(0,c(215,144,3))
anngrid2<-anngridNT2<-array(0,c(886,691,3))
for(i in 1:3)
{
  load(paste("Dailyrain_R",i,"_0708.RData",sep=""))
  anngridNT[,,i]=apply(daygridNT,c(1,2),sum)/2
  anngrid[,,i]=apply(daygrid,c(1,2),sum)/2
  
  rr=as.vector(anngrid[,,i])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  anngrid2[,,i]=rr2$z*t(Useful$mask)
  rr=as.vector(anngridNT[,,i])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  anngridNT2[,,i]=rr2$z*t(Useful$mask)
  
}
anndiff=anngridNT2-anngrid2
anndiff2=((anngridNT2/anngrid2)-1)*100

a=density(anndiff[,,1],na.rm=T,from=-800,to=400)
dens=data.frame(X=a$x,R1=a$y,R2=rep(0,length(a$y)),R3=rep(0,length(a$y)))
for(i in 2:3){
  a=density(anndiff[,,i],na.rm=T,from=-800,to=400)
  dens[,i+1]=a$y
} 
plot(NA,xlim=range(dens[,1]),ylim=range(dens[,2:4]),xlab="Annual rainfall difference, by grid point",
     ylab="Frequency",main="ESB rainfall differences in 2007-2008",cex.main=1.2)
for(i in 1:3) lines(dens[,1],dens[,i+1],col=i,lwd=2)
legend("topleft",legend=wrfv,col=1:3,lwd=2)


######## Winds

tmp=read.csv("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007/out/uwnd_gdi_3hrly.csv")
wind=cbind(tmp[tmp$Hour==0,c(1,2,3,5)],matrix(0,length(which(tmp$Hour==0)),2))
colnames(wind)=c("Year","Month","Day","R1","R2","R3")
windNT=wind

for(i in 1:3)
{
  tmp=read.csv(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",i,"_nudging_default_2007/out/uwnd_gdi_3hrly.csv",sep=""))
  wind[,3+i]=tmp[tmp$Hour==0,5]
  tmp=read.csv(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",i,"_nudging_default_2007_notopo/out/uwnd_gdi_3hrly.csv",sep=""))
  windNT[,3+i]=tmp[tmp$Hour==0,5]
}

year=c(2007,2008)
mwind=matrix(0,24,8)
colnames(mwind)=c("Year","Month","R1","R2","R3","R1nt","R2nt","R3nt")
n=1
for(y in 1:2)
  for(m in 1:12)
  {
    I=which(wind$Year==year[y] & wind$Month==m)
    mwind[n,1]=year[y]
    mwind[n,2]=m
    mwind[n,3:5]=apply(wind[I,4:6],2,mean)
    mwind[n,6:8]=apply(windNT[I,4:6],2,mean)
    n=n+1
  }

clist=c(1:3,1:3)
tlist=c(1,1,1,2,2,2)

plot(NA,xlim=c(1,24),ylim=c(-8,8),xlab="Monthly mean U-wind",
     ylab="Wind",main="ESB monthly mean 00UTC windflow",cex.main=1.2)
for(i in 1:6) lines(1:24,mwind[,i+2],col=clist[i],lwd=2,lty=tlist[i])
legend("topleft",legend=c(wrfv,"d01","d02"),col=c(1:3,1,1),lwd=2,lty=c(1,1,1,2,2),ncol=2)
abline(h=0,col="grey")

load("monthrain_0708.RData")
mrain=cbind(mwind[,1:2],monthrain,monthrainNT)

for(i in 3:5) print(cor(mwind[c(1:12,14:24),i],mrain[c(1:12,14:24),i+3]-mrain[c(1:12,14:24),i]))

### Daily wind vs. raindiff
load("dayrain_d02_0708.RData")
daydiff_d02=dayrainNT[,2:4]-dayrain[,2:4]
daydiff2_d02=(dayrainNT[,2:4]/dayrain[,2:4]-1)*100
dayrain_d02=dayrain
dayrainNT_d02=dayrainNT
load("dayrain_0708_fix.RData")
daydiff=dayrainNT[,2:4]-dayrain[,2:4]
daydiff2=(dayrainNT[,2:4]/dayrain[,2:4]-1)*100



t=5
I=which(dayrain[,2]>=t | dayrainNT[,2]>=t)
a=density(daydiff[I,1],na.rm=T,from=-50,to=50)
dens=data.frame(X=a$x,R1=a$y,R2=rep(0,length(a$y)),R3=rep(0,length(a$y)))
for(i in 2:3){
  I=which(dayrain[,i+1]>=t | dayrainNT[,i+1]>=t)
  a=density(daydiff[I,i],na.rm=T,from=-50,to=50)
  dens[,i+1]=a$y
} 

plot(NA,xlim=c(-20,20),ylim=range(dens[,2:4]),xlab="Daily rainfall difference",
     ylab="Frequency",main="ESB rainfall differences in 2007-2008",cex.main=1.2)
for(i in 1:3) lines(dens[,1],dens[,i+1],col=i,lwd=2)
legend("topleft",legend=wrfv,col=1:3,lwd=2)

## Okay, W v E days v avg
## What if thresh is 2?

stats=matrix(0,5,3)
for(i in 1:3)
{
  stats[1,i]=length(which(wind[,3+i]>5))
  stats[2,i]=length(which(wind[,3+i]<=5 & wind[,3+i]>2))
  stats[3,i]=length(which(wind[,3+i]<=2 & wind[,3+i]>=-2))
  stats[4,i]=length(which(wind[,3+i]<(-2) & wind[,3+i]>=-5))
  stats[5,i]=length(which(wind[,3+i]<(-5)))
}

## Mean raindiff as fnction of wind & rain thresholds
thresh=c(1,5,10,25)
meandiff<-catcount<-propneg<-proppos<-array(NaN,c(5,3,4))

wmin=c(-Inf,-5,-2,2,5)
wmax=c(-5,-2,2,5,Inf)

for(i in 1:3)
  for(t in 1:4)
    for(w in 1:5)
    {
      I=which(wind[,3+i]>wmin[w] & wind[,3+i]<=wmax[w] & (dayrain[,i+1]>=thresh[t] | dayrainNT[,i+1]>=thresh[t]))
      catcount[w,i,t]=length(I)
      if(length(I)>0)
      {
        meandiff[w,i,t]=mean(daydiff[I,i])
        propneg[w,i,t]=length(which(daydiff[I,i]<(-1)))/length(I)
        proppos[w,i,t]=length(which(daydiff[I,i]>1))/length(I)
      }
    }

meandiff_d02<-array(NaN,c(5,3,4))
for(i in 1:3)
  for(t in 1:4)
    for(w in 1:5)
    {
      I=which(wind[,3+i]>wmin[w] & wind[,3+i]<=wmax[w] & (dayrain_d02[,i+1]>=thresh[t] | dayrainNT_d02[,i+1]>=thresh[t]))
      catcount[w,i,t]=length(I)
      if(length(I)>0) meandiff_d02[w,i,t]=mean(daydiff_d02[I,i])
    }

plot(NA,xlim=c(0.5,5.5),ylim=c(-6,4),xlab="Mean coastal zonal wind (m/s)",
     ylab="Mean rainfall difference",main="ESB rainfall differences where R>5 vs. zonal wind",cex.main=1.2,axes=F)
axis(1,at=1:5,labels=c("< -5","< -2","-2 to 2",">2",">5"))
axis(2,at=-6:4) 
abline(h=0,col="gray")
for(i in 1:3) lines(1:5,meandiff[,i,2],col=i+1,lwd=2)
for(i in 1:3) lines(1:5,meandiff_d02[,i,2],col=i+1,lwd=2,lty=2)
#legend("topleft",legend=wrfv,col=1:3,lwd=2)
legend("topleft",legend=c(wrfv,"d01","d02"),col=c(2:4,1,1),lwd=2,lty=c(1,1,1,1,2),ncol=2)

### Okay, as a function of mean rainfall in base case

rmin=c(1,2,5,10,25,Inf)
stats=matrix(0,5,3)
for(i in 1:3)
  for(t in 1:5) stats[t,i]=length(which(dayrain[,i+1]>rmin[t] & dayrain[,i+1]<=rmin[t+1]))

meandiff2<-meandiff2_d02<-propneg<-proppos<-meandiff3<-meandiff3_d02<-array(NaN,c(5,3))
for(i in 1:3)
  for(t in 1:5)
    {
      I=which(dayrain[,i+1]>rmin[t] & dayrain[,i+1]<=rmin[t+1])
      if(length(I)>0)
      {
        meandiff2[t,i]=mean(daydiff[I,i])
        meandiff3[t,i]=mean(daydiff2[I,i])
        propneg[t,i]=length(which(daydiff[I,i]<(-1)))/length(I)
        proppos[t,i]=length(which(daydiff[I,i]>1))/length(I)
      }
      I=which(dayrain_d02[,i+1]>rmin[t] & dayrain_d02[,i+1]<=rmin[t+1])
      if(length(I)>0){
        meandiff2_d02[t,i]=mean(daydiff_d02[I,i])
        meandiff3_d02[t,i]=mean(daydiff2_d02[I,i])
      } 
    }

plot(NA,xlim=c(0.5,5.5),ylim=c(-20,2),xlab="Mean rainfall in control (mm)",
     ylab="Mean rainfall difference",main="ESB rainfall differences vs rainfall",cex.main=1.2,axes=F)
axis(1,at=1:5,labels=c("1 - 2","2 - 5","5 - 10","10 - 25","> 25"))
axis(2,at=seq(-20,2,2)) 
abline(h=0,col="gray")
for(i in 1:3) lines(1:5,meandiff2[,i],col=i+1,lwd=2)
for(i in 1:3) lines(1:5,meandiff2_d02[,i],col=i+1,lwd=2,lty=2)
#legend("topleft",legend=wrfv,col=1:3,lwd=2)
legend("bottomleft",legend=c(wrfv,"d01","d02"),col=c(2:4,1,1),lwd=2,lty=c(1,1,1,1,2),ncol=2)

plot(NA,xlim=c(0.5,5.5),ylim=c(-70,30),xlab="Mean rainfall in control (mm)",
     ylab="Mean rainfall difference (%)",main="ESB rainfall differences vs rainfall",cex.main=1.2,axes=F)
axis(1,at=1:5,labels=c("1 - 2","2 - 5","5 - 10","10 - 25","> 25"))
axis(2,at=seq(-70,30,10)) 
abline(h=0,col="gray")
for(i in 1:3) lines(1:5,meandiff3[,i],col=i+1,lwd=2)
for(i in 1:3) lines(1:5,meandiff3_d02[,i],col=i+1,lwd=2,lty=2)
#legend("topleft",legend=wrfv,col=1:3,lwd=2)
legend("bottomleft",legend=c(wrfv,"d01","d02"),col=c(2:4,1,1),lwd=2,lty=c(1,1,1,1,2),ncol=2)


a=cbind(dayrain[,2:4],dayrainNT[,2:4])
I=which(apply(a,1,max)>=25)
I2=sort(unique(c(I,I+1,I-1)))
b=a[I2,]

setwd("~/Documents/ECLs/WRFruns/0708/")
load("dayrain_0708_fix.RData")
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

dates=seq.Date(as.Date("2006-12-30"),as.Date("2009-01-02"),by=1)
AWAPrain=cbind(as.numeric(format.Date(dates,"%Y")),as.numeric(format.Date(dates,"%Y%m%d")),rep(0,length(dates)))
dir="/media/Seagate Expansion Drive/Data/daily rainfall/rainfall_2000-2009/"
for(i in 1:length(dates))
{
  fname<-paste(dir,'/rainfall-',AWAPrain[i,1],'/r',AWAPrain[i,2],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  AWAPrain[i,3]=mean(data*Useful$mask,na.rm=T)
}

save(AWAPrain,file="AWAPrain.RData")

### PDFs of daily rain

load("AWAPrain.RData")

x=seq(1,60,length.out=512)
Rdens=matrix(0,512,9)
a=density(AWAPrain[,3],from=1,to=60)
Rdens[,1]=a$y
a=density(as.matrix(dayrain[,2:4]),from=1,to=60)
Rdens[,2]=a$y
a=density(as.matrix(dayrainNT[,2:4]),from=1,to=60)
Rdens[,3]=a$y
for(i in 1:3)
{
  a=density(dayrain[,i+1],from=1,to=60)
  Rdens[,3+i]=a$y
  a=density(dayrainNT[,i+1],from=1,to=60)
  Rdens[,6+i]=a$y
}

plot(x,Rdens[,1],log="x",xlim=c(1,60),ylim=range(0,0.3),col="black",lwd=3,type="l",xlab="ESB rainfall",ylab="Frequency")
for(i in 1:3) 
  {
  lines(x,Rdens[,i+3],col=i+1,lwd=3)
  lines(x,Rdens[,i+6],col=i+1,lwd=3,lty=2)
}

for(i in 2:4) print(length(which(dayrainNT[,i]>=5))/length(which(dayrain[,i]>=5)))
for(i in 2:4) print(length(which(dayrainNT_d02[,i]>=5))/length(which(dayrain_d02[,i]>=5)))

#### Days with heavy rain

library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0

a=open.nc("/srv/ccrc/data13/z3393020/Analyses/share/geo_em_files/geo_em.d01.narclim.nc")
latW=var.get.nc(a,"XLAT_M",c(1,1,1),c(215,144,1))
lonW=var.get.nc(a,"XLONG_M",c(1,1,1),c(215,144,1))
lont=as.vector(lonW)
latt=as.vector(latW)

mask2<-mask3<-matrix(NaN,215,144)
for(i in 2:214)
  for(j in 2:143)
  {
    I=which(Useful$x>=mean(lonW[(i-1):i,j]) & Useful$x<=mean(lonW[(i+1):i,j]))
    J=which(Useful$y>=mean(latW[i,j:(j-1)]) & Useful$y<=mean(latW[i,j:(j+1)]))
    mask2[i,j]=mean(mask[I,J])
  }
I=which(mask2>=0.25)
mask3[I]=1

dayrain_max<-dayrainNT_max<-cells25<-cells25nt<-matrix(NaN,731,3)
for(i in 1:3)
{
  load(paste("Dailyrain_R",i,"_0708.RData",sep=""))
  for(j in 1:731)
  {
    dayrain_max[j,i]=max(daygrid[,,j]*mask3,na.rm=T)
    dayrainNT_max[j,i]=max(daygridNT[,,j]*mask3,na.rm=T)
    cells25[j,i]<-length(which(daygrid[,,j]*mask3>=25))
    cells25nt[j,i]<-length(which(daygridNT[,,j]*mask3>=25))
  }
}

thresh=c(1,5,10,25,50,100,250,Inf)
tab<-matrix(0,7,6)
rownames(tab)=thresh[1:7]

for(i in 1:7)
for(j in 1:3)
{
  I=which(dayrain_max[,j]>=thresh[i] & dayrain_max[,j]<thresh[i+1])
  tab[i,j]=length(I)
  I=which(dayrainNT_max[,j]>=thresh[i] & dayrainNT_max[,j]<thresh[i+1])
  tab[i,j+3]=length(I)
}

thresh=c(1,5,11,27,54,Inf)
tab2<-matrix(0,5,6)
rownames(tab2)=thresh[1:5]

for(i in 1:5)
  for(j in 1:3)
  {
    I=which(cells25[,j]>=thresh[i] & cells25[,j]<thresh[i+1])
    tab2[i,j]=length(I)
    I=which(cells25nt[,j]>=thresh[i] & cells25nt[,j]<thresh[i+1])
    tab2[i,j+3]=length(I)
    
  }


# More PDF stuff - days>25 mm by area

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
readMat('~/Documents/Data/mask_escci2.mat')->escci
mask<-t(Useful$mask)
mask[is.na(mask)]=0
mask2<-t(escci$mask.escci)
mask2[is.na(mask2)]=0

days25<-days25nt<-array(0,c(215,144,3))
for(i in 1:3)
{
  load(paste("Dailyrain_R",i,"_0708.RData",sep=""))
  days25[,,i]=apply(daygrid>=25,c(1,2),sum)
  days25nt[,,i]=apply(daygridNT>=25,c(1,2),sum)
}

days25a<-days25nta<-array(0,c(886,691,3))

a=open.nc("/srv/ccrc/data13/z3393020/Analyses/share/geo_em_files/geo_em.d01.narclim.nc")
latW=var.get.nc(a,"XLAT_M",c(1,1,1),c(215,144,1))
lonW=var.get.nc(a,"XLONG_M",c(1,1,1),c(215,144,1))
lont=as.vector(lonW)
latt=as.vector(latW)
library(abind)
days25a<-days25nta<-array(0,c(886,691,3))

for(i in 1:3)
{
  rr=as.vector(days25[,,i])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  days25a[,,i]=rr2$z*t(Useful$mask)
  rr=as.vector(days25nt[,,i])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  days25nta[,,i]=rr2$z*t(Useful$mask)
}

change=100*((days25nta/days25a)-1)

ColorBar <- function(brks,cols,vert=T,subsampleg=1,blab=brks)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = blab[seq(2, length(brks)-1, subsampleg)])
}

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(11)
cm[6]="white"
breaks=c(-1000,-50,-40,-30,-20,-10,10,20,30,40,50,1000)
layout(cbind(1,2,3,4),width=c(1,1,1,0.5))
for(i in 1:3) 
  {
  par(mar=c(0,0,3,0))
  filled.contour3(Useful$x,Useful$y,change[,,i]*t(escci$mask),lev=breaks,col=cm,
                  xlim=c(146,154),ylim=c(-40,-25),axes=F,main=paste("RCM",i),cex.main=2)
  contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F,xlim=c(146,154),ylim=c(-40,-25))
}
ColorBar(breaks,col=cm,blab=paste(breaks,"%",sep="")) 
