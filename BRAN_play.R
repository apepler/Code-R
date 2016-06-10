########### Bluelink playtime
rm(list=ls())
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
library(akima)
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
library(RNetCDF)
library(fields)

a=open.nc("/srv/ccrc/data36/z3478332/WRF/DATA/BRAN/bran3p5_sst_2007.nc")
latB=var.get.nc(a,"yt_ocean")
lonB=var.get.nc(a,"xt_ocean")
timeB=var.get.nc(a,"Time")
time=as.Date(timeB,origin="1993-01-01 00:00:00")
month=as.numeric(format(time, "%m"))
I=which(month==6)
sstB=var.get.nc(a,"temp",c(1,1,1,I[1]),c(length(lonB),length(latB),1,length(I)),unpack=T)
a=open.nc("/srv/ccrc/data36/z3478332/WRF/DATA/BRAN/bran3p5_sst_2007_noeac.nc")
sstB1=var.get.nc(a,"sst_NoEac_3",c(1,1,1,I[1]),c(length(lonB),length(latB),1,length(I)),unpack=T)
sstB2=var.get.nc(a,"sst_NoEac_7",c(1,1,1,I[1]),c(length(lonB),length(latB),1,length(I)),unpack=T)
sstB3=var.get.nc(a,"sst_NoEac_14",c(1,1,1,I[1]),c(length(lonB),length(latB),1,length(I)),unpack=T)

pdf(file="bran_sstchange_8June07_NCLa.pdf",width=15,height=15)
layout(cbind(c(1,3),c(2,4)))
image(lonB,latB,sstB[,,8],col=rich.colors(20),xlim=c(140,170),ylim=c(-50,-20),zlim=c(7,27),main="Default")
image(lonB,latB,sstB1[,,8],col=rich.colors(20),xlim=c(140,170),ylim=c(-50,-20),zlim=c(7,27),main="Remove 7-day running mean")
image(lonB,latB,sstB2[,,8],col=rich.colors(20),xlim=c(140,170),ylim=c(-50,-20),zlim=c(7,27),main="Remove 15-day running mean")
image(lonB,latB,sstB3[,,8],col=rich.colors(20),xlim=c(140,170),ylim=c(-50,-20),zlim=c(7,27),main="Remove 29-day running mean")
dev.off()




I=which(lonB>=160 & lonB<=165)
sstA=apply(sstB[I,,],c(2,3),mean,na.rm=T)
sstB3<-sstB2<-array(0,dim(sstB))
for(i in 1:length(lonB)) sstB2[i,,]=sstA
diff=apply(sstB-sstB2,c(1,2),mean)
for(i in 1:dim(sstB)[3]) sstB3[,,i]=sstB[,,i]-diff

I=which(lonB>=140 & lonB<=160)
J=which(latB>=-50 & latB<=-20)
sstB4<-sstB5<-sstB6<-sstB7<-sstB
sstB4a<-sstB5a<-sstB6a<-sstB[I,J,]
sstB2a=sstB2[I,J,]
sstB3a=sstB3[I,J,]
K=which((sstB4a-sstB2a)>0)
sstB4a[K]=sstB2a[K]
sstB5a[K]=sstB3a[K]
sstB4[I,J,]=sstB4a
sstB5[I,J,]=sstB5a

###
lat2=seq(-87.5,87.5,2.5)
lon2=seq(2.5,357.5,2.5)
sstCoarse<-array(NaN,dim=c(143,71,30))
library(fields)
for(i in 1:length(lat2))
  for(j in 1:length(lon2))
{
    I=which(latB>=(lat2[i]-1.25) & latB<=(lat2[i]+1.25))
    J=which(lonB>=(lon2[j]-1.25) & lonB>=(lon2[j]+1.25))
    if(length(I)>1 & length(J)>1) sstCoarse[j,i,]=apply(sstB[J,I,],3,mean,na.rm=T)
  }
   
for(i in 1:30)
{
#  a=interp.surface.grid(list(x=lonB,y=latB,z=sstB[,,i]),list(x=lon2,y=lat2))
#  b=interp.surface.grid(list(x=lon2,y=lat2,z=a$z),list(x=lonB,y=latB))
#  sstB6[,,i]=b$z*mask
  b=interp.surface.grid(list(x=lon2,y=lat2,z=sstCoarse[,,i]),list(x=lonB,y=latB))
  sstB7[,,i]=b$z*mask
}


pdf(file="bran_sstchange_8June07_NCL.pdf",width=15,height=15)
layout(cbind(c(1,3),c(2,4)))
image(lonB,latB,sstB[,,8],col=rich.colors(10),xlim=c(140,170),ylim=c(-50,-20),zlim=c(7,27),main="Default")
image(lonB,latB,sstB2[,,8],col=rich.colors(10),xlim=c(140,170),ylim=c(-50,-20),zlim=c(7,27),main="Constant SST")
image(lonB,latB,sstB3[,,8],col=rich.colors(10),xlim=c(140,170),ylim=c(-50,-20),zlim=c(7,27),main="Remove warm (daily) anomaly")
image(lonB,latB,sstB4a[,,8],col=rich.colors(10),xlim=c(140,170),ylim=c(-50,-20),zlim=c(7,27),main="Smooth to 2.5 degree grid")
dev.off()

### Look at variability throughout 2007
a=open.nc("/srv/ccrc/data36/z3478332/WRF/DATA/BRAN/bran3p5_sst_2007.nc")
latB=var.get.nc(a,"yt_ocean")
lonB=var.get.nc(a,"xt_ocean")
I=which(latB<=-10 & latB>=-50)
latB=latB[I]
J=which(lonB>=110 & lonB<=170)
lonB=lonB[J]
timeB=var.get.nc(a,"Time")
time=as.Date(timeB,origin="1993-01-01 00:00:00")
sstB=var.get.nc(a,"temp",c(J[1],I[1],1,1),c(length(lonB),length(latB),1,length(timeB)),unpack=T)

I=which(lonB>=160 & lonB<=165)
sstA=apply(sstB[I,,],c(2,3),mean,na.rm=T)
sstAa=cbind(sstA[,351:365],sstA,sstA[,1:15])
sstAb=array(0,dim(sstA))
for(i in 1:365) sstAb[,i]=apply(sstAa[,i:(i+30)],1,mean,na.rm=T)

#diff<-array(0,dim(sstB))
#for(i in 1:length(lonB)) diff[i,,]=sstB[i,,]-sstA
#diff2<-apply(diff,c(1,2),mean,na.rm=T)

diff1a<-array(0,dim(sstB))
for(i in 1:length(lonB)) diff1a[i,,]=sstB[i,,]-sstAb
diff2a<-apply(diff1a,c(1,2),mean,na.rm=T)

I1=which(latB<=-25 & latB>=-40)
J1=which(lonB>=150 & lonB<=155)
J2=which(lonB>=140 & lonB<=160)
I2=which(latB>=-50 & latB<=-20)
diff2b=diff2a[J2,I2]

daily=data.frame(Date=time,BoxMeanDiff=rep(0,length(time)),MeanDiff=rep(0,length(time)),
                 WarmArea=rep(0,length(time)),MeanDiffWarm=rep(0,length(time)))
for(i in 1:length(time))
{
  daily[i,2]=mean(diff1a[J1,I1,i],na.rm=T)
  
  a=diff1a[J2,I2,i]
  K=which(diff2b>0)
  daily[i,3]=mean(a[K],na.rm=T)
  K=which(a>0)
  daily[i,4]=length(K)
  daily[i,5]=mean(a[K],na.rm=T)
}

plot(daily[,1],daily[,2],type="l",lwd=3,xlab="Date",ylab="SST difference in EAC region",ylim=c(0,3))
lines(daily[,1],daily[,3],lty=2,lwd=3)
lines(daily[,1],daily[,5],lty=3,lwd=3)
legend("bottomleft",c("25-40, 1550-155 box","Area where mean diff >0","Area where daily diff >0"),col="black",lty=1:3)

month=as.numeric(format(time, "%m"))
mmean<-msd<-matrix(0,12,4)
for(i in 1:12)
{
  I=which(month==i)
  mmean[i,]=apply(daily[I,2:5],2,mean)
  msd[i,]=apply(daily[I,2:5],2,sd)
}

diff3a=array(NaN,c(dim(diff2),12))
for(i in 1:12)
{
  I=which(month==i)
  diff3a[,,i]=apply(diff1a[,,I],c(1,2),mean,na.rm=T)
}

pdf(file="bran_meanEACanomaly_2007.pdf",width=15,height=15)
layout(cbind(c(1,3),c(2,4)))
image(lonB,latB,diff3a[,,1],col=rich.colors(10),xlim=c(140,170),ylim=c(-50,-10),zlim=c(-5,5),main="January 2007")
image(lonB,latB,diff3a[,,2],col=rich.colors(10),xlim=c(140,170),ylim=c(-50,-10),zlim=c(-5,5),main="February 2007")
image(lonB,latB,diff3a[,,6],col=rich.colors(10),xlim=c(140,170),ylim=c(-50,-10),zlim=c(-5,5),main="June 2007")
image(lonB,latB,diff3a[,,7],col=rich.colors(10),xlim=c(140,170),ylim=c(-50,-10),zlim=c(-5,5),main="July 2007")
dev.off()


I=which(lonB>=160 & lonB<=165)
sstA=apply(sstB[I,,],c(2,3),mean,na.rm=T)
sstB3<-sstB2<-array(0,dim(sstB))
for(i in 1:length(lonB)) sstB2[i,,]=sstA
diff=apply(sstB-sstB2,c(1,2),mean)
for(i in 1:dim(sstB)[3]) sstB3[,,i]=sstB[,,i]-diff

I=which(lonB>=140 & lonB<=160)
J=which(latB>=-50 & latB<=-20)
sstB4<-sstB5<-sstB6<-sstB
sstB4a<-sstB5a<-sstB6a<-sstB[I,J,]
sstB2a=sstB2[I,J,]
sstB3a=sstB3[I,J,]
d2=diff[I,J]
K=which((sstB4a-sstB2a)>0)
sstB4a[K]=sstB2a[K]
sstB5a[K]=sstB3a[K]
sstB6a[K]=sstB6a[K]+d2[K]
sstB4[I,J,]=sstB4a
sstB5[I,J,]=sstB5a
sstB6[I,J,]=sstB6a

######## What about regrid to 2.5 and then back to 0.5?
## Test NCL results

f1=open.nc("/srv/ccrc/data36/z3478332/WRF/geogrid/geo_em.d02.notopo.nc")
lat1=var.get.nc(f1,"XLAT_M")
lon1=var.get.nc(f1,"XLONG_M")
lon1[lon1<0]=lon1[lon1<0]+360

test=array(NaN,c(dim(lat1),4))
a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_2007/wrflowinp_d02_2007-01")
sst2=var.get.nc(a,"SST")
test[,,1]=sst2[,,1]
a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_bran_2007/wrflowinp_d02_2007-01")
sst2=var.get.nc(a,"SST")
test[,,2]=sst2[,,1]
a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_bran_2007_noeac/wrflowinp_d02_2007-01")
sst2=var.get.nc(a,"SST")
test[,,3]=sst2[,,1]
a=open.nc("/srv/ccrc/data37/z3478332/WRF/WRF_boundary/default_bran_2007_2eac/wrflowinp_d02_2007-01")
sst2=var.get.nc(a,"SST")
test[,,4]=sst2[,,1]
I=which(test==0)
test[I]=NaN

pdf(file="bran_sstchange_1Jan2007_NCL2.pdf",width=15,height=10)
layout(cbind(c(1,3),c(2,4)))
names=c("Default","BRAN","Remove EAC","Coarse SSTs")
for(i in 1:4) image(test[,,i],zlim=c(285,305),col=rich.colors(10),axes=FALSE,main=names[i])
dev.off()


################
ColorBar <- function(brks,cols,labels=NA)
{
  par(mar = c(3, 1, 3, 4), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  if(is.na(labels[1])) labels=brks[seq(2, length(brks)-1)]
  
  axis(4, at = seq(1.5, length(brks) - 1.5), tick = TRUE, 
       labels = labels)
}

pdf(file="bran_sstchange_1June2007_NoNZ.pdf",width=7,height=9)
layout(cbind(c(1,2),c(3,3)),width=c(1,0.3))
par(mar=c(1,1,3,1))
a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_nznoland/wrflowinp_d01_2007-06")
sst=var.get.nc(a,"SST")
image(sst[,,1],zlim=c(275,305),col=rich.colors(10),axes=F,main="Default")
a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_nznoland_chris/wrflowinp_d01_2007-06")
sst2=var.get.nc(a,"SST")
image(sst2[,,1],zlim=c(275,305),col=rich.colors(10),axes=F,main="Chris's SSTs")
ColorBar(seq(275,305,3),rich.colors(10))
dev.off()


########## June 2007, change in average SSTs
library(RNetCDF)
f1=open.nc("/srv/ccrc/data36/z3478332/WRF/geogrid/geo_em.d02.notopo.nc")
lat1=var.get.nc(f1,"XLAT_M")
lon1=var.get.nc(f1,"XLONG_M")
lon1[lon1<0]=lon1[lon1<0]+360

test=array(NaN,c(dim(lat1),4))
a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_2007/wrflowinp_d02_2007-06")
sst2=var.get.nc(a,"SST")
test[,,1]=apply(sst2,c(1,2),mean)
a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_bran_2007/wrflowinp_d02_2007-06")
sst2=var.get.nc(a,"SST")
test[,,2]=apply(sst2,c(1,2),mean)
a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_bran_2007_noeac/wrflowinp_d02_2007-06")
sst2=var.get.nc(a,"SST")
test[,,3]=apply(sst2,c(1,2),mean)
a=open.nc("/srv/ccrc/data37/z3478332/WRF/WRF_boundary/default_bran_2007_2eac/wrflowinp_d02_2007-06")
sst2=var.get.nc(a,"SST")
test[,,4]=apply(sst2,c(1,2),mean)
I=which(test==0)
test[I]=NaN
test=test-273.15

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("blue","cyan","white","yellow","red"),c(20,5,5,20))
bb=c(-10,seq(0,2,0.25),10)
cm=pal(1)

tmp=test[,,4]-test[,,2]

library(fields)
image.plot(test[,,3],col=pal(17),zlim=c(10,30),main="EAC in June 2007")

pdf(file="bran_sstchange_June2007_v2.pdf",width=11,height=10)
layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.3))
par(mar=c(2,1,3,1))
names=c("Default","BRAN","No EAC","Doubled EAC")
for(i in 1:4) image(test[130:325,,i],col=pal(17),zlim=c(10,27),axes=F,main=names[i],cex.main=2)
ColorBar(10:27,pal(17),subsampleg = 2)
dev.off()




ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1.5)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = brks[seq(1, length(brks), subsampleg)])
}

cm=rich.colors(9)
bb=c(-100,5,10,15,17.5,20,22.5,25,30,100)

pdf(file="/home/nfs/z3478332/Documents/ECLs/WRFruns/SSTcomp_June07_d02.pdf",width=15,height=10)
layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.3))
names=c("Default","BRAN","Remove EAC","Double EAC")
for(i in 1:4) {
  image(test[,,i],col=cm,breaks=bb,axes=FALSE,main=names[i])
  map(add=T,lwd=2)
}
ColorBar(bb,cm)
dev.off()

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("white","yellow","red"))
bb=c(-10,seq(0,2,0.25),10)
cm=pal(10)

tmp=test[,,4]-test[,,2]

image.plot(test[,,4]-test[,,2],col=cm,breaks=bb,main="EAC in June 2007")




