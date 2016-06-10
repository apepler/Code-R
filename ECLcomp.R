rm(list=ls())
setwd('~/Documents/ECLs')
ncep1<-read.csv('ECLfixes_NCEP1.csv')
ncep1=ncep1[ncep1[,1]>=19900000 & ncep1[,1]<20090000,]
ncep2<-read.csv('ECLfixes_NCEP2.csv')
erai<-read.csv('ECLfixes_ERAI.csv',sep=";")
erai=erai[erai[,1]<20090000,]

lats=seq(-40.5,-23.5,1)
lons=seq(149.5,160.5,1)

ERA<-N2<-N1<-matrix(0,18,12)
for(i in 1:12)
  for(j in 1:18)
  {
    I=which(ncep2[,4]>=floor(lons[i]) & ncep2[,4]<ceiling(lons[i]) & ncep2[,5]>=floor(lats[j]) & ncep2[,5]<ceiling(lats[j]))
    N2[j,i]=length(I)
    I=which(ncep1[,4]>=floor(lons[i]) & ncep1[,4]<ceiling(lons[i]) & ncep1[,5]>=floor(lats[j]) & ncep1[,5]<ceiling(lats[j]))
    N1[j,i]=length(I)
    I=which(erai[,4]>=floor(lons[i]) & erai[,4]<ceiling(lons[i]) & erai[,5]>=floor(lats[j]) & erai[,5]<ceiling(lats[j]))
    ERA[j,i]=length(I)
  }

###### Now the easy way, with using programs
source('~/Documents/R/ECLextract.R')
setwd('~/output_merraNUa')
ECLdatabase2(1980,2008,'MERRANUa')
setwd('~/output_merraU_50')
ECLdatabase2(1980,2008,'MERRAU50')
setwd('~/output_merraU_50a')
ECLdatabase2(1980,2008,'MERRAU50a')
setwd('~/output_merraU_100')
ECLdatabase2(1980,2008,'MERRAU100')
setwd('~/output_merraU_150')
ECLdatabase2(1980,2008,'MERRAU150')
setwd('~/output_merraU_250')
ECLdatabase2(1980,2008,'MERRAU250')

source('~/Documents/R/ECLextract.R')
setwd('~/output_cfsr_100')
ECLdatabase2(1981,2008,'cfsr_100')
setwd('~/output_cfsr_150')
ECLdatabase2(1981,2008,'cfsr_150')
setwd('~/output_cfsr_250')
ECLdatabase2(1981,2008,'cfsr_250')

source('~/Documents/R/ECLextract.R')
setwd('~/outputUM_erai_150_topo_avrad5')
ECLdatabase2(1980,2009,'umelb_erai_topo_rad5')
setwd('~/outputUM_erai_150_topo_rad5_diff1')
ECLdatabase2(1980,2009,'umelb_erai_topo_rad5_diff1')
setwd('~/outputUM_erai_150_topo_rad5_diff0')
ECLdatabase2(1980,2009,'umelb_erai_topo_rad5_diff0')
setwd('~/outputUM_erai_150_topo_rad5_proj75')
ECLdatabase2(1980,2009,'umelb_erai_topo_rad5_proj75')

source('~/Documents/R/ECLextract.R')
setwd('~/outputUM_wrfR1_50_default')
ECLdatabase2(1980,2009,'umelb_wrfR1_default')
setwd('~/outputUM_wrfR2_50_default')
ECLdatabase2(1980,2009,'umelb_wrfR2_default')
setwd('~/outputUM_wrfR3_50_default')
ECLdatabase2(1980,2009,'umelb_wrfR3_default')




rm(list=ls())
setwd('~/Documents/ECLs')
source('~/Documents/R/ECLextract.R')

ncep1<-read.csv('ECLfixes_aus_NCEP1.csv')
ncep2<-read.csv('ECLfixes_aus_NCEP2.csv')
erai<-read.csv('ECLfixes_aus_ERAI.csv')
erai=erai[erai[,8]>=0.3,]
eraih<-read.csv('ECLfixes_aus_ERAIh.csv')
jra<-read.csv('ECLfixes_aus_JRA25.csv')


ECLfig_ann(ncep1,"NCEP1")
ECLfig_ann(ncep2,"NCEP2")
ECLfig_ann(erai,"ERAI")
ECLfig_ann(eraih,"ERAIH")
ECLfig_ann(jra,"JRA25")

seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7)
{
ECLfig_seas(ncep1,"NCEP1",snum[i,],seasons[i])
ECLfig_seas(ncep2,"NCEP2",snum[i,],seasons[i])
ECLfig_seas(erai,"ERAI",snum[i,],seasons[i])
ECLfig_seas(eraih,"ERAIH",snum[i,],seasons[i])
ECLfig_seas(jra,"JRA25",snum[i,],seasons[i])
}

##ERAI with land-sea mask

setwd('~/Documents/Data')
library("RNetCDF")
library(fields)
library("R.matlab")
readMat('Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
f1=open.nc('erai_landmask.nc')
lat=var.get.nc(f1,'latitude')
lon=var.get.nc(f1,'longitude')
lsm=var.get.nc(f1,'lsm',unpack=T)

erai<-read.csv('~/Documents/ECLs/ECLfixes_aus_ERAI.csv')
smon=c(12,2)
data=erai[erai[,8]>=0.3,]
mon=floor((data[,2] %% 10000)/100)
data=data[mon>=smon[1] | mon<=smon[2],]
all<-matrix(0,length(lat),length(lon))
for(i in 1:length(lon))
  for(j in 1:length(lat))
  {
    I=which(data[,5]>=lon[i]-0.75 & data[,5]<lon[i]+0.75 & data[,6]>=lat[j]-0.75 & data[,6]<lat[j]+0.75)
    all[j,i]=length(I)
  }
all=(all/19)
I=which(all>1.95)
all[I]=1.95
cols=gray(seq(1,0.1,-0.1))
col2=add.alpha(cols,0.8)

image(lon,lat[121:1],lsm[,121:1],col=c("lightskyblue","palegreen"),xlim=c(110,160),ylim=c(-45,-10),xlab="",ylab="")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat[121:1],t(all[121:1,]),xlab="",ylab="",breaks=seq(0,2,0.2),col=col2,zlim=c(0,2),add=T)

add.alpha <- function(COLORS, ALPHA){
  if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
  RGB <- col2rgb(COLORS, alpha=TRUE)
  RGB[4,] <- round(RGB[4,]*ALPHA)
  NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
  return(NEW.COLORS)
}

##Different version - using the heatmap approach.

smoothScatter(data[,5],data[,6],colramp=colorRampPalette(c("white", "black")),xlab="", ylab="",xlim=c(110,160),ylim=c(-45,-10),col="black")
cc=c("white","palegreen")
cc=add.alpha(cc,0.1)
image(lon,lat[121:1],lsm[,121:1],col=cc,xlim=c(110,160),ylim=c(-45,-10),xlab="",ylab="",add=T)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

##Upper level lows 
source('~/Documents/R/ECLextract.R')
setwd('~/output_erai_925')
ECLdatabase(1990,2008,'ERA925')
##Now load, & change CV thresh to 2.5
setwd('~/Documents/ECLs')
erai<-read.csv('ECLfixes_aus_ERA925.csv')
ECLfig_ann(erai,"ERAI_925_cv1.5")
erai2=erai[erai[,8]>=3,]
ECLfig_ann(erai2,"ERAI_925_cv3")
ECLfig_seas(erai,"ERAI_925_cv1.5",c(12,2),'djf')
ECLfig_seas(erai2,"ERAI_925_cv3",c(12,2),'djf')

##Spatial compare to MLDB?
rm(list=ls())
source('~/Documents/R/ECLextract.R')
setwd('~/Documents/ECLs')
ncep1<-read.csv('ECLfixes_aus_NCEP1.csv')
ncep2<-read.csv('ECLfixes_aus_NCEP2.csv')
erai<-read.csv('ECLfixes_aus_ERAI.csv')
erai2=erai[erai[,8]>=0.3,]
eraih<-read.csv('ECLfixes_aus_ERAIh.csv')
eraiU<-read.csv('ECLfixes_aus_ERA925.csv')
jra<-read.csv('ECLfixes_aus_JRA25.csv')

seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7)
{
  ECL_FA_seas(ncep1,"NCEP1",snum[i,],seasons[i])
  ECL_FA_seas(ncep2,"NCEP2",snum[i,],seasons[i])
  ECL_FA_seas(erai,"ERAI",snum[i,],seasons[i])
  ECL_FA_seas(erai2,"ERAIcv3",snum[i,],seasons[i])
  ECL_FA_seas(eraih,"ERAIH",snum[i,],seasons[i])
  ECL_FA_seas(eraiU,"ERAI925",snum[i,],seasons[i])
  ECL_FA_seas(jra,"JRA25",snum[i,],seasons[i])
}

for(i in 1:7)
{
  ECL_rain_seas(ncep1,"NCEP1",snum[i,],seasons[i])
  ECL_rain_seas(ncep2,"NCEP2",snum[i,],seasons[i])
  ECL_rain_seas(erai,"ERAI",snum[i,],seasons[i])
  ECL_rain_seas(erai2,"ERAIcv3",snum[i,],seasons[i])
  ECL_rain_seas(eraih,"ERAIH",snum[i,],seasons[i])
  ECL_rain_seas(eraiU,"ERAI925",snum[i,],seasons[i])
  ECL_rain_seas(jra,"JRA25",snum[i,],seasons[i])
}

##
setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
erai1<-read.csv('CSV/ECLfixes_aus_ERAIa.csv')
erai1a=erai1[erai1[,8]>=0.3,]
erai2<-read.csv('CSV/ECLfixes_aus_ERAIa_925.csv')
erai3<-read.csv('CSV/ECLfixes_aus_ERAIa_850.csv')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7)
{
  ECLfig_seas2(erai1,"ERAI",snum[i,],seasons[i],33)
  ECLfig_seas2(erai1a,"ERAIcv3",snum[i,],seasons[i],33)
  ECLfig_seas2(erai2,"ERAI925",snum[i,],seasons[i],33)
  ECLfig_seas2(erai3,"ERAI850",snum[i,],seasons[i],33)
}

##
rm(list=ls())
setwd('~/Documents/ECLs')
source('~/Documents/R/ECLextract.R')
wrf<-read.csv('CSV/ECLfixes_aus_NCEP1.csv')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7) ECLfig_seas2a(wrf,"NCEP",snum[i,],seasons[i])

##
cfsr<-read.csv('CSV/ECLfixes_aus_cfsr_100.csv')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7) ECLfig_seas2b(cfsr,"CFSR_100",snum[i,],seasons[i],28,10)

setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
m1<-read.csv('CSV/ECLfixes_aus_merra_250.csv')
m2<-read.csv('CSV/ECLfixes_aus_MERRA150.csv')
m3<-read.csv('CSV/ECLfixes_aus_MERRA100.csv')
m4<-read.csv('CSV/ECLfixes_aus_MERRA50.csv')
m5<-read.csv('CSV/ECLfixes_aus_MERRA50a.csv')

seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7)
{
  ECLfig_seas2b(m1,"MERRA250",snum[i,],seasons[i],29,10)
  ECLfig_seas2b(m2,"MERRA150",snum[i,],seasons[i],29,20)
  ECLfig_seas2b(m3,"MERRA100",snum[i,],seasons[i],29,40)
  ECLfig_seas2b(m4,"MERRA50",snum[i,],seasons[i],29,100)
  ECLfig_seas2b(m5,"MERRA50a",snum[i,],seasons[i],29,60)
}

length(which(m1[,5]>=155 & m1[,5]<=160 & m1[,6]>=(-40) & m1[,6]<=(-30) & m1[,8]>=0.25))/29
length(which(m1[,5]>=145 & m1[,5]<=150 & m1[,6]>=(-40) & m1[,6]<=(-30) & m1[,8]>=0.25))/29

lats=seq(-48,-10,2)
lons=seq(110,180,2)
a250<-a100<-matrix(0,length(lats),length(lons))
for(i in 1:length(lons))
  for(j in 1:length(lats))
  {
    I=which(m1[,5]>=lons[i]-1 & m1[,5]<lons[i]+1 & m1[,6]>=lats[j]-1 & m1[,6]<lats[j]+1)
    a250[j,i]=length(I)
    I=which(m3[,5]>=lons[i]-1 & m3[,5]<lons[i]+1 & m3[,6]>=lats[j]-1 & m3[,6]<lats[j]+1)
    a100[j,i]=length(I)
  }

diff=a100/a250
diff[a250<0.5]=NaN
cols=gray(seq(1,0,length.out=8))
diff2=diff
diff2[diff>105]=105
image.plot(lons,lats,t(diff2),xlab="",breaks=c(0,1,2,5,10,20,50,100,110),ylab="",col=cols,zlim=c(1,110),axis.args=list(at=c(1,5,10,20,50,100),labels=c("1x","5x","10x","20x","50x","100x")))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
m1<-read.csv('CSV/ECLfixes_aus_WRF250.csv')
m2<-read.csv('CSV/ECLfixes_aus_WRF150.csv')
m3<-read.csv('CSV/ECLfixes_aus_WRF50.csv')

seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7)
{
  ECLfig_seas2b(m1,"WRF250",snum[i,],seasons[i],29,10)
  ECLfig_seas2b(m2,"WRF150",snum[i,],seasons[i],29,20)
  ECLfig_seas2b(m3,"WRF50",snum[i,],seasons[i],29,60)
}

lats=seq(-48,-10,2)
lons=seq(110,180,2)
a250<-a150<-matrix(0,length(lats),length(lons))
for(i in 1:length(lons))
  for(j in 1:length(lats))
  {
    I=which(m1[,5]>=lons[i]-1 & m1[,5]<lons[i]+1 & m1[,6]>=lats[j]-1 & m1[,6]<lats[j]+1)
    a250[j,i]=length(I)
    I=which(m2[,5]>=lons[i]-1 & m2[,5]<lons[i]+1 & m2[,6]>=lats[j]-1 & m2[,6]<lats[j]+1)
    a150[j,i]=length(I)
  }

diff=a150/a250
diff[a250<0.5]=NaN
cols=gray(seq(1,0,length.out=9))
diff2=diff
diff2[diff>25]=25
image.plot(lons,lats,t(diff2),xlab="",breaks=c(0,1,2,4,6,8,10,15,20,30),ylab="",col=cols,zlim=c(1,25),axis.args=list(at=c(1,5,10,15,20),labels=c("1x","5x","10x","15x","20x")))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

##Average CV over land/ocean?
setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
m1<-read.csv('CSV/ECLfixes_aus_MERRA250.csv')
m2<-read.csv('CSV/ECLfixes_aus_MERRA150.csv')
m3<-read.csv('CSV/ECLfixes_aus_MERRA100.csv')
m4<-read.csv('CSV/ECLfixes_aus_MERRA50.csv')
m5<-read.csv('CSV/ECLfixes_aus_MERRA50a.csv')

lats=seq(-48,-10,2)
lons=seq(110,180,2)
a250<-a150<-a100<-a50<-a50a<-matrix(0,length(lats),length(lons))
for(i in 1:length(lons))
  for(j in 1:length(lats))
  {
    I=which(m1[,5]>=lons[i]-1 & m1[,5]<lons[i]+1 & m1[,6]>=lats[j]-1 & m1[,6]<lats[j]+1)
    a250[j,i]=mean(m1[I,8])
    I=which(m2[,5]>=lons[i]-1 & m2[,5]<lons[i]+1 & m2[,6]>=lats[j]-1 & m2[,6]<lats[j]+1)
    a150[j,i]=mean(m2[I,8])
    I=which(m3[,5]>=lons[i]-1 & m3[,5]<lons[i]+1 & m3[,6]>=lats[j]-1 & m3[,6]<lats[j]+1)
    a100[j,i]=mean(m3[I,8])
    I=which(m4[,5]>=lons[i]-1 & m4[,5]<lons[i]+1 & m4[,6]>=lats[j]-1 & m4[,6]<lats[j]+1)
    a50[j,i]=mean(m4[I,8])
    I=which(m5[,5]>=lons[i]-1 & m5[,5]<lons[i]+1 & m5[,6]>=lats[j]-1 & m5[,6]<lats[j]+1)
    a50a[j,i]=mean(m5[I,8])
  }

a50a[a50a>1.05]=1.05
image.plot(lons,lats,t(a50a),xlab="",breaks=seq(0,1.1,0.1),ylab="",col=cols,zlim=c(0,1.1))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

cols=gray(seq(1,0,-0.1))
diff=a50-a250
diff[diff>0.49]=0.49
image.plot(lons,lats,t(diff),xlab="",breaks=seq(-0.05,0.5,0.05),ylab="",col=cols,zlim=c(-0.05,0.5))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
name=c("50","50a","100","150","250")
b=c(100,60,40,20,10)
m1<-read.csv('CSV/ECLfixes_aus_MERRANUa.csv')
m1=m1[m1[,8]>=1.2,]
ECLfig_seas2b(m1,'MERRANUa',snum[1,],seasons[1],29,60)

setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
name=c("50","50a","100","150","250")
b=c(100,60,40,20,10)
for(j in 1:5)
{
  name1=paste('MERRA',name[j],sep="")
  fname1=paste('CSV/ECLfixes_aus_',name1,'.csv',sep="")
  m1<-read.csv(fname1)
  m1=m1[m1[,8]>=0.25,]
  name2=paste('MERRAU',name[j],sep="")
  fname2=paste('CSV/ECLfixes_aus_',name2,'.csv',sep="")
  m2<-read.csv(fname2)
  m2=m2[m2[,8]>=1.2,]
  ECLfig_seas2c(m1,paste(name1,'_v3',sep=""),snum[1,],seasons[1],29,25)
  ECLfig_seas2c(m2,paste(name2,'_v3',sep=""),snum[1,],seasons[1],29,25)
}

##
setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
res=c(50,100,150,250)
for(k in 1:4)
{
  data<-read.csv(paste('CSV/ECLfixes_aus_cfsr_',res[k],'.csv',sep=""))
  for(i in 1:7) ECLfig_seas2(data,paste("CFSR_",res[k],sep=""),snum[i,],seasons[i],28)
}
for(k in 1:4)
{
  data<-read.csv(paste('CSV/ECLfixes_aus_cfsr_',res[k],'.csv',sep=""))
  for(i in 1:7) ECLfig_seas2b(data,paste("CFSR_",res[k],sep=""),snum[i,],seasons[i],28,10)
}


##KS test for whether same distribtuion of duration/intensity
rm(list=ls())
setwd('~/Documents/ECLs')
ncep1<-read.csv('CSV/ECLevents_unsw_ncep1_250.csv')
ncep2<-read.csv('CSV/ECLevents_unsw_ncep2_250.csv')
erai<-read.csv('CSV/ECLevents_unsw_erai_300.csv')
jra25<-read.csv('CSV/ECLevents_unsw_jra_250.csv')
merra<-read.csv('CSV/ECLevents_unsw_merra_250.csv')
cfsr<-read.csv('CSV/ECLevents_unsw_cfsr_250.csv')

ncep1=ncep1[ncep1$Date1>=19900000 & ncep1$Date1<=20081231,]
ncep2=ncep2[ncep2$Date1>=19900000 & ncep2$Date1<=20081231,]
erai=erai[erai$Date1>=19900000 & erai$Date1<=20081231,]
jra25=jra25[jra25$Date1>=19900000 & jra25$Date1<=20081231,]
merra=merra[merra$Date1>=19900000 & merra$Date1<=20081231,]
cfsr=cfsr[cfsr$Date1>=19900000 & cfsr$Date1<=20081231,]

datasets=list(ncep1,ncep2,erai,jra25,merra,cfsr)

KS<-matrix(NaN,6,6)
cv=0.5
for(i in 1:6)
  for(j in i:6)
  {
    a=ks.test(as.numeric(datasets[[i]][datasets[[i]]$CV>=cv,]$Length),as.numeric(datasets[[j]][datasets[[j]]$CV>=cv,]$Length))
    KS[i,j]=a$p.value
  }



#Version 2 - the edits in the dgof package
library(dgof)
KSa<-KS2a<-matrix(NaN,6,6)
for(i in 1:6)
  for(j in i:6)
  {
    a=ks.test(as.numeric(datasets[[i]]$Length),as.numeric(datasets[[j]]$Length))
    KSa[i,j]=a$p.value
    a=ks.test(as.numeric(datasets[[i]][datasets[[i]]$CV>=0.25,]$Length),as.numeric(datasets[[j]][datasets[[j]]$CV>=0.25,]$Length))
    KS2a[i,j]=a$p.value
  }

######Umelb

setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
res=c("topo_rad5","topo_rad5_diff1","topo_rad5_diff0","topo_rad5_proj75")

cols=gray(seq(1,0.1,-0.1))
cols2=gray(seq(1,0,-0.2))
lats=seq(-49.5,-10.5,1)
lons=seq(110.5,179.5,1)
for(k in 1:4)
{
  data<-read.csv(paste('CSV/ECLfixes_aus_umelb_erai_',res[k],'.csv',sep=""))
  m1<-ECLcount_all(data[data$CV>=0.25,],lats,lons)/30/4
  m1a=m1/sum(m1)
  m1[m1>2]=2
  m1a[m1a>0.06]=0.06  
  tiff(file=paste("Locs_aus_UM_erai_",res[k],sep=""), height=450, width=600)
  image.plot(lons,lats,t(m1),breaks=seq(0,2,0.2),xlab="",ylab="",col=cols,zlim=c(0,2))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

lats=seq(-40,-24,2)
lons=seq(149,161,2)
for(k in 1:4)
{
  data<-read.csv(paste('CSV/ECLfixes_aus_umelb_erai_',res[k],'.csv',sep=""))
  m1<-ECLcount_all(data[data$CV>=0.25,],lats,lons)/30/4
  m1a=m1/sum(m1)
  m1[m1>2]=2
  m1a[m1a>0.06]=0.06
  
  tiff(file=paste("Locs_esb_UM_erai_",res[k],sep=""), height=400, width=400)
  image.plot(lons,lats,t(m1),breaks=seq(0,2,0.2),xlab="",ylab="",col=cols,zlim=c(0,2))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("Locs_esb_UM_erai_",res[k],"_PC",sep=""), height=400, width=400)
  image.plot(lons,lats,t(m1a),breaks=seq(0,0.06,0.01),xlab="",ylab="",col=cols2,zlim=c(0,0.06),lab.breaks=c("0","1%","2%","3%","4%","5%","6%"))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  #ECLfig_seas2(data[data$CV>=0.25,],paste("UM_erai_",res[k],sep=""),c(1,12),"ann",30)
}

###
setwd('~/Documents/ECLs/Algorithm Comparison')
warm=read.csv('events_warm.csv')
cool=read.csv('events_cool.csv')

corrs=matrix(0,5,5)
for(i in 1:5)
  for(j in 1:5)
    if(j>i) corrs[i,j]=cor(warm[,i+1],warm[,j+1],use="pairwise.complete.obs") else corrs[i,j]=cor(cool[,i+1],cool[,j+1],use="pairwise.complete.obs")


##UM NCEP vs WRF
setwd('~/Documents/ECLs')
rm(list=ls())
source('~/Documents/R/ECLextract.R')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))

data<-read.csv('CSV/ECLfixes_aus_umelb_ncep1_default.csv')
I=which(data$CV>=0.25)
ECLfig_seas2b(data[I,],"UM_NCEPa",snum[1,],seasons[1],29,20)

data<-read.csv('CSV/ECLfixes_aus_umelb_wrfR1_default.csv')
I=which(data$CV>=0.25)
for(i in 2:3) ECLfig_seas2b(data[I,],"UM_WRF_R1",snum[i,],seasons[i],29,20)


