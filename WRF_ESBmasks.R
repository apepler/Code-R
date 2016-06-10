rm(list=ls())

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

dir="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_ensemble/out/2007052700/"

fname="wrfhrly_d01_2007-06-01_00:00:00"
a=open.nc(paste(dir,fname,sep=""))
lat50=var.get.nc(a,"XLAT",c(1,1,1),c(215,144,1))
lon50=var.get.nc(a,"XLONG",c(1,1,1),c(215,144,1))

fname="wrfhrly_d02_2007-06-01_00:00:00"
a=open.nc(paste(dir,fname,sep=""))
lat10=var.get.nc(a,"XLAT",c(1,1,1),c(325,200,1))
lon10=var.get.nc(a,"XLONG",c(1,1,1),c(325,200,1))

mask50<-mask50a<-matrix(0,215,144)
mask10<-mask10a<-matrix(0,325,200)

for(i in 2:214)
  for(j in 2:143)
  {
    I=which(Useful$x>(lon50[i-1,j]+lon50[i,j])/2 & Useful$x<(lon50[i+1,j]+lon50[i,j])/2)
    J=which(Useful$y>(lat50[i,j-1]+lat50[i,j])/2 & Useful$y<(lat50[i,j+1]+lat50[i,j])/2)
    
    if(length(I)>0 & length(J)>0) mask50[i,j]=mean(mask[I,J])
  }
mask50a[mask50>=0.5]=1

for(i in 2:324)
  for(j in 2:199)
  {
    I=which(Useful$x>(lon10[i-1,j]+lon10[i,j])/2 & Useful$x<(lon10[i+1,j]+lon10[i,j])/2)
    J=which(Useful$y>(lat10[i,j-1]+lat10[i,j])/2 & Useful$y<(lat10[i,j+1]+lat10[i,j])/2)
    
    if(length(I)>0 & length(J)>0) mask10[i,j]=mean(mask[I,J])
  }
mask10a[mask10>=0.5]=1

##Oh, need to make sure has land-sea mask of WRF taken into account too.

a=open.nc("/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_ensemble/wrfinput_d02_2007052700")
land10=var.get.nc(a,"LANDMASK")
I=which(mask10>0.5 & land10==0)
mask10a[I]=0

mask10=mask10a
mask50=mask50a
I=which(mask10a==0)
mask10a[I]=NaN
I=which(mask50a==0)
mask50a[I]=NaN

#save(lat50,lon50,mask50,mask50a,lat10,lon10,mask10,mask10a,file='~/Documents/Data/Useful_WRF.RData')

library(RNetCDF)

# define dimensions
nc<-create.nc("~/WRF_d01_ESB_mask.nc")
dim.def.nc(nc,"west_east",215)
dim.def.nc(nc,"south_north",144)
var.def.nc(nc,"ESB","NC_INT",c(0,1))
var.put.nc(nc,"ESB",mask50)
close.nc(nc)

nc<-create.nc("~/WRF_d02_ESB_mask.nc")
dim.def.nc(nc,"west_east",325)
dim.def.nc(nc,"south_north",200)
var.def.nc(nc,"ESB","NC_INT",c(0,1))
var.put.nc(nc,"ESB",mask10)
close.nc(nc)

