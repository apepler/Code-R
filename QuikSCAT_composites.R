### MLDB composites
setwd("/srv/ccrc/data37/z3478332/ECL_evaluation/QuikSCAT/MLDB/")
ECLs<-read.csv("MLDB_Update_19702008.csv",header=T)
ECLs=ECLs[ECLs$Date>=19880000,]
ECLs$Year=floor(ECLs$Date/10000)
ECLs$Month=floor(ECLs$Date/100)%%100

library(RNetCDF)
library(fields)
library(akima)

a=open.nc("ECLwind_QuikSCAT_MLDB.nc")
u=var.get.nc(a,"ECL_U10")
v=var.get.nc(a,"ECL_V10")
w=var.get.nc(a,"ECL_WS10")
lat<-lon<-seq(-9.875,9.875,0.25)

source('~/Documents/R/color.palette.R')
pal1 <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
pal2 <- color.palette(c("white","cyan","blue","black"), c(20,20,20))

image.plot(apply(w,c(1,2),mean),breaks=c(seq(0,10,1),100),col=pal2(11),zlim=c(0,10))
image.plot(apply(u,c(1,2),mean),breaks=c(seq(-7,7,1),100),col=pal1(15),zlim=c(-7,7))
image.plot(apply(v,c(1,2),mean),breaks=c(seq(-7,7,1),100),col=pal1(15),zlim=c(-7,7))

ECLs$MaxWind=apply(w,3,max)
ECLs$MeanWind=apply(w[21:60,21:60,],3,mean)

I=which(ECLs$MaxWind>=16.7) ## 50 km/h
image.plot(apply(w[,,I],c(1,2),mean),breaks=c(seq(0,14,1),100),col=pal2(15),zlim=c(0,14))
image.plot(apply(u[,,I],c(1,2),mean),breaks=c(seq(-7,7,1),100),col=pal1(15),zlim=c(-7,7))
image.plot(apply(v[,,I],c(1,2),mean),breaks=c(seq(-7,7,1),100),col=pal1(15),zlim=c(-7,7))

#### Look at the mean/maximum u & v gradients

grads<-array(0,c(length(ECLs[,1]),2,11))
dimnames(grads)[[2]]=c("N-S gradient of zonal wind","E-W gradient of meridional wind")
dimnames(grads)[[3]]=c(paste(1:9,"degrees radius"),"Maximum gradient","Radius of maximum gradient")
for(i in 1:9)
{
  grads[,1,i]=apply(u[37:44,(40:41+(i*4)),]-u[37:44,(40:41-(i*4)),],3,mean)
  grads[,2,i]=apply(v[(40:41-(i*4)),37:44,]-v[(40:41+(i*4)),37:44,],3,mean)
}

grads[,,10]=apply(grads[,,1:9],c(1,2),max)
for(i in 1:length(ECLs[,1])) 
  for(j in 1:2)
    grads[i,j,11]=which(grads[i,j,1:9]==grads[i,j,10])

I=which(grads[,2,10]<=0 | grads[,1,10]<=0)
ECLs[I,]

I=which(ECLs$Date==20070608)
image.plot(u[,,I],breaks=c(seq(-12,12,2),100),col=pal1(13),zlim=c(-12,12),main="Zonal wind")
image.plot(v[,,I],breaks=c(seq(-12,12,2),100),col=pal1(13),zlim=c(-12,12),main="Meridional wind")
image.plot(w[,,I],breaks=c(seq(0,14,1),100),col=pal2(15),zlim=c(0,14))
