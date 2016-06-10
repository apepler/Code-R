rm(list=ls())
setwd('~/Documents/ECLs')
library(RNetCDF)
#load("/home/nfs/z3478332/Documents/GDI/JE_WRF/Rain_daily.RData")
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->UsefulE
maskE<-t(UsefulE$mask)
maskE[is.na(maskE)]=0
readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
library(RNetCDF)
library(akima)
a=open.nc("/srv/ccrc/data36/z3478332/WRF/output/WRF_d01_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
latt=as.vector(lat)
lont=as.vector(lon)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_AWAP_proj100_rad2cv1.RData")
ratio<-array(0,c(691,886,5,4))
dimnames(ratio)[[3]]<-c("MLDB","NCEP","ERAI","NCEP-WRF","CMIP-WRF")
dimnames(ratio)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ratio_9006<-ratio

totalrain_9009<-array(0,c(691,886,3,4))
dimnames(totalrain_9009)[[3]]<-c("AWAP","NCEP-WRF","CMIP-WRF")
dimnames(totalrain_9009)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")

ratio[,,1,]=apply(eclrain_mldb[,,1:37,],c(1,2,4),sum)/apply(rain[,,1:37,],c(1,2,4),sum)
ratio[,,2,]=apply(eclrain_ncep,c(1,2,4),sum)/apply(rain,c(1,2,4),sum)
ratio[,,3,]=apply(eclrain_erai[,,11:40,],c(1,2,4),sum)/apply(rain[,,11:40,],c(1,2,4),sum)

ratio_9006[,,1]=apply(eclrain_mldb[,,21:37,],c(1,2,4),sum)/apply(rain[,,21:37,],c(1,2,4),sum)
ratio_9006[,,2]=apply(eclrain_ncep[,,21:37,],c(1,2,4),sum)/apply(rain[,,21:37,],c(1,2,4),sum)
ratio_9006[,,3]=apply(eclrain_erai[,,21:37,],c(1,2,4),sum)/apply(rain[,,21:37,],c(1,2,4),sum)

totalrain_9009[,,1,]=apply(rain[,,21:40,],c(1,2,4),mean)

rm(c("rain","eclrain_mldb","eclrain_ncep","eclrain_erai"))

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_NCEP-WRF_proj100_rad2cv1.RData")
ratio2=apply(ECLrain37,c(1,2,4,5),sum)/apply(rain,c(1,2,4,5),sum)
ratio3=apply(ratio2,c(1,2,4),mean)
for(i in 1:4)
{
a=as.vector(ratio3[,,i])
a[which(is.na(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
ratio[,,4,i]=t(b$z)
}

ratio2=apply(ECLrain37[,,1:17,,],c(1,2,4,5),sum)/apply(rain[,,1:17,,],c(1,2,4,5),sum)
ratio3=apply(ratio2,c(1,2,4),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratio_9006[,,4,i]=t(b$z)
}

for(i in 1:4)
{
  a=as.vector(apply(rain[,,,,i],c(1,2),mean))
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  totalrain_9009[,,2,i]=t(b$z)
}


# par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
# filled.contour3(Useful$x,Useful$y,t(ratio[,,4,4]*UsefulE$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
# par(xpd = NA)
# contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
# par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
# filled.legend(Useful$x,Useful$y,t(ratio[,,4,4]*UsefulE$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "",    key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%","60%","70%"),at=seq(0,0.7,0.1))) 
# 

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_CMIP-WRF_proj100_rad2cv0.8.RData")
ratio2=apply(ECLrain37,c(1,2,4,5,6),sum)/apply(rain,c(1,2,4,5,6),sum)
ratio3=apply(ratio2,c(1,2,5),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratio[,,5,i]=t(b$z)
}
ratio2=apply(ECLrain37[,,1:17,,,],c(1,2,4,5,6),sum)/apply(rain[,,1:17,,,],c(1,2,4,5,6),sum)
ratio3=apply(ratio2,c(1,2,5),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratio_9006[,,5,i]=t(b$z)
}

for(i in 1:4)
{
  a=as.vector(apply(rain[,,,,,i],c(1,2),mean))
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  totalrain_9009[,,3,i]=t(b$z)
}

save(Useful,UsefulE,ratio,ratio_9006,totalrain_9009,file="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_ratio.RData")


##########

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_ratio.RData")
load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/WRFrain_d02.RData")

library(akima)
library(RNetCDF)
a=open.nc("/srv/ccrc/data36/z3478332/WRF/output/WRF_d02_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
latt=as.vector(lat)
lont=as.vector(lon)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->UsefulE
maskE<-t(UsefulE$mask)
maskE[is.na(maskE)]=0
readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
ColorBar <- function(brks,cols,vert=T,subsampleg=1,breaklab=brks)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = breaklab[seq(1, length(brks), subsampleg)])
}



rain_9009_d02<-array(0,c(691,886,2,4))
for(i in 1:4)
{
  a=as.vector(apply(rain_ncep[,,,,i],c(1,2),mean))
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  rain_9009_d02[,,1,i]=t(b$z)
  a=as.vector(apply(rain_cmip[,,,,,i],c(1,2),mean))
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  rain_9009_d02[,,2,i]=t(b$z)
}

library(abind)
rain_9009a=abind(rain_9009,rain_9009_d02,along=3)
names1=c("AWAP","NCEP-WRF_d01","CMIP-WRF_d01","NCEP-WRF_d02","CMIP-WRF_d02")
for(i in 1:5)
{
  pdf(file=paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Figures/Totalrain_9009_",names1[i],".pdf",sep=""),width=7,height=4)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
tmp=t(rain_9009a[,,i,1]*Useful$mask)
tmp[tmp>2000]=2050
image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,2100,100),col=rich.colors(21),frame.plot=F,
                xlim=c(135,160),ylim=c(-39,-24),xlab="",ylab="")
par(xpd = NA)
contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(ratio[,,4,1]*Useful$mask),lev=seq(0,2100,100),col=rich.colors(21),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "") 
dev.off()
}

pdf(file=paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Figures/Totalrain_9009_panel.pdf",sep=""),width=16,height=4,pointsize=12)
layout(cbind(1,2,3,4),width=c(1,1,1,0.2))
names=c("AWAP","NCEP-WRF","CMIP-WRF","NCEP-WRF","CMIP-WRF")
for(i in c(1,4,5))
{
  par(mar=c(2,3,3,1))
  tmp=t(rain_9009a[,,i,1]*Useful$mask)
  tmp[tmp>2000]=2050
  image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,2100,100),col=rich.colors(21),main=names[i],frame.plot=F,
      xlim=c(135,155),ylim=c(-39,-24),xlab="",ylab="",cex.axis=1.5)
  contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
}
ColorBar(seq(0,2100,100),rich.colors(21),T,2) 
dev.off()



for(i in 1:5)
{
  pdf(file=paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Figures/Rain25mm_9009_",names1[i],".pdf",sep=""),width=7,height=4)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  tmp=t(rain_9009a[,,i,3]*Useful$mask)
  tmp[tmp>20]=20.50
  image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,21),col=rich.colors(21),frame.plot=F,
        xlim=c(135,160),ylim=c(-39,-24),xlab="",ylab="")
  par(xpd = NA)
  contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,t(ratio[,,4,1]*Useful$mask),lev=seq(0,21),col=rich.colors(21),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "") 
  dev.off()
}

####### Now, plot ECL ratio


names2=dimnames(ratio)[[3]]
for(i in 1:5)
{
  pdf(file=paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Figures/Ratio_total_",names2[i],".pdf",sep=""),width=7,height=4)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  tmp=t(ratio[,,i,1]*Useful$mask)
  tmp[tmp>0.5]=0.525
  image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,0.55,0.05),col=rich.colors(11),frame.plot=F,
        xlim=c(135,160),ylim=c(-39,-24),xlab="",ylab="")
  par(xpd = NA)
  contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,t(ratio[,,4,1]*Useful$mask),lev=seq(0,0.55,0.05),col=rich.colors(11),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "", key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%"),at=seq(0,0.5,0.1))) 
  dev.off()
}

pdf(file=paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Figures/Ratio_9006_total_panel.pdf",sep=""),width=10,height=4,pointsize=12)
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))
for(i in c(2,4,5))
{
  par(mar=c(2,3,3,1))
  tmp=t(ratio_9006[,,i,1]*Useful$mask)
  tmp[tmp>0.5]=0.525
  image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,0.55,0.05),col=rich.colors(11),main=names2[i],frame.plot=F,
        xlim=c(145,155),ylim=c(-39,-24),xlab="",ylab="",cex.axis=1.5)
  contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
}
ColorBar(seq(0,0.55,0.05),rich.colors(11),T,2,c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%")) 
dev.off()

names2=dimnames(ratio)[[3]]
for(i in 1:5)
{
  pdf(file=paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Figures/Ratio_9009_25mm_",names2[i],".pdf",sep=""),width=7,height=4)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  tmp=t(ratio_9009[,,i,3]*Useful$mask)
  tmp[tmp>0.5]=0.525
  image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,0.55,0.05),col=rich.colors(11),frame.plot=F,
        xlim=c(135,160),ylim=c(-39,-24),xlab="",ylab="")
  par(xpd = NA)
  contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,t(ratio[,,4,1]*Useful$mask),lev=seq(0,0.55,0.05),col=rich.colors(11),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "", key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%"),at=seq(0,0.5,0.1))) 
  dev.off()
}

stat<-matrix(0,5,4)
for(i in 1:5)
  for(j in 1:4)
  {
    stat[i,j]=mean(eclrain_9009[,,i,j]*UsefulE$mask,na.rm=T)
  }

###### Comparing WRF versions
load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_CMIP-WRF_proj100_rad2cv0.8.RData")
ratio2=apply(ECLrain37,c(1,2,4,5,6),sum)/apply(rain,c(1,2,4,5,6),sum)
library(abind)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrain_NCEP-WRF_proj100_rad2cv1.RData")
ratio2=abind(ratio2,apply(ECLrain37,c(1,2,4,5),sum)/apply(rain,c(1,2,4,5),sum),along=3)
ratio3<-array(0,c(691,886,5,3,4))
dimnames(ratio3)[[3]]<-dimnames(ratio2)[[3]]
dimnames(ratio3)[[4]]<-dimnames(ratio2)[[4]]
dimnames(ratio3)[[5]]<-dimnames(ratio2)[[5]]

for(i in 1:5)
  for(j in 1:3)
    for(k in 1:4)
{
  a=as.vector(ratio2[,,i,j,k])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratio3[,,i,j,k]=t(b$z)
    }

stat<-array(0,c(5,3,4))
for(i in 1:5)
  for(j in 1:3)
    for(k in 1:4) stat[i,j,k]=mean(ratio3[,,i,j,k]*UsefulE$mask,na.rm=T)
