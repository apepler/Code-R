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
a=open.nc("/srv/ccrc/data45/z3478332/WRF/output/extracted_data/WRF_d01_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
latt=as.vector(lat)
lont=as.vector(lon)

ratioC<-array(0,c(691,886,5,4))
dimnames(ratioC)[[3]]<-c("MLDB","NCEP","ERAI","NCEP-WRF","CMIP-WRF")
dimnames(ratioC)[[4]]=c("Total","Rain>=5mm","Rain>25mm","Rain>50mm")
ratioW<-ratioC

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainC_AWAP_proj100_rad2cv1.RData")
ratioW[,,1,]=apply(eclrainW_mldb[,,21:36,],c(1,2,4),sum)/apply(rainW[,,21:36,],c(1,2,4),sum)
ratioW[,,2,]=apply(eclrainW_ncep[,,21:39,],c(1,2,4),sum)/apply(rainW[,,21:39,],c(1,2,4),sum)
ratioW[,,3,]=apply(eclrainW_erai[,,21:39,],c(1,2,4),sum)/apply(rainW[,,21:39,],c(1,2,4),sum)

rm(rainW,eclrainW_mldb,eclrainW_ncep,eclrainW_erai)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainC_NCEP-WRF_proj100_rad2cv1.RData")
ratio2=apply(ECLrain37C,c(1,2,4,5),sum)/apply(rainC,c(1,2,4,5),sum)
ratio3=apply(ratio2,c(1,2,4),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratioC[,,4,i]=t(b$z)
}
rm(ECLrainC,ECLrain37C,rainC)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainW_NCEP-WRF_proj100_rad2cv1.RData")
ratio2=apply(ECLrain37W,c(1,2,4,5),sum)/apply(rainW,c(1,2,4,5),sum)
ratio3=apply(ratio2,c(1,2,4),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratioW[,,4,i]=t(b$z)
}
rm(ECLrainW,ECLrain37W,rainW)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainC_CMIP-WRF_proj100_rad2cv0.8.RData")
ratio2=apply(ECLrain37C,c(1,2,4,5,6),sum)/apply(rainC,c(1,2,4,5,6),sum)
ratio3=apply(ratio2,c(1,2,5),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratioC[,,5,i]=t(b$z)
}
rm(ECLrainC,ECLrain37C,rainC)

load("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/outputUM/ECLrainW_CMIP-WRF_proj100_rad2cv0.8.RData")
ratio2=apply(ECLrain37W,c(1,2,4,5,6),sum)/apply(rainW,c(1,2,4,5,6),sum)
ratio3=apply(ratio2,c(1,2,5),mean)
for(i in 1:4)
{
  a=as.vector(ratio3[,,i])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ratioW[,,5,i]=t(b$z)
}
rm(ECLrainW,ECLrain37W,rainW)

names1=c("total","5mm","25mm","50mm")
names2=dimnames(ratioC)[[3]]
for(j in 1:4)
for(i in 4:5)
{
  pdf(file=paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Figures/Ratio_",names1[j],"_",names2[i],"_season.pdf",sep=""),width=10,height=4)
  par(plt = c(0.05,0.42,0.1,0.85),las = 1,cex.axis = 1)
  tmp=t(ratioC[,,i,j]*Useful$mask)
  tmp[tmp>0.5]=0.525
  image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,0.55,0.05),col=rich.colors(11),frame.plot=F,
        xlim=c(135,155),ylim=c(-39,-24),xlab="",ylab="",main="May-October")
  par(xpd = NA)
  contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
  par(new="TRUE",plt = c(0.47,0.84,0.1,0.85),las = 1,cex.axis = 1)
  tmp=t(ratioW[,,i,j]*Useful$mask)
  tmp[tmp>0.5]=0.525
  image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,0.55,0.05),col=rich.colors(11),frame.plot=F,
        xlim=c(135,155),ylim=c(-39,-24),xlab="",ylab="",main="November-April")
  par(xpd = NA)
  contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.86,0.91,0.1,0.85),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,t(ratio[,,4,1]*Useful$mask),lev=seq(0,0.55,0.05),col=rich.colors(11),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "", key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%"),at=seq(0,0.5,0.1))) 
  dev.off()
}

ColorBar <- function(brks,cols,vert=T,subsampleg=1,breaklab=brks)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = breaklab[seq(1, length(brks), subsampleg)])
}

for(j in 1:4)
{

pdf(file=paste("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Figures/Ratio_",names1[j],"_panel_warm.pdf",sep=""),width=10,height=4,pointsize=12)
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))
for(i in c(2,4,5))
{
  par(mar=c(2,3,3,1))
  tmp=t(ratioW[,,i,1]*Useful$mask)
  tmp[tmp>0.5]=0.525
  image(Useful$x[1,],Useful$y[1,],tmp,breaks=seq(0,0.55,0.05),col=rich.colors(11),main=names2[i],frame.plot=F,
        xlim=c(145,155),ylim=c(-39,-24),xlab="",ylab="",cex.axis=1.5)
  contour(Useful$x,Useful$y,maskE,add=T,drawlabels=F)
}
ColorBar(seq(0,0.55,0.05),rich.colors(11),T,2,c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%")) 
dev.off()
}

