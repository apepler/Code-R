rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')##Plotty stuff   
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
library("akima")
library("abind")
library("RNetCDF")
##Run for all years/months 1979-2009
years=seq(1979,2009)
months=seq(1,12)
load("DailyGDI.RData")
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(20)
cm[9:12]="white"
source('~/Documents/R/corr.plot.R')

load("Rain_monthly.RData")
load("AWAP_monthly.RData")
latt=as.vector(lat)
lont=as.vector(lon)

##Still want to compare to the AWAP averages 1979-2009
##Make a loop that does for every month
for(j in 1:12)
{
  WRFm<-array(0,dim=c(215,144,length(years)))
  AWAPm<-array(0,dim=c(886,691,length(years)))
  for(i in 1:length(years))
  {
    I=which(time2[,2]==years[i] & time2[,3]==j)
    WRFm[,,i]=apply(monthR[,,I],c(1,2),sum)
    AWAPm[,,i]=apply(AWAP[,,I],c(1,2),sum)
  }
  a=as.vector(apply(WRFm,c(1,2),mean))
  b=interp(lont,latt,a,Useful$x,Useful$y)
  meanWRF=b$z
  meanAWAP=apply(AWAPm,c(1,2),mean,na.rm=T)
  diff=meanWRF-meanAWAP
  diff2=diff/meanAWAP
  diff2[meanAWAP==0]=NaN
  
  fout1=paste("Plots/rain/WRF_averain_",j,".tiff",sep="")
  fout2=paste("Plots/rain/WRF_diffAWAP_",j,".tiff",sep="")
  fout3=paste("Plots/rain/WRF_diffAWAP_prop_",j,".tiff",sep="")
  
  tiff(file=fout1, height=450, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,meanWRF,levels=seq(0,500,50),col=rainbow(10))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,meanA2,levels=seq(0,500,50),col=rainbow(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
  dev.off()
  
  tiff(file=fout2, height=450, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,diff*t(Useful$mask),levels=seq(-200,200,40),col=rainbow(10))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,diff*t(Useful$mask),levels=seq(-200,200,40),col=rainbow(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
  dev.off()
  
  tiff(file=fout3, height=450, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,diff2*t(Useful$mask),levels=seq(-2,2,0.2),col=rainbow(20))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,diff2*t(Useful$mask),levels=seq(-2,2,0.2),col=rainbow(20),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
  dev.off()
}

##Look at proportion on E/W days
load("Rain_daily.RData")
latt=as.vector(lat)
lont=as.vector(lon)
rm(rain)
read.csv("DailyGDI.csv",sep=";")->gdi
I=which(time[,2]>=2 | time[,2]==12)
gdi=gdi[I,]
rain2=rain2[,,I]
wrfE<-wrfW<-wrfEa<-wrfWa<-wrfEr<-wrfWr<-matrix(0,215,144)
for(i in 1:length(gdi[,1]))
{
  if(gdi[i,4]>0) wrfE=wrfE+rain2[,,i] else wrfW=wrfW+rain2[,,i]
  if(gdi[i,5]>0) wrfEa=wrfEa+rain2[,,i] else wrfWa=wrfWa+rain2[,,i]
  if(is.na(gdi[i,8])==F) 
     {
       if(gdi[i,8]>0) wrfEr=wrfEr+rain2[,,i] else wrfWr=wrfWr+rain2[,,i]
     }
}
prop=wrfE/(wrfE+wrfW)
propa=wrfEa/(wrfEa+wrfWa)
propr=wrfEr/(wrfEr+wrfWr)
rel=(wrfE/length(which(gdi[,4]>0)))/(wrfW/length(which(gdi[,4]<=0)))
rela=(wrfEa/length(which(gdi[,5]>0)))/(wrfWa/length(which(gdi[,5]<=0)))
relr=(wrfEr/length(which(gdi[,8]>0)))/(wrfWr/length(which(gdi[,8]<=0)))

a=as.vector(prop)
a[which(is.na(a) | is.infinite(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
mat=b$z
tiff(file="Plots/rain/WRF_propE_djf.tiff", height=450, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,mat,levels=seq(0,1,0.1),col=rainbow(10))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,mat,levels=seq(0,1,0.1),col=rainbow(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
dev.off()

load('AWAP_props_djf.RData')
a=as.vector(propa)
a[which(is.na(a) | is.infinite(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
mat=b$z
diff=mat/propA
anom=mat-propA
tiff(file="Plots/rain/WRF_propEa_djf_anom.tiff", height=450, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,anom*t(Useful$mask),levels=seq(-0.25,0.25,0.05),col=rainbow(10))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,anom*t(Useful$mask),levels=seq(-0.25,0.25,0.05),col=rainbow(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
dev.off()

##AWAP GDI
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
load('DailyGDI.RData')
I=which(GDI[,2]>=9 & GDI[,2]<=11)
time=GDI[I,1:3]
yr2=floor(time[,1]/10)
rm(GDI)
read.csv("DailyGDI.csv",sep=";")->gdi
gdi=gdi[I,]
awapE<-awapW<-matrix(0,886,691)

for(i in 1:length(time[,1]))
  if(is.na(gdi[i,8])==F)     
  {
    fin=paste('/media/Seagate Expansion Drive/daily rainfall/rainfall_',yr2[i],'0-',yr2[i],'9/rainfall-',time[i,1],'/r',time[i,1],sprintf("%02d",time[i,2]),sprintf("%02d",time[i,3]),'.txt',sep="")
    read.table(fin, sep="",skip=6,nrows=691)->data
    as.matrix(data)->data
    data[data<0]=NaN
    data<-t(data[nrow(data):1,])
    if(gdi[i,8]>0) awapE=awapE+data else awapW=awapW+data
  }

#save(awapE,awapW,gdi,time,Useful,file="AWAP_EWrain.RData")

##All the initial load/organise stuff
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
source('~/Documents/R/WRF_corrs.R')
##Requires the data in a standard 31x12 matrix
load('STR.RData')
STR=array(0,dim=c(31,12,6))
yind=time2[,2]-1979
for(i in 1:length(yind))
{
  STR[yind[i],time2[i,3],1]=STRd[i,3]
  STR[yind[i],time2[i,3],2]=STRd[i,4]
  STR[yind[i],time2[i,3],3]=STRn[i,3]
  STR[yind[i],time2[i,3],4]=STRn[i,4]
  STR[yind[i],time2[i,3],5]=STRw[i,3]
  STR[yind[i],time2[i,3],6]=STRw[i,4]
}

indnames=c('STRLd','STRId','STRLn','STRIn','STRLw','STRIw')
for(i in 1:6)
{
  gen.corr(c(1,12),'ann',STR[,,i],indnames[i])
  gen.corr(c(3,5),'mam',STR[,,i],indnames[i])
  gen.corr(c(6,8),'jja',STR[,,i],indnames[i])
  gen.corr(c(9,11),'son',STR[,,i],indnames[i])
  gen.corr(c(12,2),'djf',STR[,,i],indnames[i])
  gen.corr(c(5,10),'cool',STR[,,i],indnames[i])
  gen.corr(c(11,4),'warm',STR[,,i],indnames[i])
}


##Re-do some figures for the SEA region
read.table('~/Documents/Timeseries/gdi.txt',header=T,sep="")->gdi
GDIr=gdi[86:116,2:13]
GDIm=matrix(0,length(years),length(months))
for(i in 1:length(years))
  for(j in 1:length(months))
  {
    I=which(GDI[,1]==years[i] & GDI[,2]==j)
    GDIm[i,j]=mean(GDI[I,4])    
  }

seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7)
{
  gen.corr2(snum[i,],seasons[i],GDIr,'GDIr')
  gen.corr2(snum[i,],seasons[i],GDIm,'GDIw')
}

load("AWAP_EWrain.RData")
read.csv("DailyGDI.csv",sep=";")->gdi
propA=awapE/(awapE+awapW)
relA=(awapE/length(which(gdi[,8]>0)))/(awapW/length(which(gdi[,8]<=0)))
gen.propA(propA,relA,'AWAP','E','ann')

load("Rain_daily.RData")
rm(rain)
for(i in 1:7) gen.propW(snum[i,],seasons[i],rain2)



