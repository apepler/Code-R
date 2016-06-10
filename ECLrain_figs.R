rm(list=ls())
setwd('~/Documents/ECLs')
#load("ECLrain_NCEP_8008_fix.RData")
#load("ECLrain_MLDB_8006_fix.RData")
#load("ECLrain_WRF_8008.RData")
prop=apply(ECL_all,c(1,2),sum)/apply(ECL_all+nECL_all,c(1,2),sum)
rm(ECL_all,nECL_all)
propW=apply(ECL_warm,c(1,2),sum)/apply(ECL_warm+nECL_warm,c(1,2),sum)
rm(ECL_warm,nECL_warm)
propC=apply(ECL_cool,c(1,2),sum)/apply(ECL_cool+nECL_cool,c(1,2),sum)
rm(ECL_cool,nECL_cool)

library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

##For WRF only - need to convert to same resolution
library("akima")
load("~/Documents/Data/JE_WRF/Rain_monthly.RData")
latt=as.vector(lat)
lont=as.vector(lon)
a=as.vector(prop)
a[which(is.na(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
prop=t(b$z)
a=as.vector(propW)
a[which(is.na(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
propW=t(b$z)
a=as.vector(propC)
a[which(is.na(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
propC=t(b$z)

#Best ratio = 450 by 500
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(prop*Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(prop*Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "",    key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%","60%","70%"),at=seq(0,0.7,0.1))) 

##Best ratio = 800 by 500
par(plt = c(0.06,0.47,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(propC*Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
text(115,-40,"a)",cex=2)
par(new=TRUE, plt = c(0.47,0.88,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(propW*Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),axes=F,frame.plot=F)
axis(side = 1, at = seq(146,154,2))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
text(115,-40,"a)",cex=2)
par(new = "TRUE",plt = c(0.88,0.93,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(propC*Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "", key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%","60%","70%"),at=seq(0,0.7,0.1)))


##Alternative for using days > 25 mm

#Best ratio = 450 by 500
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(prop*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(10),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(prop*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "",    key.axes = axis(4,labels=c("0%","20%","40%","60%","80%","100%"),at=seq(0,1,0.2))) 

##Best ratio = 800 by 500
par(plt = c(0.06,0.47,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(propC*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(10),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
text(115,-40,"a)",cex=2)
par(new=TRUE, plt = c(0.47,0.88,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(propW*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(10),xlim=c(145,155),ylim=c(-40,-25),axes=F,frame.plot=F)
axis(side = 1, at = seq(146,154,2))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
text(115,-40,"a)",cex=2)
par(new = "TRUE",plt = c(0.88,0.93,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(propC*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "", key.axes = axis(4,labels=c("0%","20%","40%","60%","80%","100%"),at=seq(0,1,0.2)))


###New version - including the locations of the same lows

rm(list=ls())
setwd('~/Documents/ECLs')
load("ECLrain_NCEP_8008_fix.RData")
propN=apply(ECL_all,c(1,2),sum)/apply(ECL_all+nECL_all,c(1,2),sum)
load("ECLrain_WRF_8008.RData")
propW=apply(ECL_all,c(1,2),sum)/apply(ECL_all+nECL_all,c(1,2),sum)
rm(ECL_all,nECL_all,ECL_cool,nECL_cool,ECL_warm,nECL_warm)

library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
readMat('~/Documents/Data/Useful_ECL.mat')->Useful2
mask2<-t(Useful2$mask)
mask2[is.na(mask2)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

##For WRF only - need to convert to same resolution
library("akima")
load("~/Documents/Data/JE_WRF/Rain_monthly.RData")
latt=as.vector(lat)
lont=as.vector(lon)
a=as.vector(propW)
a[which(is.na(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
propW=t(b$z)

source('~/Documents/R/ECLextract.R')
ncep<-read.csv('CSV/ECLfixes_aus_NCEP1_all.csv')
wrf<-read.csv('CSV/ECLfixes_aus_WRF25_all.csv')

lats=seq(-48,-10,2)
lons=seq(110,180,2)
locN<-locW<-matrix(0,length(lats),length(lons))
for(i in 1:length(lons))
  for(j in 1:length(lats))
  {
    I=which(ncep[,5]>=lons[i]-1 & ncep[,5]<lons[i]+1 & ncep[,6]>=lats[j]-1 & ncep[,6]<lats[j]+1)
    locN[j,i]=length(I)
    I=which(wrf[,5]>=lons[i]-1 & wrf[,5]<lons[i]+1 & wrf[,6]>=lats[j]-1 & wrf[,6]<lats[j]+1)
    locW[j,i]=length(I)
  }
locN=locN/29
locW=locW/29

I=which(locN>9.95)
locN[I]=9.95 
I=which(locW>9.95)
locW[I]=9.95 
bb=seq(0,10)

cols=gray(seq(1,0.1,-0.1))
plot.new()
image(lons,lats,t(locN),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)),xlim=c(145,180),ylim=c(-45,-20),frame.plot=F)
par(new=T)
filled.contour3(Useful$x,Useful$y,t(propN*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(n=10),xlim=c(145,180),ylim=c(-45,-20),frame.plot=F)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
plot.new()
image(lons,lats,t(locW),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)),xlim=c(145,180),ylim=c(-45,-20),frame.plot=F)
par(new=T)
filled.contour3(Useful$x,Useful$y,t(propW*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(n=10),xlim=c(145,180),ylim=c(-45,-20),frame.plot=F)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)


###############I just want to do the NCEP1 UM equivalent of my WRF

###So, I need a nice csv file with day1, day2, and csv for each dataset - I'll prep this in excel
rm(list=ls())
setwd('~/Documents/ECLs')
ECLs<-read.csv('CSV/ECLdays_umelb_ncep1_default.csv',header=T,sep=",")

##Let's just do the top 37/year = top 1110
a=order(ECLs[,3],decreasing=TRUE)
ECLs2=ECLs[sort(a[1:1073]),]

##V2 - all days with local csv >= 0.25
#ECLs2=ECLs[ECLs[,3]>=0.25,]

date1=cbind(floor(ECLs2[,1]/10000),floor((ECLs2[,1] %% 10000)/100), ECLs2[,1]%%100, floor(ECLs2[,1]/100000))
date2=cbind(floor(ECLs2[,2]/10000),floor((ECLs2[,2] %% 10000)/100), ECLs2[,2]%%100, floor(ECLs2[,2]/100000))

ECLrain<-ECLcool<-array(0,dim=c(691,886,29))
ECLwarm<-array(0,dim=c(691,886,28))

for(i in 1:length(date1[,1]))
{
  pos=date2[i,1]-1979
  fname<-paste('/media/Seagate Expansion Drive/Data/daily rainfall/rainfall_',date2[i,4],'0-',date2[i,4],'9/rainfall-',date2[i,1],'/r',date2[i,1],sprintf("%2.2i",date2[i,2]),sprintf("%2.2i",date2[i,3]),'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->rain
  as.matrix(rain)->rain
  rain[rain<0]=0
  rain<-rain[nrow(rain):1,]
  
  ##Add it to year, index is yy-1979
  ECLrain[,,pos]=ECLrain[,,pos]+rain
  
  ##Add it to warm/cool depending on month
  if(date2[i,2]>=5 & date2[i,2]<=10) { 
    ECLcool[,,pos]=ECLcool[,,pos]+rain 
  } else if(date2[i,2]>=11 & date2[i,1]<2008) {
    ECLwarm[,,pos]=ECLwarm[,,pos]+rain 
  } else if(date2[i,2]<=4 & date2[i,1]>1980) {
    ECLwarm[,,pos-1]=ECLwarm[,,pos-1]+rain }  ##i.e. add to previous year, no initial year
}

#save(date1,date2,ECLs2,ECLrain,ECLcool,ECLwarm,file="ECLrain_ncep1_cv25.RData")

###Part 2 - compare to total rain

#rm(list=ls())
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

years=seq(1980,2008)
rainA<-rainC<-array(0,dim=c(691,886,29))
rainW<-array(0,dim=c(691,886,28))
for(i in 1:length(years))
  for(j in 1:12)
{
    fname<-paste('/media/Seagate Expansion Drive/Data/monthly rainfall/',years[i],sprintf("%2.2i",j),'.grid',sep="")
    read.table(fname, sep="",skip=6,nrows=691)->rain
    as.matrix(rain)->rain
    rain[rain<0]=0
    rain<-rain[nrow(rain):1,]
    
    rainA[,,i]=rainA[,,i]+rain
    
    if(j>=5 & j<=10) { 
      rainC[,,i]=rainC[,,i]+rain 
    } else if(j>=11 & years[i]<2008) {
      rainW[,,i]=rainW[,,i]+rain 
    } else if(j<=4 & years[i]>1980) {
      rainW[,,i-1]=rainW[,,i-1]+rain  }  ##i.e. add to previous year, no initial year
}

propA=apply(ECLrain,c(1,2),sum)/apply(rainA,c(1,2),sum)
propW=apply(ECLwarm,c(1,2),sum)/apply(rainW,c(1,2),sum)
propC=apply(ECLcool,c(1,2),sum)/apply(rainC,c(1,2),sum)

##Best ratio = 450 by 500
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(propA*Useful$mask),lev=seq(0,0.5,0.05),col=rich.colors(10),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.82,0.87,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(propA*Useful$mask),lev=seq(0,0.5,0.05),col=rich.colors(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "",    key.axes = axis(4,labels=c("0%","10%","20%","30%","40%","50%","60%","70%"),at=seq(0,0.7,0.1))) 

##Best ratio = 800 by 500
par(plt = c(0.06,0.47,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(propC*Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),frame.plot=F)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
text(115,-40,"a)",cex=2)
par(new=TRUE, plt = c(0.47,0.88,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(propW*Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim=c(145,155),ylim=c(-40,-25),axes=F,frame.plot=F)
axis(side = 1, at = seq(146,154,2))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
text(115,-40,"a)",cex=2)
par(new = "TRUE",plt = c(0.88,0.93,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(propC*Useful$mask),lev=seq(0,0.7,0.05),col=rich.colors(14),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "", key.axes = axis(4,labels=c("0%","20%","40%","60%","80%","100%"),at=seq(0,1,0.2)))


