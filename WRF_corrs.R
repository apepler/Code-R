##Will only work atm if we've done all the normal load stuff. 
setwd('~/Documents/Data/JE_WRF')
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')##Plotty stuff   
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
readMat('~/Documents/Data/mask_escci2.mat')->escci
escci<-t(escci$mask)
escci[is.na(escci)]=0
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
rm(annA,annR,coolA,coolR,warmA,warmR,monthM)

##First - extract the seasonal totals
##Next, make the corresponding matrix of GDIm/GDIr
corr.plot <- function (smon,sname)
{
if(smon[2]>smon[1])
{
  srainW<-array(0,dim=c(215,144,31))
  srainA<-array(0,dim=c(886,691,31))
  for(i in 1:31)
  {
    I=which(time2[,2]==years[i] & time2[,3]>=smon[1] & time2[,3]<=smon[2])
    srainW[,,i]=apply(monthR[,,I],c(1,2),sum,na.rm=T)
    srainA[,,i]=apply(AWAP[,,I],c(1,2),sum,na.rm=T)
  }  
  GDI=matrix(0,31,2)
  GDI[,1]=rowMeans(GDIm[,smon[1]:smon[2]])
  GDI[,2]=rowMeans(GDIr[,smon[1]:smon[2]])
}
else
{
  srainW<-array(0,dim=c(215,144,30))
  srainA<-array(0,dim=c(886,691,30))
  for(i in 1:30)
  {
    I=which((time2[,2]==years[i] & time2[,3]>=smon[1]) | (time2[,2]==years[i+1] & time2[,3]<=smon[2]))
    srainW[,,i]=apply(monthR[,,I],c(1,2),sum,na.rm=T)
    srainA[,,i]=apply(AWAP[,,I],c(1,2),sum,na.rm=T)
  }
  GDI=matrix(0,30,2)
  GDI[,1]=rowMeans(cbind(GDIm[1:30,smon[1]:12],GDIm[2:31,1:smon[2]]))
  GDI[,2]=rowMeans(cbind(GDIr[1:30,smon[1]:12],GDIr[2:31,1:smon[2]]))
}

corrW<-corrR<-matrix(0,215,144)
for(i in 1:215)
  for(j in 1:144)
  {
    corrW[i,j]<-cor(srainW[i,j,],GDI[,1])
    corrR[i,j]<-cor(srainW[i,j,],GDI[,2])
  }
a=as.vector(corrW)
a[which(is.na(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
corW=b$z
a=as.vector(corrR)
a[which(is.na(a))]=0
b=interp(lont,latt,a,Useful$x,Useful$y)
corR=b$z

corA<-matrix(0,886,691)
for(i in 1:886)
  for(j in 1:691)
  {
    corA[i,j]<-cor(srainA[i,j,],GDI[,2])
  }

fout1=paste("Plots/rain/WRF_corGDIw_",sname,".tiff",sep="")
fout2=paste("Plots/rain/WRF_corGDIr_",sname,".tiff",sep="")
fout3=paste("Plots/rain/AWAP_corGDI_",sname,".tiff",sep="")

tiff(file=fout1, height=450, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,corW,levels=seq(-1,1,0.1),col=cm)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,corW,levels=seq(-1,1,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
dev.off()
tiff(file=fout2, height=450, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,corR,levels=seq(-1,1,0.1),col=cm)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,corR,levels=seq(-1,1,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
dev.off()
tiff(file=fout3, height=450, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,corA*t(Useful$mask),levels=seq(-1,1,0.1),col=cm)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,corA,levels=seq(-1,1,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
dev.off()
}

##Does both WRF and AWAP corrs with an index for a season
gen.corr <- function (smon,sname,index,iname)
{
  if(smon[2]>smon[1])
  {
    srainW<-array(0,dim=c(215,144,31))
    srainA<-array(0,dim=c(886,691,31))
    for(i in 1:31)
    {
      I=which(time2[,2]==years[i] & time2[,3]>=smon[1] & time2[,3]<=smon[2])
      srainW[,,i]=apply(monthR[,,I],c(1,2),sum,na.rm=T)
      srainA[,,i]=apply(AWAP[,,I],c(1,2),sum,na.rm=T)
    }  
    ind=rowMeans(index[,smon[1]:smon[2]])
  }
  else
  {
    srainW<-array(0,dim=c(215,144,30))
    srainA<-array(0,dim=c(886,691,30))
    for(i in 1:30)
    {
      I=which((time2[,2]==years[i] & time2[,3]>=smon[1]) | (time2[,2]==years[i+1] & time2[,3]<=smon[2]))
      srainW[,,i]=apply(monthR[,,I],c(1,2),sum,na.rm=T)
      srainA[,,i]=apply(AWAP[,,I],c(1,2),sum,na.rm=T)
    }
    ind=rowMeans(cbind(index[1:30,smon[1]:12],index[2:31,1:smon[2]]))
  }
  
  corrW<-matrix(0,215,144)
  for(i in 1:215)
    for(j in 1:144)
    {
      corrW[i,j]<-cor(srainW[i,j,],ind)
    }
  a=as.vector(corrW)
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  corW=b$z
  
  corA<-matrix(0,886,691)
  for(i in 1:886)
    for(j in 1:691)
    {
      corA[i,j]<-cor(srainA[i,j,],ind)
    }
  
  fout1=paste("Plots/rain/WRF_cor",iname,"_",sname,".tiff",sep="")
  fout2=paste("Plots/rain/AWAP_cor",iname,"_",sname,".tiff",sep="")
  
  tiff(file=fout1, height=450, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,corW,levels=seq(-1,1,0.1),col=cm)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,corW,levels=seq(-1,1,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
  dev.off()
  tiff(file=fout2, height=450, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,corA*t(Useful$mask),levels=seq(-1,1,0.1),col=cm)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,corA,levels=seq(-1,1,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
  dev.off()
}


##Restricted to SEA
gen.corr2 <- function (smon,sname,index,iname)
{
  if(smon[2]>smon[1])
  {
    srainW<-array(0,dim=c(215,144,31))
    srainA<-array(0,dim=c(886,691,31))
    for(i in 1:31)
    {
      I=which(time2[,2]==years[i] & time2[,3]>=smon[1] & time2[,3]<=smon[2])
      srainW[,,i]=apply(monthR[,,I],c(1,2),sum,na.rm=T)
      srainA[,,i]=apply(AWAP[,,I],c(1,2),sum,na.rm=T)
    }  
    ind=rowMeans(index[,smon[1]:smon[2]])
  }
  else
  {
    srainW<-array(0,dim=c(215,144,30))
    srainA<-array(0,dim=c(886,691,30))
    for(i in 1:30)
    {
      I=which((time2[,2]==years[i] & time2[,3]>=smon[1]) | (time2[,2]==years[i+1] & time2[,3]<=smon[2]))
      srainW[,,i]=apply(monthR[,,I],c(1,2),sum,na.rm=T)
      srainA[,,i]=apply(AWAP[,,I],c(1,2),sum,na.rm=T)
    }
    ind=rowMeans(cbind(index[1:30,smon[1]:12],index[2:31,1:smon[2]]))
  }
  
  corrW<-matrix(0,215,144)
  for(i in 1:215)
    for(j in 1:144)
    {
      corrW[i,j]<-cor(srainW[i,j,],ind)
    }
  a=as.vector(corrW)
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  corW=b$z
  
  corA<-matrix(0,886,691)
  for(i in 1:886)
    for(j in 1:691)
    {
      corA[i,j]<-cor(srainA[i,j,],ind)
    }
  
  fout1=paste("Plots/rain/WRF_cor",iname,"_",sname,"_sea.tiff",sep="")
  fout2=paste("Plots/rain/AWAP_cor",iname,"_",sname,"_sea.tiff",sep="")
  
  tiff(file=fout1, height=300, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,corW*t(Useful$mask),levels=seq(-1,1,0.1),col=cm,xlim=c(130,160),ylim=c(-40,-25))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  contour(Useful$x,Useful$y,escci,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,corW,levels=seq(-1,1,0.1),col=cm,xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
  dev.off()
  tiff(file=fout2, height=300, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,corA*t(Useful$mask),levels=seq(-1,1,0.1),col=cm,xlim=c(130,160),ylim=c(-40,-25))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  contour(Useful$x,Useful$y,escci,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,corA,levels=seq(-1,1,0.1),col=cm,xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
  dev.off()
}

##Restricted to SEA, just for AWAP
gen.corrA <- function (smon,sname,index,iname)
{
  if(smon[2]>smon[1])
  {
    srainA<-array(0,dim=c(886,691,31))
    for(i in 1:31)
    {
      I=which(time2[,2]==years[i] & time2[,3]>=smon[1] & time2[,3]<=smon[2])
      srainA[,,i]=apply(AWAP[,,I],c(1,2),sum,na.rm=T)
    }  
    ind=rowMeans(index[,smon[1]:smon[2]])
  }
  else
  {
    srainA<-array(0,dim=c(886,691,30))
    for(i in 1:30)
    {
      I=which((time2[,2]==years[i] & time2[,3]>=smon[1]) | (time2[,2]==years[i+1] & time2[,3]<=smon[2]))
      srainA[,,i]=apply(AWAP[,,I],c(1,2),sum,na.rm=T)
    }
    ind=rowMeans(cbind(index[1:30,smon[1]:12],index[2:31,1:smon[2]]))
  }
  
  corA<-matrix(0,886,691)
  for(i in 1:886)
    for(j in 1:691)
    {
      corA[i,j]<-cor(srainA[i,j,],ind)
    }
  
  fout2=paste("Plots/rain/AWAP_cor",iname,"_",sname,"_sea.tiff",sep="")
  tiff(file=fout2, height=300, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,corA*t(Useful$mask),levels=seq(-1,1,0.1),col=cm,xlim=c(130,160),ylim=c(-40,-25))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  contour(Useful$x,Useful$y,escci,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,corA,levels=seq(-1,1,0.1),col=cm,xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
  dev.off()
}

##PropEW, restricted to SEA - AWAP version - including all the loading of data
gen.propA <- function (smon,sname,ind,iname)
{
  if(smon[2]>smon[1]) I=which(time[,2]>=smon[1] & time[,2]<=smon[2]) 
  else I=which(time[,2]>=smon[1] | time[,2]<=smon[2])
  ind2=ind[I]
  tt=time[I,]
  dec=floor(tt[,1]/10)
  dd=tt[,1]*10000+tt[,2]*100+tt[,3]
  awapE<-awapW<-matrix(0,886,691)
  for(i in 1:length(dec))
  {
    ##Includes loading all the AWAP data
    fname<-paste('/media/Seagate Expansion Drive/daily rainfall/rainfall_',dec[i],'0-',dec[i],'9/rainfall-',tt[i,1],'/r',dd[i],'.txt',sep="")
    read.table(fname, sep="",skip=6,nrows=691)->data
    as.matrix(data)->data
    data[data<0]=0
    data<-t(data[nrow(data):1,])   
    if(ind2[i]>0) awapE=awapE+data else awapW=awapW+data
  }
  prop=awapE/(awapE+awapW)
  rel=(awapE/length(which(ind2>0)))/(awapW/length(which(ind2<=0)))
  
  fout1=paste("Plots/rain/AWAP_propE_",iname,"_",sname,"_sea.tiff",sep="")
  fout2=paste("Plots/rain/AWAP_relE_",iname,"_",sname,"_sea.tiff",sep="")
  tiff(file=fout1, height=300, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,prop*t(Useful$mask),levels=seq(0,1,0.05),col=rainbow(20),xlim=c(130,160),ylim=c(-40,-25))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  contour(Useful$x,Useful$y,escci,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,prop,levels=seq(0,1,0.05),col=rainbow(20),xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
  dev.off()
  tiff(file=fout2, height=300, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,rel*t(Useful$mask),levels=seq(0,3,0.2),col=rainbow(15),xlim=c(130,160),ylim=c(-40,-25))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  contour(Useful$x,Useful$y,escci,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,rel,levels=seq(0,3,0.2),col=rainbow(15),xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
  dev.off()
}

##PropEW, restricted to SEA - WRF version (based on daily rain), need to load first
##Now generalised that "gdi" can be any index of days of E/W
gen.propW <- function (smon,sname,rain,ind,iname)
{
  if(smon[2]>smon[1]) I=which(time[,2]>=smon[1] & time[,2]<=smon[2]) 
  else I=which(time[,2]>=smon[1] | time[,2]<=smon[2])
  rain2=rain[,,I]
  ind2=ind[I]
  wrfE<-wrfW<-matrix(0,215,144)
  for(i in 1:length(ind2))
    if(is.na(ind2)==F)
      if(ind2[i]>0) wrfE=wrfE+rain2[,,i] else wrfW=wrfW+rain2[,,i]
  prop=wrfE/(wrfE+wrfW)
  rel=(wrfE/length(which(ind2>0)))/(wrfW/length(which(ind2<=0)))
  a=as.vector(prop)
  a[which(is.na(a) | is.infinite(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  propW=b$z
  a=as.vector(rel)
  a[which(is.na(a) | is.infinite(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  relW=b$z
  
  fout1=paste("Plots/rain/WRF_propE_",iname,"_",sname,"_sea.tiff",sep="")
  fout2=paste("Plots/rain/WRF_relE_",iname,"_",sname,"_sea.tiff",sep="")
  
  tiff(file=fout1, height=300, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,propW*t(Useful$mask),levels=seq(0,1,0.05),col=rainbow(20),xlim=c(130,160),ylim=c(-40,-25))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  contour(Useful$x,Useful$y,escci,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,propW,levels=seq(0,1,0.05),col=rainbow(20),xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
  dev.off()
  tiff(file=fout2, height=300, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,relW*t(Useful$mask),levels=seq(0,3,0.2),col=rainbow(15),xlim=c(130,160),ylim=c(-40,-25))
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  contour(Useful$x,Useful$y,escci,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,relW,levels=seq(0,3,0.2),col=rainbow(15),xlim=c(130,160),ylim=c(-40,-25),xlab = "",ylab = "")
  dev.off()
}



##Takes in the months needed, the index (as a latxlonx372 file), and esb (as vec372)
esb.corr <- function (smon,sname,index,iname,rname,esb,lat,lon)
{
  d=dim(index)
  corr<-matrix(0,d[1],d[2])
  if(smon[2]>smon[1])
  {
    ind2<-array(0,dim=c(d[1],d[2],31))
    ESB2<-rep(0,31)
    for(j in 1:length(years))
    {
      I=which(time2[,2]==years[j] & time2[,3]>=smon[1] & time2[,3]<=smon[2])
      ind2[,,j]<-apply(index[,,I],c(1,2),mean,na.rm=T)
      ESB2[j]<-sum(esb[I])
    }
  } else
  {
    ind2<-array(0,dim=c(d[1],d[2],30))
    ESB2<-rep(0,30)
    for(j in 1:30)
    {
      I=which((time2[,2]==years[j] & time2[,3]>=smon[1]) | (time2[,2]==years[j+1] & time2[,3]<=smon[2]))
      ind2[,,j]<-apply(index[,,I],c(1,2),mean,na.rm=T)
      ESB2[j]<-sum(esb[I])
    }
  }
  for(i in 1:d[1]) for(j in 1:d[2]) corr[i,j]=cor(ind2[i,j,],ESB2)
  
  fout=paste("Plots/rain/",rname,"_corESBrain_",iname,"_",sname,".tiff",sep="")
  
  ##Different approaches depending on if NCEP or WRF
  
  if(is.na(dim(lat)[2])) ##Lat and Lon are vectors
  {
    a=1
    if(lat[2]<lat[1]) ##If lat is not ascending, swap
    {
      lat=lat[d[2]:1]
      corr=corr[,d[2]:1]
    }
    tiff(file=fout, height=400, width=600)
    par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
    filled.contour3(lon,lat,corr,levels=seq(-1,1,0.1),col=cm,xlim=c(110,160),ylim=c(-45,-10))
    par(xpd = NA)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
    filled.legend(lon,lat,corr,levels=seq(-1,1,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
    dev.off()
  }
  else ##Lat is an array - must be WRF
  {
    ##In that case, I want to convert to a 0.5 degree res for plotting.
    latt=as.vector(lat)
    lont=as.vector(lon)
    lat2=seq(-50,10,0.5)
    lon2=seq(100,180,0.5)
    a=as.vector(corr)
    a[which(is.na(a))]=0
    b=interp(lont,latt,a,lon2,lat2)
    corr2=b$z
    tiff(file=fout, height=400, width=600)
    par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
    filled.contour3(lon2,lat2,corr2,levels=seq(-1,1,0.1),col=cm,xlim=c(110,160),ylim=c(-45,-10))
    par(xpd = NA)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
    filled.legend(lon2,lat2,corr2,levels=seq(-1,1,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
    dev.off()
  }
}