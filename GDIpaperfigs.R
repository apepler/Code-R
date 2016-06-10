##Figure 1

setwd('~/Documents/Data')
library('RNetCDF')
f1=open.nc('ETOPO2v2c_f4.nc')
lon=var.get.nc(f1,'x')
lat=var.get.nc(f1,'y')
I=which(lat>=(-45) & lat<=(-10))
J=which(lon>=110 & lon<=155)
lat2=lat[I]
lon2=lon[J]
topo=var.get.nc(f1,'z',c(J[1],I[1]),c(length(J),length(I)))

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
readMat('~/Documents/Data/mask_escci2.mat')->escci
mask<-t(Useful$mask)
mask[is.na(mask)]=0
mask2<-t(escci$mask.escci)
mask2[is.na(mask2)]=0
cols=gray(c(1,0.6,0.4,0.2))  																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																														


plot.new()
filled.contour3(lon2,lat2,topo,levels=seq(0,2000,500),col=cols,axes=F)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
#rect(138,-39.3,150,-35,lwd=3)
## Add Sydney, Gayndah, Deniliquin
points(151.21,-33.86,pch=16,cex=1)
points(151.61,-25.63,pch=16,cex=1)
points(144.95,-35.56,pch=16,cex=1)

##Figure 2

rm(list=ls())
setwd('~/Documents/Data')
load('AWAPrain.RData')
read.table('~/Documents/Timeseries/tri.txt',header=T)->tri
TRI<-matrix(NaN,113,3)
TRI[,1]<-tri[,1]
TRI[,2]<-rowSums(tri[,7:11])
TRI[1:112,3]<-rowSums(cbind(tri[1:112,12:13],tri[2:113,2:4]))

corC<-corW<-matrix(0,691,886)
for(i in 1:691)
  for(j in 1:886)
  {
    a=cor.test(warm[i,j,],TRI[1:112,3])
    if(is.na(a$p.value)) corW[i,j]=0 else
       if(a$p.value<=0.05) corW[i,j]=a$estimate else corW[i,j]=0
    a=cor.test(cool[i,j,],TRI[,2])
    if(is.na(a$p.value)) corC[i,j]=0 else
      if(a$p.value<=0.05) corC[i,j]=a$estimate else corC[i,j]=0
  }

##Plot
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("white","black"))
cm=pal(7)
library("R.matlab")
readMat('Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

##Best aspect is 1200 by 500
plot.new()
par(plt = c(0.05,0.45,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(corC*Useful$mask),levels=seq(0,0.7,0.1),col=cm)
text(115,-40,"a)",cex=2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new=TRUE, plt = c(0.5,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(corW*Useful$mask),levels=seq(0,0.7,0.1),col=cm)
text(115,-40,"b)",cex=2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.92,0.95,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(corC*Useful$mask),levels=seq(0,0.7,0.1),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")

##Figure 3

rm(list=ls())
setwd('~/Documents/Data')
read.table('~/Documents/Timeseries/gdi.txt',header=T)->gdi
GDIm=colMeans(gdi[,2:13])
GDIsd=rep(0,12)
for(i in 1:12) GDIsd[i]=sd(GDI[,i+1])
GDIu=GDIm+GDIsd
GDIl=GDIm-GDIsd
library(hydroGOF)
plot.new()
plot(seq(1,12),GDIm,type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-4,5))
plotbandsonly(GDIl,GDIu, bands.col=rgb(0.9,0.9,0.9),ylim=c(-4,5),alpha=0.5)
lines(seq(1,12),GDIm,col="black",lwd=2)
axis(1, at = seq(1,12,1), labels = names(gdi[2:13]))
abline(0,0,col=rgb(0.6,0.6,0.6))

##Or, version 2, with EN/LN

gdi2=matrix(NaN,113,17)
gdi2[,1:11]=as.matrix(gdi[,3:13])
gdi2[1:112,12:17]=as.matrix(gdi[2:113,2:7])
gdi3=matrix(NaN,113,15)
for(i in 1:15) gdi3[,i]=rowMeans(gdi2[,i:(i+2)])
GDI=matrix(0,5,15)
I=which(enso[,2]=="A")
GDI[1,]=colMeans(gdi3[I,],na.rm=T)
GDIsd=rep(0,15)
for(i in 1:15) GDIsd[i]=sd(gdi3[I,i],na.rm=T)
GDI[2,]=GDI[1,]+GDIsd
GDI[3,]=GDI[1,]-GDIsd
read.table('~/Documents/Timeseries/enso.txt',header=T)->enso
I=which(enso[,2]=="EN")
GDI[4,]=colMeans(gdi3[I,])
I=which(enso[,2]=="LN")
GDI[5,]=colMeans(gdi3[I,])

plot.new()
plot(seq(1,15),GDI[1,],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-4,5))
plotbandsonly(GDI[3,],GDI[2,], bands.col=rgb(0.9,0.9,0.9),ylim=c(-4,5),alpha=0.5)
lines(seq(1,15),GDI[1,],col="black",lwd=2)
axis(1, at = seq(1,15,2), labels = c("Mar","May","Jul","Sep","Nov","Jan","Mar","May"))
abline(0,0,col=rgb(0.6,0.6,0.6))
lines(seq(1,15),GDI[4,],col="black",lwd=2,lty=2)
lines(seq(1,15),GDI[5,],col="black",lwd=2,lty=3)
legend("bottomright",c('El Nino','Neutral','La Nina'),col=rep("black",3),lwd=c(2,2,2),lty=c(1,2,3))

##Figure 4a
rm(list=ls())
setwd('~/Documents/Data')
read.csv('GDIlist_cool.csv',sep=";")->GDI

E<-W<-array(0,dim=c(691,886))

for(i in 1:length(GDI[,1]))
{
  if(GDI[i,3]>=2010) fname<-paste('/media/Seagate Expansion Drive/daily rainfall/rainfall_2010/rainfall-',GDI[i,3],'/r',GDI[i,1],'.txt',sep="") 
  else fname<-paste('/media/Seagate Expansion Drive/daily rainfall/rainfall_',GDI[i,4],'0-',GDI[i,4],'9/rainfall-',GDI[i,3],'/r',GDI[i,1],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  if(GDI[i,8]=="E") E=E+data else W=W+data
  rm(data)
}

library("R.matlab")
readMat('Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

prop=E/(E+W)
plot.new()
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(prop*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(10))
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
text(115,-40,"a)",cex=2)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(prop*Useful$mask),lev=seq(0,1,0.1),col=rich.colors(10),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "",    key.axes = axis(4,labels=c("0%","20%","40%","60%","80%","100%"),at=seq(0,1,0.2))) 


##Figure 4b

rm(list=ls())
setwd('~/Documents/Data')
load('AWAPrain.RData')
rm(warm)
source('~/Documents/R/pcor.R')
read.table('~/Documents/Timeseries/gdi.txt',header=T)->gdi
G<-rowMeans(gdi[,7:11])
corC<-matrix(0,691,886)
for(i in 1:691)
  for(j in 1:886)
  {
    a=cor.test(cool[i,j,],G)
    if(is.na(a$p.value)) corC[i,j]=0 else
      if(a$p.value<=0.05) corC[i,j]=a$estimate else corC[i,j]=0
  }

##Plot
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(14)
cm[7:8]="white"
library("R.matlab")
readMat('Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
##Best aspect is 1200 by 500
plot.new()
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(corC*Useful$mask),levels=seq(-0.7,0.7,0.1),col=cm)
text(115,-40,"b)",cex=2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(corC*Useful$mask),levels=seq(-0.7,0.7,0.1),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")


##Figure 5

rm(list=ls())
setwd('~/Documents/Data')
read.table('~/Documents/Timeseries/gdi.txt',header=T)->gdi
read.table('~/Documents/Timeseries/n34.txt',header=T)->n34
read.table('~/Documents/Timeseries/dmi.txt',header=T)->dmi
read.table('~/Documents/Timeseries/soi.txt',header=T)->soi

mCorrs3<-function(x,y)
{
  corrs<-rep(0,12)
  corrs[1]<-cor(rowMeans(cbind(x[1:(length(x[,1])-1),13],x[2:length(x[,1]),2:3])),rowMeans(cbind(y[1:(length(y[,1])-1),13],y[2:length(y[,1]),2:3])))
  corrs[12]<-cor(rowMeans(cbind(x[1:(length(x[,1])-1),12:13],x[2:length(x[,1]),2])),rowMeans(cbind(y[1:(length(y[,1])-1),12:13],y[2:length(y[,1]),2])))
  for(i in 2:11) corrs[i]<-cor(rowMeans(x[1:(length(x[,1])-1),i:(i+2)]),rowMeans(y[1:(length(y[,1])-1),i:(i+2)]))
  return(corrs)
}

mCorrs3b<-function(x,y1,y2)
{
  corrs<-rep(0,12)
  months=1:12
  source('~/Documents/R/pcor.R')
  a<-pcor.test(rowMeans(cbind(x[1:(length(x[,1])-1),13],x[2:length(x[,1]),2:3])),rowMeans(cbind(y1[1:(length(y1[,1])-1),13],y1[2:length(y1[,1]),2:3])),rowMeans(cbind(y2[1:(length(y2[,1])-1),13],y2[2:length(y2[,1]),2:3])))
  corrs[1]<-as.numeric(a[1])
  a<-pcor.test(rowMeans(cbind(x[1:(length(x[,1])-1),12:13],x[2:length(x[,1]),2])),rowMeans(cbind(y1[1:(length(y1[,1])-1),12:13],y1[2:length(y1[,1]),2])),rowMeans(cbind(y2[1:(length(y2[,1])-1),12:13],y2[2:length(y2[,1]),2])))
  corrs[12]<-as.numeric(a[1])
  for(i in 2:11)
  {
    a<-pcor.test(rowMeans(x[1:(length(x[,1])-1),i:(i+2)]),rowMeans(y1[1:(length(y1[,1])-1),i:(i+2)]),rowMeans(y2[1:(length(y2[,1])-1),i:(i+2)]))
    corrs[i]<-as.numeric(a[1])
  }  
  return(corrs)
}

corrS=mCorrs3(-soi,gdi)
corrN_D=mCorrs3b(gdi,n34,dmi)
corrD_N=mCorrs3b(gdi,dmi,n34)

plot.new()
par(mar=c(3,3,2,2))
plot(seq(1,12),corrS,type="o",col="black",xaxt="n", xlab="", ylab="",pch=18,lwd=2,ylim=range(-0.4,0.4))
rect(xleft=-1000,xright=1000,ybottom=-0.186, ytop=0.185, density=100, col=rgb(0.9,0.9,0.9))
lines(seq(1,12),corrS,type="o",col="black",xaxt="n", xlab="", ylab="",pch=18,lwd=2)
lines(seq(1,12),corrN_D,type="o",col=rgb(0.5,0.5,0.5),xaxt="n", xlab="", ylab="",pch=4,lwd=2)
axis(1, at = seq(1,12,1), labels = names(gdi[2:13]))
lines(seq(1,12),corrD_N,type="o",col=rgb(0.5,0.5,0.5),lty=2,xaxt="n", xlab="", ylab="",pch=16,lwd=2)
legend("topleft",c('-SOI','N34 (no DMI)','DMI (no N34)'),col=c("black",rgb(0.5,0.5,0.5),rgb(0.5,0.5,0.5)),lwd=c(2,2,2),lty=c(1,1,2),pch=c(18,4,16))

##Figure 8
rm(list=ls())
setwd('~/Documents/Data')
load('AWAPrain.RData')
rm(warm)
years=seq(1900,2012)
I=which(years>=1958)
cool=cool[,,I]
read.csv('~/Documents/Timeseries/state.csv',sep=";",header=T)->enso
enso=enso[enso[,1]>=1958,]
medR<-apply(cool,c(1,2),median)

EN<-which(enso[,2]=="E")
LN<-which(enso[,2]=="L")
ENn<-which(enso[,2]=="E" & enso[,3]!="P")
ENy<-which(enso[,2]=="E" & enso[,3]=="P")
LNn<-which(enso[,2]=="L" & enso[,3]!="N")
LNy<-which(enso[,2]=="L" & enso[,3]=="N")
Pn<-which(enso[,2]!="E" & enso[,3]=="P")
Nn<-which(enso[,2]!="L" & enso[,3]=="N")

state=list(EN,LN,ENn,LNn,ENy,LNy,Pn,Nn)
labels=c("a) El Nino","b) La Nina","c) El Nino, not pIOD", "d) La Nina, not nIOD","e) Positive IOD + El Nino","f) Negative IOD + La Nina","g) Positive IOD, not El Nino","h) Negative IOD, not La Nina")
fname=c("ENall","LNall","ENn","LNn","Dry","Wet","pIODn","nIODn")

library("R.matlab")
readMat('Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

for (l in 1:8)
{
  I=unlist(state[l])
  abovemed<-matrix(0,691,886)
  for(i in 1:length(I))
  {
    J<-which(cool[,,I[i]]>medR)
    abovemed[J]<-abovemed[J]+1
  }
  if(length(I)<6)
  {
    cm<-rich.colors(length(I)+2)
    cm2<-rich.colors(length(I)+1)
    for(i in 1:(length(I)+1)) cm2[i]=cm[(length(I)+3-i)]
  } else {
    cm<-cm2<-rich.colors(length(I)+1)
    for(i in 1:(length(I)+1)) cm2[i]=cm[(length(I)+2-i)]
  }
  tiff(file=paste("~/Documents/GDI paper/Cool_abovemed2_",fname[l],".tiff",sep=""), height=500, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 2)
  filled.contour3(Useful$x,Useful$y,t(abovemed*Useful$mask),lev=seq(-0.5,length(I)+0.5),col=cm2,main=paste(labels[l]),cex.main=2)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 2)
  filled.legend(Useful$x,Useful$y,t(abovemed*Useful$mask),lev=seq(-0.5,length(I)+0.5),col=cm2,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "") 
  dev.off()
}




##Figure 10

rm(list=ls())
setwd('~/Documents/Data')
load('AWAPrain.RData')
source('~/Documents/R/pcor.R')
read.table('~/Documents/Timeseries/n34.txt',header=T)->n34
read.table('~/Documents/Timeseries/dmi.txt',header=T)->dmi
E<-D<-matrix(NaN,113,3)
E[,1]<-D[,1]<-n34[,1]
E[,2]<-rowSums(n34[,7:11])
E[1:112,3]<-rowSums(cbind(n34[1:112,12:13],n34[2:113,2:4]))
D[,2]<-rowSums(dmi[,7:11])
D[1:112,3]<-rowSums(cbind(dmi[1:112,12:13],dmi[2:113,2:4]))

corC<-corW<-matrix(0,691,886)
for(i in 1:691)
  for(j in 1:886)
  {
    a=cor.test(warm[i,j,],E[1:112,3])
    if(is.na(a$p.value)) corW[i,j]=0 else
      if(a$p.value<=0.05) corW[i,j]=a$estimate else corW[i,j]=0
    a=pcor.test(cool[i,j,],E[,2],D[,2])
    if(is.na(a[2])) corC[i,j]=0 else
      if(a[2]<=0.05) corC[i,j]=as.numeric(a[1]) else corC[i,j]=0
  }

##Plot
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("white","black"))
cm=pal(7)
library("R.matlab")
readMat('Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
cm2=cm;
for(i in 1:7) cm2[i]=cm[8-i]
##Best aspect is 1200 by 500
plot.new()
par(plt = c(0.05,0.45,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(corC*Useful$mask),levels=seq(-0.7,0,0.1),col=cm2)
text(115,-40,"a)",cex=2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new=TRUE, plt = c(0.5,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(corW*Useful$mask),levels=seq(-0.7,0,0.1),col=cm2)
text(115,-40,"b)",cex=2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.92,0.95,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(corC*Useful$mask),levels=seq(-0.7,0,0.1),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm2,xlab = "",ylab = "")


####Ummenhofer
rm(list=ls())
setwd('~/Documents/Data')
load('AWAPrain.RData')
rm(warm)
years=seq(1900,2012)
I=which(years<=2006)
cool=cool[,,I]
read.csv('~/Documents/Timeseries/meyers.csv',header=T)->enso
medR<-apply(cool,c(1,2),median)

EN<-which(enso[,2]=="E")
LN<-which(enso[,2]=="L")
ENn<-which(enso[,2]=="E" & enso[,3]!="P")
ENy<-which(enso[,2]=="E" & enso[,3]=="P")
LNn<-which(enso[,2]=="L" & enso[,3]!="N")
LNy<-which(enso[,2]=="L" & enso[,3]=="N")
Pn<-which(enso[,2]!="E" & enso[,3]=="P")
Nn<-which(enso[,2]!="L" & enso[,3]=="N")

state=list(EN,LN,ENn,LNn,ENy,LNy,Pn,Nn)
labels=c("a) El Nino","b) La Nina","c) El Nino, not pIOD", "d) La Nina, not nIOD","e) Positive IOD + El Nino","f) Negative IOD + La Nina","g) Positive IOD, not El Nino","h) Negative IOD, not La Nina")
fname=c("ENall","LNall","ENn","LNn","Dry","Wet","pIODn","nIODn")

library("R.matlab")
readMat('Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

for (l in 1:8)
{
  I=unlist(state[l])
  abovemed<-matrix(0,691,886)
  for(i in 1:length(I))
  {
    J<-which(cool[,,I[i]]>medR)
    abovemed[J]<-abovemed[J]+1
  }
  if(length(I)<6)
  {
    cm<-rich.colors(length(I)+2)
    cm2<-rich.colors(length(I)+1)
    for(i in 1:(length(I)+1)) cm2[i]=cm[(length(I)+3-i)]
  } else {
    cm<-cm2<-rich.colors(length(I)+1)
    for(i in 1:(length(I)+1)) cm2[i]=cm[(length(I)+2-i)]
  }
  tiff(file=paste("~/GDI paper/Meyers_Cool_abovemed_",fname[l],".tiff",sep=""), height=500, width=600)
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(Useful$x,Useful$y,t(abovemed*Useful$mask),lev=seq(-0.5,length(I)+0.5),col=cm2,main=paste(labels[l]),cex.main=2)
  par(xpd = NA)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(Useful$x,Useful$y,t(abovemed*Useful$mask),lev=seq(-0.5,length(I)+0.5),col=cm2,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "") 
  dev.off()
}






