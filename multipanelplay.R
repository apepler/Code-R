##Multi-panel plots

##Data

rm(list=ls())
setwd('G:/monthly-rainfall')
source('~/R/load.AWAP.season.ff.R')
source('~/R/load.AWAP.wseason.ff.R')
library(ff)
m1<-c("06","07","08","09","10")
m2<-seq(6,10,1)
name<-"JJASO"
rain<-getRainff(m1)
m1a<-c("11","12")
m1b<-c("01","02","03")
m2a<-seq(11,12,1)
m2b<-seq(1,3,1)
name<-"NDJFM"
rain2<-getRainWff(m1a,m1b)

##Correlations

read.table('n34.txt',header=T,sep="")->n34
read.table('dmi.txt',header=T,sep="")->dmi
read.table('tri.txt',header=T,sep="")->tri
read.table('soi.txt',header=T,sep="")->soi
read.table('gdi.txt',header=T,sep="")->gdi

source('~/R/pcor.R')
corrC<-corrW<-matrix(0,691,886)
for(i in 1:691)
  for(j in 1:886)
  {
    a<-pcor.test(rowMeans(n34[31:142,7:11]),rain[i,j,],rowMeans(dmi[31:142,7:11]))
    b<-as.numeric(a[1])
    corrC[i,j]<-b
    #corrC[i,j]<-cor(rowMeans(tri[31:142,7:11]),rain[i,j,])
    corrW[i,j]<-cor(rowMeans(cbind(n34[31:141,12:13],n34[32:142,2:4])),rain2[i,j,])
  }

##Plot
source('~/R/filled.contour3.R')
source('~/R/filled.legend.R')
source('~/R/color.palette.R')
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
filled.contour3(Useful$x,Useful$y,t(corrC*Useful$mask),levels=seq(-0.7,0,0.1),col=cm2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new=TRUE, plt = c(0.5,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(corrW*Useful$mask),levels=seq(-0.7,0,0.1),col=cm2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.92,0.95,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(corrC*Useful$mask),levels=seq(-0.7,0,0.1),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm2,xlab = "",ylab = "")


##Fig3
read.table('~/GDI/Figures/fig3.csv',header=T,sep=",")->fig
tiff(file='~/GDI/Figures/paper_fig3.tif', height=6, width=10, units='cm', res=600, compression='lzw', pointsize=8)
par(mar=c(3,3,2,2))
plot(seq(1,15),fig[,2],type="l",col="black",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=16,lwd=2,ylim=range(-2,3))
axis(1, at = seq(1,15,1), labels = fig[,1],cex.axis=0.75)
lines(seq(1,15),fig[,3],type="l",lty=2,col=rgb(0.7,0.7,0.7),xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=16,lwd=2)
lines(seq(1,15),fig[,4],type="l",col=rgb(0.4,0.4,0.4),xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=16,lwd=2)
abline(0,0,col=rgb(0.7,0.7,0.7))
legend("bottomright",c('El Nino','Neutral','La Nina'),col=c("black",rgb(0.7,0.7,0.7),rgb(0.4, 0.4,0.4)),lwd=c(2,2,2),lty=c(1,2,1))
dev.off()

##Fig4
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
  source('~/R/pcor.R')
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
corrT=mCorrs3(n34[n34[,1]>=1900,],gdi[gdi[,1]>=1900,])
corrN_D=mCorrs3b(gdi[gdi[,1]>=1900,],n34[n34[,1]>=1900,],dmi[dmi[,1]>=1900,])
corrD_N=mCorrs3b(gdi[gdi[,1]>=1900,],dmi[dmi[,1]>=1900,],n34[n34[,1]>=1900,])

tiff(file='~/GDI/Figures/paper_fig4a.tif', height=6, width=10, units='cm', res=600, compression='lzw', pointsize=8)
par(mar=c(3,3,2,2))
plot(seq(1,12),corrT,type="o",col="black",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=18,lwd=2,ylim=range(-0.4,0.4))
lines(seq(1,12),corrN_D,type="o",col=rgb(0.7,0.7,0.7),xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=4,lwd=2)
axis(1, at = seq(1,12,1), labels = names(gdi[2:13]),cex.axis=0.75)
lines(seq(1,12),corrD_N,type="o",col=rgb(0.7,0.7,0.7),lty=2,xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=16,lwd=2)
legend("topleft",c('N34','N34 (no DMI)','DMI (no N34)'),col=c("black",rgb(0.7,0.7,0.7),rgb(0.7,0.7,0.7)),lwd=c(2,2,2),lty=c(1,1,2),pch=c(18,4,16))
color <- rgb(190, 190, 190, alpha=160, maxColorValue=255)
rect(xleft=-1000,xright=1000,ybottom=-0.186, ytop=0.185, density=100, col=color)
dev.off()

corrN=mCorrs3(n34[n34[,1]>=1900,],gdi[gdi[,1]>=1900,])
corrS=mCorrs3(soi[soi[,1]>=1900,],gdi[gdi[,1]>=1900,])
tiff(file='~/GDI/Figures/paper_fig4d.tif', height=6, width=10, units='cm', res=600, compression='lzw', pointsize=8)
par(mar=c(3,3,2,2))
plot(seq(1,12),corrS,type="o",col="black",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=18,lwd=2,ylim=range(-0.4,0.4))
lines(seq(1,12),corrN,type="o",col=rgb(0.7,0.7,0.7),xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=4,lwd=2)
axis(1, at = seq(1,12,1), labels = names(gdi[2:13]),cex.axis=0.75)
legend("topleft",c('SOI','N34'),col=c("black",rgb(0.7,0.7,0.7)),lwd=c(2,2,2),lty=c(1,1),pch=c(18,4))
color <- rgb(190, 190, 190, alpha=160, maxColorValue=255)
rect(xleft=-1000,xright=1000,ybottom=-0.186, ytop=0.185, density=100, col=color)
dev.off()


##Re-do figure 8. Black & white too. 

rm(list=ls())
setwd('G:/monthly-rainfall')
source('~/R/load.AWAP.season.ff.R')
source('~/R/load.AWAP.wseason.ff.R')
library(ff)
m1<-c("06","07","08","09","10")
m2<-seq(6,10,1)
name<-"JJASO"
rain<-getRainff(m1)

source('~/R/filled.contour3.R')
source('~/R/filled.legend.R')
source('~/R/color.palette.R')
pal <- color.palette(c("white","black"))
library("R.matlab")
readMat('Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

cm=pal(7)
cm2=cm;
for(i in 1:7) cm2[i]=cm[8-i]

read.table('enso.csv',header=T,sep=",")->enso
I<-which(enso[,1]<1958)
enso<-enso[-I,]
#rain<-rain[,,58:112]
state=rep(0,length(enso[,1]))

I<-which(enso[,2]!="EN" & enso[,3]=="P")
state[I]=1
I<-which(enso[,2]=="EN" & enso[,3]=="P")
state[I]=3
I<-which(enso[,2]=="EN" & enso[,3]!="P")
state[I]=5
I<-which(enso[,2]!="LN" & enso[,3]=="N")
state[I]=2
I<-which(enso[,2]=="LN" & enso[,3]=="N")
state[I]=4
I<-which(enso[,2]=="LN" & enso[,3]!="N")
state[I]=6

medR<-matrix(0,691,886)
for(i in 1:691)
{
  medR[i,]<-apply(rain[i,,],1,median)
}
labels=c("Positive IOD, not El Nino","Negative IOD, not La Nina","Positive IOD + El Nino","Negative IOD + La Nina","El Nino, not pIOD", "La Nina, not nIOD")

loc<-matrix(0,6,4)
loc[1,]=c(0.1,0.4,0.1,0.3)
loc[2,]=c(0.5,0.8,0.1,0.3)
loc[3,]=c(0.1,0.4,0.4,0.6)
loc[4,]=c(0.5,0.8,0.4,0.6)
loc[5,]=c(0.1,0.4,0.7,0.9)
loc[6,]=c(0.5,0.8,0.7,0.9)


cm=pal(10)
cm2=cm;
for(i in 1:10) cm2[i]=cm[10-i]

tiff(file='~/GDI/Figures/paper_fig8.tif', height=15, width=10, units='cm', res=600, compression='lzw', pointsize=8)
for (l in 1:6)
{
  I=which(state==l)
  abovemed<-matrix(0,691,886)
  for(i in 1:length(I))
  {
    J<-which(rain[,,I[i]]>medR)
    abovemed[J]<-abovemed[J]+1
  }
  abovemed=abovemed/length(I)

  par(new="TRUE",plt=loc[l,],las = 1,cex.axis = 0.75,cex.main=0.75)
  filled.contour3(Useful$x,Useful$y,t(abovemed*Useful$mask),lev=seq(0,1,0.1),col=cm2,main=paste(labels[l]))
}
par(new = "TRUE",plt = c(0.85,0.9,0.3,0.7),las = 1,cex.axis = 1)
filled.legend(Useful$x,Useful$y,t(abovemed*Useful$mask),levels=seq(0,1,0.1),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm2,xlab = "",ylab = "")
dev.off()


plot.new()
par(new="TRUE",plt=loc[l,],las = 1,cex.axis = 1)
filled.contour3(Useful$x,Useful$y,t(abovemed*Useful$mask),lev=seq(0,length(I),1),col=cm2,main=paste(labels[l]))
