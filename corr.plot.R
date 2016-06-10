##Plots the filled contour plot with a set colourmap and the mask outline

##Recommended colormap: 
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(20)
cm[9:12]="white"

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
corr.plot <- function (x, y, z, mask, cm=rich.colors(20),filename) 
  
{
  if(missing(filename)) plot.new() else pdf(file=filename, height=5, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)  
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(x,y,z,levels=seq(-1,1,0.1),col=cm)
  par(xpd = NA)
  contour(x,y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(x,y,z,levels=seq(-1,1,0.1),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")
  if(!missing(filename)) dev.off()
}

mCorrs<-function(x,y,names,fname)
{
  corrs=rep(0,12)
  for(i in 1:12) corrs[i]<-cor(x[,i+1],y[,i+1])
  months=1:12
  if(missing(fname)) plot.new() else pdf(file=fname, height=4.5, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
  plot(months,corrs,type="o",col="black",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=4,lwd=2)
  axis(1, at = seq(1,12,1), labels = names,cex.axis=0.75)
  color <- rgb(190, 190, 190, alpha=120, maxColorValue=255)
  rect(xleft=-1000,xright=1000,ybottom=-0.16, ytop=0.16, density=100, col=color)
  if(!missing(fname)) dev.off()
}


mCorrs3<-function(x,y,names,fname)
{
  corrs<-rep(0,12)
  corrs[1]<-cor(rowMeans(cbind(x[1:(length(x[,1])-1),13],x[2:length(x[,1]),2:3])),rowMeans(cbind(y[1:(length(y[,1])-1),13],y[2:length(y[,1]),2:3])))
  corrs[12]<-cor(rowMeans(cbind(x[1:(length(x[,1])-1),12:13],x[2:length(x[,1]),2])),rowMeans(cbind(y[1:(length(y[,1])-1),12:13],y[2:length(y[,1]),2])))
  for(i in 2:11) corrs[i]<-cor(rowMeans(x[1:(length(x[,1])-1),i:(i+2)]),rowMeans(y[1:(length(y[,1])-1),i:(i+2)]))
  
  months=1:12
  if(missing(fname)) plot.new() else pdf(file=fname, height=4.5, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
  plot(months,corrs,type="o",col="black",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=4,lwd=2)
  axis(1, at = seq(1,12,1), labels = names,cex.axis=0.75)
  color <- rgb(190, 190, 190, alpha=120, maxColorValue=255)
  rect(xleft=-1000,xright=1000,ybottom=-0.16, ytop=0.16, density=100, col=color)
  if(!missing(fname)) dev.off()
}

mCorrs3a<-function(x,y1,y2,names,fname,name1,name2)
{
  corrs1<-corrs2<-rep(0,12)
  corrs1[1]<-cor(rowMeans(cbind(x[1:(length(x[,1])-1),13],x[2:length(x[,1]),2:3])),rowMeans(cbind(y1[1:(length(y1[,1])-1),13],y1[2:length(y1[,1]),2:3])))
  corrs1[12]<-cor(rowMeans(cbind(x[1:(length(x[,1])-1),12:13],x[2:length(x[,1]),2])),rowMeans(cbind(y1[1:(length(y1[,1])-1),12:13],y1[2:length(y1[,1]),2])))
  for(i in 2:11) corrs1[i]<-cor(rowMeans(x[1:(length(x[,1])-1),i:(i+2)]),rowMeans(y1[1:(length(y1[,1])-1),i:(i+2)]))
  corrs2[1]<-cor(rowMeans(cbind(x[1:(length(x[,1])-1),13],x[2:length(x[,1]),2:3])),rowMeans(cbind(y2[1:(length(y2[,1])-1),13],y2[2:length(y2[,1]),2:3])))
  corrs2[12]<-cor(rowMeans(cbind(x[1:(length(x[,1])-1),12:13],x[2:length(x[,1]),2])),rowMeans(cbind(y2[1:(length(y2[,1])-1),12:13],y2[2:length(y2[,1]),2])))
  for(i in 2:11) corrs2[i]<-cor(rowMeans(x[1:(length(x[,1])-1),i:(i+2)]),rowMeans(y2[1:(length(y2[,1])-1),i:(i+2)]))
  months=1:12
  
  if(missing(fname)) plot.new() else pdf(file=fname, height=4.5, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
  plot(months,corrs1,type="o",col="blue",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=4,lwd=2,ylim=range(corrs1,corrs2))
  axis(1, at = seq(1,12,1), labels = names(soi[,2:13]),cex.axis=0.75)
  lines(months,corrs2,type="o",col="red",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=4,lwd=2)
  legend("right",c(name1,name2),fill=c("blue","red"))
  color <- rgb(190, 190, 190, alpha=120, maxColorValue=255)
  rect(xleft=-1000,xright=1000,ybottom=-0.16, ytop=0.16, density=100, col=color)
  if(!missing(fname)) dev.off()
}

##Partial correlations
mCorrs3b<-function(x,y1,y2,names,fname,name1,name2)
{
  corr1<-corr2<-rep(0,12)
  months=1:12
  source('~/R/pcor.R')
  
  a<-pcor.test(rowMeans(cbind(x[1:(length(x[,1])-1),13],x[2:length(x[,1]),2:3])),rowMeans(cbind(y1[1:(length(y1[,1])-1),13],y1[2:length(y1[,1]),2:3])),rowMeans(cbind(y2[1:(length(y2[,1])-1),13],y2[2:length(y2[,1]),2:3])))
  corr1[1]<-as.numeric(a[1])
  a<-pcor.test(rowMeans(cbind(x[1:(length(x[,1])-1),12:13],x[2:length(x[,1]),2])),rowMeans(cbind(y1[1:(length(y1[,1])-1),12:13],y1[2:length(y1[,1]),2])),rowMeans(cbind(y2[1:(length(y2[,1])-1),12:13],y2[2:length(y2[,1]),2])))
  corr1[12]<-as.numeric(a[1])
  a<-pcor.test(rowMeans(cbind(x[1:(length(x[,1])-1),13],x[2:length(x[,1]),2:3])),rowMeans(cbind(y2[1:(length(y2[,1])-1),13],y2[2:length(y2[,1]),2:3])),rowMeans(cbind(y1[1:(length(y1[,1])-1),13],y1[2:length(y1[,1]),2:3])))
  corr2[1]<-as.numeric(a[1])
  a<-pcor.test(rowMeans(cbind(x[1:(length(x[,1])-1),12:13],x[2:length(x[,1]),2])),rowMeans(cbind(y2[1:(length(y2[,1])-1),12:13],y2[2:length(y2[,1]),2])),rowMeans(cbind(y1[1:(length(y1[,1])-1),12:13],y1[2:length(y1[,1]),2])))
  corr2[12]<-as.numeric(a[1])  
  
  for(i in 2:11)
    {
    a<-pcor.test(rowMeans(x[1:(length(x[,1])-1),i:(i+2)]),rowMeans(y1[1:(length(y1[,1])-1),i:(i+2)]),rowMeans(y2[1:(length(y2[,1])-1),i:(i+2)]))
    corr1[i]<-as.numeric(a[1])
    a<-pcor.test(rowMeans(x[1:(length(x[,1])-1),i:(i+2)]),rowMeans(y2[1:(length(y2[,1])-1),i:(i+2)]),rowMeans(y1[1:(length(y1[,1])-1),i:(i+2)]))
    corr2[i]<-as.numeric(a[1])
  }

  if(missing(fname)) plot.new() else pdf(file=fname, height=4.5, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
  plot(months,corr1,type="o",col="blue",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=4,lwd=2,ylim=range(corr1,corr2))
  axis(1, at = seq(1,12,1), labels = names(soi[,2:13]),cex.axis=0.75)
  lines(months,corr2,type="o",col="red",xaxt="n", xlab="", ylab="",cex.axis=0.75,pch=4,lwd=2)
  legend("topleft",c(name1,name2),fill=c("blue","red"))
  color <- rgb(190, 190, 190, alpha=120, maxColorValue=255)
  rect(xleft=-1000,xright=1000,ybottom=-0.16, ytop=0.16, density=100, col=color)
  if(!missing(fname)) dev.off()
}



gen.plot <- function (x, y, z, mask, lev=seq(-1,1,0.1), cm=rich.colors(20),filename)
  
{
  if(missing(filename)) plot.new() else pdf(file=filename, height=5, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)  
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(x,y,z,levels=lev,col=cm)
  par(xpd = NA)
  contour(x,y,mask,add=T,drawlabels=F)
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(x,y,z,levels=lev,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")
  if(!missing(filename)) dev.off()
}

gen.plot.nsw <- function (x, y, z, mask, lev=seq(-1,1,0.1), cm=rich.colors(20),filename)
  
{
  library("R.matlab")
  readMat('Useful3.mat')->Useful
  coastline<-Useful$coastline
  
  if(missing(filename)) plot.new() else pdf(file=filename, height=5, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)  
  par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
  filled.contour3(x,y,z,levels=lev,col=cm,xlim=c(140,155),ylim=c(-40,-25))
  par(xpd = NA)
  lines(coastline[,1],coastline[,2])
  par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
  filled.legend(x,y,z,levels=lev,,xlim=c(140,155),ylim=c(-40,-25),col=cm,xlab = "",ylab = "")
  if(!missing(filename)) dev.off()
}
