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
filled.contour3(Useful$x,Useful$y,mask2,col=c("white",gray(0.8)),xlim=c(140,165),ylim=c(-45,-20))
contour(Useful$x,Useful$y,mask,drawlabels=F,add=T)
par(xpd=NA)
lines(c(149,161,161,152,152,149,149),c(-41,-41,-24,-24,-31,-38,-41),lwd=3,col="black")
#Add wind locations - Amberley, Coffs, Williamtown, Sydney, Nowra, East Sale
points(c(152.7111,153.1187,151.8359,151.1731,150.5353,147.1322),c(-27.6297,-30.3107,-32.7932,-33.9465,-34.9469,-38.1156),pch=16,cex=2) #Circles
#Add wave locations - Byron, Coffs, Crowdy Head, Sydney, Pt Kembla, Batemans, Eden, 
points(c(153.7166666667,153.2666666667,152.85,151.4166666667,151.0166666667,150.3333333333,150.1833333333),c(-28.8166666667,-30.35,-31.8166666667,-33.7666666667,-34.4666666667,-35.7,-37.3),pch=17,cex=2,col=gray(0.4)) #Triangles
#Add Sydney and Brisbane - add 1 deg Long so can see
points(152.211111,-33.859972,pch="S",cex=2) 
points(154.027778,-27.467917,pch="B",cex=2) 
rect(141,-36,156,-25,lwd=3,col=NA,border="black",lty=2)
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																								

