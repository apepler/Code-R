setwd('~/Documents/ECLs/Algorithm Comparison')
colors=c("blue","red", "orange", "darkgreen","purple")

#Table 2
warm=read.csv('events_warm.csv')
cool=read.csv('events_cool.csv')
corrs=matrix(0,5,5)
for(i in 1:5)
  for(j in 1:5)
    if(j>i) corrs[i,j]=cor(warm[,i+1],warm[,j+1],use="pairwise.complete.obs") else corrs[i,j]=cor(cool[,i+1],cool[,j+1],use="pairwise.complete.obs")

##Figure 1
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

##Figure 2
data=read.csv('duration2.csv')
plot(data[,1],data[,2],col=colors[1],type="l",lwd=4,xlab="Duration (days)",ylab="Proportion of events",cex=2,yaxt="n",ylim=c(0,1))
axis(2,at=seq(0,1,0.2),labels=c('0%','20%','40%','60%','80%','100%'))
for(i in 2:5) lines(data[,1],data[,i+1],col=colors[i],type="l",lwd=4)
legend("topright",lwd=4,col=colors,legend=c("MLD","LAPB","LAPM","PG","GV"))

#Figure 3
data=read.csv('season.csv')
plot(seq(1,12),data[,2],col=colors[1],type="l",lwd=4,xlab="Month",ylab="Proportion of events",cex=2,yaxt="n",xaxt="n",ylim=c(0,0.15))
axis(1,at=seq(1,12),labels=data[,1])
axis(2,at=seq(0,0.15,0.05),labels=c('0%','5%','10%','15%'))
for(i in 2:5) lines(seq(1,12),data[,i+1],col=colors[i],type="l",lwd=4)
legend("topleft",lwd=4,col=colors,legend=c("MLD","LAPB","LAPM","PG","GV"))

#Figure 4
data=read.csv('seasonB.csv')
plot(seq(1,12),data[,2],col=colors[1],type="l",lwd=4,xlab="Month",ylab="Proportion of events",cex=2,yaxt="n",xaxt="n",ylim=c(0,0.25))
axis(1,at=seq(1,12),labels=data[,1])
axis(2,at=seq(0,0.25,0.05),labels=c('0%','5%','10%','15%','20%','25%'))
for(i in 2:4) lines(seq(1,12),data[,i+1],col=colors[i],type="l",lwd=4)
legend("topleft",lwd=4,col=colors,legend=c("MLD","LAPB","LAPM","PG"))

##Figure5
data=read.csv('counts_btest.csv')
tests=data.frame(Variable=character(36),Method=character(36),HR=rep(0,36),CIlower=rep(0,36),CIupper=rep(0,36),stringsAsFactors=FALSE)
for(i in 1:9)
  for(j in 1:4)
  {
    tests[(i-1)*4+j,1]=as.character(data[i,1])
    tests[(i-1)*4+j,2]=names(data)[j+2]
    a=prop.test(data[i,j+2],data[i,2])
    tests[(i-1)*4+j,3]=a$estimate
    tests[(i-1)*4+j,4:5]=a$conf.int[1:2]      
  }
tests$Variable=factor(tests$Variable,levels=as.character(data[,1]))
tests$Method=factor(tests$Method,levels=c("LAPB","LAPM","PG","GV"))
#make a percentage
tests[,3:5]=tests[,3:5]*100

library(ggplot2)
tiff(file=paste("Fig5.tiff",sep=""), height=400, width=800,pointsize=20)
zp1 <- ggplot(tests, aes(colour = Method)) + scale_colour_manual(values=c("red", "orange", "darkgreen","purple"))
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = HR, ymin = CIlower,
                                 ymax = CIupper),
                             lwd = 1, position = position_dodge(width = 1/2),
                             shape = 21, fill = "WHITE")
zp1 <- zp1 + theme_bw() + theme(text = element_text(size=20)) + theme(axis.title.x = element_blank()) +  ylab("Hit Rate (%)")   
print(zp1) # The trick to these is position_dodge().
dev.off()

##Figure 6
read.csv('Mine_cv4.csv')->mine
read.csv('MLDB.csv')->mldb
read.csv('UM_cv35.csv')->UM
read.csv('Ale_v12.csv')->Ale
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

tiff(file="Fig6.tiff",height=600,width=500)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),cex.axis=1.5)
I=which(mine[,3]==20070628)
a=unique(mine[I,1])
for(i in 1:length(a)) lines(mine[(mine[,1]==a[i]),6],mine[(mine[,1]==a[i]),7],lwd=3,col="red",pch=4, type="b")
I=which(UM[,4]==20070628)
a=unique(UM[I,1])
for(i in 1:length(a)) lines(UM[(UM[,1]==a[i]),7],UM[(UM[,1]==a[i]),8],lwd=3,col="orange",pch=4, type="b")
I=which(Ale[,4]==20070628)
a=unique(Ale[I,1])
for(i in 1:length(a)) lines(Ale[(Ale[,1]==a[i]),6],Ale[(Ale[,1]==a[i]),7],lwd=3,col="darkgreen",pch=4, type="b")
legend("bottomright",legend=c("LAPB","LAPM","PG"),pch=4,col=c("red","orange","darkgreen"),cex=1.5)
dev.off()

#Figure 7
data1=read.csv('HR_days.csv')
data2=read.csv('HR_mslp.csv')
plot.new()
par(plt = c(0.12,0.52,0.2,0.9),las = 1)
plot(data1[,1],data1[,3],col=colors[2],type="l",lwd=4,xlab="Duration (days)",ylab="Hit rate",cex=2,yaxt="n",xaxt="n",ylim=c(0,1))
axis(1,at=seq(1,4),labels=c('0','1','2','3+'))
axis(2,at=seq(0,1,0.2),labels=c('0%','20%','40%','60%','80%','100%'))
for(i in 3:5) lines(data1[,1],data1[,i+1],col=colors[i],type="l",lwd=4)
par(xpd = NA)
par(new = "TRUE",plt = c(0.55,0.95,0.2,0.9),las = 1)
plot(data2[,1],data2[,3],col=colors[2],type="l",lwd=4,xlab="Minimum central pressure (hPa)",ylab="",cex=2,yaxt="n",xaxt="n",ylim=c(0,1))
axis(1,at=seq(985,1015,5),labels=c('<990','990','995','1000','1005','1010','>1015'))
for(i in 3:5) lines(data2[,1],data2[,i+1],col=colors[i],type="l",lwd=4)
legend("bottomleft",lwd=4,col=colors[2:5],legend=c("LAPB","LAPM","PG","GV"))
