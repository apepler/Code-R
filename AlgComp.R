library("RNetCDF")
setwd("~/Documents/ECLs/Algorithm Comparison/")
#lats=seq(-39,-23,2)
#lons=seq(141,159,2)
lats=seq(-39,-23,2)
lons=seq(149,159,2)
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
cols=gray(seq(1,0,-0.1))
source('~/Documents/R/ECLextract.R')

setwd("~/Documents/ECLs/Algorithm Comparison/")
read.csv('Mine_cv4.csv')->mine
read.csv('Fei.csv')->Fei
read.csv('Ale_v12.csv')->Ale
read.csv('MLDB.csv')->MLDB
read.csv('UM_cv35.csv')->UM

p1<-ECLcount_seas(mine[mine$CV>=0.4,],lats,lons,c(1,12))
p1=apply(p1,c(1,2),sum)/30
a1<-ECLcount_seas(Ale,lats,lons,c(1,12))
a1=apply(a1,c(1,2),sum)/30

u1<-ECLcount_seas(UM[UM$Location==1,],lats,lons,c(11,4))
u1=apply(u1,c(1,2),sum)/30


m1<-ECLcount_seas(MLDB,lats,lons,c(1,12))
m1=apply(m1,c(1,2),sum)/27
m1[m1>2]=2
tiff(file=paste("Locs_MLDB_v1a.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(m1),breaks=seq(0,2,0.2),xlab="",ylab="",col=cols,zlim=c(0,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

m1a=m1/sum(m1)
m1a[m1a>0.06]=0.06
cols2=gray(seq(1,0,-0.2))
tiff(file=paste("LocsPC_MLDB_v1a.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(m1a),breaks=seq(0,0.06,0.01),xlab="",ylab="",col=cols2,zlim=c(0,0.06),lab.breaks=c("0","1%","2%","3%","4%","5%","6%"))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

p1<-ECLcount_seas(mine[mine$CV>=0.4 & mine$Time==0,],lats,lons,c(1,12))
p1=apply(p1,c(1,2),sum)/30
p1[p1>2]=2
tiff(file=paste("Locs_Pepler_cv4_v1a.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(p1),breaks=seq(0,2,0.2),xlab="",ylab="",col=cols,zlim=c(0,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

p1<-ECLcount_seas(mine[mine$CV>=0.4,],lats,lons,c(1,12))
p1=apply(p1,c(1,2),sum)/30/4
p1[p1>2]=2
tiff(file=paste("Locs_Pepler_cv4_v2a.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(p1),breaks=seq(0,2,0.2),xlab="",ylab="",col=cols,zlim=c(0,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

p1<-ECLcount_seas(mine[mine$CV>=0.4,],lats,lons,c(1,12))
p1=apply(p1,c(1,2),sum)/30/4
p1a=p1/sum(p1)
p1a[p1a>0.06]=0.06
tiff(file=paste("LocsPC_Pepler_cv4_v2a.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(p1a),breaks=seq(0,0.06,0.01),xlab="",ylab="",col=cols2,zlim=c(0,0.06),lab.breaks=c("0","1%","2%","3%","4%","5%","6%"))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

p1<-ECLcount_seas(Ale[Ale$Closed==-1,],lats,lons,c(1,12))
p1=apply(p1,c(1,2),sum)/30/4
p1a=p1/sum(p1)
p1a[p1a>0.06]=0.06
tiff(file=paste("LocsPC_Alejandro_v2a.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(p1a),breaks=seq(0,0.06,0.01),xlab="",ylab="",col=cols2,zlim=c(0,0.06),lab.breaks=c("0","1%","2%","3%","4%","5%","6%"))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

lats=seq(-39,-23,1)
lons=seq(141,159,1)

f1<-ECLcount_seas(Fei,lats,lons,c(1,12))
f1=apply(f1,c(1,2),sum)/30
tiff(file=paste("LocsPC_Fei_3.tiff",sep=""), height=500, width=600)
image.plot(lons,lats,t(f1),xlab="",ylab="",col=cols)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()



names=c("Pepler_cv25","Pepler_cv4","Pepler_cv4e","Fei","Alejandro")
for(k in 1:5)
{
  hr<-array(0,dim=c(length(lats),length(lons)))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
      {
        I=which(MLDB[,6]>=lons[i]-1 & MLDB[,6]<lons[i]+1 & MLDB[,7]>=lats[j]-1 & MLDB[,7]<lats[j]+1)
        hr[j,i]=sum(MLDB[I,10+k]/length(I))
      }
  tiff(file=paste("HR_",names[k],".tiff",sep=""), height=600, width=500)
  image.plot(lons,lats,t(hr),breaks=seq(0,1,0.1),xlab="",ylab="",col=cols,zlim=c(0,1))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off() 
}

##Locations of MLDB lows where detected by others
library("RNetCDF")
setwd("~/Documents/ECLs/Algorithm Comparison/")
lats=seq(-39,-23,2)
lons=seq(149,159,2)
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
cols=gray(seq(1,0,-0.1))
source('~/Documents/R/ECLextract.R')

setwd("~/Documents/ECLs/Algorithm Comparison/")
read.csv('Mine_cv4.csv')->mine
read.csv('Fei.csv')->Fei
read.csv('Ale.csv')->Ale
read.csv('MLDB.csv')->MLDB

HR=matrix(0,length(lats),5)
HR[,1]=lats
for(i in 1:length(lats))
{
  I=which(MLDB[,7]>=lats[i]-1 & MLDB[,7]<lats[i]+1)
  HR[i,2]=length(I)
  HR[i,3:5]=apply(MLDB[I,13:15],2,sum)
}

plot(HR[,2],HR[,1],type="l",col="blue",xlab="",ylab="",main="ECLs vs. latitude",lwd=2)
cc=c("red","green","orange")
for(i in 3:5) lines(HR[,i],HR[,1],type="l",col=cc[i-2],lwd=2)
legend("topright",lwd=2,col=c("blue","red","orange","green"),legend=c("Speer","Acacia","Alejandro","Fei"))

HR=matrix(0,length(lats),5)
HR[,1]=lats
for(i in 1:length(lats))
{
  I=which(MLDB[,7]>=lats[i]-1 & MLDB[,7]<lats[i]+1 & MLDB[,6]<=160 & MLDB[,6]>=150)
  HR[i,2]=length(I)
  I=which(mine[,7]>=lats[i]-1 & mine[,7]<lats[i]+1 & mine[,6]<=160 & mine[,6]>=150)
  HR[i,3]=length(I)
  I=which(Fei[,7]>=lats[i]-1 & Fei[,7]<lats[i]+1 & Fei[,6]<=160 & Fei[,6]>=150)
  HR[i,4]=length(I)
  I=which(Ale[,7]>=lats[i]-1 & Ale[,7]<lats[i]+1 & Ale[,6]<=160 & Ale[,6]>=150 & Ale$Closed==-1)
  HR[i,5]=length(I)
}

for(i in 2:5) HR[,i]=HR[,i]/sum(HR[,i])
plot(HR[,2],HR[,1],type="l",col="blue",xlab="",ylab="",main="ECLs vs. latitude",lwd=2)
cc=c("red","green","orange")
for(i in c(3,5)) lines(HR[,i],HR[,1],type="l",col=cc[i-2],lwd=2)
legend("topright",lwd=2,col=c("blue","red","orange"),legend=c("Speer","Acacia","Alejandro"))



############################

m1<-ECLcount_seas(MLDB,lats,lons,c(11,4))
m1=apply(m1,c(1,2),sum)/30
m1a=m1/sum(m1)
m1a[m1a>0.06]=0.06
cols2=gray(seq(1,0,-0.2))
tiff(file=paste("LocsPC_MLDB_v1a_warm.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(m1a),breaks=seq(0,0.06,0.01),xlab="",ylab="",col=cols2,zlim=c(0,0.06),lab.breaks=c("0","1%","2%","3%","4%","5%","6%"))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

p1<-ECLcount_seas(mine[mine$CV>=0.4,],lats,lons,c(11,4))
p1=apply(p1,c(1,2),sum)/30/4
p1a=p1/sum(p1)
p1a[p1a>0.06]=0.06
tiff(file=paste("LocsPC_Pepler_cv4_v2a_warm.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(p1a),breaks=seq(0,0.06,0.01),xlab="",ylab="",col=cols2,zlim=c(0,0.06),lab.breaks=c("0","1%","2%","3%","4%","5%","6%"))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

a1<-ECLcount_seas(Ale[Ale$Closed==-1,],lats,lons,c(11,4))
a1=apply(a1,c(1,2),sum)/30/4
a1a=a1/sum(a1)
a1a[a1a>0.06]=0.06
tiff(file=paste("LocsPC_Alejandro_v2a_warm.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(a1a),breaks=seq(0,0.06,0.01),xlab="",ylab="",col=cols2,zlim=c(0,0.06),lab.breaks=c("0","1%","2%","3%","4%","5%","6%"))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

####################
setwd("~/Documents/ECLs/CSV/")
library("RNetCDF")
#lats=seq(-39,-23,2)
#lons=seq(141,159,2)
lats=seq(-39,-23,2)
lons=seq(149,159,2)
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
cols=gray(seq(1,0,-0.1))
source('~/Documents/R/ECLextract.R')

read.csv('ECLfixes_wrf25_cv25.csv')->wrf
read.csv('ECLfixes_ncep25_cv25.csv')->mine
n1<-ECLcount_seas(mine,lats,lons,c(1,12))
n1=apply(n1,c(1,2),sum)/29/4
w1<-ECLcount_seas(wrf,lats,lons,c(1,12))
w1=apply(w1,c(1,2),sum)/29/4

cols=gray(seq(1,0.1,-0.12))
tiff(file=paste("../Locs_NCEP1_cv25events_2.tiff",sep=""), height=600, width=500,pointsize=20)
image.plot(lons,lats,t(n1),breaks=seq(0,2,0.25),xlab="",ylab="",col=cols,zlim=c(0,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=paste("../Locs_WRF25_cv25events_2.tiff",sep=""), height=600, width=500,pointsize=20)
image.plot(lons,lats,t(w1),breaks=seq(0,2,0.25),xlab="",ylab="",col=cols,zlim=c(0,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

read.csv('../Algorithm Comparison/MLDB.csv')->MLDB
m1<-ECLcount_seas(MLDB,lats,lons,c(1,12))
m1=apply(m1,c(1,2),sum)/27
m1[m1>2]=2
tiff(file=paste("../Locs_MLDB.tiff",sep=""), height=600,width=500,pointsize=20)
image.plot(lons,lats,t(m1),breaks=seq(0,2,0.25),xlab="",ylab="",col=cols,zlim=c(0,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

n1<-ECLcount_seas(UM_T,lats,lons,c(1,12))
n1=apply(n1,c(1,2),sum)/29/4
w1<-ECLcount_seas(UM_noT,lats,lons,c(1,12))
w1=apply(w1,c(1,2),sum)/29/4

cols=gray(seq(1,0.1,-0.12))
tiff(file=paste("../Locs_NCEP1_cv25events_2.tiff",sep=""), height=600, width=500,pointsize=20)
image.plot(lons,lats,t(n1),breaks=seq(0,2,0.25),xlab="",ylab="",col=cols,zlim=c(0,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=paste("../Locs_WRF25_cv25events_2.tiff",sep=""), height=600, width=500,pointsize=20)
image.plot(lons,lats,t(w1),breaks=seq(0,2,0.25),xlab="",ylab="",col=cols,zlim=c(0,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()



###########
rm(list=ls())
library("RNetCDF")
setwd("~/Documents/ECLs/Algorithm Comparison/")
#lats=seq(-39,-23,2)
#lons=seq(141,159,2)
lats=seq(-39,-23,2)
lons=seq(149,159,2)
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
cols=gray(seq(1,0,-0.1))
source('~/Documents/R/ECLextract.R')
read.csv('Dates.csv')->all

comp=matrix(0,5,5)
for(i in 1:5)
  for(j in 1:5)
  {
    I=which(all[,i+9]==1 & all[,j+9]==1)
    J=which(all[,i+9]==1 | all[,j+9]==1)
    comp[i,j]=length(I)/length(J)
  }

read.csv('Mine_cv4.csv')->mine
read.csv('Fei.csv')->Fei
read.csv('Ale.csv')->Ale
read.csv('MLDB.csv')->MLDB
read.csv('UM_cv35.csv')->UM


###Playing with a binomial test for hit rat
setwd("~/Documents/ECLs/Algorithm Comparison/")
data=read.csv('counts_btest.csv')

k=1
diff=matrix(0,4,4)

for(k in 1:9)
{
for(i in 1:4)
for(j in 1:4)
{
a=binom.test(data[k,i+2],data[k,2],data[k,j+2]/data[k,2],alternative="two.sided")
diff[i,j]=a$p.value
}
I=which(diff<=.05,arr.ind=T)
data[k,1]
I[(I[,1]!=I[,2]),]

for(i in 1:4)
  for(j in 1:4)
  {
    a=prop.test(c(data[k,i+2],data[k,j+2]),c(data[k,2],data[k,2]))
    diff[i,j]=a$p.value
  }
I=which(diff<=.05,arr.ind=T)
print(data[k,1])
print(I[(I[,1]!=I[,2]),])
}

##Trying to make a pretty plot
##This is a matrix of the CIs for the hit rates!
tests=data.frame(Variable=character(36),Method=character(36),HR=rep(0,36),CIlower=rep(0,36),CIupper=rep(0,36),stringsAsFactors=FALSE)
for(i in 1:9)
  for(j in 1:4)
  {
    tests[(i-1)*4+j,1]=as.character(data[i,1])
    tests[(i-1)*4+j,2]=names(data)[j+2]
    a=binom.test(data[i,j+2],data[i,2],0,alternative="two.sided")
    tests[(i-1)*4+j,3]=a$estimate
    tests[(i-1)*4+j,4:5]=a$conf.int[1:2]    
  }
tests$Variable=factor(tests$Variable,levels=as.character(data[,1]))
tests$Method=factor(tests$Method,levels=c("LAPB","LAPM","PG","GV"))
#make a percentage
tests[,3:5]=tests[,3:5]*100

tests2=tests
for(i in 1:9)
  for(j in 1:4)
  {
    a=prop.test(data[i,j+2],data[i,2])
    tests2[(i-1)*4+j,3]=a$estimate
    tests2[(i-1)*4+j,4:5]=a$conf.int[1:2]    
  }
#make a percentage
tests2[,3:5]=tests2[,3:5]*100

library(ggplot2)

tiff(file=paste("HitRates_Methods.tiff",sep=""), height=400, width=800,pointsize=20)
zp1 <- ggplot(tests, aes(colour = Method)) + scale_colour_manual(values=c("red", "orange", "darkgreen","purple"))
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = HR, ymin = CIlower,
                                 ymax = CIupper),
                             lwd = 1, position = position_dodge(width = 1/2),
                             shape = 21, fill = "WHITE")
zp1 <- zp1 + theme_bw() + theme(text = element_text(size=20)) + theme(axis.title.x = element_blank()) +  ylab("Hit Rate (%)")   
print(zp1) # The trick to these is position_dodge().
dev.off()

tiff(file=paste("HitRates_Methods2.tiff",sep=""), height=400, width=800,pointsize=20)
zp1 <- ggplot(tests2, aes(colour = Method)) + scale_colour_manual(values=c("red", "orange", "darkgreen","purple"))
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = HR, ymin = CIlower,
                                 ymax = CIupper),
                             lwd = 1, position = position_dodge(width = 1/2),
                             shape = 21, fill = "WHITE")
zp1 <- zp1 + theme_bw() + theme(text = element_text(size=20)) + theme(axis.title.x = element_blank()) +  ylab("Hit Rate (%)")   
print(zp1) # The trick to these is position_dodge().
dev.off()


