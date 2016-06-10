setwd('~/Documents/ECLs/Algorithm Comparison')
colors=c("blue","red", "orange", "darkgreen","purple","purple4")

#Table 2
warm=read.csv('events_warm_v3.csv')
cool=read.csv('events_cool_v3.csv')
#warm=warm[,-6]
#cool=cool[,-6]

corrs=matrix(0,5,5)
for(i in 1:5)
  for(j in 1:5)
    if(j>i) corrs[i,j]=cor(warm[,i+1],warm[,j+1],use="pairwise.complete.obs") else corrs[i,j]=cor(cool[,i+1],cool[,j+1],use="pairwise.complete.obs")

write.csv(corrs,file='table2_v3.csv')

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
#par(xpd=NA)
lines(c(149,161,161,152,152,149,149),c(-41,-41,-24,-24,-31,-38,-41),lwd=3,col="black")
#Add wind locations - Amberley, Coffs, Williamtown, Sydney, Nowra, East Sale
points(c(152.7111,153.1187,151.8359,151.1731,150.5353,147.1322),c(-27.6297,-30.3107,-32.7932,-33.9465,-34.9469,-38.1156),pch=16,cex=2) #Circles
#Add wave locations - Byron, Coffs, Crowdy Head, Sydney, Pt Kembla, Batemans, Eden, 
points(c(153.7166666667,153.2666666667,152.85,151.4166666667,151.0166666667,150.3333333333,150.1833333333),c(-28.8166666667,-30.35,-31.8166666667,-33.7666666667,-34.4666666667,-35.7,-37.3),pch=17,cex=2,col=gray(0.4)) #Triangles
#Add Sydney and Brisbane - add 1 deg Long so can see
points(152.211111,-33.859972,pch="S",cex=2) 
points(154.027778,-27.467917,pch="B",cex=2) 
rect(145.5,-34.25,160.5,-23.75,lwd=3,col=NA,border="black",lty=2)

##Figure 2
data=read.csv('duration2_v3.csv')
#data=data[,-6]
plot(data[,1],data[,2],col=colors[1],type="l",lwd=4,xlab="Duration (days)",ylab="Proportion of events",cex=2,yaxt="n",ylim=c(0,1))
axis(2,at=seq(0,1,0.2),labels=c('0%','20%','40%','60%','80%','100%'))
for(i in 2:6) lines(data[,1],data[,i+1],col=colors[i],type="l",lwd=4)
legend("topright",lwd=4,col=colors,legend=c("MLD","LAPv1","LAPv2","PG","ULGV 48hr","ULGV 12hr"),ncol=2)

#Figure 3
data=read.csv('season_v3.csv')
#data=data[,-6]
pdf('Figure3.pdf',height=10,width=15,pointsize=32)
par(mar=c(5,4,2,2)+0.1)
plot(seq(1,12),data[,2],col=colors[1],type="l",lwd=8,xlab="Month",ylab="Proportion of events",yaxt="n",xaxt="n",ylim=c(0,0.15))
axis(1,at=seq(1,12),labels=data[,1])
axis(2,at=seq(0,0.15,0.05),labels=c('0%','5%','10%','15%'))
for(i in 2:5) lines(seq(1,12),data[,i+1],col=colors[i],type="l",lwd=8)
legend("bottomright",lwd=8,col=colors,legend=c("MLD","LAPv1","LAPv2","PG","ULGV"),ncol=2,bty="n")
dev.off()

#Figure 4
data=read.csv('seasonB_v3.csv')
plot(seq(1,12),data[,2],col=colors[1],type="l",lwd=4,xlab="Month",ylab="Proportion of events",cex=2,yaxt="n",xaxt="n",ylim=c(0,0.25))
axis(1,at=seq(1,12,2),labels=data[seq(1,12,2),1])
axis(2,at=seq(0,0.25,0.05),labels=c('0%','5%','10%','15%','20%','25%'))
for(i in 2:4) lines(seq(1,12),data[,i+1],col=colors[i],type="l",lwd=4)
legend("topleft",lwd=4,col=colors,legend=c("MLD","LAPv1","LAPv2","PG"))


##Figure5
data=read.csv('counts_btest.csv')
data=data[,-6]
names(data)[6]="ULGV"
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
tests$Method=factor(tests$Method,levels=c("LAPv1","LAPv2","PG","ULGV"))
#make a percentage
tests[,3:5]=tests[,3:5]*100

library(ggplot2)
tiff(file=paste("Fig5_v4.tiff",sep=""), height=400, width=800,pointsize=20)
zp1 <- ggplot(tests, aes(colour = Method)) + scale_colour_manual(values=c("red", "orange", "darkgreen","purple"))
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = HR, ymin = CIlower,
                                 ymax = CIupper),
                             lwd = 1, position = position_dodge(width = 1/2),
                             shape = 21, fill = "WHITE")
zp1 <- zp1 + theme_bw() + theme(text = element_text(size=20)) + theme(axis.title.x = element_blank()) +  ylab("Hit Rate (%)")   
print(zp1) # The trick to these is position_dodge().
dev.off()

##Version 2 - using a bootstrap
data=read.csv('MLDBevents_hits.csv')
tests=data.frame(Variable=character(36),Method=character(36),HR=rep(0,36),CIlower=rep(0,36),CIupper=rep(0,36),stringsAsFactors=FALSE)

names1=c("All","E","W","MayOct","NovApr","Rain","Wind","Wave","Expl")
names2=c("LAPv1","LAPv2","PG","ULGV")

for(i in 1:9)  
  for(j in 1:4)
    {
    tests[(i-1)*4+j,1]=names1[i]
    tests[(i-1)*4+j,2]=names2[j]
  }

hr.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)

for(i in 1:4)
{
  tests[i,3]=sum(data[,9+i])/length(data[,1])
  a=boot(data[,9+i],hr.fun,R=10000)
  b=boot.ci(a,type="bca")
  tests[i,4:5]=b$bca[4:5]
}

for(i in 1:8)
{
  I=which(data[,i+1]==1)
  data2=data[I,]
  for(j in 1:4)
  {
    tests[i*4+j,3]=sum(data2[,9+j])/length(data2[,1])
    a=boot(data2[,9+j],hr.fun,R=10000)
    b=boot.ci(a,type="bca")
    tests[i*4+j,4:5]=b$bca[4:5]
  }
}



#Figure 6
data1=read.csv('HR_days_v3.csv')
data2=read.csv('HR_mslp_v3.csv')
#data1=data1[,-6]
#data2=data2[,-6]
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
legend("bottomleft",lwd=4,col=colors[2:5],legend=c("LAPv1","LAPv2","PG","ULGV"))

##Figure 7
read.csv('Mine_rad2.csv')->mine
read.csv('MLDB.csv')->mldb
read.csv('UM_rad2_p100.csv')->UM
read.csv('Ale_v15CTL.csv')->Ale
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

tiff(file="Fig7_v4.tiff",height=600,width=500)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),cex.axis=1.5)
I=which(mine[,3]==20070628)
a=unique(mine[I,1])
for(i in 1:length(a)) lines(mine[(mine[,1]==a[i]),6],mine[(mine[,1]==a[i]),7],lwd=3,col="red",pch=4, type="b")
I=which(UM[,4]==20070628)
a=unique(UM[I,2])
for(i in 1:length(a)) lines(UM[(UM[,2]==a[i]),7],UM[(UM[,2]==a[i]),8],lwd=3,col="orange",pch=4, type="b")
I=which(Ale[,5]==20070628)
a=unique(Ale[I,1])
for(i in 1:length(a)) lines(Ale[(Ale[,1]==a[i]),8],Ale[(Ale[,1]==a[i]),9],lwd=3,col="darkgreen",pch=4, type="b")
legend("bottomright",legend=c("LAPv1","LAPv2","PG"),pch=4,col=c("red","orange","darkgreen"),cex=1.5)
dev.off()

##Figure 8
library(RNetCDF)
a<-open.nc("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_1994-01_1998-12.nc")
lonP=var.get.nc(a,'longitude')
latP=var.get.nc(a,'latitude')
time=var.get.nc(a,'time')
time=as.Date(time/24,origin="1900-01-01")
yyyymm=as.numeric(substr(time,1,4))*100+as.numeric(substr(time,6,7))
I=which(yyyymm==199808)
MSLP=var.get.nc(a,"msl",c(1,1,I[1]),c(length(lonP),length(latP),length(I)),unpack=T)/100

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
readMat('~/Documents/Data/mask_escci2.mat')->escci
mask<-t(Useful$mask)
mask[is.na(mask)]=0
mask2<-t(escci$mask.escci)
mask2[is.na(mask2)]=0																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																								

contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-20))
contour(lonP,latP[length(latP):1],MSLP[,length(latP):1,29],levels=seq(992,1050,4),add=T) #19980808 = 29 (1 + 7*4)

setwd("~/Documents/ECLs/Algorithm Comparison/")
read.csv('Mine_rad2.csv')->mine
read.csv('Fei_Dowdy.csv')->Fei
read.csv('Ale_v15CTL.csv')->Ale
read.csv('MLDB.csv')->MLDB
read.csv('UM_rad2_p100.csv')->UM
mine$Time=(as.numeric(mine$Time)-1)*6 #Convert to integer
UM$Time=(as.numeric(UM$Time)-1)*6 #Convert to integer

I=which(MLDB$Date==19980808)
points(MLDB$Lon[I],MLDB$Lat[I],col="blue",pch=4,cex=2,lwd=4)
I=which(mine$Date==19980808 & mine$Time==0 & mine$Location==1)
points(mine$Lon[I],mine$Lat[I],col="red",pch=4,cex=2,lwd=4)
I=which(UM$Date==19980808 & UM$Time==0 & UM$Location==1)
points(UM$Lon[I],UM$Lat[I],col="orange",pch=4,cex=2,lwd=4)
I=which(Ale$Date==19980808 & Ale$Time==0 &Ale$Location==1)
points(Ale$Lon[I],Ale$Lat[I],col="darkgreen",pch=4,cex=2,lwd=4)
legend("bottomleft",col=colors[1:4],legend=c("MLD","LAPv1","LAPv2","PG"),pch=rep(4,4),pt.lwd=rep(4,4))


##Well, this didn't seem to work
#a<-open.nc("/srv/ccrc/data04/z3346206/Acacia/GV_ERA-I/VorZ_6hrs_pl_1998_08.nc")
a<-open.nc("/srv/ccrc/data04/z3346206/IPV+VOR/ERA-Interim/GV_cordex/VorZ_6hrs_pl_1998_08.nc")
lonGV=var.get.nc(a,'lon')
latGV=var.get.nc(a,'lat')
GV=var.get.nc(a,'lap_new',unpack=T)
a=apply(GV[,,21:29],c(1,2),mean)
contour(lonGV,latGV,a,add=T,lty=2,levels=seq(-8,20,2)) #19980808 = 29 (1 + 7*4)

a<-open.nc("/srv/ccrc/data34/z3478332/ERAI/erai_500hPa.nc")
GH=var.get.nc(a,"z",unpack=T)
V0=var.get.nc(a,"vo",unpack=T)
contour(lonP,latP[length(latP):1],V0[,length(latP):1,29],lty=2,add=T) #19980808 = 29 (1 + 7*4)



##Figure 7
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

list=list(read.csv('Ale_v15.csv'),read.csv('Ale_v15L.csv'),read.csv('Ale_v15CTL.csv'),read.csv('Ale_v15LHR.csv'))
names=c("50km B","50km L","10km B","10km L")
cm=c("red","orange","green","cyan")

pdf(file="Ale_comp_20070828.pdf",height=6,width=5)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,165),ylim=c(-45,-25),cex.axis=1.5)
for(i in 1:4)
{
  b=list[[i]]
  I=which(b$Date==20070628)
  a=unique(b$ID2[I])
  for(j in 1:length(a)) lines(b$Lon[b$ID2==a[j]],b$Lat[b$ID2==a[j]],lwd=3,col=cm[i],pch=4, type="b")
}
legend("bottomright",legend=names,pch=4,col=cm,cex=1.5)
dev.off()


##########For poster

##Tracks for ECLs
library(RNetCDF)
f1<-open.nc("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_1994-01_1998-12.nc")
lonP=var.get.nc(f1,'longitude')
latP=var.get.nc(f1,'latitude')
time=var.get.nc(f1,'time')
hh=time%%24
time=as.Date(time/24,origin="1900-01-01")
date=as.numeric(substr(time,1,4))*10000+as.numeric(substr(time,6,7))*100+
  as.numeric(substr(time,9,10))

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
readMat('~/Documents/Data/mask_escci2.mat')->escci
mask<-t(Useful$mask)
mask[is.na(mask)]=0
mask2<-t(escci$mask.escci)
mask2[is.na(mask2)]=0  																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

setwd("~/Documents/ECLs/Algorithm Comparison/")
read.csv('Mine_rad2.csv')->mine
read.csv('Ale_v15CTL.csv')->Ale
read.csv('MLDB.csv')->MLDB
read.csv('UM_rad2_p100.csv')->UM
mine$Time=(as.numeric(mine$Time)-1)*6 #Convert to integer
UM$Time=(as.numeric(UM$Time)-1)*6 #Convert to integer

list=list(MLDB,mine,UM,Ale)
names=c("MLD","LAPv1","LAPv2","PG")
cm=c("blue","red","orange","darkgreen")

keydates=c(19980623,19980808)

for(t in 1:2)
{
  I=which(date==keydates[t] & hh==0)
  MSLP=var.get.nc(f1,"msl",c(1,1,I),c(length(lonP),length(latP),1),unpack=T)/100
  
  pdf(file=paste("Fig8_forposter_",keydates[t],".pdf",sep=""),width=4.5,height=5)
  filled.contour3(Useful$x,Useful$y,mask,col=c("white","lightgreen"),xlim=c(145,170),ylim=c(-45,-20),cex.axes=0.8)
  
  polygon(x=c(149,161,161,152,152,149,149),y=c(-41,-41,-24,-24,-31,-38,-41),col="gray",density=20,border=NA)
  contour(lonP,latP[length(latP):1],MSLP[,length(latP):1],levels=seq(992,1050,4),add=T,lwd=2,labcex=1) #19980808 = 29 (1 + 7*4)
  for(i in 1:4)
  {
    b=list[[i]]
    I=which(b$Date==keydates[t] & b$Time==0)
    points(b$Lon[I],b$Lat[I],col=cm[i],pch=4,cex=2,lwd=3)
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,lwd=3,col=cm[i])
    }
  }
  legend("topright",legend=names,lwd=3,col=cm,bg="white",cex=0.8)
  dev.off()
}

##Tracks for ECLs - June 2007
library(RNetCDF)
f1<-open.nc("/srv/ccrc/data34/z3478332/ERAI/ERAI_mslp_2005-01_2010-12.nc")
lonP=var.get.nc(f1,'longitude')
latP=var.get.nc(f1,'latitude')
time=var.get.nc(f1,'time')
hh=time%%24
time=as.Date(time/24,origin="1900-01-01")
date=as.numeric(substr(time,1,4))*10000+as.numeric(substr(time,6,7))*100+
  as.numeric(substr(time,9,10))

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
readMat('~/Documents/Data/mask_escci2.mat')->escci
mask<-t(Useful$mask)
mask[is.na(mask)]=0
mask2<-t(escci$mask.escci)
mask2[is.na(mask2)]=0    																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																						

setwd("~/Documents/ECLs/Algorithm Comparison/")
read.csv('Mine_rad2.csv')->mine
read.csv('Ale_v15CTL.csv')->Ale
read.csv('MLDB.csv')->MLDB
read.csv('UM_rad2_p100.csv')->UM
mine$Time=(as.numeric(mine$Time)-1)*6 #Convert to integer
UM$Time=(as.numeric(UM$Time)-1)*6 #Convert to integer

list=list(mine,UM,Ale)
names=c("LAPv1","LAPv2","PG")
cm=c("red","orange","darkgreen")

keydates=c(20070608,20070619,20070627)

for(t in 3)
{
  I=which(date==keydates[t] & hh==0)
  MSLP=var.get.nc(f1,"msl",c(1,1,I),c(length(lonP),length(latP),1),unpack=T)/100
  
  pdf(file=paste("Fig8_forposter_",keydates[t],".pdf",sep=""),width=9,height=10,pointsize=36)
  filled.contour3(Useful$x,Useful$y,mask,col=c("white","lightgreen"),xlim=c(145,170),ylim=c(-45,-20),cex.axes=0.8)
  
  polygon(x=c(149,161,161,152,152,149,149),y=c(-41,-41,-24,-24,-31,-38,-41),col="gray",density=20,border=NA)
  contour(lonP,latP[length(latP):1],MSLP[,length(latP):1],levels=seq(992,1050,4),add=T,lwd=2,labcex=1) #19980808 = 29 (1 + 7*4)
  for(i in 1:3)
  {
    b=list[[i]]
    I=which(b$Date==keydates[t] & b$Time==0)
    points(b$Lon[I],b$Lat[I],col=cm[i],pch=4,cex=2,lwd=3)
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,lwd=3,col=cm[i])
    }
  }
  legend("topright",legend=names,lwd=3,col=cm,bg="white",cex=0.8)
  dev.off()
}


#### ECL average locations vs. EAC separation above/below ave (cool)

library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

read.csv('~/Documents/Timeseries/eacsep.csv')->eacsep
read.csv('~/Documents/ECLs/Algorithm Comparison/UM_rad2_p100.csv')->UM

month=floor(UM$Date/100)%%100
year=floor(UM$Date/10000)

eacsep=eacsep[c(-1,-19),] ##Remove 1992 (half-year) and 2010
I=which((year>=1993 & year <=2009) & (month>=5 & month<=10))
UM=UM[I,]
year=year[I]

UM$Loc2=rep(0,length(UM$Location))
I<-which(UM$Lon>=149 & UM$Lon<=155 & UM$Lat<(-37) & UM$Lat>=-41)
UM$Loc2[I]<-1
I<-which(UM$Lon>=(149+(37+UM$Lat)/2) & UM$Lon<=(155+(37+UM$Lat)/2) & UM$Lat<(-31) & UM$Lat>=-37)
UM$Loc2[I]<-1
I<-which(UM$Lon>=152 & UM$Lon<=158 & UM$Lat<=(-24) & UM$Lat>=-31)
UM$Loc2[I]<-1

ylist=rep(0,17)
for(n in 1:length(eacsep$Year))
{
  I=which(UM$Loc2==1 & year==eacsep$Year[n])
  J=unique(UM$ID[I])
  ylist[n]=length(J)
}


lats=seq(-45,-10,2.5)
lons=seq(145,170,2.5)
locN<-locS<-matrix(0,length(lats),length(lons))
for(n in 1:length(eacsep$Year))
for(i in 1:length(lons))
  for(j in 1:length(lats))
  {
    I=which(UM$Lon>=lons[i]-1.25 & UM$Lon<lons[i]+1.25 & UM$Lat>=lats[j]-1.25 & UM$Lat<lats[j]+1.25 & year==eacsep$Year[n])
if(length(I)>0)
  if(eacsep$Cool[n] >= median(eacsep$Cool)) locN[j,i]=locN[j,i]+length(I) else locS[j,i]=locS[j,i]+length(I)
  }

locN=locN/9
locS=locS/8
bb=c(-0.5,0,0.5,1,2,3,5,7)
cols=gray(seq(1,0.1,-0.15))
ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}

layout(cbind(1,2,3),c(1,1,0.3))
image(lons,lats,t(locN),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,7),xlim=c(145,160),ylim=c(-45,-25),main="North of average",cex.axis=1.5,cex.main=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lons,lats,t(locS),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,7),xlim=c(145,160),ylim=c(-45,-25),main="South of average",cex.axis=1.5,cex.main=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb,cols)

