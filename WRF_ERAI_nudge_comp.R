rm(list=ls())
setwd('~/Documents/ECLs/ERAI-WRF')

erai_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
erai_E$Year=floor(erai_E$Date1/10000)
erai_E$Month=floor(erai_E$Date1/100) %% 100
erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")

wrfR2_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ERA-nudge/p100_rad2cv1/ECLevents_ERA-nudge_rad2_p100_cv1.0.csv",sep=";")
wrfR2_E$Year=floor(wrfR2_E$Date1/10000)
wrfR2_E$Month=floor(wrfR2_E$Date1/100) %% 100
wrfR2=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ERA-nudge/p100_rad2cv1/ECLfixes_ERA-nudge_rad2_p100_cv1.0.csv",sep=";")

wrfR2N_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ERA-nonudge/p100_rad2cv1/ECLevents_ERA-nonudge_rad2_p100_cv1.0.csv",sep=";")
wrfR2N_E$Year=floor(wrfR2N_E$Date1/10000)
wrfR2N_E$Month=floor(wrfR2N_E$Date1/100) %% 100
wrfR2N=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_ERA-nonudge/p100_rad2cv1/ECLfixes_ERA-nonudge_rad2_p100_cv1.0.csv",sep=";")

ann=cbind(1990:2009,matrix(0,20,3))
colnames(ann)=c("Year","ERAI","WRFR2 nudge","WRFR2 nonudge")
for(i in 1:20)
{
  ann[i,2]=length(which(erai_E$Year==ann[i,1]))
  ann[i,3]=length(which(wrfR2_E$Year==ann[i,1]))
  ann[i,4]=length(which(wrfR2N_E$Year==ann[i,1]))
}
mon=cbind(1:12,matrix(0,12,3))
colnames(mon)=c("Year","ERAI","WRFR2 nudge","WRFR2 nonudge")
for(i in 1:12)
{
  mon[i,2]=length(which(erai_E$Month==i & erai_E$Year>=1990 & erai_E$Year<=2009))/20
  mon[i,3]=length(which(wrfR2_E$Month==i))/20
  mon[i,4]=length(which(wrfR2N_E$Month==i))/20
}

plot(mon[,1],mon[,2],type="l",lwd=3,ylim=c(0,5))
for(i in 3:4) lines(mon[,1],mon[,i],lwd=3,col=i)

stats=rbind(apply(erai_E[erai_E$Year>=1990 & erai_E$Year<=2009,c(3,6:10)],2,mean),
            apply(wrfR2_E[,c(3,6:10)],2,mean),apply(wrfR2N_E[,c(3,6:10)],2,mean))

dates=seq(as.Date("1990/1/1"), as.Date("2009/12/31"),by="day",format="%Y%m%d")
dates2=cbind(as.numeric(format.Date(dates,"%Y%m%d")),matrix(0,length(dates),3))
tmp=sort(unique(erai$Date[erai$Location==1]))
I=which(dates2[,1] %in% tmp)
dates2[I,2]=1
tmp=sort(unique(wrfR2$Date[wrfR2$Location==1]))
I=which(dates2[,1] %in% tmp)
dates2[I,3]=1
tmp=sort(unique(wrfR2N$Date[wrfR2N$Location==1]))
I=which(dates2[,1] %in% tmp)
dates2[I,4]=1

eracomp=matrix(0,4,2)
rownames(eracomp)=c("Cor","HR","FAR","CSI")
colnames(eracomp)=c("Nudge","No-nudge")
eracomp[1,1]=cor(ann[,3],ann[,2])
eracomp[1,2]=cor(ann[,4],ann[,2])

for(i in 1:2)
{
  eracomp[2,i]=length(which(dates2[,2]==1 & dates2[,i+2]==1))/length(which(dates2[,i+2]==1))
  eracomp[3,i]=length(which(dates2[,2]==1 & dates2[,i+2]==0))/length(which(dates2[,2]==1))
  eracomp[4,i]=length(which(dates2[,2]==1 & dates2[,i+2]==1))/length(which(dates2[,i+2]==1 | dates2[,2]==1))
}

###Okay, compare locations

 lat=seq(-40,-25,2.5)
 lon=seq(145,160,2.5)
 loc<-array(0,c(7,7,3))


for(i in 1:7)
  for(j in 1:7)
    {
      I=which(erai$Lat>=lat[i]-1.25 & erai$Lat<lat[i]+1.25 & erai$Lon>=lon[j]-1.25 & erai$Lon<lon[j]+1.25 & erai$Date>=19900000 & erai$Date<=200100000)
      loc[i,j,1]=length(I)
      I=which(wrfR2$Lat>=lat[i]-1.25 & wrfR2$Lat<lat[i]+1.25 & wrfR2$Lon>=lon[j]-1.25 & wrfR2$Lon<lon[j]+1.25)
      loc[i,j,2]=length(I)
      I=which(wrfR2N$Lat>=lat[i]-1.25 & wrfR2N$Lat<lat[i]+1.25 & wrfR2N$Lon>=lon[j]-1.25 & wrfR2N$Lon<lon[j]+1.25)
      loc[i,j,3]=length(I)
  }
 
loc=loc/20

cols=gray(seq(1,0.1,-0.15))
ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}

bb=c(-0.5,0,1,2,3,4,5,100)
names=c("ERAI","WRFR2 nudge","WRFR2 nonudge")
library("R.matlab")
readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

pdf(file=paste("ECL_locations_vERAI.pdf",sep=""),width=11,height=4,pointsize=12)
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))
for(i in 1:3)
{
image(lon,lat,t(loc[,,i]),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=names[i],cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb,cols)
dev.off()

#### Now look at CV

cvthresh=seq(1,4,0.25)
dens<-matrix(0,12,3)
  for(i in 1:12)
  {
    I=which(erai_E$CV2>=cvthresh[i] & erai_E$CV2<cvthresh[i+1] & erai_E$Year>=1990 & erai_E$Year<=2009)
    dens[i,1]=length(I)
    I=which(wrfR2_E$CV2>=cvthresh[i] & wrfR2_E$CV2<cvthresh[i+1])
    dens[i,2]=length(I)
    I=which(wrfR2N_E$CV2>=cvthresh[i] & wrfR2N_E$CV2<cvthresh[i+1])
    dens[i,3]=length(I)
  }

pdf(file="ECL_CVpdf_d01_vERAI.pdf",width=6,height=4,pointsize=12)
plot(cvthresh[1:12],dens[,1],type="l",lwd=3,xlab="Intensity",ylab="Frequency",ylim=range(dens))
for(i in 2:3) lines(cvthresh[1:12],dens[,i],lwd=3,col=i)
legend("topright",names,lwd=3,col=1:3)
dev.off()
