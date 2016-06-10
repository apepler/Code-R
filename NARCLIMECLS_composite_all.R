rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)
cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")
namelist=rep("aaa",12)

now<-array(NaN,c(12,9,2,4))
dimnames(now)[[1]]<-rep("ncep",12)
dimnames(now)[[2]]=c("Count","Mean(meanR)","Mean(maxR)","Mean(meanW)","Mean(maxW)","MeanRain>=6","MaxRain>=50","MeanWind>=50km/h","MaxWind>=80km/h")
dimnames(now)[[3]]=c("Cool","Warm")
dimnames(now)[[4]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now


#### Four methods
#### 22 ECLs p.a, by season
#### Count, Mean(meanR), mean(maxR), MeanR>=5, MaxR>=50


thresh=22
intcol=c(10,10,13,13) # CVmax/PG200max
impcol=cbind(c(14:17),c(14:17),c(15:18),c(15:18)) # MaxMeanR, MaxMaxR, MaxMeanW, MaxMaxW
#impcol=cbind(c(19:22),c(19:22),c(20:23),c(20:23)) # MaxMeanR, MaxMaxR, MaxMeanW, MaxMaxW
xthresh=c(6,50,13.9,22.2)


n=1
for(i in 1:4)
  for(j in 1:3)
  {
    namelist[n]<-paste(cmip[i],wrf[j])
    
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.2_impacts_v2.csv",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelist1[ff])
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      data=data[data[,intcol[ff]]>=thresh2,]
      mm=floor((data$Date1%%10000)/100)
      data2=read.csv(filelist2[ff])
      data2=data2[data2[,intcol[ff]]>=thresh2,]
      mm2=floor((data2$Date1%%10000)/100)
      data[data==-Inf | data==Inf]=NaN
      data2[data2==-Inf | data2==Inf]=NaN
      
      ### Count & mean
      I=which(mm>=5 & mm<=10)
      now[n,1,1,ff]=length(I)
      for(x in 1:4) 
        {  
        now[n,x+1,1,ff]=mean(data[I,impcol[x,ff]],na.rm=T)
        now[n,x+5,1,ff]=length(which(data[I,impcol[x,ff]]>=xthresh[x]))
      }
      
      I=which(mm2>=5 & mm2<=10)
      future[n,1,1,ff]=length(I)
      for(x in 1:4) 
      {  
        future[n,x+1,1,ff]=mean(data2[I,impcol[x,ff]],na.rm=T)
        future[n,x+5,1,ff]=length(which(data2[I,impcol[x,ff]]>=xthresh[x]))
      }

      I=which(mm>=11 | mm<=4)
      now[n,1,2,ff]=length(I)
      for(x in 1:4) 
      {  
        now[n,x+1,2,ff]=mean(data[I,impcol[x,ff]],na.rm=T)
        now[n,x+5,2,ff]=length(which(data[I,impcol[x,ff]]>=xthresh[x]))
      }
      
      I=which(mm2>=11 | mm2<=4)
      future[n,1,2,ff]=length(I)
      for(x in 1:4) 
      {  
        future[n,x+1,2,ff]=mean(data2[I,impcol[x,ff]],na.rm=T)
        future[n,x+5,2,ff]=length(which(data2[I,impcol[x,ff]]>=xthresh[x]))
      }
    }
    
    n=n+1  
  }

change=100*((future/now)-1)
names(dimnames(change))=c("Source","Statistic","Season","Method")

stat=c("Count","MeanRmean","MeanRmax","MeanWmean","MeanWmax",
       "MeanR6","MaxR50","MeanW50","MaxW80")

for(i in 1:9)
{
  pdf(file=paste("ECLevents_ALL_thresh",thresh,"_change_byseasonmethod_",stat[i],"_500km.pdf",sep=""),width=8,height=4,pointsize=10)
  tmp=melt(change[,i,,])
  print(ggplot(tmp, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-500, 500, 25)) +
          theme_bw() + ylab("Percentage change") + xlab("") +  geom_hline(yintercept = 0))
  dev.off()
  pdf(file=paste("ECLevents_ALL_thresh",thresh,"_change_byseasonmethod_",stat[i],"_500km_trunc.pdf",sep=""),width=8,height=4,pointsize=10)
  print(ggplot(tmp, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 100, 25),limits = c(-50, 100)) +
          theme_bw() + ylab("Percentage change") + xlab("") +  geom_hline(yintercept = 0))
  
  dev.off()
}


######
##  Awkward version for multipanel

library(abind)

now2=abind(array(NaN,c(12,9,2)),now,along=4)
future2=abind(array(NaN,c(12,9,2)),future,along=4)
dimnames(now2)[[4]]=c("ULGV","LAP 150km","LAP 50km","PG 150km","PG 50km")
dimnames(future2)[[4]]=c("ULGV","LAP 150km","LAP 50km","PG 150km","PG 50km")

n=1
for(i in 1:4)
  for(j in 1:3)
  {
    data=read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""))
    b=order(data$MaxGV,decreasing=T)
    thresh2=data$MaxGV[b[20*thresh]]
    if(is.na(thresh2)) thresh2=min(data$MaxGV,na.rm=T)
    data=data[data$MaxGV>=thresh2,]
    mm=floor((data$Date1%%10000)/100)
    data2=read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""))
    data2=data2[data2$MaxGV>=thresh2,]
    mm2=floor((data2$Date1%%10000)/100)
    ### Count & mean
    I=which(mm>=5 & mm<=10)
    now2[n,1,1,1]=length(I)
    now2[n,1,2,1]=length(mm)-length(I)
    I=which(mm2>=5 & mm2<=10)
    future2[n,1,1,1]=length(I)
    future2[n,1,2,1]=length(mm2)-length(I)
    n=n+1  
  }

change=100*((future2/now2)-1)
names(dimnames(change))=c("Source","Statistic","Season","Method")
stat=c("Count","MeanRmean","MeanRmax","MeanWmean","MeanWmax",
       "MeanR10","MaxR50","MeanW50","MaxW80")

tmp=melt(change[,1,,])
plot1<-(ggplot(tmp, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("a) All ECLs") + #stat_boxplot(coef=4) +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-500, 500, 25)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Change in number of lows (%)") + xlab("") +  geom_hline(yintercept = 0))
tmp=melt(change[,6,,])
plot2<-(ggplot(tmp, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("b) Heavy rain (areal average)") + #stat_boxplot(coef=4) +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-500, 500, 50)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Change in number of lows (%)") + xlab("") +  geom_hline(yintercept = 0))
tmp=melt(change[,7,,])
plot3<-(ggplot(tmp, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("c) Heavy rain (grid point)") +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + #stat_boxplot(coef=4) + 
          scale_y_continuous(breaks=seq(-500, 500, 50)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Change in number of lows (%)") + xlab("") +  geom_hline(yintercept = 0))
tmp=melt(change[,8,,])
plot4<-(ggplot(tmp, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("d) Strong wind (areal average)") + #stat_boxplot(coef=4) +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-500, 500, 50)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Change in number of lows (%)") + xlab("") +  geom_hline(yintercept = 0))

library(grid)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf(file=paste("ECLevents_ALL_thresh",thresh,"_change_byseasonmethod_500km_panel_v8.pdf",sep=""),width=5,height=11,pointsize=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,1)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(2, 1))
print(plot3, vp = vplayout(3, 1))
print(plot4, vp = vplayout(4, 1))
dev.off()


##### Change in PDF of rain associated with ECLs
##### Mean or max?
##### What range, 0-100mm?

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)
cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")
namelist=rep("aaa",12)

now<-array(NaN,c(12,512,2,4))
dimnames(now)[[1]]<-rep("ncep",12)
dimnames(now)[[3]]=c("Cool","Warm")
dimnames(now)[[4]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now

thresh=22
intcol=c(10,10,13,13) # CVmax/PG200max
intcolF=c(10,10,10,10)
denscol=c(14,14,16,16)

n=1
for(i in 1:4)
  for(j in 1:3)
  {
    namelist[n]<-paste(cmip[i],wrf[j])
    filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.2_impacts_v2.csv",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelistE[ff])
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      
      data=read.csv(filelist1[ff])
      data=data[(data[,intcolF[ff]]>=thresh2 & data$Location==1),]
      mm=floor((data$Date%%10000)/100)
      data2=read.csv(filelist2[ff])
      data2=data2[(data2[,intcolF[ff]]>=thresh2 & data2$Location==1),]
      mm2=floor((data2$Date%%10000)/100)
      data[data==-Inf | data==Inf]=NaN
      data2[data2==-Inf | data2==Inf]=NaN
      
      ### Density
      I=which(mm>=5 & mm<=10)
      a=density(data[I,denscol[ff]],from=0,to=50,na.rm=T)
      now[n,,1,ff]=a$y
      a=density(data[-I,denscol[ff]],from=0,to=50,na.rm=T)
      now[n,,2,ff]=a$y      
      
      I=which(mm2>=5 & mm2<=10)
      a=density(data2[I,denscol[ff]],from=0,to=50,na.rm=T)
      future[n,,1,ff]=a$y
      a=density(data2[-I,denscol[ff]],from=0,to=50,na.rm=T)
      future[n,,2,ff]=a$y      
    }
    n=n+1  
  }

change=future-now
change2=apply(change,c(2,3,4),median)

plot(1,NA,xlim=c(0,20),ylim=range(0,change),xlab="Mean rainfall (mm)",ylab="Density")
abline(h=0,col="black")
for(i in 1:4) 
  for(j in 1:12) lines(a$x,change[j,,1,i],lwd=3,col=i+1)
legend("topright",c("LAP 150km","LAP 50km","PG 150km","PG 50km"),lwd=3,col=2:5)

change3=apply(change,c(2,3),quantile,c(0.025,0.25,0.5,0.75,0.95))
plot(1,NA,xlim=c(0,20),ylim=range(0,change3[,,1]),xlab="Mean rainfall (mm)",ylab="Density")
polygon(c(a$x, rev(a$x)), c(change3[1,,1], rev(change3[5,,1])), col=rgb(0.8,0.8,0.8,0.3), border=NA)
polygon(c(a$x, rev(a$x)), c(change3[2,,1], rev(change3[4,,1])), col=rgb(0.4,0.4,0.4,0.3), border=NA)
abline(h=0,col="black",lty=2)
lines(a$x,change3[3,,1],lwd=3,col=rgb(0.3,0.3,0.3,1))

plot(1,NA,xlim=c(0,20),ylim=range(0,change3[,,1]),xlab="Mean rainfall (mm)",ylab="Density")
polygon(c(a$x, rev(a$x)), c(change3[1,,2], rev(change3[5,,2])), col=rgb(0.8,0.8,0.8,0.3), border=NA)
polygon(c(a$x, rev(a$x)), c(change3[2,,2], rev(change3[4,,2])), col=rgb(0.4,0.4,0.4,0.3), border=NA)
abline(h=0,col="black",lty=2)
lines(a$x,change3[3,,2],lwd=3,col=rgb(0.3,0.3,0.3,1))


cdfN<-cdfF<-array(0,dim(now))
for(i in 2:512) cdfN[,i,,]=apply(now[,1:i,,],c(1,3,4),sum)
for(i in 2:512) cdfF[,i,,]=apply(future[,1:i,,],c(1,3,4),sum)
change4=apply(cdfF-cdfN,c(2,3),quantile,c(0.025,0.25,0.5,0.75,0.95))
plot(1,NA,xlim=c(0,50),ylim=range(0,change4),xlab="Mean rainfall (mm)",ylab="Density")
polygon(c(a$x, rev(a$x)), c(change4[1,,2], rev(change4[5,,2])), col=rgb(0.8,0.8,0.8,0.3), border=NA)
polygon(c(a$x, rev(a$x)), c(change4[2,,2], rev(change4[4,,2])), col=rgb(0.4,0.4,0.4,0.3), border=NA)
abline(h=0,col="black",lty=2)
lines(a$x,change4[3,,2],lwd=3,col=rgb(0.3,0.3,0.3,1))


############ Adding the spatial plot to this script so I know where it lives.

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)
cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")
namelist=rep("aaa",16)
library(abind)
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

ColorBar <- function(brks,cols,labels=NA)
{
  par(mar = c(3, 1, 3, 4), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  if(is.na(labels[1])) labels=brks[seq(2, length(brks)-1)]
  
  axis(4, at = seq(1.5, length(brks) - 1.5), tick = TRUE, 
       labels = labels)
}
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(10)
bb=c(-10000,-50,-25,-10,-5,0,5,10,25,50,10000)
blab=c("-50%","-25%","-10%","-5%","","+5%","+10%","+25%","+50%")

thresh=22
intcol=c(10,10,13,13) # CVmax/PG200max
intcolF=c(10,10,10,10)
denscol=c(14,14,16,16)

res=c(2.5,2,1.5,1)
latlist=list(seq(-40,-25,2.5),seq(-40,-24,2),seq(-40.5,-24,1.5),seq(-41,-24,1))
lonlist=list(seq(147.5,160,2.5),seq(148,160,2),seq(147,160.5,1.5),seq(148,161,1))

for(r in 2)
{
  lat=latlist[[r]]
  lon=lonlist[[r]]
  now<-array(NaN,c(length(lat),length(lon),20,16,2,4))
  dimnames(now)[[4]]<-rep("ncep",16)
  dimnames(now)[[5]]=c("Cool","Warm")
  dimnames(now)[[6]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
  names(dimnames(now))<-c("Lat","Lon","Year","Source","Month","Method")
  future<-now
  
  
  n=1
for(i in 1:4)
  for(j in 1:3)
  {
    namelist[n]<-paste(cmip[i],wrf[j])
    filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.2_impacts_v2.csv",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelistE[ff])
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      
      data=read.csv(filelist1[ff])
      data=data[(data[,intcolF[ff]]>=thresh2 & data$Location==1),]
      yy=floor(data$Date/10000)
      data2=read.csv(filelist2[ff])
      data2=data2[(data2[,intcolF[ff]]>=thresh2 & data2$Location==1),]
      yy2=floor(data2$Date/10000)
      data[data==-Inf | data==Inf]=NaN
      data2[data2==-Inf | data2==Inf]=NaN
      
      years1=unique(yy)
      years2=unique(yy2)
      
      ### Density
      for(y in 1:length(lat))
        for(x in 1:length(lon))
          for(t in 1:20)
          {
            I=which(data$Lat>=(lat[y]-res[r]/2) & data$Lat<=(lat[y]+res[r]/2) & data$Lon<=(lon[x]+res[r]/2) & data$Lon>=(lon[x]-res[r]/2) & yy==years1[t])
            mm=floor((data[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            now[y,x,t,n,1,ff]=length(I)
            now[y,x,t,n,2,ff]=length(mm)-length(I)     
            
            I=which(data2$Lat>=(lat[y]-res[r]/2) & data2$Lat<=(lat[y]+res[r]/2) & data2$Lon<=(lon[x]+res[r]/2) & data2$Lon>=(lon[x]-res[r]/2) & yy2==years2[t])
            mm=floor((data2[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            future[y,x,t,n,1,ff]=length(I)
            future[y,x,t,n,2,ff]=length(mm)-length(I)                 
          }
    }
    n=n+1  
  }

now2=apply(now[,,,1:12,,],c(1,2,4,5,6),sum)
future2=apply(future[,,,1:12,,],c(1,2,4,5,6),sum)

pval<-array(NaN,dim(now2))

for(i in 1:length(lat))
  for(j in 1:length(lon))
    for(k in 1:12)
      for(l in 1:2)
        for(m in 1:4)
        {
          a=t.test(now[i,j,,k,l,m],future[i,j,,k,l,m])
          pval[i,j,k,l,m]=a$p.value
        }

change=100*((future2/now2)-1)
changeM=apply(100*((future2/now2)-1),c(1,2,4),median)
pval2=abind(pval[,,,,1],pval[,,,,2],pval[,,,,3],pval[,,,,4],along=3)
change2=abind(change[,,,,1],change[,,,,2],change[,,,,3],change[,,,,4],along=3)

sigM=array(NaN,c(length(lat),length(lon),2,4))
for(i in 1:length(lat))
  for(j in 1:length(lon))
    for(k in 1:2)
    {
      I=which(pval2[i,j,,k]<0.05)
      sigM[i,j,k,1]=length(I)/48
      if(!is.na(changeM[i,j,k])) sigM[i,j,k,2]=length(which(change2[i,j,,k]<0))/48
      if(!is.na(changeM[i,j,k])) sigM[i,j,k,3]=length(which(change2[i,j,,k]>0))/48
      if(length(I)>0) sigM[i,j,k,4]=length(which(change2[i,j,I,k]<0))/length(I)
    }

pdf(file=paste("ECLlocation_ALL_thresh",thresh,"_change_season_prop75_res",res[r],".pdf",sep=""),width=7,height=4.2,pointsize=10)
layout(cbind(1,2,3),width=c(1,1,0.4))
image(lon,lat,t(changeM[,,1]),xlab="",ylab="",breaks=bb,col=cm,
      zlim=c(-1000,1000),xlim=c(149,161),ylim=c(-41,-24),cex.axis=1.5,cex.main=2,axes=F)
mtext("a) Cool", side=3, adj=0, cex=1.5) 
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
axis(1,at=seq(145,165,5),cex.axis=1.5,labels=c("",expression(150~degree),expression(155~degree),expression(160~degree),""))
axis(2,at=seq(-45,-20,5),cex.axis=1.5,labels=c("",expression(40~degree),expression(35~degree),expression(30~degree),expression(25~degree),""))
sigmask=which(sigM[,,1,2]>=0.75 | sigM[,,1,3]>=0.75,arr.ind=T)
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=2,lwd=2)

image(lon,lat,t(changeM[,,2]),xlab="",ylab="",breaks=bb,col=cm,
      zlim=c(-1000,1000),xlim=c(149,161),ylim=c(-41,-24),cex.axis=1.5,cex.main=2,axes=F)
mtext("b) Warm", side=3, adj=0, cex=1.5) 
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
axis(1,at=seq(145,165,5),cex.axis=1.5,labels=c("",expression(150~degree),expression(155~degree),expression(160~degree),""))
axis(2,at=seq(-45,-20,5),cex.axis=1.5,labels=c("",expression(40~degree),expression(35~degree),expression(30~degree),expression(25~degree),""))
sigmask=which(sigM[,,2,2]>=0.75 | sigM[,,2,3]>=0.75,arr.ind=T)
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=2,lwd=2)
ColorBar(bb,cm,labels=blab)
dev.off()

changeM2=abind(apply(100*((future2[,,,,1:2]/now2[,,,,1:2])-1),c(1,2,4),median),
               apply(100*((future2[,,,,3:4]/now2[,,,,3:4])-1),c(1,2,4),median),along=3)

sigM=array(NaN,c(length(lat),length(lon),4,2))
for(i in 1:length(lat))
  for(j in 1:length(lon))
    for(k in 1:2)
    {
      if(!is.na(changeM2[i,j,k])) sigM[i,j,k,1]=length(which(change2[i,j,1:24,k]<0))/24
      if(!is.na(changeM2[i,j,k])) sigM[i,j,k,2]=length(which(change2[i,j,1:24,k]>0))/24
      if(!is.na(changeM2[i,j,k+2])) sigM[i,j,k+2,1]=length(which(change2[i,j,25:48,k]<0))/24
      if(!is.na(changeM2[i,j,k+2])) sigM[i,j,k+2,2]=length(which(change2[i,j,25:48,k]>0))/24
    }

names=c("a) LAP (Cool)","c) LAP (Warm)","b) PG (Cool)","d) PG (Warm)")
pdf(file=paste("ECLlocation_ALL_thresh",thresh,"_change_season_prop75_res",res[r],"_bymethod2.pdf",sep=""),width=7,height=8.4,pointsize=10)
layout(rbind(cbind(1,2,5),c(3,4,5)),width=c(1,1,0.4))
for(i in 1:4)
{
image(lon,lat,t(changeM2[,,i]),xlab="",ylab="",breaks=bb,col=cm,
      zlim=c(-1000,1000),xlim=c(149,161),ylim=c(-41,-24),cex.axis=1.5,cex.main=2,axes=F)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
axis(1,at=seq(145,165,5),cex.axis=1.5,labels=c("",expression(150~degree),expression(155~degree),expression(160~degree),""))
axis(2,at=seq(-45,-20,5),cex.axis=1.5,labels=c("",expression(40~degree),expression(35~degree),expression(30~degree),expression(25~degree),""))
mtext(names[i], side=3, adj=0, cex=1.5) 
sigmask=which(sigM[,,i,1]>=0.75 | sigM[,,i,2]>=0.75,arr.ind=T)
points(lon[sigmask[,2]],lat[sigmask[,1]],col="black",pch=4,cex=2,lwd=2)
}
ColorBar(bb,cm,labels=blab)
dev.off()


}

#### Now I need to plot/validate the biases in ECL frequency

i=5
  for(j in 1:3)
  {
    namelist[n]<-paste(cmip[i],wrf[j])
    filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_final.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_final.csv",sep=""))
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_final.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_final.csv",sep=""))
    for(ff in 1:4)
    {
      data=read.csv(filelistE[ff])
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      
      data=read.csv(filelist1[ff])
      data=data[(data[,intcolF[ff]]>=thresh2 & data$Location==1),]
      yy=floor(data$Date/10000)
      data[data==-Inf | data==Inf]=NaN
      years1=unique(yy)
      
      ### Density
      for(y in 1:18)
        for(x in 1:14)
          for(t in 1:20)
          {
            I=which(data$Lat>=(lat[y]-0.5) & data$Lat<=(lat[y]+0.5) & data$Lon<=(lon[x]+0.5) & data$Lon>=(lon[x]-0.5) & yy==years1[t])
            mm=floor((data[I,]$Date%%10000)/100)
            I=which(mm>=5 & mm<=10)
            now[y,x,t,n,1,ff]=length(I)
            now[y,x,t,n,2,ff]=length(mm)-length(I)     
          }
    }
    n=n+1  
  }
intcolE=c(12,12,13,13)
intcolF=c(10,10,10,10)

dimnames(nowC)[[3]][n]<-"ERAI"
filelistE=c("~/Documents/ECLs/Algorithm Comparison/events_UM_rad2_p100.csv",NaN,
            paste("Alejandro/final/ECLevents_Alejandro_ERAI_res150_9009_pg0.7_final.csv",sep=""))
filelist1=c("~/Documents/ECLs/Algorithm Comparison/UM_rad2_p100.csv",NaN,
            paste("Alejandro/final/ECLfixes_Alejandro_ERAI_res150_9009_pg0.7_final.csv",sep=""))

for(ff in c(1,3))
{
  data=read.csv(filelistE[ff])
  data=data[data$Date1>=19900000 & data$Date2<=20100000,]
  b=order(data[,intcolE[ff]],decreasing=T)
  thresh2=data[b[20*thresh],intcolE[ff]]
  if(is.na(thresh2)) thresh2=min(data[,intcolE[ff]],na.rm=T)
  
  data=read.csv(filelist1[ff])
  data=data[(data[,intcolF[ff]]>=thresh2 & data$Location==1),]
  yy=floor(data$Date/10000)
  data[data==-Inf | data==Inf]=NaN
  years1=unique(yy)

  
  ### Density
  for(y in 1:18)
    for(x in 1:14)
      for(t in 1:20)
      {
        I=which(data$Lat>=(lat[y]-0.5) & data$Lat<=(lat[y]+0.5) & data$Lon<=(lon[x]+0.5) & data$Lon>=(lon[x]-0.5) & yy==years1[t])
        mm=floor((data[I,]$Date%%10000)/100)
        I=which(mm>=5 & mm<=10)
        now[y,x,t,n,1,ff]=length(I)
        now[y,x,t,n,2,ff]=length(mm)-length(I)     
      }
}

now2=apply(now,c(1,2,4,5,6),sum)
now3=abind(apply(now2[,,16,,c(1,3)],c(1,2,4),mean),
           apply(now2[,,13:15,,c(1,3)],c(1,2,4),mean),
           apply(now2[,,1:12,,c(1,3)],c(1,2,4),mean),along=4)
now3=now3/20

now4=apply(now,c(1,2,4,6),sum)
now5=abind(apply(now4[,,16,c(1,3)],c(1,2),mean),
           apply(now4[,,13:15,],c(1,2),mean),
           apply(now4[,,1:12,],c(1,2),mean),along=3)
now5=now5/20

cm2=gray(seq(1,0.1,-0.15))
#bb2=cbind(c(-0.5,0,0.2,0.5,1,1.5,2,100),c(-0.5,0,0.5,1,1.5,2,3,100),c(-0.5,0,1,2,3,4,5,100),c(-0.5,0,1,2.5,5,7.5,10,100),c(-0.5,0,1,2.5,5,7.5,10,100))
#bb2=c(-0.5,0,0.5,1,1.5,2,3,100)
bb2=c(-0.5,0,1,2,3,4,5,100)
bb2=c(-0.5,0,0.2,0.4,0.6,0.8,1,100)
pdf(file=paste("ECLlocation_ALL_thresh",thresh,"_mean_v2_fixed.pdf",sep=""),width=10,height=4.2,pointsize=10)
layout(cbind(1,2,3,4),width=c(1,1,1,0.4))
names=c("ERAI","NCEP-WRF","CMIP-WRF")
for(i in 1:3)
{
image(lon,lat,t(now5[,,i]),xlab="",ylab="",breaks=bb2,col=cm2,main=names[i],
      zlim=c(-1000,1000),xlim=c(149,161),ylim=c(-41,-24),cex.axis=1.5,cex.main=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb2,cm2)
dev.off()



##### Actual rainfall patterns

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)
cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")

nowC<-array(NaN,c(21,21,16,4))
dimnames(nowC)[[3]]<-rep("ncep",16)
dimnames(nowC)[[4]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(nowC))<-c("Lon","Lat","Source","Method")

nowW<-nowC
#### Four methods
#### 22 ECLs p.a, by season
#### Count, Mean(meanR), mean(maxR), MeanR>=5, MaxR>=50


thresh=22
intcolE=c(10,10,13,13) # CVmax/PG200max
intcolF=c(10,10,10,10) # CVmax/PG200max

n=1

for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(nowC)[[3]][n]<-dimnames(nowW)[[3]][n]<-paste(cmip[i],wrf[j])
    
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelistC=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_compositerain_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.nc",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_compositerain_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.nc",sep=""),
                paste("Alejandro/final/ECLfixes_compositerain_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7.nc",sep=""),
                paste("Alejandro/final/ECLfixes_compositerain_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2.nc",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelist1[ff])
      b=order(data[,intcolE[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcolE[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcolE[ff]],na.rm=T)
      
      data=read.csv(filelist2[ff])
      mm=floor((data$Date%%10000)/100)
      a=open.nc(filelistC[ff])
      tmp=var.get.nc(a,"ECLrain")
      tmp[tmp>=600]=NaN      
      
      I=which(data$Location==1 & data[,intcolF[ff]]>=thresh2 & (mm>=11 | mm<=4))
      nowW[,,n,ff]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
      I=which(data$Location==1 & data[,intcolF[ff]]>=thresh2 & (mm>=5 & mm<=10))
      nowC[,,n,ff]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
    }
    
    n=n+1  
  }
i=5
for(j in 1:3)
{
  dimnames(nowC)[[3]][n]<-dimnames(nowW)[[3]][n]<-paste(cmip[i],wrf[j])
  
  filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
              paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
              paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_final.csv",sep=""),
              paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_final.csv",sep=""))
  filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
              paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
              paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_final.csv",sep=""),
              paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_final.csv",sep=""))
  filelistC=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_compositerain_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.nc",sep=""),
              paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_compositerain_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.nc",sep=""),
              paste("Alejandro/final/ECLfixes_compositerain_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7.nc",sep=""),
              paste("Alejandro/final/ECLfixes_compositerain_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2.nc",sep=""))
  
  for(ff in 1:4)
  {
    data=read.csv(filelist1[ff])
    b=order(data[,intcolE[ff]],decreasing=T)
    thresh2=data[b[20*thresh],intcolE[ff]]
    if(is.na(thresh2)) thresh2=min(data[,intcolE[ff]],na.rm=T)
    
    data=read.csv(filelist2[ff])
    mm=floor((data$Date%%10000)/100)
    a=open.nc(filelistC[ff])
    tmp=var.get.nc(a,"ECLrain")
    tmp[tmp>=600]=NaN      
    
    I=which(data$Location==1 & data[,intcolF[ff]]>=thresh2 & (mm>=11 | mm<=4))
    nowW[,,n,ff]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
    I=which(data$Location==1 & data[,intcolF[ff]]>=thresh2 & (mm>=5 & mm<=10))
    nowC[,,n,ff]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
  }
  n=n+1  
}

library(abind)
nowM1<-abind(apply(nowC[,,1:12,],c(1,2,4),mean),apply(nowW[,,1:12,],c(1,2,4),mean),along=4)
nowM2<-abind(apply(nowC,c(1,2,3),mean),apply(nowW,c(1,2,3),mean),along=4)

source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
library(fields)
ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = brks[seq(1, length(brks), subsampleg)])
}


namelist<-c("LAP 150km","LAP 50km","PG 150km","PG 50km")
pdf(file="ECL_compositerain_mean_bymethod_warm.pdf",width=6.5,height=6)
layout(cbind(c(1,2),c(3,4),c(5,5)),width=c(1,1,0.4))
par(mar=c(1,1,3,1))
for(i in 1:4) 
{
  image(nowM1[,,i,2],breaks=c(seq(0,14,1),100),col=tim.colors(15),main=namelist[i],axes=F,cex.main=2)
  points(0.5,0.5,pch=4,col="black",lwd=3,cex=3)
}
ColorBar(c(seq(0,14,1),100),tim.colors(15),subsampleg=2)
dev.off()

namelist<-dimnames(nowC)[[3]]
pdf(file="ECL_compositerain_mean_bymodel_warm_v2.pdf",width=12,height=7)
layout(cbind(1:3,4:6,7:9,10:12,13:15,c(16,16,16)),width=c(1,1,1,1,1,0.6))
par(mar=c(1,1,3,1))
for(i in 1:15) 
{
  image(nowM2[,,i,2],breaks=c(seq(0,14,1),100),col=tim.colors(15),main=namelist[i],axes=F,cex.main=2)
  points(0.5,0.5,pch=4,col="black",lwd=3,cex=3)
}
ColorBar(c(seq(0,14,1),100),tim.colors(15),subsampleg=2)
dev.off()


######### Re-do the location validation.





### Re-do the location validation
intcolE=c(12,12,13,13)
intcolF=c(10,10,10,10)

dimnames(nowC)[[3]][n]<-"ERAI"
filelist1=c("~/Documents/ECLs/Algorithm Comparison/events_UM_rad2_p100.csv",NaN,
            paste("Alejandro/final/ECLevents_Alejandro_ERAI_res150_9009_pg0.7_final.csv",sep=""))
filelist2=c("~/Documents/ECLs/Algorithm Comparison/UM_rad2_p100.csv",NaN,
            paste("Alejandro/final/ECLfixes_Alejandro_ERAI_res150_9009_pg0.7_final.csv",sep=""))

for(ff in c(1,3))
{
  data=read.csv(filelist1[ff])
  data=data[data$Date1>=19900000 & data$Date2<=20100000,]
  b=order(data[,intcolE[ff]],decreasing=T)
  thresh2=data[b[20*thresh],intcolE[ff]]
  if(is.na(thresh2)) thresh2=min(data[,intcolE[ff]],na.rm=T)
  
  I=which(data$Location==1 & data[,intcolF[ff]]>=thresh2 & (mm>=11 | mm<=4))
  nowW[,,n,ff]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
  I=which(data$Location==1 & data[,intcolF[ff]]>=thresh2 & (mm>=5 & mm<=10))
  nowC[,,n,ff]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
}
n=n+1  

