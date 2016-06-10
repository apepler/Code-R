rm(list=ls())
setwd('/home/nfs/z3478332/output/outputUM_wrf_2007_all/')

makePDF = function(data1,data2,xlabel,labloc="topleft") {
  a=density(data1,na.rm=T)
  b=density(data2,na.rm=T)
  
  lims=range(data1,data2,na.rm=T)
  if((lims[2]-lims[1])<10)
  {
    lims[1]=floor(lims[1])
    lims[2]=ceiling(lims[2])
  } else {
    lims[1]=floor(lims[1]/5)*5
    lims[2]=ceiling(lims[2]/5)*5
  }
  
  plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
       xlab=xlabel,ylab="Frequency",cex.main=1.2,
       main="2007-2008 ECL statistics")
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=c("Control","NoTopo"),
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
}

############

year=c(2007,2008)
domain=c("d01","d02")
cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
count<-countNT<-count_d02<-countNT_d02<-matrix(0,4,3)

for(c in 1:4)
for(r in 1:3)
  for(year in 2007:2008)
{

count[c,r]=count[c,r]+length(read.csv(paste("ECLevents_d01_",year,"_R",r,"_",cat[c],".csv",sep=""))[,1])
countNT[c,r]=countNT[c,r]+length(read.csv(paste("ECLevents_d01_",year,"_R",r,"_notopo_",cat[c],".csv",sep=""))[,1])
count_d02[c,r]=count_d02[c,r]+length(read.csv(paste("ECLevents_d02_",year,"_R",r,"_",cat[c],".csv",sep=""))[,1])
countNT_d02[c,r]=countNT_d02[c,r]+length(read.csv(paste("ECLevents_d02_",year,"_R",r,"_notopo_",cat[c],".csv",sep=""))[,1])
}

prop=(rbind((countNT/count)-1,(countNT_d02/count_d02)-1))*100
colnames(prop)=c("R1","R2","R3")
boxplot(prop,xlab="WRF version",ylab="% change when topography removed")
abline(h=0,col="red")

prop2=cbind(as.vector(prop[c(1,5),]),as.vector(prop[c(3,7),]),as.vector(prop[c(2,6),]),as.vector(prop[c(4,8),]))
colnames(prop2)=c("1.5deg 500km","1.5 deg 200km","0.5deg 500km","0.5 deg 200km")
boxplot(prop2,xlab="Settings",ylab="% change when topography removed")
abline(h=0,col="red")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
dom=c("d01","d02")
cvthresh=c(40,20,10,5,2)

cv<-cvNT<-array(0,c(8,3,5))

for(r in 1:3)
{
  n=1
  for(dom in c("d01","d02"))
    for(c in 1:4)
    {
      data=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                 read.csv(paste("ECLevents_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
      data2=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_notopo_",cat[c],".csv",sep="")),
                  read.csv(paste("ECLevents_",dom,"_2008_R",r,"_notopo_",cat[c],".csv",sep="")))
      
      a=order(data$CV2,decreasing=T)
      b=order(data2$CV2,decreasing=T)
      
      for(t in 1:5)
      {
        cv[n,r,t]=data$CV2[a[cvthresh[t]*2]]
        cvNT[n,r,t]=data2$CV2[b[cvthresh[t]*2]]
      }
      n=n+1 
    }
  
}

apply(cvNT-cv,3,median,na.rm=T)
diff=rbind((cvNT-cv)[,1,],(cvNT-cv)[,2,],(cvNT-cv)[,3,])
colnames(diff)=cvthresh
boxplot(diff,xlab="Events p.a.",ylab="CV difference when topography removed")
abline(h=0,col="red")

#### Comparing to ERAI - p100, rad2

rm(list=ls())
library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
setwd('~/Documents/ECLs/WRFruns/0708/')

wrfv=c("R1","R2","R3")
dom=c("d01","d02")
monthECLs<-matrix(0,24,9)

eventsE=read.csv("~/output/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
eventsE$Year=floor(eventsE$Date1/10000)
eventsE$Month=floor(eventsE$Date1/100)%%100

year=c(2007,2008)
n=1
for(y in 1:2)
  for(m in 1:12)
  {
    I=which(eventsE$Year==year[y] & eventsE$Month==m)
    monthECLs[n,1]=length(I)
    n=n+1
  }

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3
dir='~/output/outputUM_wrf_2007_all/'
events<-eventsNT<-list()

n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    events[[n]]=rbind(read.csv(paste(dir,"ECLevents_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                      read.csv(paste(dir,"ECLevents_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    
    eventsNT[[n]]=rbind(read.csv(paste(dir,"ECLevents_",dom,"_2007_R",r,"_notopo_",cat[c],".csv",sep="")),
                        read.csv(paste(dir,"ECLevents_",dom,"_2008_R",r,"_notopo_",cat[c],".csv",sep="")))
    eventsNT[[n]]$Year=floor(eventsNT[[n]]$Date1/10000)
    eventsNT[[n]]$Month=floor(eventsNT[[n]]$Date1/100)%%100
    
#     k=1
#     for(y in 1:2)
#       for(m in 1:12)
#       {
#         I=which(events[[n]]$Year==year[y] & events[[n]]$Month==m)
#         monthECLs[k,n+1]=length(I)
#         k=k+1
#       }
    
    n=n+1     
  }
monthECLs[,8]=apply(monthECLs[,2:4],1,mean)
monthECLs[,9]=apply(monthECLs[,5:7],1,mean)

monthECLsNT=monthECLs
for(n in 1:6)
{
k=1
for(y in 1:2)
  for(m in 1:12)
  {
    I=which(eventsNT[[n]]$Year==year[y] & eventsNT[[n]]$Month==m)
    monthECLsNT[k,n+1]=length(I)
    k=k+1
  }
}  
monthECLsNT[,8]=apply(monthECLsNT[,2:4],1,mean)
monthECLsNT[,9]=apply(monthECLsNT[,5:7],1,mean)

plot(1:24,monthECLs[,1],ylim=c(0,6),axes=F,lwd=3,type="l",col="black",xlab="Month",ylab="Number of ECLs")
lines(1:24,monthECLs[,8],lwd=3,col="blue")
lines(1:24,monthECLsNT[,8],lwd=3,col="red")
axis(1,at=c(1,7,13,19,25),labels=c("Jan 07","Jul 07","Jan 08","Jul 08","Jan 09"))
axis(2,at=0:6)
legend("topright",c("ERAI","WRF 50km","WRF 50km NoTopo"),lwd=3,col=c("black","blue","red"),bty="n")

### ECL hours

rm(list=ls())
library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

wrfv=c("R1","R2","R3")
dom=c("d01","d02")
monthECLs<-monthECLsNT<-matrix(0,24,9)

fixesE=read.csv("~/output/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")
fixesE$Year=floor(fixesE$Date/10000)
fixesE$Month=floor(fixesE$Date/100)%%100

year=c(2007,2008)
n=1
for(y in 1:2)
  for(m in 1:12)
  {
    I=which(data$Year==year[y] & data$Month==m & data$Location==1)
    b=unique(data$Date[I]+(as.numeric(data$Time[I])-1)/24)
    monthECLs[n,1]<-monthECLsNT[n,1]<-length(b)*6
    n=n+1
  }

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3
dir='~/output/outputUM_wrf_2007_all/'
fixes<-fixesNT<-list()

n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    fixes[[n]]=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                      read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    
    fixesNT[[n]]=rbind(read.csv(paste(dir,"ECLfixes_",dom,"_2007_R",r,"_notopo_",cat[c],".csv",sep="")),
                        read.csv(paste(dir,"ECLfixes_",dom,"_2008_R",r,"_notopo_",cat[c],".csv",sep="")))
    fixesNT[[n]]$Year=floor(fixesNT[[n]]$Date/10000)
    fixesNT[[n]]$Month=floor(fixesNT[[n]]$Date/100)%%100
    
#     k=1
#     for(y in 1:2)
#       for(m in 1:12)
#       {
#         I=which(fixes[[n]]$Year==year[y] & fixes[[n]]$Month==m & fixes[[n]]$Location==1)
#         b=unique(fixes[[n]]$Date[I]+(as.numeric(fixes[[n]]$Time[I])-1)/24)
#         monthECLs[k,n+1]<-length(b)*6
#         I=which(fixesNT[[n]]$Year==year[y] & fixesNT[[n]]$Month==m & fixesNT[[n]]$Location==1)
#         b=unique(fixesNT[[n]]$Date[I]+(as.numeric(fixesNT[[n]]$Time[I])-1)/24)
#         monthECLsNT[k,n+1]<-length(b)*6        
# 
#         k=k+1
#       }
    
    n=n+1     
  }



k=1
for(y in 1:2)
  for(m in 1:12)
  {
    I=which(data$Year==year[y] & data$Month==m & data$Location==1 & data$CV)
    b=unique(data$Date[I]+(as.numeric(data$Time[I])-1)/24)
    monthECLs[k,1]<-monthECLsNT[k,1]<-length(b)*6
    
    for(n in 1:6)
    {
      I=which(fixes[[n]]$Year==year[y] & fixes[[n]]$Month==m & fixes[[n]]$Location==1 & fixes[[n]]$CV)
      b=unique(fixes[[n]]$Date[I]+(as.numeric(fixes[[n]]$Time[I])-1)/24)
      monthECLs[k,n+1]<-length(b)*6
      I=which(fixesNT[[n]]$Year==year[y] & fixesNT[[n]]$Month==m & fixesNT[[n]]$Location==1 & fixesNT[[n]]$CV)
      b=unique(fixesNT[[n]]$Date[I]+(as.numeric(fixesNT[[n]]$Time[I])-1)/24)
      monthECLsNT[k,n+1]<-length(b)*6              
    }    
    k=k+1
  }
monthECLs[,8]=apply(monthECLs[,2:4],1,mean)
monthECLs[,9]=apply(monthECLs[,5:7],1,mean)
monthECLsNT[,8]=apply(monthECLsNT[,2:4],1,mean)
monthECLsNT[,9]=apply(monthECLsNT[,5:7],1,mean)

plot(1:24,monthECLs[,1],ylim=c(0,250),axes=F,lwd=3,type="l",col="black",xlab="Month",ylab="Number of hours with ECLs present")
lines(1:24,monthECLs[,8],lwd=3,col="blue")
lines(1:24,monthECLs[,9],lwd=3,col="blue",lty=2)
axis(1,at=c(1,7,13,19,25),labels=c("Jan 07","Jul 07","Jan 08","Jul 08","Jan 09"))
axis(2,at=seq(0,50,250))
legend("topright",c("ERAI","WRF 50km","WRF 10km"),lwd=3,col=c("black","blue","blue"),lty=c(1,1,2),bty="n")

###### Okay, now on to ECL locations

fixes2=rbind(fixes[[1]],fixes[[2]],fixes[[3]])
fixes_notopo=rbind(fixesNT[[1]],fixesNT[[2]],fixesNT[[3]])

lat=seq(-40,-25,2.5)
lon=seq(145,160,2.5)
loc<-locNT<-locERAI<-matrix(0,7,7)

for(i in 1:7)
  for(j in 1:7)
  {
    I=which(fixesE$Lat>=lat[i]-1.25 & fixesE$Lat<lat[i]+1.25 & fixesE$Lon>=lon[j]-1.25 & fixesE$Lon<lon[j]+1.25 & fixesE$Year>=2007 & fixesE$Year<=2008)
    locERAI[i,j]=length(I)
    
    I=which(fixes2$Lat>=lat[i]-1.25 & fixes2$Lat<lat[i]+1.25 & fixes2$Lon>=lon[j]-1.25 & fixes2$Lon<lon[j]+1.25)
    loc[i,j]=length(I)/3
    
    I=which(fixes_notopo$Lat>=lat[i]-1.25 & fixes_notopo$Lat<lat[i]+1.25 & fixes_notopo$Lon>=lon[j]-1.25 & fixes_notopo$Lon<lon[j]+1.25)
    locNT[i,j]=length(I)/3
  }


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

bb=c(-0.5,0,1,2,5,10,15,100)
pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_",cat[c],"_poly.pdf",sep=""),width=12,height=4,pointsize=12)
layout(cbind(1,2,3,4),c(1,1,1,0.3))
image(lon,lat,t(locERAI),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="ERAI",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
polygon(c(149,149,152,152,157,157,154,154,149),c(-41,-37,-31,-24,-24,-31,-37,-41,-41),border="red",lwd=3)
image(lon,lat,t(loc),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="WRF 50km - Control",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,t(locNT),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="WRF 50km - No Topography",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb,cols)
dev.off()

locdiff=locNT-loc

bb2=c(-100,-5,-2,-1,1,2,5,100)
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
cm=pal(7)
cm[4]="white"

pdf(file=paste("ECL_locations_0708_",cat[c],"_change_d02.pdf",sep=""),width=4.5,height=4,pointsize=12)
layout(cbind(1,2),c(1,0.3))
image(lon,lat,t(locdiff),xlab="",ylab="",breaks=bb,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="Change in ECL frequency")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb2,cm)
dev.off() 

bb2=c(-100,-5,-2,-1,0,1,2,5,100)
cm=pal(8)
locdiff[locdiff==0]=NaN
pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_",cat[c],"_poly_change2.pdf",sep=""),width=13,height=4,pointsize=12)
layout(cbind(1,2,4,3,5),c(1,1,0.3,1,0.3))
par(mar=c(3,3,3,1)+0.1)
image(lon,lat,t(locERAI),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="ERAI",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
polygon(c(149,149,152,152,157,157,154,154,149),c(-41,-37,-31,-24,-24,-31,-37,-41,-41),border="red",lwd=3)
image(lon,lat,t(loc),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="WRF 50km - Control",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,t(locdiff),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="Change in ECL frequency",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb,cols)
ColorBar(bb2,cm)
dev.off()


lat=seq(-40,-25,2.5)
lon=seq(145,160,2.5)
locERAI<-matrix(0,7,7)
loc<-locNT<-array(0,c(7,7,6))

thresh=1.5

for(i in 1:7)
  for(j in 1:7)
  {
    I=which(fixesE$Lat>=lat[i]-1.25 & fixesE$Lat<lat[i]+1.25 & fixesE$Lon>=lon[j]-1.25 & fixesE$Lon<lon[j]+1.25)
    locERAI[i,j]=length(I)
    
    for(k in 1:6)
    {
    I=which(fixes[[k]]$Lat>=lat[i]-1.25 & fixes[[k]]$Lat<lat[i]+1.25 & fixes[[k]]$Lon>=lon[j]-1.25 & fixes[[k]]$Lon<lon[j]+1.25)
    loc[i,j,k]=length(I)
    
    I=which(fixesNT[[k]]$Lat>=lat[i]-1.25 & fixesNT[[k]]$Lat<lat[i]+1.25 & fixesNT[[k]]$Lon>=lon[j]-1.25 & fixesNT[[k]]$Lon<lon[j]+1.25)
    locNT[i,j,k]=length(I)
    }
  }

bb=c(-0.5,0,1,2,5,10,15,100)
pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_",cat[c],"_wrfv_d02.pdf",sep=""),width=8,height=12,pointsize=12)
layout(cbind(c(1,3,5),c(2,4,6),c(7,7,7)),height=c(1,1,1),width=c(1,1,0.3))
for(i in 1:3)
{
image(lon,lat,t(loc[,,i+3]),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=paste("R",i,sep=""),cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,t(locNT[,,i+3]),xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=paste("R",i," NoTopo",sep=""),cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb,cols)
dev.off()

#####Change by WRF version

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
bb2=c(-100,-5,-2,-1,0,1,2,5,100)
cm=pal(8)

lat=seq(-40,-25,5)
lon=seq(145,160,5)
loc<-locNT<-array(0,c(4,4,3))

for(r in 1:3)
{
  month=floor(fixes[[r]]$Date/100)%%100
  monthNT=floor(fixesNT[[r]]$Date/100)%%100
  
  for(i in 1:4)
    for(j in 1:4)
    {
      I=which(fixes[[r]]$Lat>=lat[i]-2.5 & fixes[[r]]$Lat<lat[i]+2.5 & fixes[[r]]$Lon>=lon[j]-2.5 & fixes[[r]]$Lon<lon[j]+2.5)
      loc[i,j,r]=length(I)
      
      I=which(fixesNT[[r]]$Lat>=lat[i]-2.5 & fixesNT[[r]]$Lat<lat[i]+2.5 & fixesNT[[r]]$Lon>=lon[j]-2.5 & fixesNT[[r]]$Lon<lon[j]+2.5)
      locNT[i,j,r]=length(I)
    }
}

locdiff=locNT-loc
locdiff[loc==0]=NaN

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_",cat[c],"_change_vwrf_v3.pdf",sep=""),width=12,height=4,pointsize=12)
layout(cbind(1,2,3,4),c(1,1,1,0.3))
par(mar=c(3,3,3,1)+0.1)
for(i in 1:3)
{
  image(lon,lat,t(locdiff[,,i]/2),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=paste("RCM",r),cex.axis=1.5,cex.main=1.5)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb2,cm)
dev.off()





### Distribution of intensity
events=rbind(events[[1]],events[[2]],events[[3]])
events_notopo=rbind(eventsNT[[1]],eventsNT[[2]],eventsNT[[3]])

a=density(data$CV2[data$Year>=2007 & data$Year<=2008],na.rm=T)
b=density(events$CV2,na.rm=T)
c=density(events_notopo$CV2,na.rm=T)

plot(a,col="black",xlim=c(0,4),ylim=range(0,1.5),
     xlab="ECL intensity",ylab="Frequency",cex.main=1.2,lwd=3,main="")
polygon(b,col=rgb(0,0,1,1/4),density=-1)
polygon(c,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("ERAI","Control","NoTopo"),
       col=c("black",rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   


####################
## ECL R2_events[[r]]$CV2s?
#######################

fixesE$Location2<-0
I<-which(fixesE[,7]>=149 & fixesE[,7]<=154 & fixesE[,8]<(-37) & fixesE[,8]>=-41)
fixesE$Location2[I]<-1
I<-which(fixesE[,7]>=(149+(37+fixesE[,8])/2) & fixesE[,7]<=(154+(37+fixesE[,8])/2) & fixesE[,8]<(-31) & fixesE[,8]>=-37)
fixesE$Location2[I]<-1
I<-which(fixesE[,7]>=152 & fixesE[,7]<=157 & fixesE[,8]<=(-24) & fixesE[,8]>=-31)
fixesE$Location2[I]<-1

for(n in 1:6)
{
  fixes[[n]]$Location2<-0
  I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
  fixes[[n]]$Location2[I]<-1
  I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
  fixes[[n]]$Location2[I]<-1
  I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
  fixes[[n]]$Location2[I]<-1
  
  fixesNT[[n]]$Location2<-0
  I<-which(fixesNT[[n]][,7]>=149 & fixesNT[[n]][,7]<=154 & fixesNT[[n]][,8]<(-37) & fixesNT[[n]][,8]>=-41)
  fixesNT[[n]]$Location2[I]<-1
  I<-which(fixesNT[[n]][,7]>=(149+(37+fixesNT[[n]][,8])/2) & fixesNT[[n]][,7]<=(154+(37+fixesNT[[n]][,8])/2) & fixesNT[[n]][,8]<(-31) & fixesNT[[n]][,8]>=-37)
  fixesNT[[n]]$Location2[I]<-1
  I<-which(fixesNT[[n]][,7]>=152 & fixesNT[[n]][,7]<=157 & fixesNT[[n]][,8]<=(-24) & fixesNT[[n]][,8]>=-31)
  fixesNT[[n]]$Location2[I]<-1
}

load("~/Documents/ECLs/WRFruns/0708/dayrain_0708_fix.RData")
ECLrainchange[t,]=(apply(ECLday2NT*dayrainNT[,2:4],2,mean)/apply(ECLday2*dayrain[,2:4],2,mean)-1)*100
thresh=seq(1,2,0.1)
ECLchange<-ECLrainchange<-nECLrainchange<-matrix(0,11,3)
for(t in 1:11)
{
  date=as.numeric(format.Date(dayrain[,1],"%Y%m%d"))
  ECLday<-ECLday2<-ECLdayNT<-ECLday2NT<-matrix(0,length(date),3)
  for(n in 1:3)
    for(i in 1:length(date))
    {
      I=which(fixes[[n]]$Date==date[i] & fixes[[n]]$Location==1 & fixes[[n]]$CV>=thresh[t])
      if(length(I)>0) ECLday2[i,n]<-1
      I=which(fixesNT[[n]]$Date==date[i] & fixesNT[[n]]$Location==1 & fixesNT[[n]]$CV>=thresh[t])
      if(length(I)>0) ECLday2NT[i,n]<-1
    }
  
  ECLchange[t,]=(apply(ECLday2NT,2,sum)/apply(ECLday2,2,sum)-1)*100
  ECLrainchange[t,]=(apply(ECLday2NT*dayrainNT[,2:4],2,mean)/apply(ECLday2*dayrain[,2:4],2,mean)-1)*100
  nECLrainchange[t,]=(apply((1-ECLday2NT)*dayrainNT[,2:4],2,mean)/apply((1-ECLday2)*dayrain[,2:4],2,mean)-1)*100
}

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECLmeanrain_",cat[c],"_vCV_vwrf_loc1_2day.pdf",sep=""),width=6,height=4)
plot(NA,xlim=c(1,2),ylim=c(-75,50),xlab="Intensity threshold",ylab="% change")
for(i in 1:3) lines(thresh,ECLrainchange[,i],lwd=3,col=i+1)
for(i in 1:3) lines(thresh,nECLrainchange[,i],lwd=2,col=i+1,lty=2)
abline(h=0,col="gray",lwd=2,lty=2)
legend("topright",legend=c("Mean rainfall on ECL days","Mean rainfall on non-ECL days","","R1","R2","R3"),
       col=c(1,1,NA,2:4),lty=c(1,2,1,1,1,1),lwd=3,bty="n",ncol=2,inset=c(-0.25,0),cex=0.8)
dev.off()

### Distribution of intensity
fixes2=rbind(fixes[[4]],fixes[[5]],fixes[[6]])
fixes_notopo=rbind(fixesNT[[4]],fixesNT[[5]],fixesNT[[6]])

a=density(fixesE$CV[fixesE$Date>=20070000 & fixesE$Date<=20089999 & fixesE$Location==1],na.rm=T)
b=density(fixes2$CV[fixes2$Location==1],na.rm=T)
c=density(fixes_notopo$CV[fixes_notopo$Location==1],na.rm=T)

plot(a,col="black",xlim=c(0,4),ylim=range(0,1.5),
     xlab="ECL intensity",ylab="Frequency",cex.main=1.2,lwd=3,main="")
polygon(b,col=rgb(0,0,1,1/4),density=-1)
polygon(c,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("ERAI","Control","NoTopo"),
       col=c("black",rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   

######## Somehow, the PDF?
###

Rthresh=c(1,2,5,10,25,Inf)
ECLrain=ECLday2*dayrain[,2:4]
ECLrainNT=ECLday2NT*dayrainNT[,2:4]

Rcount=array(0,c(5,4,3))
for(r in 1:3)
  for(t in 1:5)
  {
    I=which(ECLrain[,r]>=Rthresh[t] & ECLrain[,r]<Rthresh[t+1])
    Rcount[t,1,r]=length(I)
    I=which(ECLrainNT[,r]>=Rthresh[t] & ECLrainNT[,r]<Rthresh[t+1])
    Rcount[t,2,r]=length(I)
    I=which(dayrain[,r+1]>=Rthresh[t] & dayrain[,r+1]<Rthresh[t+1])
    Rcount[t,3,r]=length(I)
    I=which(dayrainNT[,r+1]>=Rthresh[t] & dayrainNT[,r+1]<Rthresh[t+1])
    Rcount[t,4,r]=length(I)
  }

Rcount2=apply(Rcount,c(1,2),mean)

plot(1,NA,xlim=c(1,5),ylim=c(0,40),xlab="Mean ESB rainfall",ylab="Frequency",axes=F)
axis(1,at=1:5,Rthresh[1:5])
axis(2,at=seq(0,40,5))
for(r in 1:3)
{
  lines(1:5,Rcount[,1,r],col=r+1,lwd=3,lty=1)
  lines(1:5,Rcount[,2,r],col=r+1,lwd=3,lty=2)
}




ERAIday<-matrix(0,length(date),2)
  for(i in 1:length(date))
  {
    I=which(fixesE$Date==date[i] & fixesE$Location==1)
    if(length(I)>0) ERAIday[i,1]=1
    I=which(fixesE$Date==date[i] & fixesE$Location2==1)
    if(length(I)>0) ERAIday[i,2]=1
  }

#AWAP dayrain? 
t=1
a=apply(ECLday2*dayrain[,2:4],2,sum)/apply(dayrain[,2:4],2,sum)
b=apply(ECLday2NT*dayrainNT[,2:4],2,sum)/apply(dayrainNT[,2:4],2,sum)



for(i in 2:4) print(length(which(dayrainNT[,i]>=20))/length(which(dayrain[,i]>=20)))


######## Change in ECL PDF

fixes2=rbind(fixes[[1]],fixes[[2]],fixes[[3]])
fixes_notopo=rbind(fixesNT[[1]],fixesNT[[2]],fixesNT[[3]])

pdfF<-pdfFnt<-data.frame(X=seq(0,5,length.out=512),R1=rep(0,512),R2=rep(0,512),R3=rep(0,512),All=rep(0,512))
pval<-matrix(0,4,2)
for(r in 1:3)
{ 
  a=density(fixes[[r]]$CV[fixes[[r]]$Location==1],na.rm=T,from=0,to=5)
  pdfF[,r+1]=a$y
  a=density(fixesNT[[r]]$CV[fixesNT[[r]]$Location==1],na.rm=T,from=0,to=5)
  pdfFnt[,r+1]=a$y
  
  a=ks.test(fixes[[r]]$CV[fixes[[r]]$Location==1],fixesNT[[r]]$CV[fixesNT[[r]]$Location==1])
  pval[r,1]=a$p.value
}
a=ks.test(fixes2$CV[fixes2$Location==1],fixes_notopo$CV[fixes_notopo$Location==1])
pval[4,1]=a$p.value

a=density(fixes2$CV[fixes2$Location==1],na.rm=T,from=0,to=5)
pdfF[,5]=a$y
a=density(fixes_notopo$CV[fixes_notopo$Location==1],na.rm=T,from=0,to=5)
pdfFnt[,5]=a$y
pdfF2=cbind(pdfF[,1],pdfFnt[,2:5]-pdfF[,2:5])


yy=c(-signif(max(abs(pdfF2[,2:5]))+0.1,digits=1),signif(max(abs(pdfF2[,2:5]))+0.1,digits=1))
pdf(file=paste("~/Documents/ECLs/WRFruns/0708/FixCV_PDFchange_",cat[c],"_vwrf.pdf",sep=""),width=6,height=4)
plot(pdfF2[,1],pdfF2[,5],lwd=4,col="black",type="l",
     xlim=c(0,5),ylim=yy,xlab="Intensity",ylab="Change in PDF")
for(i in 2:4) lines(pdfF2[,1],pdfF2[,i],col=i,lwd=2)
abline(h=0,col="gray",lwd=2,lty=2)
legend("topright",legend=c("All","R1","R2","R3"),col=1:4,lwd=c(4,2,2,2),bty="n")
dev.off()


for(r in 1:3) print(length(which(fixes[[r]]$Location2==1))/length(which(fixes[[r]]$Location==1)))
for(r in 1:3) print(length(which(fixesNT[[r]]$Location2==1))/length(which(fixesNT[[r]]$Location==1)))



######### Copy the intensity figure

c=3
n=1
dates=seq(20070101,20081231.75,0.25)
CV<-CVnt<-matrix(NaN,length(dates),3)
for(r in 1:3)
    {
      fixes[[r]]$Date2=fixes[[r]]$Date+(as.numeric(fixes[[r]]$Time)-1)/4
      fixesNT[[r]]$Date2=fixesNT[[r]]$Date+(as.numeric(fixesNT[[r]]$Time)-1)/4
      
      for(i in 1:length(dates))
      {
        I=which(fixes[[r]]$Date2==dates[i] & fixes[[r]]$Location==1)
        if(length(I)>0) CV[i,r]=mean(fixes[[r]]$CV[I])          
        I=which(fixesNT[[r]]$Date2==dates[i] & fixesNT[[r]]$Location==1)
        if(length(I)>0) CVnt[i,r]=mean(fixesNT[[r]]$CV[I])
      }
}

yy=range(CV,CVnt,na.rm=T)
pdf(file=paste("~/Documents/ECLs/WRFruns/0708/DailyCV_scatter_",cat[c],"_vwrf.pdf",sep=""),width=4,height=4)
plot(1:ceiling(yy[2]),1:ceiling(yy[2]),type="l",col="black",lwd=3,xlab="Default case",ylab="No topography")
for(i in 1:3) points(CV[,i],CVnt[,i],col=i+1,pch=4,lwd=2)
legend("topleft",legend=c("R1","R2","R3"),col=2:4,pch=4,pt.lwd=2,bty="n")
dev.off()

thrsh=seq(1,4,0.25)
CVln=matrix(NaN,length(thrsh),3)
for(r in 1:3)
  for(t in 1:length(thrsh)-1)
  {
    I=which(CV[,r]>=thrsh[t] & CV[,r]<=thrsh[t+1])
    if(length(I)>1) CVln[t,r]=mean(CVnt[I,r]-CV[I,r],na.rm=T)
  }

CVln2=cbind(thrsh,CVln)

yy=range(CV,CVnt,na.rm=T)
pdf(file=paste("~/Documents/ECLs/WRFruns/0708/DailyCV_line_",cat[c],"_vwrf.pdf",sep=""),width=4,height=4)
plot(1,NA,xlim=c(1,3),ylim=c(-1,0.5),lwd=3,xlab="Default case",ylab="Mean change in intensity")
abline(h=0,col="gray",lty=2)
for(i in 1:3) lines(CVln2[,1],CVln2[,i+1],col=i+1,pch=4,lwd=2)
legend("bottomleft",legend=c("R1","R2","R3"),col=2:4,pch=4,pt.lwd=2,bty="n",ncol=3)
dev.off()



#####
  date=as.numeric(format.Date(dayrain[,1],"%Y%m%d"))
  ECLday2<-ECLday2NT<-matrix(0,length(date),3)
  for(n in 1:3)
    for(i in 1:length(date))
    {
      I=which(fixes[[n]]$Date==date[i] & fixes[[n]]$Location==1)
      if(length(I)>0) ECLday2[i,n]<-ECLday2[i-1,n]<-1
      I=which(fixesNT[[n]]$Date==date[i] & fixesNT[[n]]$Location==1)
      if(length(I)>0) ECLday2NT[i,n]<-ECLday2NT[i-1,n]<-1
    }
  
  thresh=c(1,5,10,25,50,100,250,Inf)
  tabE<-matrix(0,7,6)
  rownames(tabE)=thresh[1:7]
  
  for(i in 1:7)
    for(j in 1:3)
    {
      I=which(dayrain_max[,j]>=thresh[i] & dayrain_max[,j]<thresh[i+1] & ECLday2[,j]==1)
      tabE[i,j]=length(I)
      I=which(dayrainNT_max[,j]>=thresh[i] & dayrainNT_max[,j]<thresh[i+1] & ECLday2NT[,j]==1)
      tabE[i,j+3]=length(I)
    }
  
  


pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECLmeanrain_",cat[c],"_vCV_vwrf_loc1_2day.pdf",sep=""),width=6,height=4)
plot(NA,xlim=c(1,2),ylim=c(-75,50),xlab="Intensity threshold",ylab="% change")
for(i in 1:3) lines(thresh,ECLrainchange[,i],lwd=3,col=i+1)
for(i in 1:3) lines(thresh,nECLrainchange[,i],lwd=2,col=i+1,lty=2)
abline(h=0,col="gray",lwd=2,lty=2)
legend("topright",legend=c("Mean rainfall on ECL days","Mean rainfall on non-ECL days","","R1","R2","R3"),
       col=c(1,1,NA,2:4),lty=c(1,2,1,1,1,1),lwd=3,bty="n",ncol=2,inset=c(-0.25,0),cex=0.8)
dev.off()


## Testing ECL days
rm(list=ls())
setwd('/home/nfs/z3478332/output/outputUM_wrf_2007_all/')
year=c(2007,2008)
domain=c("d01","d02")
type=c("_","_notopo_")
cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
days<-array(0,c(4,3,4))

for(c in 1:4)
  for(r in 1:3)
    {
      n=1
      for(d in 1:2)
        for(t in 1:2)
      
    {
      
      a=rbind(read.csv(paste("ECLfixes_",domain[d],"_2007_R",r,type[t],cat[c],".csv",sep="")),
              read.csv(paste("ECLfixes_",domain[d],"_2008_R",r,type[t],cat[c],".csv",sep="")))
      days[c,r,n]=length(unique(a$Date[a$Location==1]))
      n=n+1
    }
    }

