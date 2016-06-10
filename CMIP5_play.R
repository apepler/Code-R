lat=seq(-90,90,2.5)
lon=seq(0,357.5,2.5)
maskS<-maskN<-matrix(NaN,144,73)
reglats=seq(-25,-32.5,-2.5)
for(i in 1:4)
{
  I=which(lon==152.5)
  J=which(lat==reglats[i])
  maskN[I,J]=1
}
# 
# ##v2 - based on 2.5 degrees
# 
# reglats=c(rep(-32.5,4),rep(-35,5),rep(-37.5,3))
# reglons=c(seq(135,142.5,2.5),seq(137.5,147.5,2.5),seq(140,145,2.5))
# for(i in 1:length(reglats))
# {
#   I=which(lat==reglats[i])
#   J=which(lon==reglons[i])
#   maskS[J,I]=1
# }

## v3 - complex polygon
reglats=c(rep(-32.5,4),rep(-35,4),rep(-37.5,3))
reglons=c(seq(135,142.5,2.5),seq(137.5,145,2.5),seq(140,145,2.5))
for(i in 1:length(reglats))
{
  I=which(lat==reglats[i])
  J=which(lon==reglons[i])
  maskS[J,I]=1
}

# polygon(x=c(131.25,141.25,146.75,146.75,138.75,131.25),y=c(-31.25,-31.25,-36.25,-38.75,-38.75,-31.25),lwd=4,border="gray")
# J=which(lat>=-38.75 & lat<=-31.25)
# for(i in 1:length(J))
# {
#   l1=131.25-(10/7.5)*(31.25+lat[J[i]])
#   if(lat[J[i]])<(-36.25) l2=146.75 else l2=138.75-(8/5)*(31.25+lat[J[i]])
#   I=which(lon<=l2 & lon>=l1)
#   maskS[I,J[i]]=1
# }


CMIP_wind<-function(name,date1,date2,loc)
{ 
library("RNetCDF")

if(loc=="S")
  {
  a=21
  b=49
  } else {
  a=23
  b=47
  }

Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
Ufile=paste("/srv/ccrc/data02/z3394369/CMIP5_uw850_historical/ua_850hPa_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
fig1=paste("~/Documents/Data/CMIP5/Uwind_",name,".tiff",sep="")
fig2=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_",name,"_cool.tiff",sep="")
fig3=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_",name,"_warm.tiff",sep="")

f1=open.nc(Rfile)
time=var.get.nc(f1,"time")
ori=att.get.nc(f1,"time","units")
ori=unlist(strsplit(ori, split=" "))[3]

time2=as.Date(time,origin=ori)
time=data.frame(Year=as.numeric(format(time2,"%Y")),Month=as.numeric(format(time2,"%m")))
lat=var.get.nc(f1,"lat")
lon=var.get.nc(f1,"lon")
rain=var.get.nc(f1,"pr")
I=which(time[,2] %in% c(1,3,5,7,8,10,12))
rain[,,I]=rain[,,I]*31*60*60*24
I=which(time[,2] %in% c(4,6,9,11))
rain[,,I]=rain[,,I]*30*60*60*24
I=which(time[,2]==2)
rain[,,I]=rain[,,I]*28*60*60*24

f1=open.nc(Ufile)
u=var.get.nc(f1,"ua")
##What was the region we decided on for NCEP?
##Oh, woops, opposite hemisphere in order!
Uave=apply(u[61,a:(a+4),],2,mean)

##Need to have NCEP for comparison! 
# f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
# timeN=var.get.nc(f1,'time')
# timeN=as.Date((timeN/24-2),origin="0001-01-01")
# 
# timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
# uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(144,73,1,length(timeN[,1])),unpack=T)
# NCEP=apply(var.get.nc(f1,'uwnd',c(61,b,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)
# 
# ##So, for colocated years 1948-2005, want the mean wind for each month & plot
# Ucomp=matrix(0,12,3)
# for(i in 1:12)
# {
#   Ucomp[i,1]=i
#   I=which(timeN[,1]<=2005 & timeN[,2]==i)
#   Ucomp[i,2]=mean(NCEP[I])
#   I=which(time[,1]>=1948 & time[,1]<=2005 & time[,2]==i)
#   Ucomp[i,3]=mean(Uave[I])
#   }
# 
# mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# yl=c(floor(min(Ucomp[,2:3])),ceiling(max(Ucomp[,2:3])))
# tiff(file=fig1, height=500, width=800)
# plot(seq(1:12),Ucomp[,2],ty="l",lwd=2,col="blue",xlab="",ylab="",axes=F,ylim=yl)
# axis(1, at = seq(1:12), labels = c(mnames),cex.axis=1.5) 
# axis(2, at = seq(yl[1],yl[2]),cex.axis=1.5)
# title(main="Average Uwind (m/s)",cex.main=2)
# lines(seq(1:12),Ucomp[,3],ty="l",lwd=2,col="red")
# abline(h=0,col="gray",lty=2)
# legend("topleft",legend=c("NCEP",name),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
# dev.off()

#Finally, let's do the two correlation plots
years=seq(1950,2005)
years=cbind(seq(1950,2005),matrix(0,56,2))
coolR<-warmR<-array(0,c(144,73,56))
for(i in 1:length(years[,1]))
{
  I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
  years[i,2]=mean(Uave[I])
  coolR[,,i]=apply(rain[,,I],c(1,2),sum)
  I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
  years[i,3]=mean(Uave[I])
  warmR[,,i]=apply(rain[,,I],c(1,2),sum)
}

coolC<-warmC<-matrix(0,144,73)
for(i in 1:144)
  for(j in 1:73)
  {
    coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
    warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
  }

library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(14)
cm[7:8]="white"
coolC[coolC>=0.69]=0.69
coolC[coolC<=(-0.69)]=-0.69
warmC[warmC>=0.69]=0.69
warmC[warmC<=(-0.69)]=-0.69
coolC[abs(coolC)<=0.27]=NaN
warmC[abs(warmC)<=0.27]=NaN

tiff(file=fig2, height=500, width=800)
image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=800)
image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

return(list(coolC,warmC))
}

CMIP_data<-function(name,date1,date2,loc) ##Loc is N vs. South
{ 
  library("RNetCDF")
  if(loc=="S") a=21 else a=23
  
  Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  Ufile=paste("/srv/ccrc/data02/z3394369/CMIP5_uw850_historical/ua_850hPa_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  
  f1=open.nc(Rfile)
  time=var.get.nc(f1,"time")
  ori=att.get.nc(f1,"time","units")
  ori=unlist(strsplit(ori, split=" "))[3]
  
  time2=as.Date(time,origin=ori)
  time=data.frame(Year=as.numeric(format(time2,"%Y")),Month=as.numeric(format(time2,"%m")))
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  rain=var.get.nc(f1,"pr")
  I=which(time[,2] %in% c(1,3,5,7,8,10,12))
  rain[,,I]=rain[,,I]*31*60*60*24
  I=which(time[,2] %in% c(4,6,9,11))
  rain[,,I]=rain[,,I]*30*60*60*24
  I=which(time[,2]==2)
  rain[,,I]=rain[,,I]*28*60*60*24
  f1=open.nc(Ufile)
  u=var.get.nc(f1,"ua")
  ##What was the region we decided on for NCEP?
  ##Oh, woops, opposite hemisphere in order!
  Uave=apply(u[61,a:(a+4),],2,mean)
  U2=rep(0,12)
  
  for(i in 1:12)
  {
    I=which(time[,1]>=1948 & time[,1]<=2005 & time[,2]==i)
    U2[i]=mean(Uave[I])
  }
  
  #Finally, let's do the two correlation plots
  years=seq(1950,2005)
  years=cbind(seq(1950,2005),matrix(0,56,2))
  coolR<-warmR<-array(0,c(144,73,56))
  for(i in 1:length(years[,1]))
  {
    I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
    years[i,2]=mean(Uave[I])
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
    years[i,3]=mean(Uave[I])
    warmR[,,i]=apply(rain[,,I],c(1,2),sum)
   }
  
  coolC<-warmC<-matrix(0,144,73)
  for(i in 1:144)
    for(j in 1:73)
    {
      coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
      warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
    }

 return(list(U2,coolC,warmC))
}

CMIP_rain<-function(name,date1,date2)
{ 
  library("RNetCDF")
  
  Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  fig2=paste("~/Documents/Data/CMIP5/AveRain/Rain_",name,"_cool.tiff",sep="")
  fig3=paste("~/Documents/Data/CMIP5/AveRain/Rain_",name,"_warm.tiff",sep="")
  
  f1=open.nc(Rfile)
  time=var.get.nc(f1,"time")
  ori=att.get.nc(f1,"time","units")
  ori=unlist(strsplit(ori, split=" "))[3]
  
  time2=as.Date(time,origin=ori)
  time=data.frame(Year=as.numeric(format(time2,"%Y")),Month=as.numeric(format(time2,"%m")))
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  rain=var.get.nc(f1,"pr")
  
  I=which(time[,2] %in% c(1,3,5,7,8,10,12))
  rain[,,I]=rain[,,I]*31*60*60*24
  I=which(time[,2] %in% c(4,6,9,11))
  rain[,,I]=rain[,,I]*30*60*60*24
  I=which(time[,2]==2)
  rain[,,I]=rain[,,I]*28*60*60*24
  
  #Finally, let's do the averages
  years=seq(1950,2005)
  years=cbind(seq(1950,2005),matrix(0,56,2))
  coolR<-warmR<-array(0,c(144,73,56))
  for(i in 1:length(years[,1]))
  {
    I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
    warmR[,,i]=apply(rain[,,I],c(1,2),sum)
  }
  
  coolC=apply(coolR,c(1,2),mean,na.rm=T)
  warmC=apply(warmR,c(1,2),mean,na.rm=T)
#     
#   library(fields)
#   library("R.matlab")
#   readMat('~/Documents/Data/Useful.mat')->Useful
#   mask<-t(Useful$mask)
#   mask[is.na(mask)]=0
#   
#   source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
#   tiff(file=fig2, height=500, width=800)
#   par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
#   image(lon,lat,coolC,xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
#   image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
#   dev.off()
#   tiff(file=fig3, height=500, width=800)
#   par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
#   image(lon,lat,warmC,xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
#   image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
#   dev.off()
  
  return(list(coolC,warmC))
}

CMIP_v2<-function(name,date1,date2,loc)
{ 
  library("RNetCDF")
  
  if(loc=="S")
  {
    a=21
    b=49
  } else {
    a=23
    b=47
  }
  
yy=seq(1850,2012)
time=matrix(0,length(yy)*12,3)
n=1
for(i in 1:163)
  for(j in 1:12)
  {
    time[n,1]=yy[i]
    time[n,2]=j
    time[n,3]=yy[i]*100+j
    n=n+1
  }
  
  Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  Ufile=paste("/srv/ccrc/data02/z3394369/CMIP5_uw850_historical/ua_850hPa_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  fig1=paste("~/Documents/Data/CMIP5/Uwind_",name,".tiff",sep="")
  fig2=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_",name,"_cool.tiff",sep="")
  fig3=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_",name,"_warm.tiff",sep="")
  fig2a=paste("~/Documents/Data/CMIP5/Rain_corrUwind_",name,"_cool.tiff",sep="")
  fig3a=paste("~/Documents/Data/CMIP5/Rain_corrUwind_",name,"_warm.tiff",sep="")
  
  f1=open.nc(Rfile)
  time=time[(time[,3]>=date1 & time[,3]<=date2),1:2]
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  rain=var.get.nc(f1,"pr")
  I=which(time[,2] %in% c(1,3,5,7,8,10,12))
  rain[,,I]=rain[,,I]*31*60*60*24
  I=which(time[,2] %in% c(4,6,9,11))
  rain[,,I]=rain[,,I]*30*60*60*24
  I=which(time[,2]==2)
  rain[,,I]=rain[,,I]*28*60*60*24
  
  f1=open.nc(Ufile)
  u=var.get.nc(f1,"ua")
  ##What was the region we decided on for NCEP?
  ##Oh, woops, opposite hemisphere in order!
  Uave=apply(u[61,a:(a+4),],2,mean)
  U2=rep(0,12)  
  for(i in 1:12)
  {
    I=which(time[,1]>=1948 & time[,1]<=2005 & time[,2]==i)
    U2[i]=mean(Uave[I])
  }  

  #Finally, let's do the two correlation plots
  years=seq(1950,2005)
  years=cbind(seq(1950,2005),matrix(0,56,2))
  coolR<-warmR<-array(0,c(144,73,56))
  for(i in 1:length(years[,1]))
  {
    I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
    years[i,2]=mean(Uave[I])
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
    years[i,3]=mean(Uave[I])
    warmR[,,i]=apply(rain[,,I],c(1,2),sum)
  }
  
  coolC<-warmC<-matrix(0,144,73)
  for(i in 1:144)
    for(j in 1:73)
    {
      coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
      warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
    }
  
#   library(fields)
#   library("R.matlab")
#   readMat('~/Documents/Data/Useful.mat')->Useful
#   mask<-t(Useful$mask)
#   mask[is.na(mask)]=0
#   source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
#   source('~/Documents/R/color.palette.R')
#   pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
#   cm=pal(14)
#   cm[7:8]="white"
#   coolC[coolC>=0.69]=0.69
#   coolC[coolC<=(-0.69)]=-0.69
#   warmC[warmC>=0.69]=0.69
#   warmC[warmC<=(-0.69)]=-0.69
# tiff(file=fig2a, height=500, width=800)
# image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# tiff(file=fig3a, height=500, width=800)
# image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
# contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
# dev.off()
# 
# 
#   coolC[abs(coolC)<=0.27]=NaN
#   warmC[abs(warmC)<=0.27]=NaN
#   tiff(file=fig2, height=500, width=800)
#   image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   dev.off()
#   tiff(file=fig3, height=500, width=800)
#   image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   dev.off()
  
  return(list(U2,coolC,warmC))
}

CMIP_rain_v2<-function(name,date1,date2,loc)
{ 
  library("RNetCDF")
  
  if(loc=="S")
  {
    a=21
    b=49
  } else {
    a=23
    b=47
  }
  
  yy=seq(1850,2012)
  time=matrix(0,length(yy)*12,3)
  n=1
  for(i in 1:163)
    for(j in 1:12)
    {
      time[n,1]=yy[i]
      time[n,2]=j
      time[n,3]=yy[i]*100+j
      n=n+1
    }
  
  Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  Ufile=paste("/srv/ccrc/data02/z3394369/CMIP5_uw850_historical/ua_850hPa_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  fig1=paste("~/Documents/Data/CMIP5/Rain_",name,".tiff",sep="")
  fig2=paste("~/Documents/Data/CMIP5/Rain_",name,"_cool.tiff",sep="")
  fig3=paste("~/Documents/Data/CMIP5/Rain_",name,"_warm.tiff",sep="")
  
  f1=open.nc(Rfile)
  time=time[(time[,3]>=date1 & time[,3]<=date2),1:2]
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  rain=var.get.nc(f1,"pr")
  I=which(time[,2] %in% c(1,3,5,7,8,10,12))
  rain[,,I]=rain[,,I]*31*60*60*24
  I=which(time[,2] %in% c(4,6,9,11))
  rain[,,I]=rain[,,I]*30*60*60*24
  I=which(time[,2]==2)
  rain[,,I]=rain[,,I]*28*60*60*24
  

  #Finally, let's do the two correlation plots
  years=seq(1950,2005)
  years=cbind(seq(1950,2005),matrix(0,56,2))
  coolR<-annR<-array(0,c(144,73,56))
  warmR<-array(0,c(144,73,55))
  for(i in 1:length(years[,1]))
  {
    I=which(time[,1]==years[i,1])
    annR[,,i]=apply(rain[,,I],c(1,2),sum)
    I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    
    if(i<length(years[,1]))
    {
    I=which((time[,1]==years[i+1,1] & time[,2]<=4) | (time[,1]==years[i,1] & time[,2]>=11))
    warmR[,,i]=apply(rain[,,I],c(1,2),sum)
    }
  }
  
  coolC=apply(coolR,c(1,2),mean,na.rm=T)
  warmC=apply(warmR,c(1,2),mean,na.rm=T)
  annC=apply(annR,c(1,2),mean,na.rm=T)
  
  library(fields)
  library("R.matlab")
  readMat('~/Documents/Data/Useful.mat')->Useful
  mask<-t(Useful$mask)
  mask[is.na(mask)]=0
  
#   source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
#   tiff(file=fig1, height=500, width=800)
#   par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
#   image(lon,lat,annC,xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
#   image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
#   dev.off()
#   tiff(file=fig2, height=500, width=800)
#   par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
#   image(lon,lat,coolC,xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
#   image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
#   dev.off()
#   tiff(file=fig3, height=500, width=800)
#   par(plt = c(0.1,0.85,0.15,0.85),las = 1,cex.axis = 1)
#   image(lon,lat,warmC,xlab="",ylab="",breaks=c(0,10,25,50,100,200,300,400,600,800,1200,1600,2000),col=rich.colors(12),zlim=c(0,2000),xlim=c(100,180),ylim=c(-50,-10))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   par(new = "TRUE",plt = c(0.85,0.9,0.15,0.85),las = 1,cex.axis = 1)
#   image.plot(legend.only=TRUE,breaks=seq(0,12),lab.breaks=c("0","10","25","50","100","200","300","400","600","800","1200","1600","2000"),col=rich.colors(12),zlim=c(0,12))
#   dev.off()

  
  return(list(annC,coolC,warmC))
}

CMIP_wind_v2<-function(name,date1,date2,loc)
{ 
  library("RNetCDF")
  
  if(loc=="S")
  {
    a=21
    b=49
  } else {
    a=23
    b=47
  }
  
  yy=seq(1850,2012)
  time=matrix(0,length(yy)*12,3)
  n=1
  for(i in 1:163)
    for(j in 1:12)
    {
      time[n,1]=yy[i]
      time[n,2]=j
      time[n,3]=yy[i]*100+j
      n=n+1
    }
  
  Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  Ufile=paste("/srv/ccrc/data02/z3394369/CMIP5_uw850_historical/ua_850hPa_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  fig1=paste("~/Documents/Data/CMIP5/Uwind_",name,".tiff",sep="")
  fig2=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_",name,"_cool.tiff",sep="")
  fig3=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_",name,"_warm.tiff",sep="")
  fig2a=paste("~/Documents/Data/CMIP5/Rain_corrUwind_",name,"_cool.tiff",sep="")
  fig3a=paste("~/Documents/Data/CMIP5/Rain_corrUwind_",name,"_warm.tiff",sep="")
  
  f1=open.nc(Rfile)
  time=time[(time[,3]>=date1 & time[,3]<=date2),1:2]
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  
  f1=open.nc(Ufile)
  u=var.get.nc(f1,"ua")
  ##What was the region we decided on for NCEP?
  ##Oh, woops, opposite hemisphere in order!
  Uave=apply(u[61,a:(a+4),],2,mean)
  U2=rep(0,12)  
  for(i in 1:12)
  {
    I=which(time[,1]>=1948 & time[,1]<=2005 & time[,2]==i)
    U2[i]=mean(Uave[I])
  }  
  
  return(U2)
}

CMIP_corrs_v2<-function(name,date1,date2,loc)
{ 
  library("RNetCDF")
  out=rep(NaN,4)
  
  if(loc=="S")
  {
    a=21
    b=49
  } else {
    a=23
    b=47
  }
  
  yy=seq(1850,2012)
  time=matrix(0,length(yy)*12,3)
  n=1
  for(i in 1:163)
    for(j in 1:12)
    {
      time[n,1]=yy[i]
      time[n,2]=j
      time[n,3]=yy[i]*100+j
      n=n+1
    }
  
  Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  Ufile=paste("/srv/ccrc/data02/z3394369/CMIP5_uw850_historical/ua_850hPa_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")
  fig1=paste("~/Documents/Data/CMIP5/Rain_",name,".tiff",sep="")
  fig2=paste("~/Documents/Data/CMIP5/Rain_",name,"_cool.tiff",sep="")
  fig3=paste("~/Documents/Data/CMIP5/Rain_",name,"_warm.tiff",sep="")
  
  f1=open.nc(Rfile)
  time=time[(time[,3]>=date1 & time[,3]<=date2),1:2]
  latA=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  rain=var.get.nc(f1,"pr")
  I=which(time[,2] %in% c(1,3,5,7,8,10,12))
  rain[,,I]=rain[,,I]*31*60*60*24
  I=which(time[,2] %in% c(4,6,9,11))
  rain[,,I]=rain[,,I]*30*60*60*24
  I=which(time[,2]==2)
  rain[,,I]=rain[,,I]*28*60*60*24
  
  f1=open.nc(Ufile)
  u=var.get.nc(f1,"ua")
  ##What was the region we decided on for NCEP?
  ##Oh, woops, opposite hemisphere in order!
  Uave=apply(u[61,a:(a+4),],2,mean)
  U2=rep(0,12)  
  for(i in 1:12)
  {
    I=which(time[,1]>=1948 & time[,1]<=2005 & time[,2]==i)
    U2[i]=mean(Uave[I])
  }  
  
  #Finally, let's do the two correlation plots
  years<-coolM<-cbind(seq(1950,2005),matrix(0,56,2))
  coolR<-array(0,c(144,73,56))
  for(i in 1:length(years[,1]))
  {
    I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
    years[i,2]=mean(Uave[I])
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    coolM[i,1]=mean(coolR[,,i]*maskS,na.rm=T)
    coolM[i,2]=mean(coolR[,,i]*maskN,na.rm=T)
  } 
  
  coolC<-warmC<-matrix(0,144,73)
  for(i in 1:144) for(j in 1:73) coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")

  out[1]=mean(coolC*maskS,na.rm=T)
  out[2]=mean(coolC*maskN,na.rm=T)
  out[3]=cor(coolM[,1],years[,2],use="pairwise.complete.obs")
  out[4]=cor(coolM[,2],years[,2],use="pairwise.complete.obs")
  return(out)
}

str<-function(slp,lat,lon)
{
  ##Assumes slp is 3-d - lon, lat, time
  ##Calcutes subtropical ridge
  
  I=which(lat>=-45 & lat<=-10)
  lat=lat[I]
  J=which(lon>=145 & lon<=150)
  slp2=apply(slp[J,I,],c(2,3),mean) ##Profile of mean 145-150 pressure
  d=dim(slp)
  
  str=data.frame(STRI=rep(0,d[3]),STRP=rep(0,d[3]))
  for(i in 1:length(d3))
  {
    str[i,1]=max(slp2[,i])
    I=which(slp2[,i]==str[i,1])
    str[i,2]=max(lat[I])
  }
  
  return(str)
}


CMIP_mean_rain<-function(name,date1,date2,loc)
{ 
  library("RNetCDF")
  out=rep(NaN,2)
  
  if(loc=="S")
  {
    a=21
    b=49
  } else {
    a=23
    b=47
  }
  
  yy=seq(1850,2012)
  time=matrix(0,length(yy)*12,3)
  n=1
  for(i in 1:163)
    for(j in 1:12)
    {
      time[n,1]=yy[i]
      time[n,2]=j
      time[n,3]=yy[i]*100+j
      n=n+1
    }
  
  Rfile=paste("/srv/ccrc/data02/z3394369/CMIP5_precip_Historical/Regrid_2.5_cdo/pr_Amon_",name,"_historical_r1i1p1_",date1,"-",date2,"_2.5deg.nc",sep="")

  
  f1=open.nc(Rfile)
  time=time[(time[,3]>=date1 & time[,3]<=date2),1:2]
  latA=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  rain=var.get.nc(f1,"pr")
  I=which(time[,2] %in% c(1,3,5,7,8,10,12))
  rain[,,I]=rain[,,I]*31*60*60*24
  I=which(time[,2] %in% c(4,6,9,11))
  rain[,,I]=rain[,,I]*30*60*60*24
  I=which(time[,2]==2)
  rain[,,I]=rain[,,I]*28*60*60*24
  
  #Finally, let's do the two correlation plots
  coolM<-cbind(seq(1950,2005),matrix(0,56,2))
  for(i in 1:length(coolM[,1]))
  {
    I=which(time[,1]==coolM[i,1] & time[,2]>=5 & time[,2]<=10)
    coolR=apply(rain[,,I],c(1,2),sum)
    coolM[i,1]=mean(coolR*maskS,na.rm=T)
    coolM[i,2]=mean(coolR*maskN,na.rm=T)
  } 
  
  out[1]=mean(coolM[,1],na.rm=T)
  out[2]=mean(coolM[,2],na.rm=T)
  return(out)
}