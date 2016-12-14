rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->Useful
readMat('~/Documents/GDI/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

events<-fixes<-comp<-list()
comp_rain<-comp_wind<-array(0,c(101,101,3,5))
tlist=c("","_notopo","_BRAN","_notopo_","_2eac")
dir1=c(36,36,37,37,45)

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")
c=3 ## Which version do I want? 3 is the default
dom="d01"

##########
## Step 1: load all the data for my category of choice

events<-events_notopo<-fixes<-fixes_notopo<-list()

n=1
  for(r in 1:3)
  {
    events[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    events[[n]][events[[n]]==-Inf]=NaN
    
    events_notopo[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_notopo_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    events_notopo[[n]]$Year=floor(events_notopo[[n]]$Date1/10000)
    events_notopo[[n]]$Month=floor(events_notopo[[n]]$Date1/100)%%100
    events_notopo[[n]][events_notopo[[n]]==-Inf]=NaN
    
    fixes[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    fixes[[n]]$Location2<-0
    I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    fixes_notopo[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_notopo_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixes_notopo[[n]]$Year=floor(fixes_notopo[[n]]$Date/10000)
    fixes_notopo[[n]]$Month=floor(fixes_notopo[[n]]$Date/100)%%100
    fixes_notopo[[n]]$Location2<-0
    I<-which(fixes_notopo[[n]][,7]>=149 & fixes_notopo[[n]][,7]<=154 & fixes_notopo[[n]][,8]<(-37) & fixes_notopo[[n]][,8]>=-41)
    fixes_notopo[[n]]$Location2[I]<-1
    I<-which(fixes_notopo[[n]][,7]>=(149+(37+fixes_notopo[[n]][,8])/2) & fixes_notopo[[n]][,7]<=(154+(37+fixes_notopo[[n]][,8])/2) & fixes_notopo[[n]][,8]<(-31) & fixes_notopo[[n]][,8]>=-37)
    fixes_notopo[[n]]$Location2[I]<-1
    I<-which(fixes_notopo[[n]][,7]>=152 & fixes_notopo[[n]][,7]<=157 & fixes_notopo[[n]][,8]<=(-24) & fixes_notopo[[n]][,8]>=-31)
    fixes_notopo[[n]]$Location2[I]<-1
    fixes_notopo[[n]]$Date2=as.POSIXct(paste(as.character(fixes_notopo[[n]]$Date),substr(fixes_notopo[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    n=n+1     
  }

slp<-slp_notopo<-gv<-gv_notopo<-list()
rain<-rain_notopo<-wind<-wind_notopo<-list()

n=1
for(r in 1:3)
{
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLslp_0708_",cat[c],".nc",sep=""))
  slp[[n]]=var.get.nc(a,"ECL_slp")
  gv[[n]]=var.get.nc(a,"ECL_gv")
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLrain_0708_",cat[c],"_centred.nc",sep=""))
  rain[[n]]=var.get.nc(a,"ECLrain")
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLwind_d01_0708_",cat[c],".nc",sep=""))
  wind[[n]]=var.get.nc(a,"ECL_WS10")
  
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_notopo/out/ECLslp_0708_",cat[c],".nc",sep=""))
  slp_notopo[[n]]=var.get.nc(a,"ECL_slp")
  gv_notopo[[n]]=var.get.nc(a,"ECL_gv")
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_notopo/out/ECLrain_0708_",cat[c],"_centred.nc",sep=""))
  rain_notopo[[n]]=var.get.nc(a,"ECLrain")
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_notopo/out/ECLwind_d01_0708_",cat[c],".nc",sep=""))
  wind_notopo[[n]]=var.get.nc(a,"ECL_WS10")
  
  n=n+1     
}

######## Step 1 - average composite for all events

all_slp<-array(0,c(21,21,3,2))
dimnames(all_slp)[[3]]=c("R1","R2","R3")
dimnames(all_slp)[[4]]=c("Control","NoTopo")
all_gv<-all_slp

for(n in 1:3)
{
  I=which(fixes[[n]]$Location==1 & fixes[[n]]$CV>=2)
  all_slp[,,n,1]=apply(slp[[n]][,,I],c(1,2),mean)
  all_gv[,,n,1]=apply(gv[[n]][,,I],c(1,2),mean)
  
  I=which(fixes_notopo[[n]]$Location==1 & fixes_notopo[[n]]$CV>=2)
  all_slp[,,n,2]=apply(slp_notopo[[n]][,,I],c(1,2),mean)
  all_gv[,,n,2]=apply(gv_notopo[[n]][,,I],c(1,2),mean)
}

pdf(file=paste(figdir,"ECL_slpgv_all_d01_",cat[c],"_NoTopo_location_cv2.pdf",sep=""),width=7,height=4,pointsize=12)
layout(cbind(1,2))
par(mar=c(3,3,4,1))
contour(apply(all_slp[,,,1],c(1,2),mean),main="a) Control",axes=F,cex.main=1.5,lwd=2)
contour(apply(all_gv[,,,1],c(1,2),mean),levels=0:10,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(all_slp[,,,2],c(1,2),mean),main="b) NoTopo",axes=F,cex.main=1.5,lwd=2)
contour(apply(all_gv[,,,2],c(1,2),mean),levels=0:10,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()

### How to identify first point in region

for(n in 1:3)
{
  fixes[[n]]$LocFirst=0
  ids=unique(fixes[[n]]$ID)
  for(i in 1:length(ids))
  {
    I=which(fixes[[n]]$ID==ids[i] & fixes[[n]]$Location==1)
    fixes[[n]]$LocFirst[I[1]]=1
  }
  
  fixes_notopo[[n]]$LocFirst=0
  ids=unique(fixes_notopo[[n]]$ID)
  for(i in 1:length(ids))
  {
    I=which(fixes_notopo[[n]]$ID==ids[i] & fixes_notopo[[n]]$Location==1)
    fixes_notopo[[n]]$LocFirst[I[1]]=1
  }
}


########## Now, need to do matching - keeping stats for the strongest point of the event (in terms of cv within ECL domain)

match<-list()
match_slp<-array(0,c(21,21,3,3))
dimnames(match_slp)[[3]]=c("R1","R2","R3")
dimnames(match_slp)[[4]]=c("Unmatched events","Matched events","Corresponding event")
match_gv<-match_wind<-match_rain<-match_slp

for(n in 1:3)
{
  match[[n]]=rep(NaN,length(events[[n]]$ID))
  
  tmp_slp<-tmp_gv<-tmp_rain<-tmp_wind<-array(NaN,c(21,21,length(events[[n]]$ID),2))
  
  for(i in 1:length(events[[n]]$ID))
  {
    tmp=fixes[[n]][(fixes[[n]]$ID==events[[n]]$ID[i] & fixes[[n]]$Location==1),]
    rn=range(tmp$Date2)
    
    ###### Find the strongest point & identify the corresponding grids
    I=which(fixes[[n]]$ID==events[[n]]$ID[i] & fixes[[n]]$Location==1 & fixes[[n]]$CV==max(tmp$CV))
    if(length(I)==1)
      {
      tmp_slp[,,i,1]<-slp[[n]][,,I]
      tmp_gv[,,i,1]<-gv[[n]][,,I]
      tmp_rain[,,i,1]<-rain[[n]][,,I]
      tmp_wind[,,i,1]<-wind[[n]][,,I]
    } else {
      tmp_slp[,,i,1]<-apply(slp[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_gv[,,i,1]<-apply(rain[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_rain[,,i,1]<-apply(wind[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_wind[,,i,1]<-apply(wind[[n]][,,I],c(1,2),mean,na.rm=T)
      }
      
    I=which(fixes_notopo[[n]]$Date2<=rn[2]+(60*60*6) & fixes_notopo[[n]]$Date2>=rn[1]-(60*60*6) & fixes_notopo[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_notopo[[n]]$ID[I])
      match[[n]][i]=length(J) #All events that match
      
      K=which(fixes_notopo[[n]]$ID%in%J & fixes_notopo[[n]]$Location==1)
      tmp2=fixes_notopo[[n]][K,]
      L=which(fixes_notopo[[n]]$ID%in%J & fixes_notopo[[n]]$Location==1 & fixes_notopo[[n]]$CV==max(tmp2$CV))
      
      if(length(L)==1)
      {
        tmp_slp[,,i,2]<-slp_notopo[[n]][,,L]
        tmp_gv[,,i,2]<-gv_notopo[[n]][,,L]
        tmp_rain[,,i,2]<-rain_notopo[[n]][,,L]
        tmp_wind[,,i,2]<-wind_notopo[[n]][,,L]
      } else {
        tmp_slp[,,i,2]<-apply(slp_notopo[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_gv[,,i,2]<-apply(rain_notopo[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_rain[,,i,2]<-apply(wind_notopo[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_wind[,,i,2]<-apply(wind_notopo[[n]][,,L],c(1,2),mean,na.rm=T)
      }
      
    } else match[[n]][i]=0
  }
  
  I=which(match[[n]]>0)
  match_slp[,,n,1]=apply(tmp_slp[,,-I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,2]=apply(tmp_slp[,,I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,3]=apply(tmp_slp[,,I,2],c(1,2),mean,na.rm=T)
  match_gv[,,n,1]=apply(tmp_gv[,,-I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,2]=apply(tmp_gv[,,I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,3]=apply(tmp_gv[,,I,2],c(1,2),mean,na.rm=T)
  match_rain[,,n,1]=apply(tmp_rain[,,-I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,2]=apply(tmp_rain[,,I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,3]=apply(tmp_rain[,,I,2],c(1,2),mean,na.rm=T)
  match_wind[,,n,1]=apply(tmp_wind[,,-I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,2]=apply(tmp_wind[,,I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,3]=apply(tmp_wind[,,I,2],c(1,2),mean,na.rm=T)

}

pdf(file=paste(figdir,"ECL_slpgv_matched_d01_",cat[c],"_NoTopo_location.pdf",sep=""),width=10.5,height=4,pointsize=12)
layout(cbind(1,2,3))
par(mar=c(3,3,4,1))
contour(apply(match_slp[,,,1],c(1,2),mean),main="a) Unmatched Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,1],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,2],c(1,2),mean),main="b) Matched Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,3],c(1,2),mean),main="c) Corresponding Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,3],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()

### Version 2 - only look at events where CV is > 2 in the control

########## Now, need to do matching - keeping stats for the strongest point of the event (in terms of cv within ECL domain)

match<-list()
match_slp<-array(0,c(21,21,3,3))
dimnames(match_slp)[[3]]=c("R1","R2","R3")
dimnames(match_slp)[[4]]=c("Unmatched events","Matched events","Corresponding event")
match_gv<-match_wind<-match_rain<-match_slp
thresh=2

for(n in 1:3)
{
  events2=events[[n]][events[[n]]$CV2>=thresh,]
  
  match[[n]]=rep(NaN,length(events[[n]]$ID))
  
  tmp_slp<-tmp_gv<-tmp_rain<-tmp_wind<-array(NaN,c(21,21,length(events[[n]]$ID),2))
  
  for(i in 1:length(events2$ID))
  {
    tmp=fixes[[n]][(fixes[[n]]$ID==events2$ID[i] & fixes[[n]]$Location==1),]
    rn=range(tmp$Date2)
    
    ###### Find the strongest point & identify the corresponding grids
    I=which(fixes[[n]]$ID==events2$ID[i] & fixes[[n]]$Location==1 & fixes[[n]]$CV==max(tmp$CV))
    if(length(I)==1)
    {
      tmp_slp[,,i,1]<-slp[[n]][,,I]
      tmp_gv[,,i,1]<-gv[[n]][,,I]
      tmp_rain[,,i,1]<-rain[[n]][,,I]
      tmp_wind[,,i,1]<-wind[[n]][,,I]
    } else {
      tmp_slp[,,i,1]<-apply(slp[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_gv[,,i,1]<-apply(rain[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_rain[,,i,1]<-apply(wind[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_wind[,,i,1]<-apply(wind[[n]][,,I],c(1,2),mean,na.rm=T)
    }
    
    I=which(fixes_notopo[[n]]$Date2<=rn[2]+(60*60*6) & fixes_notopo[[n]]$Date2>=rn[1]-(60*60*6) & fixes_notopo[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_notopo[[n]]$ID[I])
      match[[n]][i]=length(J) #All events that match
      
      K=which(fixes_notopo[[n]]$ID%in%J & fixes_notopo[[n]]$Location==1)
      tmp2=fixes_notopo[[n]][K,]
      L=which(fixes_notopo[[n]]$ID%in%J & fixes_notopo[[n]]$Location==1 & fixes_notopo[[n]]$CV==max(tmp2$CV))
      
      if(length(L)==1)
      {
        tmp_slp[,,i,2]<-slp_notopo[[n]][,,L]
        tmp_gv[,,i,2]<-gv_notopo[[n]][,,L]
        tmp_rain[,,i,2]<-rain_notopo[[n]][,,L]
        tmp_wind[,,i,2]<-wind_notopo[[n]][,,L]
      } else {
        tmp_slp[,,i,2]<-apply(slp_notopo[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_gv[,,i,2]<-apply(rain_notopo[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_rain[,,i,2]<-apply(wind_notopo[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_wind[,,i,2]<-apply(wind_notopo[[n]][,,L],c(1,2),mean,na.rm=T)
      }
      
    } else match[[n]][i]=0
  }
  
  I=which(match[[n]]>0)
  match_slp[,,n,1]=apply(tmp_slp[,,-I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,2]=apply(tmp_slp[,,I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,3]=apply(tmp_slp[,,I,2],c(1,2),mean,na.rm=T)
  match_gv[,,n,1]=apply(tmp_gv[,,-I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,2]=apply(tmp_gv[,,I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,3]=apply(tmp_gv[,,I,2],c(1,2),mean,na.rm=T)
  match_rain[,,n,1]=apply(tmp_rain[,,-I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,2]=apply(tmp_rain[,,I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,3]=apply(tmp_rain[,,I,2],c(1,2),mean,na.rm=T)
  match_wind[,,n,1]=apply(tmp_wind[,,-I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,2]=apply(tmp_wind[,,I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,3]=apply(tmp_wind[,,I,2],c(1,2),mean,na.rm=T)
  
}

pdf(file=paste(figdir,"ECL_slpgv_matched_d01_",cat[c],"_NoTopo_location_cv",thresh,".pdf",sep=""),width=7,height=4,pointsize=12)
layout(cbind(1,2))
par(mar=c(3,3,4,1))
contour(apply(match_slp[,,,2],c(1,2),mean),main="a) Control Events",axes=F,cex.main=1.5,lwd=2)
contour(apply(match_gv[,,,2],c(1,2),mean),levels=0:12,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,3],c(1,2),mean),main="b) Corresponding NoTopo Events",axes=F,cex.main=1.5,lwd=2)
contour(apply(match_gv[,,,3],c(1,2),mean),levels=0:12,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()
