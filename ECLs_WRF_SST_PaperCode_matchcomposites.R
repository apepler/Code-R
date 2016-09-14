rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing')
figdir="~/Documents/ECLs/WRFruns/0708/EACpaper/"

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
tlist=c("","_notopo","_BRAN","_BRAN_noeac","_BRAN_2eac")
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

eventsBRAN<-events_noeac<-events_2eac<-fixesBRAN<-fixes_noeac<-fixes_2eac<-list()

n=1
  for(r in 1:3)
  {
    eventsBRAN[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    eventsBRAN[[n]]$Year=floor(eventsBRAN[[n]]$Date1/10000)
    eventsBRAN[[n]]$Month=floor(eventsBRAN[[n]]$Date1/100)%%100
    eventsBRAN[[n]][eventsBRAN[[n]]==-Inf]=NaN
    
    events_noeac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    events_noeac[[n]]$Year=floor(events_noeac[[n]]$Date1/10000)
    events_noeac[[n]]$Month=floor(events_noeac[[n]]$Date1/100)%%100
    events_noeac[[n]][events_noeac[[n]]==-Inf]=NaN
    
    fixesBRAN[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixesBRAN[[n]]$Year=floor(fixesBRAN[[n]]$Date/10000)
    fixesBRAN[[n]]$Month=floor(fixesBRAN[[n]]$Date/100)%%100
    fixesBRAN[[n]]$Location2<-0
    I<-which(fixesBRAN[[n]][,7]>=149 & fixesBRAN[[n]][,7]<=154 & fixesBRAN[[n]][,8]<(-37) & fixesBRAN[[n]][,8]>=-41)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=(149+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,7]<=(154+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,8]<(-31) & fixesBRAN[[n]][,8]>=-37)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=152 & fixesBRAN[[n]][,7]<=157 & fixesBRAN[[n]][,8]<=(-24) & fixesBRAN[[n]][,8]>=-31)
    fixesBRAN[[n]]$Location2[I]<-1
    fixesBRAN[[n]]$Date2=as.POSIXct(paste(as.character(fixesBRAN[[n]]$Date),substr(fixesBRAN[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    fixes_noeac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_noeac_",cat[c],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixes_noeac[[n]]$Year=floor(fixes_noeac[[n]]$Date/10000)
    fixes_noeac[[n]]$Month=floor(fixes_noeac[[n]]$Date/100)%%100
    fixes_noeac[[n]]$Location2<-0
    I<-which(fixes_noeac[[n]][,7]>=149 & fixes_noeac[[n]][,7]<=154 & fixes_noeac[[n]][,8]<(-37) & fixes_noeac[[n]][,8]>=-41)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=(149+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,7]<=(154+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,8]<(-31) & fixes_noeac[[n]][,8]>=-37)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=152 & fixes_noeac[[n]][,7]<=157 & fixes_noeac[[n]][,8]<=(-24) & fixes_noeac[[n]][,8]>=-31)
    fixes_noeac[[n]]$Location2[I]<-1
    fixes_noeac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_noeac[[n]]$Date),substr(fixes_noeac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    fixes_2eac[[n]]=read.csv(paste("ECLfixes_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixes_2eac[[n]]$Year=floor(fixes_2eac[[n]]$Date/10000)
    fixes_2eac[[n]]$Month=floor(fixes_2eac[[n]]$Date/100)%%100
    fixes_2eac[[n]]$Location2<-0
    I<-which(fixes_2eac[[n]][,7]>=149 & fixes_2eac[[n]][,7]<=154 & fixes_2eac[[n]][,8]<(-37) & fixes_2eac[[n]][,8]>=-41)
    fixes_2eac[[n]]$Location2[I]<-1
    I<-which(fixes_2eac[[n]][,7]>=(149+(37+fixes_2eac[[n]][,8])/2) & fixes_2eac[[n]][,7]<=(154+(37+fixes_2eac[[n]][,8])/2) & fixes_2eac[[n]][,8]<(-31) & fixes_2eac[[n]][,8]>=-37)
    fixes_2eac[[n]]$Location2[I]<-1
    I<-which(fixes_2eac[[n]][,7]>=152 & fixes_2eac[[n]][,7]<=157 & fixes_2eac[[n]][,8]<=(-24) & fixes_2eac[[n]][,8]>=-31)
    fixes_2eac[[n]]$Location2[I]<-1
    fixes_2eac[[n]]$Date2=as.POSIXct(paste(as.character(fixes_2eac[[n]]$Date),substr(fixes_2eac[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    
    
    events_2eac[[n]]=read.csv(paste("ECLevents_",dom,"_0708_R",r,"_BRAN_2eac_",cat[c],"_v2_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    events_2eac[[n]]$Year=floor(events_2eac[[n]]$Date1/10000)
    events_2eac[[n]]$Month=floor(events_2eac[[n]]$Date1/100)%%100  
    events_2eac[[n]][events_2eac[[n]]==-Inf]=NaN
    n=n+1     
  }

slpBRAN<-slp_noeac<-slp_2eac<-gvBRAN<-gv_noeac<-gv_2eac<-list()
rainBRAN<-rain_noeac<-rain_2eac<-windBRAN<-wind_noeac<-wind_2eac<-list()

n=1
for(r in 1:3)
{
  a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN/out/ECLslp_0708_",cat[c],".nc",sep=""))
  slpBRAN[[n]]=var.get.nc(a,"ECL_slp")
  gvBRAN[[n]]=var.get.nc(a,"ECL_gv")
  a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN/out/ECLrain_0708_",cat[c],"_centred.nc",sep=""))
  rainBRAN[[n]]=var.get.nc(a,"ECLrain")
  a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN/out/ECLwind_d01_0708_",cat[c],".nc",sep=""))
  windBRAN[[n]]=var.get.nc(a,"ECL_WS10")
  
  a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN_noeac/out/ECLslp_0708_",cat[c],".nc",sep=""))
  slp_noeac[[n]]=var.get.nc(a,"ECL_slp")
  gv_noeac[[n]]=var.get.nc(a,"ECL_gv")
  a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN_noeac/out/ECLrain_0708_",cat[c],"_centred.nc",sep=""))
  rain_noeac[[n]]=var.get.nc(a,"ECLrain")
  a=open.nc(paste("/srv/ccrc/data37/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN_noeac/out/ECLwind_d01_0708_",cat[c],".nc",sep=""))
  wind_noeac[[n]]=var.get.nc(a,"ECL_WS10")
  
  a=open.nc(paste("/srv/ccrc/data45/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN_2eac/out/Analysis/ECLslp_0708_",cat[c],"_v2.nc",sep=""))
  slp_2eac[[n]]=var.get.nc(a,"ECL_slp")
  gv_2eac[[n]]=var.get.nc(a,"ECL_gv")
  a=open.nc(paste("/srv/ccrc/data45/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN_2eac/out/Analysis/ECLrain_0708_",cat[c],"_v2_centred.nc",sep=""))
  rain_2eac[[n]]=var.get.nc(a,"ECLrain")
  a=open.nc(paste("/srv/ccrc/data45/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_BRAN_2eac/out/Analysis/ECLwind_d01_0708_",cat[c],"_v2.nc",sep=""))
  wind_2eac[[n]]=var.get.nc(a,"ECL_WS10")
  
  n=n+1     
}

########## Now, need to do matching - keeping stats for the strongest point of the event (in terms of cv within ECL domain)

match<-list()
match_slp<-array(0,c(21,21,3,2,3))
dimnames(match_slp)[[3]]=c("R1","R2","R3")
dimnames(match_slp)[[4]]=c("NoEAC","2EAC")
dimnames(match_slp)[[5]]=c("Unmatched events","Matched events","Corresponding event")
match_gv<-match_wind<-match_rain<-match_slp

for(n in 1:3)
{
  match[[n]]=array(NaN,c(length(eventsBRAN[[n]]$ID),2))
  dimnames(match[[n]])[[2]]=c("NoEac","2EAC")
  
  tmp_slp<-tmp_gv<-tmp_rain<-tmp_wind<-array(NaN,c(21,21,length(eventsBRAN[[n]]$ID),3))
  
  for(i in 1:length(eventsBRAN[[n]]$ID))
  {
    tmp=fixesBRAN[[n]][(fixesBRAN[[n]]$ID==eventsBRAN[[n]]$ID[i] & fixesBRAN[[n]]$Location==1),]
    rn=range(tmp$Date2)
    
    ###### Find the strongest point & identify the corresponding grids
    I=which(fixesBRAN[[n]]$ID==eventsBRAN[[n]]$ID[i] & fixesBRAN[[n]]$Location==1 & fixesBRAN[[n]]$CV==max(tmp$CV))
    if(length(I)==1)
      {
      tmp_slp[,,i,1]<-slpBRAN[[n]][,,I]
      tmp_gv[,,i,1]<-gvBRAN[[n]][,,I]
      tmp_rain[,,i,1]<-rainBRAN[[n]][,,I]
      tmp_wind[,,i,1]<-windBRAN[[n]][,,I]
    } else {
      tmp_slp[,,i,1]<-apply(slpBRAN[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_gv[,,i,1]<-apply(rainBRAN[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_rain[,,i,1]<-apply(windBRAN[[n]][,,I],c(1,2),mean,na.rm=T)
      tmp_wind[,,i,1]<-apply(windBRAN[[n]][,,I],c(1,2),mean,na.rm=T)
      }
      
    I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_noeac[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_noeac[[n]]$ID[I])
      match[[n]][i,1]=length(J) #All events that match
      
      K=which(fixes_noeac[[n]]$ID%in%J & fixes_noeac[[n]]$Location==1)
      tmp2=fixes_noeac[[n]][K,]
      L=which(fixes_noeac[[n]]$ID%in%J & fixes_noeac[[n]]$Location==1 & fixes_noeac[[n]]$CV==max(tmp2$CV))
      
      if(length(L)==1)
      {
        tmp_slp[,,i,2]<-slp_noeac[[n]][,,L]
        tmp_gv[,,i,2]<-gv_noeac[[n]][,,L]
        tmp_rain[,,i,2]<-rain_noeac[[n]][,,L]
        tmp_wind[,,i,2]<-wind_noeac[[n]][,,L]
      } else {
        tmp_slp[,,i,2]<-apply(slp_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_gv[,,i,2]<-apply(rain_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_rain[,,i,2]<-apply(wind_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_wind[,,i,2]<-apply(wind_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
      }
      
    } else match[[n]][i,1]=0
    
    
    I=which(fixes_2eac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_2eac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_2eac[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_2eac[[n]]$ID[I])
      match[[n]][i,2]=length(J) #All events that match
      
      K=which(fixes_2eac[[n]]$ID%in%J & fixes_2eac[[n]]$Location==1)
      tmp2=fixes_2eac[[n]][K,]
      L=which(fixes_2eac[[n]]$ID%in%J & fixes_2eac[[n]]$Location==1 & fixes_2eac[[n]]$CV==max(tmp2$CV))
      
      if(length(L)==1)
      {
        tmp_slp[,,i,3]<-slp_2eac[[n]][,,L]
        tmp_gv[,,i,3]<-gv_2eac[[n]][,,L]
        tmp_rain[,,i,3]<-rain_2eac[[n]][,,L]
        tmp_wind[,,i,3]<-wind_2eac[[n]][,,L]
      } else {
        tmp_slp[,,i,3]<-apply(slp_2eac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_gv[,,i,3]<-apply(rain_2eac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_rain[,,i,3]<-apply(wind_2eac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_wind[,,i,3]<-apply(wind_2eac[[n]][,,L],c(1,2),mean,na.rm=T)
      }
    } else match[[n]][i,2]=0
    
  }
  
  I=which(match[[n]][,1]>0)
  match_slp[,,n,1,1]=apply(tmp_slp[,,-I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,1,2]=apply(tmp_slp[,,I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,1,3]=apply(tmp_slp[,,I,2],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,1]=apply(tmp_gv[,,-I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,2]=apply(tmp_gv[,,I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,3]=apply(tmp_gv[,,I,2],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,1]=apply(tmp_rain[,,-I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,2]=apply(tmp_rain[,,I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,3]=apply(tmp_rain[,,I,2],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,1]=apply(tmp_wind[,,-I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,2]=apply(tmp_wind[,,I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,3]=apply(tmp_wind[,,I,2],c(1,2),mean,na.rm=T)
  
  I=which(match[[n]][,2]>0)
  match_slp[,,n,2,1]=apply(tmp_slp[,,-I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,2,2]=apply(tmp_slp[,,I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,2,3]=apply(tmp_slp[,,I,3],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,1]=apply(tmp_gv[,,-I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,2]=apply(tmp_gv[,,I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,3]=apply(tmp_gv[,,I,3],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,1]=apply(tmp_rain[,,-I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,2]=apply(tmp_rain[,,I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,3]=apply(tmp_rain[,,I,3],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,1]=apply(tmp_wind[,,-I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,2]=apply(tmp_wind[,,I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,3]=apply(tmp_wind[,,I,3],c(1,2),mean,na.rm=T)

}


contour(apply(match_gv[,,,1,1],c(1,2),mean))
contour(apply(match_gv[,,,1,2],c(1,2),mean))
contour(apply(match_gv[,,,1,3],c(1,2),mean))

pdf(file=paste(figdir,"ECL_slpgv_matched_d01_",cat[c],"_NoEAC_location.pdf",sep=""),width=10.5,height=4,pointsize=12)
layout(cbind(1,2,3))
par(mar=c(3,3,4,1))
contour(apply(match_slp[,,,1,1],c(1,2),mean),main="a) Unmatched Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,1,1],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,1,2],c(1,2),mean),main="b) Matched Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,1,2],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,1,3],c(1,2),mean),main="c) Corresponding Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,1,3],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()


######### Version 2 - keep the match for the initial point where >1 & in region

match<-list()
match_slp<-array(0,c(21,21,3,2,3))
dimnames(match_slp)[[3]]=c("R1","R2","R3")
dimnames(match_slp)[[4]]=c("NoEAC","2EAC")
dimnames(match_slp)[[5]]=c("Unmatched events","Matched events","Corresponding event")
match_gv<-match_wind<-match_rain<-match_slp

for(n in 1:3)
{
  match[[n]]=array(NaN,c(length(eventsBRAN[[n]]$ID),2))
  dimnames(match[[n]])[[2]]=c("NoEac","2EAC")
  
  tmp_slp<-tmp_gv<-tmp_rain<-tmp_wind<-array(NaN,c(21,21,length(eventsBRAN[[n]]$ID),3))
  
  for(i in 1:length(eventsBRAN[[n]]$ID))
  {
    ###### Find the strongest point & identify the corresponding grids
    I=which(fixesBRAN[[n]]$ID==eventsBRAN[[n]]$ID[i] & fixesBRAN[[n]]$Location==1)
    tmp=fixesBRAN[[n]][I,]
    rn=range(tmp$Date2)
    tmp_slp[,,i,1]<-slpBRAN[[n]][,,I[1]]
      tmp_gv[,,i,1]<-gvBRAN[[n]][,,I[1]]
      tmp_rain[,,i,1]<-rainBRAN[[n]][,,I[1]]
      tmp_wind[,,i,1]<-windBRAN[[n]][,,I[1]]
    
    I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_noeac[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_noeac[[n]]$ID[I])
      match[[n]][i,1]=length(J) #All events that match
      
      K=which(fixes_noeac[[n]]$ID%in%J & fixes_noeac[[n]]$Location==1)
      tmp2=fixes_noeac[[n]][K,]
      L=which(fixes_noeac[[n]]$ID%in%J & fixes_noeac[[n]]$Location==1 & fixes_noeac[[n]]$Date2==min(tmp2$Date2))
      
      if(length(L)==1)
      {
        tmp_slp[,,i,2]<-slp_noeac[[n]][,,L]
        tmp_gv[,,i,2]<-gv_noeac[[n]][,,L]
        tmp_rain[,,i,2]<-rain_noeac[[n]][,,L]
        tmp_wind[,,i,2]<-wind_noeac[[n]][,,L]
      } else {
        tmp_slp[,,i,2]<-apply(slp_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_gv[,,i,2]<-apply(rain_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_rain[,,i,2]<-apply(wind_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_wind[,,i,2]<-apply(wind_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
      }
      
    } else match[[n]][i,1]=0
    
    
    I=which(fixes_2eac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_2eac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_2eac[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_2eac[[n]]$ID[I])
      match[[n]][i,2]=length(J) #All events that match
      
      K=which(fixes_2eac[[n]]$ID%in%J & fixes_2eac[[n]]$Location==1)
      tmp2=fixes_2eac[[n]][K,]
      L=which(fixes_2eac[[n]]$ID%in%J & fixes_2eac[[n]]$Location==1 & fixes_2eac[[n]]$Date2==min(tmp2$Date2))
      
      if(length(L)==1)
      {
        tmp_slp[,,i,3]<-slp_2eac[[n]][,,L]
        tmp_gv[,,i,3]<-gv_2eac[[n]][,,L]
        tmp_rain[,,i,3]<-rain_2eac[[n]][,,L]
        tmp_wind[,,i,3]<-wind_2eac[[n]][,,L]
      } else {
        tmp_slp[,,i,3]<-apply(slp_2eac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_gv[,,i,3]<-apply(rain_2eac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_rain[,,i,3]<-apply(wind_2eac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_wind[,,i,3]<-apply(wind_2eac[[n]][,,L],c(1,2),mean,na.rm=T)
      }
    } else match[[n]][i,2]=0
    
  }
  
  I=which(match[[n]][,1]>0)
  match_slp[,,n,1,1]=apply(tmp_slp[,,-I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,1,2]=apply(tmp_slp[,,I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,1,3]=apply(tmp_slp[,,I,2],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,1]=apply(tmp_gv[,,-I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,2]=apply(tmp_gv[,,I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,3]=apply(tmp_gv[,,I,2],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,1]=apply(tmp_rain[,,-I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,2]=apply(tmp_rain[,,I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,3]=apply(tmp_rain[,,I,2],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,1]=apply(tmp_wind[,,-I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,2]=apply(tmp_wind[,,I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,3]=apply(tmp_wind[,,I,2],c(1,2),mean,na.rm=T)
  
  I=which(match[[n]][,2]>0)
  match_slp[,,n,2,1]=apply(tmp_slp[,,-I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,2,2]=apply(tmp_slp[,,I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,2,3]=apply(tmp_slp[,,I,3],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,1]=apply(tmp_gv[,,-I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,2]=apply(tmp_gv[,,I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,3]=apply(tmp_gv[,,I,3],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,1]=apply(tmp_rain[,,-I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,2]=apply(tmp_rain[,,I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,3]=apply(tmp_rain[,,I,3],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,1]=apply(tmp_wind[,,-I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,2]=apply(tmp_wind[,,I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,3]=apply(tmp_wind[,,I,3],c(1,2),mean,na.rm=T)
  
}


contour(apply(match_gv[,,,1,1],c(1,2),mean))
contour(apply(match_gv[,,,1,2],c(1,2),mean))
contour(apply(match_gv[,,,1,3],c(1,2),mean))

pdf(file=paste(figdir,"ECL_slpgv_matched_d01_",cat[c],"_2EAC_location_initial.pdf",sep=""),width=10.5,height=4,pointsize=12)
layout(cbind(1,2,3))
par(mar=c(3,3,4,1))
contour(apply(match_slp[,,,2,1],c(1,2),mean),main="a) Unmatched Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2,1],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,2,2],c(1,2),mean),main="b) Matched Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2,2],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,2,3],c(1,2),mean),main="c) Corresponding Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2,3],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()


###Figure6
pdf(file=paste(figdir,"Figure6.pdf",sep=""),width=7,height=4,pointsize=12)
layout(cbind(1,2))
par(mar=c(3,3,4,1))
contour(apply(match_slp[,,,1,1],c(1,2),mean),main="",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,1,1],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,1,2],c(1,2),mean),main="",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,1,2],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()

################ Part 3 - Taking 2EAC as standard, and matching control/NoEAC


match<-list()
match_slp<-array(0,c(21,21,3,2,3))
dimnames(match_slp)[[3]]=c("R1","R2","R3")
dimnames(match_slp)[[4]]=c("NoEAC","BRAN")
dimnames(match_slp)[[5]]=c("Unmatched events","Matched events","Corresponding event")
match_gv<-match_wind<-match_rain<-match_slp

for(n in 1:3)
{
  match[[n]]=array(NaN,c(length(events_2eac[[n]]$ID),2))
  dimnames(match[[n]])[[2]]=c("NoEac","2EAC")
  
  tmp_slp<-tmp_gv<-tmp_rain<-tmp_wind<-array(NaN,c(21,21,length(events_2eac[[n]]$ID),3))
  
  for(i in 1:length(events_2eac[[n]]$ID))
  {
    ###### Find the strongest point & identify the corresponding grids
    I=which(fixes_2eac[[n]]$ID==events_2eac[[n]]$ID[i] & fixes_2eac[[n]]$Location==1)
    tmp=fixes_2eac[[n]][I,]
    rn=range(tmp$Date2)
    tmp_slp[,,i,1]<-slp_2eac[[n]][,,I[1]]
    tmp_gv[,,i,1]<-gv_2eac[[n]][,,I[1]]
    tmp_rain[,,i,1]<-rain_2eac[[n]][,,I[1]]
    tmp_wind[,,i,1]<-wind_2eac[[n]][,,I[1]]
    
    I=which(fixes_noeac[[n]]$Date2<=rn[2]+(60*60*6) & fixes_noeac[[n]]$Date2>=rn[1]-(60*60*6) & fixes_noeac[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes_noeac[[n]]$ID[I])
      match[[n]][i,1]=length(J) #All events that match
      
      K=which(fixes_noeac[[n]]$ID%in%J & fixes_noeac[[n]]$Location==1)
      tmp2=fixes_noeac[[n]][K,]
      L=which(fixes_noeac[[n]]$ID%in%J & fixes_noeac[[n]]$Location==1 & fixes_noeac[[n]]$Date2==min(tmp2$Date2))
      
      if(length(L)==1)
      {
        tmp_slp[,,i,2]<-slp_noeac[[n]][,,L]
        tmp_gv[,,i,2]<-gv_noeac[[n]][,,L]
        tmp_rain[,,i,2]<-rain_noeac[[n]][,,L]
        tmp_wind[,,i,2]<-wind_noeac[[n]][,,L]
      } else {
        tmp_slp[,,i,2]<-apply(slp_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_gv[,,i,2]<-apply(rain_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_rain[,,i,2]<-apply(wind_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_wind[,,i,2]<-apply(wind_noeac[[n]][,,L],c(1,2),mean,na.rm=T)
      }
      
    } else match[[n]][i,1]=0
    
    
    I=which(fixesBRAN[[n]]$Date2<=rn[2]+(60*60*6) & fixesBRAN[[n]]$Date2>=rn[1]-(60*60*6) & fixesBRAN[[n]]$Location==1)
    if(length(I)>0)
    {
      J=unique(fixesBRAN[[n]]$ID[I])
      match[[n]][i,2]=length(J) #All events that match
      
      K=which(fixesBRAN[[n]]$ID%in%J & fixesBRAN[[n]]$Location==1)
      tmp2=fixesBRAN[[n]][K,]
      L=which(fixesBRAN[[n]]$ID%in%J & fixesBRAN[[n]]$Location==1 & fixesBRAN[[n]]$Date2==min(tmp2$Date2))
      
      if(length(L)==1)
      {
        tmp_slp[,,i,3]<-slpBRAN[[n]][,,L]
        tmp_gv[,,i,3]<-gvBRAN[[n]][,,L]
        tmp_rain[,,i,3]<-rainBRAN[[n]][,,L]
        tmp_wind[,,i,3]<-windBRAN[[n]][,,L]
      } else {
        tmp_slp[,,i,3]<-apply(slpBRAN[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_gv[,,i,3]<-apply(rainBRAN[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_rain[,,i,3]<-apply(windBRAN[[n]][,,L],c(1,2),mean,na.rm=T)
        tmp_wind[,,i,3]<-apply(windBRAN[[n]][,,L],c(1,2),mean,na.rm=T)
      }
    } else match[[n]][i,2]=0
    
  }
  
  I=which(match[[n]][,1]>0)
  match_slp[,,n,1,1]=apply(tmp_slp[,,-I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,1,2]=apply(tmp_slp[,,I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,1,3]=apply(tmp_slp[,,I,2],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,1]=apply(tmp_gv[,,-I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,2]=apply(tmp_gv[,,I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,1,3]=apply(tmp_gv[,,I,2],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,1]=apply(tmp_rain[,,-I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,2]=apply(tmp_rain[,,I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,1,3]=apply(tmp_rain[,,I,2],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,1]=apply(tmp_wind[,,-I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,2]=apply(tmp_wind[,,I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,1,3]=apply(tmp_wind[,,I,2],c(1,2),mean,na.rm=T)
  
  I=which(match[[n]][,2]>0)
  match_slp[,,n,2,1]=apply(tmp_slp[,,-I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,2,2]=apply(tmp_slp[,,I,1],c(1,2),mean,na.rm=T)
  match_slp[,,n,2,3]=apply(tmp_slp[,,I,3],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,1]=apply(tmp_gv[,,-I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,2]=apply(tmp_gv[,,I,1],c(1,2),mean,na.rm=T)
  match_gv[,,n,2,3]=apply(tmp_gv[,,I,3],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,1]=apply(tmp_rain[,,-I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,2]=apply(tmp_rain[,,I,1],c(1,2),mean,na.rm=T)
  match_rain[,,n,2,3]=apply(tmp_rain[,,I,3],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,1]=apply(tmp_wind[,,-I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,2]=apply(tmp_wind[,,I,1],c(1,2),mean,na.rm=T)
  match_wind[,,n,2,3]=apply(tmp_wind[,,I,3],c(1,2),mean,na.rm=T)
  
}


contour(apply(match_gv[,,,1,1],c(1,2),mean))
contour(apply(match_gv[,,,1,2],c(1,2),mean))
contour(apply(match_gv[,,,1,3],c(1,2),mean))

pdf(file=paste(figdir,"ECL_slpgv_matched_d01_",cat[c],"_Controlvs2EAC_location_initial.pdf",sep=""),width=10.5,height=4,pointsize=12)
layout(cbind(1,2,3))
par(mar=c(3,3,4,1))
contour(apply(match_slp[,,,2,1],c(1,2),mean),main="a) Unmatched Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2,1],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,2,2],c(1,2),mean),main="b) Matched Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2,2],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,2,3],c(1,2),mean),main="c) Corresponding Events",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2,3],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()

pdf(file=paste(figdir,"Figure7.pdf",sep=""),width=7,height=4,pointsize=12)
layout(cbind(1,2))
par(mar=c(3,3,4,1))
contour(apply(match_slp[,,,2,1],c(1,2),mean),main="",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2,1],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
contour(apply(match_slp[,,,2,2],c(1,2),mean),main="",axes=F,cex.main=2,lwd=2)
contour(apply(match_gv[,,,2,2],c(1,2),mean),levels=0:8,add=T,lwd=2,lty=2,col="darkgrey")
points(0.5,0.5,pch=4,col="black",lwd=2,cex=2)
axis(1,at=seq(0,1,0.25),seq(-500,500,250))
axis(2,at=seq(0,1,0.25),seq(-500,500,250))
dev.off()