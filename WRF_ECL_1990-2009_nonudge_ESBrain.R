rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
library(maps)

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->Useful
readMat('~/Documents/GDI/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')

cat="p100_rad2cv1"
cat2="rad2_p100_cv1.0"

tlist=c("ERA-nonudge","ERA-nonudge_notopo")
dirs=c("/srv/ccrc/data34/z3478332/WRF/ERA-nonudge/","/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nonudging_notopo/out/impact/")

a=open.nc(paste("/srv/ccrc/data45/z3478332/WRF/output/extracted_data/WRF_d01_LH_PRCP_BRAN.nc",sep=""))
lat1=var.get.nc(a,"lat")
lon1=var.get.nc(a,"lon")
a=open.nc(paste("/srv/ccrc/data34/z3478332/WRF_d01_ESB_mask.nc",sep=""))
mask1=var.get.nc(a,"ESB")
mask1[mask1==0]=NaN
a=open.nc(paste("/srv/ccrc/data45/z3478332/WRF/output/extracted_data/WRF_d01_LH_PRCP_BRAN.nc",sep=""))
lat2=var.get.nc(a,"lat")
lon2=var.get.nc(a,"lon")
a=open.nc(paste("/srv/ccrc/data34/z3478332/WRF_d02_ESB_mask.nc",sep=""))
mask2=var.get.nc(a,"ESB")
mask2[mask2==0]=NaN

allrain<-ECLrain<-ECLrain_loc<-events<-fixes<-list()

ESBrain=array(0,c(20,4,3,4))
dimnames(ESBrain)[[2]]=c("All","6mm","12mm","24mm")
dimnames(ESBrain)[[3]]=c("All","ECL","ECL_loc")
dimnames(ESBrain)[[4]]=c(paste(tlist,"d01"),paste(tlist,"d02"))

end2=c("_2","")
n=1
for(dom in c("d01","d02"))
  for(r in 1:2)
  {
    events[[n]]<-read.csv(paste("outputUM_",tlist[r],"/",cat,end2[r],"/",dom,"/ECLevents_",tlist[r],"_",dom,"_",cat2,"_typing_impactsC2.csv",sep=""))
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    events[[n]][events[[n]]==-Inf]=NaN
    
    fixes[[n]]<-read.csv(paste("outputUM_",tlist[r],"/",cat,end2[r],"/",dom,"/ECLfixes_",tlist[r],"_",dom,"_",cat2,"_typing_impactsC2.csv",sep=""))
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    fixes[[n]]$Location2<-0
    #       I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    #       fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    
    ### Rain stuff
    if(dom=="d02") {
      mask=mask2
      dom2="d02_" 
      } else {
        dom2=""
        mask=mask1
      }

    a=open.nc(paste(dirs[r],"ECLrain_annual_",dom2,tlist[r],"_",cat,end2[r],".nc",sep=""))
    ECLrain[[n]]<-var.get.nc(a,"ECLrain")
    allrain[[n]]<-var.get.nc(a,"allrain")
    ECLrain_loc[[n]]<-var.get.nc(a,"ECLrain_loc")
    
    for(i in 1:20)
      for(j in 1:4)
      {
        ESBrain[i,j,1,n]=mean(allrain[[n]][,,j,i]*mask,na.rm=T)
        ESBrain[i,j,2,n]=mean(ECLrain[[n]][,,j,i]*mask,na.rm=T)
        ESBrain[i,j,3,n]=mean(ECLrain_loc[[n]][,,j,i]*mask,na.rm=T)
      }

    n=n+1     
  }

