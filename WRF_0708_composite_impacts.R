rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)

lat=seq(-500,500,10)
lon=seq(-500,500,10)
library(sp)

dist=matrix(0,101,101)
for(i in 1:101)
  for(j in 1:101)
    dist[i,j]=sqrt(lat[i]^2 + lon[j]^2)

dist2<-dist3<-matrix(NaN,101,101)
dist2[dist<=500]=1
dist3[dist<=250]=1

dom="d02"
cat="rad2_p100"

wdirs=c("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007_notopo/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN_noeac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN_2eac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_2eac/out/",
          "/srv/ccrc/data45/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN_2eac/out/")

wnames=c("R1","R2","R3",
         "R1_notopo","R2_notopo","R3_notopo",
         "R1_BRAN","R2_BRAN","R3_BRAN",
         "R1_BRAN_noeac","R2_BRAN_noeac","R3_BRAN_noeac",
         "R1_BRAN_2eac","R2_BRAN_2eac","R3_BRAN_2eac")

for(w in 1)
      {
  data=read.csv(paste("ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_typing.csv",sep=""),stringsAsFactors=F)
  events=read.csv(paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing.csv",sep=""),stringsAsFactors=F)
  data=data[,-1]
  events=events[,-1]
  
  filelistC=paste(wdirs[w],"ECLrain_d02_0708_",cat,"_centred.nc",sep="")
  filelistW=paste(wdirs[w],"ECLwind_d02_0708_",cat,".nc",sep="")
  
  a=open.nc(filelistC)
        tmp=var.get.nc(a,"ECLrain")
        tmp[tmp>=500]=NaN
        close.nc(a)
        a=dim(tmp)
        for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist2
        data$MeanRain500=apply(tmp,3,mean,na.rm=T)
        data$MaxRain500=apply(tmp,3,max,na.rm=T)
        for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist3
        data$MeanRain250=apply(tmp,3,mean,na.rm=T)
        data$MaxRain250=apply(tmp,3,max,na.rm=T)
        
        a=open.nc(filelistW)
        tmp=var.get.nc(a,"ECL_WS10")
        close.nc(a)
        a=dim(tmp)
        for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist2
        data$MeanWind500=apply(tmp,3,mean,na.rm=T)
        data$MaxWind500=apply(tmp,3,max,na.rm=T)
        for(x in 1:a[3]) tmp[,,x]=tmp[,,x]*dist3
        data$MeanWind250=apply(tmp,3,mean,na.rm=T)
        data$MaxWind250=apply(tmp,3,max,na.rm=T)
        
        events$MaxPointWind500<-events$MaxMeanWind500<-events$MaxPointRain500<-events$MaxMeanRain500<-events$TotalRain500<-0
        events$MaxPointWind250<-events$MaxMeanWind250<-events$MaxPointRain250<-events$MaxMeanRain250<-events$TotalRain250<-0
        for(k in 1:length(events$ID))
        {
          I=which(data$ID==events$ID[k] & data$Location==1)
          events$TotalRain500[k]=sum(data$MeanRain500[I],na.rm=T)
          events$MaxMeanRain500[k]=max(data$MeanRain500[I],na.rm=T)
          events$MaxPointRain500[k]=max(data$MaxRain500[I],na.rm=T)
          events$MaxMeanWind500[k]=max(data$MeanWind500[I],na.rm=T)
          events$MaxPointWind500[k]=max(data$MaxWind500[I],na.rm=T)
          events$TotalRain250[k]=sum(data$MeanRain250[I],na.rm=T)
          events$MaxMeanRain250[k]=max(data$MeanRain250[I],na.rm=T)
          events$MaxPointRain250[k]=max(data$MaxRain250[I],na.rm=T)
          events$MaxMeanWind250[k]=max(data$MeanWind250[I],na.rm=T)
          events$MaxPointWind250[k]=max(data$MaxWind250[I],na.rm=T)
        }
        write.csv(data,paste("ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_typing_impactsC.csv",sep=""))
        write.csv(events,paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing_impactsC.csv",sep=""))
}


####### Now, analyse!

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)

dom="d01"
cat="rad2_p100"

wnames=c("R1","R2","R3",
         "R1_notopo","R2_notopo","R3_notopo",
         "R1_BRAN","R2_BRAN","R3_BRAN",
         "R1_BRAN_noeac","R2_BRAN_noeac","R3_BRAN_noeac",
         "R2_BRAN_2eac","R2_BRAN_2eac","R3_BRAN_2eac")

events<-fixes<-list()

for(w in 1:length(wnames))
{
  fixes[[w]]=read.csv(paste("ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_typing_impacts.csv",sep=""),stringsAsFactors=F)
  fixes[[w]]$Date2=as.POSIXct(paste(as.character(fixes[[w]]$Date),substr(fixes[[w]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  events[[w]]=read.csv(paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing_impacts.csv",sep=""))
}


### Average of each statistic by dataset

statave=matrix(0,15,9)
rownames(statave)=wnames
colnames(statave)=c("Count","Mean(meanR)","Mean(maxR)","Mean(meanW)","Mean(maxW)","MeanRain>=6","MaxRain>=50","MeanWind>=50km/h","MaxWind>=80km/h")

impcol=17:20 ## For 500, 22:25 for 250
xthresh=c(6,50,13.9,22.2)

for(w in 1:length(wnames))
{
  tmp=events[[w]]
  statave[w,1]=length(tmp[,1])
  statave[w,2:5]=apply(tmp[,impcol],2,mean,na.rm=T)
  for(x in 1:4) statave[w,x+5]=length(which(tmp[,impcol[x]]>=xthresh[x]))
}

## What about by type, across the  different versions 
#types=c("EC","SC","TC","Mixed")
types=c("ET","IT","SSL","CL")
statave=array(0,c(15,4,9))
dimnames(statave)[[1]]=wnames
dimnames(statave)[[2]]=types
dimnames(statave)[[3]]=c("Count","Mean(meanR)","Mean(maxR)","Mean(meanW)","Mean(maxW)","MeanRain>=6","MaxRain>=50","MeanWind>=50km/h","MaxWind>=80km/h")

impcol=17:20 ## For 500, 22:25 for 250
xthresh=c(6,50,13.9,22.2)

for(w in 1:length(wnames))
  for(t in 1:length(types))
{
  I=which(events[[w]]$TypeSB==types[t])
    tmp=events[[w]][I,]
  statave[w,t,1]=length(tmp[,1])
  statave[w,t,2:5]=apply(tmp[,impcol],2,mean,na.rm=T)
  for(x in 1:4) statave[w,t,x+5]=length(which(tmp[,impcol[x]]>=xthresh[x]))
  }

###
statave2=matrix(0,15,10)
colnames(statave2)<-c("Count","Length2","MSLP2","CV2","CV>2","Bomb","Deepening rate","Formed","Entered","Intensified")

for(w in 1:length(wnames))
{
  events[[w]]$MaxNDR=0
  for(i in 1:length(events[[w]]$ID))
  {
    I=which(fixes[[w]]$ID==events[[w]]$ID[i] & fixes[[w]]$Location==1 & !is.na(fixes[[w]]$NDR))
    if(length(I)>0) events[[w]]$MaxNDR[i]=max(fixes[[w]]$NDR[I],na.rm=T)
  }
  
  tmp=events[[w]]
  statave2[w,1]=length(tmp[,1])
  statave2[w,2:4]=apply(tmp[,8:10],2,mean,na.rm=T)
  statave2[w,5]=length(which(tmp$CV2>=2))
  statave2[w,6]=sum(tmp$Bomb)
  statave2[w,7]=mean(tmp$MaxNDR,na.rm=T)
  
  for(i in 1:3) statave2[w,7+i]=length(which(tmp$EnteredFormed==i))

}