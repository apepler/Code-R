setwd("/srv/ccrc/data45/z3478332/WRF/output/ECLrain")
library(RNetCDF)

P<-array(0,c(4,2,3,2,2))
dimnames(P)[[1]]<-c("All","6 mm","12 mm","24 mm")
dimnames(P)[[2]]<-types
dimnames(P)[[3]]<-paste("R",1:3)
dimnames(P)[[4]]<-c("Total","Ratio")
dimnames(P)[[5]]<-c("50km","10km")

P_ESB=P

dom=c("d01","d02")
dom2=c("","_d02")
dims=cbind(c(215,144),c(325,200))
types=c("","_notopo")
r=1:3
cat="rad2_p100"

for(d in 1:2)
{
  
  a=open.nc(paste("/srv/ccrc/data45/z3478332/WRF/output/extracted_data/WRF_",dom[d],"_LH_PRCP_BRAN.nc",sep=""))
  lat=var.get.nc(a,"lat")
  lon=var.get.nc(a,"lon")
  
  mask=matrix(NaN,dims[1,d],dims[2,d])
  mask[(lat>=-40 & lat<=-24 & lon<=160 & lon>=150)]=1
  
  a=open.nc(paste("/srv/ccrc/data34/z3478332/WRF_",dom[d],"_ESB_mask.nc",sep=""))
  maskE=var.get.nc(a,"ESB")
  maskE[maskE==0]=NaN
  AllRain<-ECLrain<-NoECLrain<-array(0,c(dims[1,d],dims[2,d],4,2,3))
  
  for(t in 1:2)
    for(r in 1:3)
    {
      a=open.nc(paste("ECLrain_0708_extremes",dom2[d],"_R",r,types[t],"_",cat,"_v3.nc",sep=""))
      AllRain[,,,t,r]<-var.get.nc(a,"allrain")
      ECLrain[,,,t,r]<-var.get.nc(a,"ECLrain")
    }
  
  dimnames(AllRain)[[4]]=types
  dimnames(AllRain)[[5]]=paste("R",1:3)
  Ratio=ECLrain/AllRain
  
  for(t in 1:2)
    for(r in 1:3)
      for(x in 1:4)
      {
        P[x,t,r,1,d]=mean(AllRain[,,x,t,r]*mask,na.rm=T)
        P_ESB[x,t,r,1,d]=mean(AllRain[,,x,t,r]*maskE,na.rm=T)
        P[x,t,r,2,d]=mean(Ratio[,,x,t,r]*mask,na.rm=T)
        P_ESB[x,t,r,2,d]=mean(Ratio[,,x,t,r]*maskE,na.rm=T)
      }
  
}

P2<-array(0,c(4,2,4,2))
dimnames(P2)[[1]]<-c("All","6 mm","12 mm","24 mm")
dimnames(P2)[[2]]<-types
dimnames(P2)[[3]]<-c("Total","ECL","NoECL","Ratio")
dimnames(P2)[[4]]<-c("50km","10km")

P2_ESB=P2

dom=c("d01","d02")
dom2=c("","d02_")
dims=cbind(c(215,144),c(325,200))
types=c("ERA-nonudge","ERA-nonudge_notopo")
dirs=c("/srv/ccrc/data34/z3478332/WRF/ERA-nonudge/","/srv/ccrc/data45/z3478332/WRF/output/ERAI_R2_nonudging_notopo/out/impact/")
end2=c("_2","")
r=1:3
cat="p100_rad2cv1"

for(d in 1:2)
{
  
  a=open.nc(paste("/srv/ccrc/data45/z3478332/WRF/output/extracted_data/WRF_",dom[d],"_LH_PRCP_BRAN.nc",sep=""))
  lat=var.get.nc(a,"lat")
  lon=var.get.nc(a,"lon")
  
  mask=matrix(NaN,dims[1,d],dims[2,d])
  mask[(lat>=-40 & lat<=-24 & lon<=160 & lon>=150)]=1
  
  a=open.nc(paste("/srv/ccrc/data34/z3478332/WRF_",dom[d],"_ESB_mask.nc",sep=""))
  maskE=var.get.nc(a,"ESB")
  maskE[maskE==0]=NaN
  AllRain<-ECLrain<-NoECLrain<-array(0,c(dims[1,d],dims[2,d],4,2))
  
  for(t in 1:2)
  {
    a=open.nc(paste(dirs[t],"ECLrain_annual_",dom2[d],types[t],"_",cat,end2[t],".nc",sep=""))
    AllRain[,,,t]<-apply(var.get.nc(a,"allrain"),c(1,2,3),mean)
    ECLrain[,,,t]<-apply(var.get.nc(a,"ECLrain"),c(1,2,3),mean)
  }
  
  NoECLrain=AllRain-ECLrain
  
  dimnames(AllRain)[[4]]=types
  Ratio=ECLrain/AllRain
  
  for(t in 1:2)
    for(x in 1:4)
    {
      P2[x,t,1,d]=mean(AllRain[,,x,t]*mask,na.rm=T)
      P2_ESB[x,t,1,d]=mean(AllRain[,,x,t]*maskE,na.rm=T)
      P2[x,t,2,d]=mean(ECLrain[,,x,t]*mask,na.rm=T)
      P2_ESB[x,t,2,d]=mean(ECLrain[,,x,t]*maskE,na.rm=T)
      P2[x,t,3,d]=mean(NoECLrain[,,x,t]*mask,na.rm=T)
      P2_ESB[x,t,3,d]=mean(NoECLrain[,,x,t]*maskE,na.rm=T)
      P2[x,t,4,d]=mean(Ratio[,,x,t]*mask,na.rm=T)
      P2_ESB[x,t,4,d]=mean(Ratio[,,x,t]*maskE,na.rm=T)
    }
  
}

