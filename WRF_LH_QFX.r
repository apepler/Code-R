setwd("/srv/ccrc/data45/z3478332/WRF/output/extracted_data")
library(RNetCDF)

a=open.nc("WRF_d01_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
LH=var.get.nc(a,"LH_d01")
P=var.get.nc(a,"PRCP_d01")

a=open.nc("WRF_d01_QFX_BRAN.nc")
QFX=var.get.nc(a,"QFX_d01")

a=open.nc("WRF_d01_LH_PRCP_BRAN_noeac.nc")
LH2=var.get.nc(a,"LH_d01")
P2=var.get.nc(a,"PRCP_d01")

a=open.nc("WRF_d01_QFX_BRAN_noeac.nc")
QFX2=var.get.nc(a,"QFX_d01")

#a=open.nc("WRF_d01_LH_PRCP_BRAN_2eac.nc")
a=open.nc("WRF_d01_LH_PRCP_BRAN_2eac_v2.nc")
LH3=var.get.nc(a,"LH_d01")
P3=var.get.nc(a,"PRCP_d01")



LH_diff=LH2-LH
P_diff=P2/P
Q_diff=QFX2-QFX
P_diff2=P3/P
LH_diff2=LH3-LH

mask=matrix(NaN,215,144)
mask[(lat>=-40 & lat<=-24 & lon<=160 & lon>=150)]=1

diff=matrix(0,3,5)
for(i in 1:3)
{
diff[i,1]=mean(LH_diff[i,,]*mask,na.rm=T)
diff[i,2]=mean(P_diff[i,,]*mask,na.rm=T)
diff[i,3]=mean(Q_diff[i,,]*mask,na.rm=T)
diff[i,4]=mean(LH_diff2[i,,]*mask,na.rm=T)
diff[i,5]=mean(P_diff2[i,,]*mask,na.rm=T)
}
apply(diff,2,mean)

clist=c("red","blue","orange")
plot(as.vector(Q_diff[1,,]*mask),as.vector(LH_diff[1,,]*mask),pch=4,lwd=2,type="p",col="red",
     xlab="Change in moisture flux (Kg.m2)",ylab="Change in precipitation (mm)",
     ylim=range(LH_diff,na.rm=T),xlim=range(Q_diff,na.rm=T))
for(i in 2:3) points(as.vector(Q_diff[i,,]*mask),as.vector(LH_diff[i,,]*mask),pch=4,lwd=2,col=clist[i])
for(i in 1:3) abline(lm(as.vector(LH_diff[i,,]*mask)~as.vector(Q_diff[i,,]*mask)),lwd=6,col=clist[i])
legend("topleft",c("R1","R2","R3"),pch=4,pt.lwd=2,col=clist)

for(i in 1:3) print(cor(as.vector(Q_diff[i,,]*mask),as.vector(LH_diff[i,,]*mask),use="pairwise.complete.obs"))

########## Mean precipitation change in ESB region.
########## Need the ESB that I regridded to the WRF resiolution.


a=open.nc("/srv/ccrc/data34/z3478332/WRF_d01_ESB_mask.nc")
mask=var.get.nc(a,"ESB")
mask[mask==0]=NaN

diff=matrix(0,3,3)
for(i in 1:3)
{
  diff[i,1]=mean(P[i,,]*mask,na.rm=T)
  diff[i,2]=mean(P2[i,,]*mask,na.rm=T)
  diff[i,3]=mean(P3[i,,]*mask,na.rm=T)
}


######### SST play

library(RNetCDF)
year=2007:2008
month=1:12

SSTdiff=matrix(0,24,3)
colnames(SSTdiff)=c("Average","Maximum","Cells>=1")
CONTDir="/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_bran_2007/"
DATADir="/srv/ccrc/data37/z3478332/WRF/WRF_boundary/default_bran_2007_2eac/"

n=1
for(y in 1:2)
  for(m in 1:12)
  {
    a=open.nc(paste(CONTDir,"wrflowinp_d01_",year[y],"-",sprintf("%2.2d",m),sep=""))
    b=open.nc(paste(DATADir,"wrflowinp_d01_",year[y],"-",sprintf("%2.2d",m),sep=""))
    sst=var.get.nc(b,"SST")-var.get.nc(a,"SST")
    c=dim(sst)
    sst2=apply(sst[,,1:(c[3]-1)],c(1,2),mean,na.rm=T)
    a=apply(is.na(sst),c(1,2),sum)
    sst2[a>6]=NaN 
    
    SSTdiff[n,1]=mean(sst2,na.rm=T)
    SSTdiff[n,2]=max(sst2,na.rm=T)
    SSTdiff[n,3]=length(which(sst2>=1))
    n=n+1
  }

library(RNetCDF)
year=2007:2008
month=1:12

SSTdiff=matrix(0,24,3)
colnames(SSTdiff)=c("Average","Maximum","Cells>=1")
CONTDir="/srv/ccrc/data36/z3478332/WRF/WRF_boundary/default_bran_2007/"
DATADir="/srv/ccrc/data37/z3478332/WRF/WRF_boundary/default_bran_2007_2eac/"

sst=array(0,c(215,144,24))
n=1
for(y in 1:2)
  for(m in 1:12)
  {
    a=open.nc(paste(CONTDir,"wrflowinp_d01_",year[y],"-",sprintf("%2.2d",m),sep=""))
    b=open.nc(paste(DATADir,"wrflowinp_d01_",year[y],"-",sprintf("%2.2d",m),sep=""))
    tmp=var.get.nc(b,"SST")-var.get.nc(a,"SST")
    c=dim(tmp)
    sst[,,n]=apply(tmp[,,1:(c[3]-1)],c(1,2),mean,na.rm=T)
    a=apply(is.na(tmp),c(1,2),sum)
    sst2[a>6]=NaN 
    
    SSTdiff[n,1]=mean(sst[,,n],na.rm=T)
    SSTdiff[n,2]=max(sst[,,n],na.rm=T)
    SSTdiff[n,3]=length(which(sst2>=1))
    n=n+1
  }

image(apply(sst,c(1,2),mean))
mask<-matrix(NaN,215,144)
I=which(apply(sst,c(1,2),mean)>=0.1)
mask[I]=1

n=1
for(y in 1:2)
  for(m in 1:12)
  {
    SSTdiff[n,1]=mean(sst[,,n]*mask,na.rm=T)
    SSTdiff[n,2]=max(sst[,,n],na.rm=T)
    SSTdiff[n,3]=length(which(sst[,,n]>=1))
    n=n+1
  }

########### Okay, need to play with ECL statistics for each time period.
EAC_mask_d01=mask
EAC_mask_d02=mask
save(EAC_mask_d01,EAC_mask_d02,file="EAC_mask.RData")

load(file="EAC_mask.RData")

n=1
for(y in 1:2)
  for(m in 1:12)
  {
    print(n)
    a=open.nc(paste(CONTDir,"wrflowinp_d01_",year[y],"-",sprintf("%2.2d",m),sep=""))
    b=open.nc(paste(DATADir,"wrflowinp_d01_",year[y],"-",sprintf("%2.2d",m),sep=""))
    tmp1=var.get.nc(b,"SST")-var.get.nc(a,"SST")
    
    a=open.nc(paste(CONTDir,"wrflowinp_d02_",year[y],"-",sprintf("%2.2d",m),sep=""))
    b=open.nc(paste(DATADir,"wrflowinp_d02_",year[y],"-",sprintf("%2.2d",m),sep=""))
    tmp2=var.get.nc(b,"SST")-var.get.nc(a,"SST")
    
    c=dim(tmp1)
    sst1<-matrix(0,c[3]-1,4)
    colnames(sst1)<-c("Average","Maximum","Cells>=1","Cells>=2")
    sst2<-sst1
    
    for(k in 1:(c[3]-1))
    {
      sst1[k,1]=mean(tmp1[,,k]*EAC_mask_d01,na.rm=T)
      sst1[k,2]=max(tmp1[,,k])
      sst1[k,3]=length(which(tmp1[,,k]>=1))
      sst1[k,4]=length(which(tmp1[,,k]>=2))
      sst2[k,1]=mean(tmp2[,,k]*EAC_mask_d02,na.rm=T)
      sst2[k,2]=max(tmp2[,,k])
      sst2[k,3]=length(which(tmp2[,,k]>=1))
      sst2[k,4]=length(which(tmp2[,,k]>=2))
    }
    
    if(n==1)
    {
      EAC_d01_all<-sst1
      EAC_d02_all<-sst2
    } else
    {
      EAC_d01_all<-rbind(EAC_d01_all,sst1)
      EAC_d02_all<-rbind(EAC_d02_all,sst2)
    }
    n=n+1
  }

save(EAC_mask_d01,EAC_mask_d02,EAC_d01_all,EAC_d02_all,file="EAC_mask.RData")

#####################
## ECL component of rainfall

#setwd("/srv/ccrc/data36/z3478332/WRF/output/")
library(RNetCDF)
a=open.nc("WRF_d02_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
mask=matrix(NaN,325,200)
mask[(lat>=-40 & lat<=-24 & lon<=160 & lon>=150)]=1

a=open.nc("/srv/ccrc/data34/z3478332/WRF_d02_ESB_mask.nc")
maskE=var.get.nc(a,"ESB")
maskE[maskE==0]=NaN

setwd("/srv/ccrc/data45/z3478332/WRF/output/ECLrain")
types=c("BRAN","BRAN_noeac","BRAN_2eac2")
r=1:3
cat="rad2_p100"
AllRain<-ECLrain<-NoECLrain<-array(0,c(325,200,3,3))

for(t in 1:3)
  for(r in 1:3)
  {
    a=open.nc(paste("ECLrain_0708_d02_R",r,"_",types[t],"_",cat,"_v2.nc",sep=""))
    AllRain[,,t,r]<-var.get.nc(a,"allrain")
    ECLrain[,,t,r]<-var.get.nc(a,"ECLrain")
  }

dimnames(AllRain)[[3]]=types
dimnames(AllRain)[[4]]=paste("R",1:3)

NoECLrain=AllRain-ECLrain

P<-array(0,c(3,3,3))
dimnames(P)[[1]]<-types
dimnames(P)[[2]]<-paste("R",1:3)
dimnames(P)[[3]]<-c("All","ECL","NoECL")

P_ESB=P

for(t in 1:3)
  for(r in 1:3)
  {
    P[t,r,1]=mean(AllRain[,,t,r]*mask,na.rm=T)
    P_ESB[t,r,1]=mean(AllRain[,,t,r]*maskE,na.rm=T)
    
    P[t,r,2]=mean(ECLrain[,,t,r]*mask,na.rm=T)
    P_ESB[t,r,2]=mean(ECLrain[,,t,r]*maskE,na.rm=T)
    
    P[t,r,3]=mean(NoECLrain[,,t,r]*mask,na.rm=T)
    P_ESB[t,r,3]=mean(NoECLrain[,,t,r]*maskE,na.rm=T)
  }

apply(P,c(1,3),mean)

tmp=ECLrain/AllRain
Ratio=matrix(0,3,3)
for(t in 1:3)
  for(r in 1:3)
  {
    Ratio[t,r]=mean(tmp[,,t,r]*mask,na.rm=T)
  }

########### What about look at all cases with 6-hourly accumulations above 6-12-24 mm/6hr

setwd("/srv/ccrc/data36/z3478332/WRF/output/")
library(RNetCDF)
a=open.nc("WRF_d02_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
mask=matrix(NaN,325,200)
mask[(lat>=-40 & lat<=-24 & lon<=160 & lon>=150)]=1

a=open.nc("/srv/ccrc/data34/z3478332/WRF_d02_ESB_mask.nc")
maskE=var.get.nc(a,"ESB")
maskE[maskE==0]=NaN

setwd("/srv/ccrc/data45/z3478332/WRF/output/")
types=c("BRAN","BRAN_noeac","BRAN_2eac2")
r=1:3
cat="rad2_p100"
AllRain<-ECLrain<-NoECLrain<-array(0,c(325,200,3,3,3))

for(t in 1:3)
  for(r in 1:3)
  {
    a=open.nc(paste("ECLrain_0708_extremes_d02_R",r,"_",types[t],"_",cat,"_v2.nc",sep=""))
    AllRain[,,,t,r]<-var.get.nc(a,"allrain")
    ECLrain[,,,t,r]<-var.get.nc(a,"ECLrain")
  }

dimnames(AllRain)[[4]]=types
dimnames(AllRain)[[5]]=paste("R",1:3)

NoECLrain=AllRain-ECLrain

Ratio=ECLrain/AllRain

P<-array(0,c(3,3,3))
dimnames(P)[[1]]<-c("6 mm","12 mm","24 mm")
dimnames(P)[[2]]<-types
dimnames(P)[[3]]<-paste("R",1:3)

P_ESB=P

for(t in 1:3)
  for(r in 1:3)
    for(x in 1:3)
  {
    P[x,t,r]=mean(Ratio[,,x,t,r]*mask,na.rm=T)
    P_ESB[x,t,r]=mean(Ratio[,,x,t,r]*maskE,na.rm=T)
    }



############# NCEP-WRF Stuff for chapter 3

setwd("/srv/ccrc/data36/z3478332/WRF/output/")
library(RNetCDF)
a=open.nc("WRF_d01_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
mask=matrix(NaN,215,144)
mask[(lat>=-40 & lat<=-24 & lon<=160 & lon>=150)]=1

a=open.nc("/srv/ccrc/data34/z3478332/WRF_d01_ESB_mask.nc")
maskE=var.get.nc(a,"ESB")
maskE[maskE==0]=NaN

r=1:3
cat="proj100_rad2cv1"
AllRain<-ECLrain<-NoECLrain<-array(0,c(215,144,4,3))

for(r in 1:3)
  {
    a=open.nc(paste("/srv/ccrc/data34/z3478332/WRF/R",r,"/ECLrain_extremes2_ncep_wrfR",r,"_",cat,"_9009.nc",sep=""))
    AllRain[,,,r]<-var.get.nc(a,"allrain")
    ECLrain[,,,r]<-var.get.nc(a,"ECLrain")
  }

NoECLrain=AllRain-ECLrain

P<-array(0,c(4,3,3))
dimnames(P)[[1]]<-c("Total","6mm","12mm","24mm")
dimnames(P)[[2]]<-paste("R",1:3)
dimnames(P)[[3]]<-c("All","ECL","NoECL")

P_ESB=P
for(t in 1:4)
for(r in 1:3)
  {
    P[t,r,1]=mean(AllRain[,,t,r]*mask,na.rm=T)
    P_ESB[t,r,1]=mean(AllRain[,,t,r]*maskE,na.rm=T)
    
    P[t,r,2]=mean(ECLrain[,,t,r]*mask,na.rm=T)
    P_ESB[t,r,2]=mean(ECLrain[,,t,r]*maskE,na.rm=T)
    
    P[t,r,3]=mean(NoECLrain[,,t,r]*mask,na.rm=T)
    P_ESB[t,r,3]=mean(NoECLrain[,,t,r]*maskE,na.rm=T)
  }

image(apply(ECLrain/AllRain,c(1,2),mean))
maskE2=matrix(0,215,144)
maskE2[maskE==1]=1
contour(maskE2,add=T)

tmp=ECLrain/AllRain
Ratio=matrix(0,4,3)
for(t in 1:4)
  for(r in 1:3)
  {
    Ratio[t,r]=mean(tmp[,,t,r]*maskE,na.rm=T)
  }

##### What about ECL contribution to daily rain
rm(list=ls())
library(RNetCDF)
a=open.nc("WRF_d01_LH_PRCP_BRAN.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
mask=matrix(NaN,215,144)
mask[(lat>=-40 & lat<=-24 & lon<=160 & lon>=150)]=1

a=open.nc("/srv/ccrc/data34/z3478332/WRF_d01_ESB_mask.nc")
maskE=var.get.nc(a,"ESB")
maskE[maskE==0]=NaN

wrfv=c("R1","R2","R3")
cmip="NNRP"
runyears="1950-2010"
ystart=c(1990,1995,2000,2005)

#dayrain=array(0,c(215,144,7305,3))

for(w in 1:3)
{
  n=1
  for(s in 1:4)
  {
    dir=paste("/srv/ccrc/data30/z3393020/NARCliM/filtered/",cmip,"/",wrfv[w],"/",runyears,"/d01/",sep="")
    fname=paste(dir,"CCRC_NARCliM_DAY_",ystart[s],"-",ystart[s]+4,"_pracc_fl.nc",sep="")
    a=open.nc(fname)
    tmp=var.get.nc(a,"pracc_fl")
    dayrain[,,n:(n+dim(tmp)[3]-1),w]=tmp
    n=dim(tmp)[3]
  }
}



c=0
do while(c.lt.dimsizes(cmip))
  w=0
do while(w.lt.dimsizes(wrfv))
  dayrain=new((/7305,144,215/),"float")
dayrain!1="south_north"
dayrain!2="east_west"
dayrain!0="threshold"
ECLrain&threshold=(/6,12,24/)
ECLrain@description="Number of hours with at least X mm/6hr accumulation within 500km radius of low centre"

dir="/srv/ccrc/data30/z3393020/NARCliM/filtered/"+cmip(c)+"/"+wrfv(w)+"/"+runyears(c)+"/d01/"


#### Need average ESB rain in 2007-2008 from AWAP, d01 and d02 - for reviewer comments
#### Calculate both bias in average rain, & average bias
#### Regridding both d01 & d02 to AWAP grid?

library(RNetCDF)
a=open.nc("/srv/ccrc/data02/z3236814/data/AWAP/MONTHLY/monthly_means_grid_0.05.dir.V3/pre_1200.nc")
time=var.get.nc(a,"time")
rain=var.get.nc(a,"pre",start=c(1,1,97),count=c(NA,NA,24))
AWAP=apply(rain,c(1,2),sum)/2
WRF<-array(0,c(886,691,3,2))

library(R.matlab)
readMat('~/Documents/GDI/Useful_ECL.mat')->UsefulE
maskE<-t(UsefulE$mask)

setwd("/srv/ccrc/data45/z3478332/WRF/output/extracted_data")
library(RNetCDF)
a=open.nc("WRF_d01_LH_PRCP.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
latt=as.vector(lat)
lont=as.vector(lon)
P=var.get.nc(a,"PRCP_d01")

library(akima)
for(i in 1:3)
{
  a=as.vector(P[i,,])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,UsefulE$x,UsefulE$y)
  WRF[,,i,1]=b$z
}

library(RNetCDF)
a=open.nc("WRF_d02_LH_PRCP.nc")
lat=var.get.nc(a,"lat")
lon=var.get.nc(a,"lon")
latt=as.vector(lat)
lont=as.vector(lon)
P=var.get.nc(a,"PRCP_d02")

library(akima)
for(i in 1:3)
{
  a=as.vector(P[i,,])
  a[which(is.na(a))]=0
  b=interp(lont,latt,a,UsefulE$x,UsefulE$y)
  WRF[,,i,2]=b$z
}
WRF=WRF/2


bias<-array(0,c(3,2,2))

for(i in 1:3)
  for(j in 1:2)
  {
    bias[i,j,1]=mean(WRF[,,i,j]*maskE,na.rm=T)-mean(AWAP*maskE,na.rm=T)
    bias[i,j,2]=mean((WRF[,,i,j]/AWAP)*maskE,na.rm=T)
  }





