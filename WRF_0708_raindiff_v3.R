rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful
readMat('~/Documents/Data/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')


dayrain<-array(NaN,c(731,4,4,3))
dimnames(dayrain)[[2]]=c("Mean rain","Max rain","Cells >= 5mm","Cells >= 25mm")
dimnames(dayrain)[[3]]=c("Control","NoTopo","BRAN","NoEAC")
dimnames(dayrain)[[4]]=c("R1","R2","R3")
dayrain_d02=dayrain

wrfv=c("R1","R2","R3")
case=c("","_notopo","_BRAN","_BRAN_noeac")
ddir=c(36,36,37,37)

for(i in 1:4)
  for(j in 1:3)
  {
    dir=paste("/srv/ccrc/data",ddir[i],"/z3478332/WRF/output/ERAI_",wrfv[j],"_nudging_default_2007",case[i],"/out/",sep="")
    a=read.delim(paste(dir,"dayrain_d01.txt",sep=""),sep=";",header=F,stringsAsFactors=F)[[1]]
    b=read.delim(paste(dir,"dayrain_d02.txt",sep=""),sep=";",header=F,stringsAsFactors=F)[[1]]
    
    for(k in 1:4)
    {
      dayrain[,k,i,j]=as.numeric(substr(a,((k-1)*7+1),k*7))
      dayrain_d02[,k,i,j]=as.numeric(substr(b,((k-1)*7+1),k*7))
    }
  }


count5a=matrix(0,4,3)
for(i in 1:4)
  for(j in 1:3)
    count5a[i,j]=length(which(dayrain_d02[,1,i,j]>=5))

change_topo=matrix(0,4,4)
change_topo[1,1:3]=count5[2,]/count5[1,]
change_topo[2,1:3]=count5a[2,]/count5a[1,]
change_topo[3,1:3]=count50[2,]/count50[1,]
change_topo[4,1:3]=count50a[2,]/count50a[1,]
change_topo[,4]=apply(change_topo[,1:3],1,mean,na.rm=T)

change_eac=matrix(0,4,4)
change_eac[1,1:3]=count5[4,]/count5[3,]
change_eac[2,1:3]=count5a[4,]/count5a[3,]
change_eac[3,1:3]=count50[4,]/count50[3,]
change_eac[4,1:3]=count50a[4,]/count50a[3,]
change_eac[,4]=apply(change_eac[,1:3],1,mean,na.rm=T)
