rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

nowW<-array(0,c(12,20))
dimnames(nowW)[[1]]<-rep("aaa",12)
dimnames(nowW)[[2]]<-seq(1990,2009)
futureW<-futureC<-nowC<-nowW
dimnames(futureW)[[2]]<-dimnames(futureC)[[2]]<-seq(2060,2079)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      dimnames(nowC)[[1]][n]<-dimnames(futureC)[[1]][n]<-dimnames(nowW)[[1]][n]<-dimnames(futureW)[[1]][n]<-paste(cmip[i],wrf[j])
      
      filelist1=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep="")
      filelist2=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep="")
      
      data=read.csv(filelist1)
      yy=floor(data$Date1/10000)
      mm=floor((data$Date1%%10000)/100)
      years=unique(yy)
      for(t in 1:length(years))
      {
        I=which(yy==years[t] & mm>=5 & mm <=10)
        nowC[n,t]=max(data[I,11])
        I=which(yy==years[t] & (mm>=11 | mm <=4))
        nowW[n,t]=max(data[I,11])
      }
      
      data=read.csv(filelist2)
      yy=floor(data$Date1/10000)
      mm=floor((data$Date1%%10000)/100)
      years=unique(yy)
      for(t in 1:length(years))
      {
        I=which(yy==years[t] & mm>=5 & mm <=10)
      futureC[n,t]=max(data[I,11])
      I=which(yy==years[t] & (mm>=11 | mm <=4))
      futureW[n,t]=max(data[I,11])
      }      

      n=n+1  
    }

### Okay, have annual max CV, how to do ARI?

library(extRemes)

rlevs=c(2,5,10,20,50,100)
nowC.ari<-matrix(0,13,6)
rownames(nowC.ari)<-c(dimnames(nowC)[[1]],"All")
colnames(nowC.ari)<-rlevs
futureW.ari<-nowW.ari<-futureC.ari<-nowC.ari

for(i in 1:12)
{
  a=fevd(nowW[i,])
  nowW.ari[i,]=return.level(a,rlevs)
  a=fevd(futureW[i,])
  futureW.ari[i,]=return.level(a,rlevs)
  a=fevd(nowC[i,])
  nowC.ari[i,]=return.level(a,rlevs)
  a=fevd(futureC[i,])
  futureC.ari[i,]=return.level(a,rlevs)
}
a=fevd(as.vector(nowW))
nowW.ari[13,]=return.level(a,rlevs)
a=fevd(as.vector(futureW))
futureW.ari[13,]=return.level(a,rlevs)
a=fevd(as.vector(nowC))
nowC.ari[13,]=return.level(a,rlevs)
a=fevd(as.vector(futureC))
futureC.ari[13,]=return.level(a,rlevs)

# plot(rlevs,rep(NaN,6),xlab=c("Return level"),ylab="Maximum event curvature",ylim=c(2,6),xlim=c(0,25))
# lines(rlevs,nowC.ari[13,],col="blue",lwd=3)
# lines(rlevs,futureC.ari[13,],col="blue",lwd=3,lty=2)
# lines(rlevs,nowW.ari[13,],col="red",lwd=3)
# lines(rlevs,futureW.ari[13,],col="red",lwd=3,lty=2)
# legend("topleft",c("Cool","Warm"),col=c("blue","red"),lwd=3)

boxplot(futureW.ari[1:12,]-nowW.ari[1:12,])
abline(h=0,col="red")
boxplot(futureC.ari[1:12,]-nowC.ari[1:12,])
abline(h=0,col="red")

###### Redo with a 100-member bootstrap

library(boot)
rlevs=c(2,5,10,20,50,100)
retlev <- function(x,inds,rlev)
{
  a=fevd(x[inds])
  return.level(a,rlev)
}

for(i in 1:12)
{
  a=boot(as.vector(nowW[i,]),retlev,100,rlev=rlevs)
  nowW.ari[i,]=apply(a$t,2,median)
  a=boot(as.vector(futureW[i,]),retlev,100,rlev=rlevs)
  futureW.ari[i,]=apply(a$t,2,median)
  a=boot(as.vector(nowC[i,]),retlev,100,rlev=rlevs)
  nowC.ari[i,]=apply(a$t,2,median)
  a=boot(as.vector(futureC[i,]),retlev,100,rlev=rlevs)
  futureC.ari[i,]=apply(a$t,2,median)
}

a=boot(as.vector(nowW),retlev,100,rlev=rlevs)
nowW.ari[13,]=apply(a$t,2,median)
a=boot(as.vector(futureW),retlev,100,rlev=rlevs)
futureW.ari[13,]=apply(a$t,2,median)
a=boot(as.vector(nowC),retlev,100,rlev=rlevs)
nowC.ari[13,]=apply(a$t,2,median)
a=boot(as.vector(futureC),retlev,100,rlev=rlevs)
futureC.ari[13,]=apply(a$t,2,median)

boxplot(futureC.ari[1:12,]-nowC.ari[1:12,])
abline(h=0,col="red")
boxplot(futureW.ari[1:12,]-nowW.ari[1:12,])
abline(h=0,col="red")

### What about PoT
### Seems to need evenly spaced data - so I need to convert to daily data?

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
daylist=read.csv("daylist.csv")

ECLdaysN<-ECLdaysF<-matrix(0,length(daylist[,1]),12)
rownames(ECLdaysN)=daylist[,1]
rownames(ECLdaysF)=daylist[,2]
colnames(ECLdaysN)<-colnames(ECLdaysF)<-rep("aaa",12)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    colnames(ECLdaysN)[n]<-colnames(ECLdaysF)[n]<-paste(cmip[i],wrf[j])
    
    filelist1=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep="")
    filelist2=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep="")
    
    data=read.csv(filelist1)
    for(t in 1:length(daylist[,1]))
    {
      I=which(data$Date==daylist[t,1] & data$Location==1)
      if(length(I)>0) ECLdaysN[t,n]=max(data$CV[I])
    }
    
    data=read.csv(filelist2)
    for(t in 1:length(daylist[,2]))
    {
      I=which(data$Date==daylist[t,2] & data$Location==1)
      if(length(I)>0) ECLdaysF[t,n]=max(data$CV[I])
    } 
    
    n=n+1  
  }

library(extRemes)

mm=(floor(daylist/100)%%100)
I=which(mm[,1]>=5 & mm[,1]<=10)
nowC=ECLdaysN[I,]
I=which(mm[,1]>=11 | mm[,1]<=4)
nowW=ECLdaysN[I,]
I=which(mm[,2]>=5 & mm[,2]<=10)
futureC=ECLdaysF[I,]
I=which(mm[,2]>=11 | mm[,2]<=4)
futureW=ECLdaysF[I,]

rlevs=c(2,5,10,20,50,100)
nowC.ari<-matrix(0,13,6)
rownames(nowC.ari)<-c(dimnames(nowC)[[2]],"All")
colnames(nowC.ari)<-rlevs
futureW.ari<-nowW.ari<-futureC.ari<-nowC.ari

for(i in 1:12)
{
  a=fevd(nowC[,i],threshold=quantile(nowC[,i],0.95),type="GP",time.units="184/year")
  nowC.ari[i,]=return.level(a,rlevs)
  a=fevd(nowW[,i],threshold=quantile(nowW[,i],0.95),type="GP",time.units="181.25/year")
  nowW.ari[i,]=return.level(a,rlevs)
  a=fevd(futureC[,i],threshold=quantile(futureC[,i],0.95),type="GP",time.units="184/year")
  futureC.ari[i,]=return.level(a,rlevs)
  a=fevd(futureW[,i],threshold=quantile(futureW[,i],0.95),type="GP",time.units="181.25/year")
  futureW.ari[i,]=return.level(a,rlevs)
}

boxplot(futureC.ari[1:12,]-nowC.ari[1:12,])
abline(h=0,col="red")
boxplot(futureW.ari[1:12,]-nowW.ari[1:12,])
abline(h=0,col="red")

########### Okay, GEV for everyone!
rm(list=ls())
rlevs=c(2,5,10,20,50,100)
nowC.ari<-array(0,c(12,6,5))
dimnames(nowC.ari)[[1]]<-rep("aaa",12)
dimnames(nowC.ari)[[2]]<-rlevs
dimnames(nowC.ari)[[3]]<-c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
futureW.ari<-nowW.ari<-futureC.ari<-nowC.ari

daylist=read.csv("daylist.csv")
mm=(floor(daylist/100)%%100)
cmip=c("echam5","csiromk3","miroc","cccma")
cmipF=c("ECHAM5","MK30","MIROC","CCCMA")
wrf=c("R1","R2","R3")

intcol=c(5,10,10,9,9)

n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(nowC.ari)[[1]][n]<-dimnames(futureC.ari)[[1]][n]<-dimnames(nowW.ari)[[1]][n]<-dimnames(futureW.ari)[[1]][n]<-paste(cmip[i],wrf[j])
    
    filelist1=c(paste("Fei/vor_",cmipF[i],"_",wrf[j],"_1990-2010_Andrew.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
    filelist2=c(paste("Fei/vor_",cmipF[i],"_",wrf[j],"_2060-2080_Andrew.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
    
    for(ff in 1:5)
    {
      data=read.csv(filelist1[ff])
      days=rep(0,length(daylist[,1]))
      for(t in 1:length(daylist[,1]))
      {
        if(ff==1) I=which(data$Date==daylist[t,1]) else I=which(data$Date==daylist[t,1] & data$Location==1)
        if(length(I)>0) days[t]=max(data[I,intcol[ff]])
      }
      
      I=which(mm[,1]>=5 & mm[,1]<=10)
      daysC=days[I]
      I=which(mm[,1]>=11 | mm[,1]<=4)
      daysW=days[I]
      
      a=fevd(daysC,threshold=quantile(daysC,0.95),type="GP",time.units="184/year")
      nowC.ari[n,,ff]=return.level(a,rlevs)
      a=fevd(daysW,threshold=quantile(daysW,0.95),type="GP",time.units="184/year")
      nowW.ari[n,,ff]=return.level(a,rlevs)
      
      data=read.csv(filelist2[ff])
      days=rep(0,length(daylist[,1]))
      for(t in 1:length(daylist[,1]))
      {
        if(ff==1) I=which(data$Date==daylist[t,2]) else I=which(data$Date==daylist[t,2] & data$Location==1)
        if(length(I)>0) days[t]=max(data[I,intcol[ff]])
      }
      
      I=which(mm[,1]>=5 & mm[,1]<=10)
      daysC=days[I]
      I=which(mm[,1]>=11 | mm[,1]<=4)
      daysW=days[I]
      
      a=fevd(daysC,threshold=quantile(daysC,0.95),type="GP",time.units="184/year")
      futureC.ari[n,,ff]=return.level(a,rlevs)
      a=fevd(daysW,threshold=quantile(daysW,0.95),type="GP",time.units="184/year")
      futureW.ari[n,,ff]=return.level(a,rlevs)
      
      
    }
    
    n=n+1  
  }

changeC=(futureC.ari[1:12,,1]-nowC.ari[1:12,,1])/nowC.ari[1:12,,1]
for(i in 2:5) changeC=rbind(changeC,(futureC.ari[1:12,,i]-nowC.ari[1:12,,i])/nowC.ari[1:12,,i])
changeW=(futureW.ari[1:12,,1]-nowW.ari[1:12,,1])/nowW.ari[1:12,,1]
for(i in 2:5) changeW=rbind(changeW,(futureW.ari[1:12,,i]-nowW.ari[1:12,,i])/nowW.ari[1:12,,1])

boxplot(100*(futureC.ari[,6,]-nowC.ari[,6,])/nowC.ari[,6,],ylab="% change",main="Change in 100-year return level (cool)")
abline(h=0,col="red")
boxplot(100*(futureW.ari[,6,]-nowW.ari[,6,])/nowW.ari[,6,],ylab="% change",main="Change in 100-year return level (warm)")
abline(h=0,col="red")
