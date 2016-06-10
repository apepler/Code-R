rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

dates_9009=as.numeric(format.Date(seq.Date(as.Date("1990-01-01"),as.Date("2009-12-31"),by=1),"%Y%m%d"))
dates_6079=as.numeric(format.Date(seq.Date(as.Date("2060-01-01"),as.Date("2079-12-31"),by=1),"%Y%m%d"))
month_9009<-floor(dates_9009/100)%%100
month_6079<-floor(dates_6079/100)%%100

now<-array(NaN,c(length(dates_9009),15,5))
future<-array(NaN,c(length(dates_9009),12,5))
dimnames(now)[[3]]<-dimnames(future)[[3]]<-c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
dimnames(now)[[2]]<-rep("aaa",15)
dimnames(future)[[2]]<-rep("aaa",12)

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
cmip2=c("NNRP","ECHAM5","MK30","MIROC","CCCMA")
wrf=c("R1","R2","R3")
intcol=c(5,10,10,9,9)
loccol=c(NA,11,11,10,10)
datecol=c(3,4,4,4,4)

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    dimnames(now)[[2]][n]<-paste(cmip[i],wrf[j])
    
    filelist=c(paste("Fei/vor_",cmip2[i],"_",wrf[j],"_1990-2010_Andrew.csv",sep=""),
               paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_9009.csv",sep=""),
               paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_9009.csv",sep=""),
               paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.5.csv",sep=""),
               paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg0.95.csv",sep=""))
    
    data=read.csv(filelist[1])
    for(m in 1:length(dates_9009))
    {
      I=which(data[,datecol[1]]==dates_9009[m])
      if(length(I)>0) now[m,n,1]=max(data[I,intcol[1]])
    }
    
    for(ff in 2:5)
    {
      data=read.csv(filelist[ff])
      data$Location2=0
      I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
      data$Location2[I]<-1
      I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
      data$Location2[I]<-1
      I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
      data$Location2[I]<-1
      
      for(m in 1:length(dates_9009))
      {
        I=which(data[,datecol[ff]]==dates_9009[m] & data$Location==1)
        if(length(I)>0) now[m,n,ff]=max(data[I,intcol[ff]])
      }
    }
    n=n+1  
  }
a=apply(!is.na(now),c(2,3),sum)/20

n=1
for(i in 2:5)
  for(j in 1:3)
  {
    dimnames(future)[[2]][n]<-paste(cmip[i],wrf[j])
    
    filelist=c(paste("Fei/vor_",cmip2[i],"_",wrf[j],"_2060-2080_Andrew.csv",sep=""),
               paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_6079.csv",sep=""),
               paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_6079.csv",sep=""),
               paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.5.csv",sep=""),
               paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg0.95.csv",sep=""))
    
    data=read.csv(filelist[1])
    for(m in 1:length(dates_6079))
    {
      I=which(data[,datecol[1]]==dates_6079[m])
      if(length(I)>0) future[m,n,1]=max(data[I,intcol[1]])
    }
    
    for(ff in 2:5)
    {
      data=read.csv(filelist[ff])
      data$Location2=0
      I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
      data$Location2[I]<-1
      I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
      data$Location2[I]<-1
      I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
      data$Location2[I]<-1
      for(m in 1:length(dates_6079))
      {
        I=which(data[,datecol[ff]]==dates_6079[m] & data$Location==1)
        if(length(I)>0) future[m,n,ff]=max(data[I,intcol[ff]])
      }
    }
    n=n+1  
  }


### Change in intensity threshold that occurs on X days per year

thresh=c(15,10,5,2,1,0.2)

nowint<-array(0,c(15,5,length(thresh),2))
dimnames(nowint)[[1]]<-dimnames(now)[[2]]
dimnames(nowint)[[2]]<-dimnames(now)[[3]]
dimnames(nowint)[[3]]<-paste(thresh,"pa")
dimnames(nowint)[[4]]<-c("Cool","Warm")
futureint<-nowint[4:15,,,]

C1<-which(month_9009>=5 & month_9009<=10)
C2<-which(month_6079>=5 & month_6079<=10)

for(i in 1:15)
  for(j in 1:5)
  {
    c=sort(now[C1,i,j],decreasing=T)
    w=sort(now[-C1,i,j],decreasing=T)
    
    for(t in 1:length(thresh)) {
      nowint[i,j,t,2]=w[20*thresh[t]]
      nowint[i,j,t,1]=c[20*thresh[t]]
    }
    
    if(i<13)
    {
      c=sort(future[C2,i,j],decreasing=T)
      w=sort(future[-C2,i,j],decreasing=T)
      
      for(t in 1:length(thresh)) {
        futureint[i,j,t,2]=w[20*thresh[t]]
        futureint[i,j,t,1]=c[20*thresh[t]]
      }
    }
  }


change=100*((futureint/nowint[4:15,,,])-1)
names(dimnames(change))<-c("Source","Method","Threshold","Season")

library(ggplot2)
library(reshape2)
data=melt(change[,,,])  
pdf(file=paste("Figures/ECLdays_maxCV_change_location2_ALL.pdf",sep=""),width=10,height=5)
ggplot(data, aes(x = Threshold, y = value, fill = Season)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  #scale_y_continuous(breaks=seq(-50, 50, 25),limits = c(-25, 50)) +
  theme_bw() + ylab("% change in intensity threshold") + xlab("")
dev.off()


########## Testing if pdf is different for top 37 days p.a. (equiv to Speer et al.)
## Or, alternately, using the threshold

tests1<-tests2<-matrix(0,13,5)
thresh=37

for(i in 1:12)
  for(j in 1:5)
  {
    a=sort(now[,i+3,j],decreasing=T)
    b=sort(future[,i,j],decreasing=T)
    
    c=ks.test(a[1:(thresh*20)],b[1:(thresh*20)])
    tests1[i,j]=c$p.value
    
    t=a[thresh*20]
    c=ks.test(a[a>=t],b[b>=t])
    tests2[i,j]=c$p.value   
  }

a=sort(now[,4:15,j],decreasing=T)
b=sort(future[,,j],decreasing=T)

############# If I want to look at change in pdf, best to multiply PG by 100
now[,,4:5]=now[,,4:5]*100
future[,,4:5]=future[,,4:5]*100

nowD<-array(NaN,c(512,15,5))
dimnames(nowD)[[2]]<-dimnames(now)[[2]]
dimnames(nowD)[[3]]<-dimnames(now)[[3]]
names(dimnames(nowD))<-c("X","Source","Method")
futureD=nowD[,4:15,]
count<-array(0,c(15,5,2))

for(i in 1:15)
  for(j in 1:5)
  {
    a=sort(now[,i,j],decreasing=T)
    t=a[thresh*20]
    b=density(a[a>=t],from=0,to=10)
    nowD[,i,j]=b$y
    count[i,j,1]=length(which(a>=t))
    
    if(i>=4)
    {
      a=sort(future[,i-3,j],decreasing=T)
      b=density(a[a>=t],from=0,to=10)
      futureD[,i-3,j]=b$y
      count[i,j,2]=length(which(a>=t))
    }
  }


changeD=futureD-nowD[,4:15,]
changeDm=apply(changeD,c(1,3),median,na.rm=T)
plot(NA,xlim=c(0,10),ylim=c(-0.11,0.11),xlab="Intensity",ylab="Change in frequency")
abline(h=0,col="grey",lwd=2)
for(i in 1:5) lines(b$x,changeDm[,i],col=i,lwd=2)
legend("topright",legend=c("GV","UM 150km","UM 50km","PG 150km","PG 50km"),lwd=2,col=1:5)


changeD2=array(NaN,dim(changeD))
for(i in 1:12)
  for(j in 1:5)
    changeD2[,i,j]=(futureD[,i,j]*count[i+3,j,2]-nowD[,i+3,j]*count[i+3,j,1])/20

changeDm2=apply(changeD2,c(1,3),median,na.rm=T)
plot(NA,xlim=c(0,10),ylim=c(-10,10),xlab="Intensity",ylab="Change in frequency")
abline(h=0,col="grey",lwd=2)
for(i in 1:5) lines(b$x,changeDm2[,i],col=i,lwd=2)
legend("topright",legend=c("GV","UM 150km","UM 50km","PG 150km","PG 50km"),lwd=2,col=1:5)

###########
## Okay, quantile matching.
## First step: divide into deciles

thresh=37

quantN<-quantF<-array(NaN,c(10,12,5))

for(i in 1:12)
  for(j in 1:5)
  {
    a=sort(now[,i+3,j],decreasing=T)
    t=a[thresh*20]
    q=c(quantile(a[a>=t],seq(0,0.9,0.1)),Inf)
    
    for(k in 1:10)
    {
      quantN[k,i,j]=length(which(a>=q[k] & a<q[k+1]))
      quantF[k,i,j]=length(which(future[,i,j]>=q[k] & future[,i,j]<q[k+1]))
    }
  }

changeQ=quantF-quantN
changeQ2=apply(changeQ,c(1,3),median)
plot(NA,xlim=c(1,10),ylim=c(-25,25),xlab="Quantile",ylab="Change in frequency")
abline(h=0,col="grey",lwd=2)
for(i in 1:5) lines(1:10,changeQ2[,i],col=i,lwd=2)
legend("topright",legend=c("GV","UM 150km","UM 50km","PG 150km","PG 50km"),lwd=2,col=1:5)

##### Re-do with seasons - but annual quantiles

thresh=37

quantN<-quantF<-array(NaN,c(10,12,5,2))
C1<-which(month_9009>=5 & month_9009<=10)
C2<-which(month_6079>=5 & month_6079<=10)

for(i in 1:12)
  for(j in 1:5)
  {
    a=sort(now[,i+3,j],decreasing=T)
    t=a[thresh*20]
    q=c(quantile(a[a>=t],seq(0,0.9,0.1)),Inf)
    
    for(k in 1:10)
    {
      quantN[k,i,j,1]=length(which(now[C1,i+3,j]>=q[k] & now[C1,i+3,j]<q[k+1]))
      quantF[k,i,j,1]=length(which(future[C2,i,j]>=q[k] & future[C2,i,j]<q[k+1]))
      quantN[k,i,j,2]=length(which(now[-C1,i+3,j]>=q[k] & now[-C1,i+3,j]<q[k+1]))
      quantF[k,i,j,2]=length(which(future[-C2,i,j]>=q[k] & future[-C2,i,j]<q[k+1]))
    }
  }

changeQ=100*((quantF/quantN)-1)
dimnames(changeQ)[[1]]<-paste("Q",1:10)
dimnames(changeQ)[[2]]<-dimnames(future)[[2]]
dimnames(changeQ)[[3]]<-dimnames(future)[[3]]
dimnames(changeQ)[[4]]<-c("Cool","Warm")
names(dimnames(changeQ))<-c("Quantile","Source","Method","Season")

changeQ2=apply(changeQ,c(1,3,4),median)
plot(NA,xlim=c(1,10),ylim=c(-50,25),xlab="Quantile",ylab="% Change in frequency")
abline(h=0,col="grey",lwd=2)
for(i in 1:5) lines(1:10,changeQ2[,i,2],col=i,lwd=2)
legend("topright",legend=c("GV","UM 150km","UM 50km","PG 150km","PG 50km"),lwd=2,col=1:5)

data=melt(changeQ)
ggplot(data[data$Season=="Cool",], aes(x = Quantile, y = value, fill = Method)) +
  geom_boxplot() +
  #scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  scale_y_continuous(breaks=seq(-50, 100, 25),limits = c(-50, 100)) +
  theme_bw() + ylab("% change in frequency") + xlab("")