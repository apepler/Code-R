setwd("~/Documents/ECLs/CSV")
listcsv <- dir(pattern = "*rad2cv1.csv")
listcsv=listcsv[substr(listcsv,1,4)=="ECLe"]

names=data.frame(model=rep("aaa",15),wrf=rep("aaa",15),stringsAsFactors=F)
a=strsplit(listcsv,"_")
for(i in 1:15){
  names[i,1]=a[[i]][3]
  names[i,2]=substr(a[[i]][4],4,5)
}

dataFiles <- lapply(listcsv, read.csv)

annual=cbind(seq(1990,2009),matrix(0,20,15))
monthly=cbind(seq(1,12),matrix(0,12,15))
MSLP=cbind(c(0,seq(980,1015,5)),matrix(0,9,15))
CV=cbind(seq(1,6,0.5),matrix(0,11,15))
DD=cbind(seq(0,16,4),matrix(0,5,15))
RAD=cbind(c(0,seq(1,6,0.5)),matrix(0,12,15))

for(i in 1:15)
{
  a=dataFiles[[i]]
  yy=floor(a$Date1/10000)
  mm=floor(a$Date1/100) %% 100
  
  for(j in 1:length(annual[,1]))
  {
    I=which(yy==annual[j,1])
    annual[j,i+1]=length(I)
  }
  
  for(j in 1:length(monthly[,1]))
  {
    I=which(mm==monthly[j,1])
    monthly[j,i+1]=length(I)
  }
  
  for(j in 1:length(MSLP[,1]))
  {
    if(j==length(MSLP[,1])) I=which(a$MSLP2>=MSLP[j,1]) else I=which(a$MSLP2>=MSLP[j,1] & a$MSLP2<MSLP[j+1,1]) 
    MSLP[j,i+1]=length(I)
  }
  
  for(j in 1:length(CV[,1]))
  {
    if(j==length(CV[,1])) I=which(a$CV2>=CV[j,1]) else I=which(a$CV2>=CV[j,1] & a$CV2<CV[j+1,1]) 
    CV[j,i+1]=length(I)
  }
  
  for(j in 1:length(DD[,1]))
  {
    if(j==length(DD[,1])) I=which(a$Length2>=DD[j,1]) else I=which(a$Length2>=DD[j,1] & a$Length2<DD[j+1,1]) 
    DD[j,i+1]=length(I)
  }
  
  for(j in 1:length(RAD[,1]))
  {
    if(j==length(RAD[,1])) I=which(a$Rad2>=RAD[j,1]) else I=which(a$Rad2>=RAD[j,1] & a$Rad2<RAD[j+1,1]) 
    RAD[j,i+1]=length(I)
  }
}

ave=apply(annual[,2:16],2,mean)
coolP=apply(monthly[5:10,2:16],2,sum)/apply(monthly[,2:16],2,sum)

names2=data.frame(Model=as.factor(names[,1]),WRF=as.factor(names[,2]))
plot(ave,coolP,col=as.integer(names2[,1]),pch=as.integer(names2[,2]),lwd=2,type="p",
     xlim=c(10,40),ylim=c(0.35,0.65),xlab="ECLs p.a.",ylab="Prop in May-Oct")
points(23.6,0.6257,pch=4,col=8,lwd=2)
legend("bottomright",c(levels(names2[,2]),levels(names2[,1]),"MERRA"),ncol=3,
       pch=c(unique(as.integer(names2[,2])),rep(4,6)),col=c(rep(6,3),unique(as.integer(names2[,1])),8),pt.lwd=2)

stats=data.frame(names2,aveDays=rep(0,15),aveDD=rep(0,15),aveMSLP=rep(0,15),aveCV=rep(0,15),aveRAD=rep(0,15),
                 cool=rep(0,15),DD2=rep(0,15),MSLP1000=rep(0,15),CV2=rep(0,15),RAD4=rep(0,15))

stats$aveDays=apply(annual[,2:16],2,mean)
stats$cool=apply(monthly[5:10,2:16],2,sum)/apply(monthly[,2:16],2,sum)
stats$DD2=apply(DD[3:5,2:16],2,sum)/apply(DD[,2:16],2,sum)
stats$MSLP1000=apply(MSLP[1:5,2:16],2,sum)/apply(MSLP[,2:16],2,sum)
stats$CV2=apply(CV[3:11,2:16],2,sum)/apply(CV[,2:16],2,sum)
stats$RAD4=apply(RAD[1:7,2:16],2,sum)/apply(RAD[,2:16],2,sum)

for(i in 1:15)
{
  a=dataFiles[[i]]
  stats$aveDD[i]=mean(a$Length2)
  stats$aveMSLP[i]=mean(a$MSLP2)
  stats$aveCV[i]=mean(a$CV2)
  stats$aveRAD[i]=mean(a$Rad2)
}

wrf=data.frame(levels(names2[,2]),aveDays=rep(0,3),aveDD=rep(0,3),aveMSLP=rep(0,3),aveCV=rep(0,3),aveRAD=rep(0,3),
               cool=rep(0,3),DD2=rep(0,3),MSLP1000=rep(0,3),CV2=rep(0,3),RAD4=rep(0,3))
cmip=data.frame(levels(names2[,1]),aveDays=rep(0,5),aveDD=rep(0,5),aveMSLP=rep(0,5),aveCV=rep(0,5),aveRAD=rep(0,5),
               cool=rep(0,5),DD2=rep(0,5),MSLP1000=rep(0,5),CV2=rep(0,5),RAD4=rep(0,5))
for(i in 2:11){
  a=aggregate(stats[,i+1],by=list(names2[,2]),mean)
  wrf[,i]=a[,2]
  a=aggregate(stats[,i+1],by=list(names2[,1]),mean)
  cmip[,i]=a[,2]
} 

plot(stats$aveCV,stats$aveMSLP,col=as.integer(names2[,1]),pch=as.integer(names2[,2]),lwd=2,type="p",
     xlim=c(1.5,1.7),ylim=c(995,1005),xlab="Mean curvature",ylab="Mean MSLP")
points(1.55,1000.2,pch=4,col=8,lwd=2)
legend("topright",c(levels(names2[,2]),levels(names2[,1]),"MERRA"),ncol=3,
       pch=c(unique(as.integer(names2[,2])),rep(4,6)),col=c(rep(6,3),unique(as.integer(names2[,1])),8),pt.lwd=2)


##Different approach - only looking at strongest 10/year = top 200

stats2=data.frame(names2,aveDD=rep(0,15),aveMSLP=rep(0,15),aveCV=rep(0,15),aveRAD=rep(0,15),
                 cool=rep(0,15),DD2=rep(0,15),MSLP1000=rep(0,15),CV2=rep(0,15),RAD4=rep(0,15))

for(i in 1:15)
{
  a=dataFiles[[i]]
  b=order(a$CV,a$Length2,decreasing=T)
  a=a[b[1:200],]
  mm=floor(a$Date1/100) %% 100
  
  stats2$aveDD[i]=mean(a$Length2)
  stats2$aveMSLP[i]=mean(a$MSLP2)
  stats2$aveCV[i]=mean(a$CV2)
  stats2$aveRAD[i]=mean(a$Rad2)
  
  stats2$cool[i]=length(which(mm>=5 & mm<=10))/200
  stats2$DD2[i]=length(which(a$Length2>=8))/200
  stats2$CV2[i]=length(which(a$CV2>=2))/200
  stats2$MSLP1000[i]=length(which(a$MSLP2<1000))/200
  stats2$RAD4[i]=length(which(a$Rad2<=4))/200
}

wrf2=data.frame(levels(names2[,2]),aveDD=rep(0,3),aveMSLP=rep(0,3),aveCV=rep(0,3),aveRAD=rep(0,3),
               cool=rep(0,3),DD2=rep(0,3),MSLP1000=rep(0,3),CV2=rep(0,3),RAD4=rep(0,3))
cmip2=data.frame(levels(names2[,1]),aveDays=rep(0,5),aveDD=rep(0,5),aveMSLP=rep(0,5),aveCV=rep(0,5),aveRAD=rep(0,5),
                cool=rep(0,5),DD2=rep(0,5),MSLP1000=rep(0,5),CV2=rep(0,5),RAD4=rep(0,5))
for(i in 2:10){
  a=aggregate(stats2[,i+1],by=list(names2[,2]),mean)
  wrf2[,i]=a[,2]
  a=aggregate(stats2[,i+1],by=list(names2[,1]),mean)
  cmip2[,i]=a[,2]
} 

