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
      if(length(I)>0) now[m,n,1]=mean(data[I,intcol[1]])
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
        I=which(data[,datecol[ff]]==dates_9009[m] & data$Location2==1)
        if(length(I)>0) now[m,n,ff]=mean(data[I,intcol[ff]])
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
      if(length(I)>0) future[m,n,1]=mean(data[I,intcol[1]])
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
        I=which(data[,datecol[ff]]==dates_6079[m] & data$Location2==1)
        if(length(I)>0) future[m,n,ff]=mean(data[I,intcol[ff]])
      }
    }
    n=n+1  
  }

dayrain_now<-read.table("DailyRain/NARCLIM_dayrain_d01_9009.dat")

dayspa=22

##Quick test for matches
dayspa=22
now2=array(0,dim(now))
for(i in 1:15)
  for(j in 1:5)
  {
    b=order(now[,i,j],decreasing=T)
    thresh=now[b[20*dayspa],i,j]
    if(is.na(thresh)) thresh=min(now[,i,j],na.rm=T)
    I=which(now[,i,j]>=thresh)
    now2[I,i,j]=1
  }

matchdays=array(NaN,c(15,5,5))
for(i in 1:15)
  for(j in 1:5)
    for(k in 1:5)
      matchdays[i,j,k]=length(which(now2[,i,j]>=1 & now2[,i,k]>=1))/length(which(now2[,i,j]>=1))


## Back to regular programming

ECLrain<-nECLrain<-ECLdays<-array(NaN,c(15,5,4))
dimnames(ECLrain)[[1]]<-dimnames(nECLrain)[[1]]<-dimnames(ECLdays)[[1]]<-dimnames(now)[[2]]
dimnames(ECLrain)[[2]]<-dimnames(nECLrain)[[2]]<-dimnames(ECLdays)[[1]]<-c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
dimnames(ECLrain)[[3]]<-dimnames(nECLrain)[[3]]<-dimnames(ECLdays)[[1]]<-c("1990-2009","2060-2079")

for(i in 1:15)
  for(j in 1:5)
  {
    b=order(now[,i,j],decreasing=T)
    thresh=now[b[20*dayspa],i,j]
    if(is.na(thresh)) thresh=min(now[,i,j],na.rm=T)
    
    I=unique(which(now[,i,j]>=thresh))
    ECLdays[i,j,3]=length(I)/20
    I=unique(c(which(now[,i,j]>=thresh),which(now[,i,j]>=thresh)-1))
    ECLdays[i,j,1]=length(I)/20
    ECLrain[i,j,1]=sum(dayrain_now[I,i],na.rm=T)
    nECLrain[i,j,1]=sum(dayrain_now[-I,i],na.rm=T)
    ECLrain[i,j,3]=mean(dayrain_now[I,i],na.rm=T)
    nECLrain[i,j,3]=mean(dayrain_now[-I,i],na.rm=T)
    
    
    if(i>3)
    {
      I=unique(which(future[,i-3,j]>=thresh))
      ECLdays[i,j,4]=length(I)/20
      I=unique(c(which(future[,i-3,j]>=thresh),which(future[,i-3,j]>=thresh)-1))
      ECLdays[i,j,2]=length(I)/20
      ECLrain[i,j,2]=sum(dayrain_future[I,i-3],na.rm=T)
      nECLrain[i,j,2]=sum(dayrain_future[-I,i-3],na.rm=T)
      ECLrain[i,j,4]=mean(dayrain_future[I,i-3],na.rm=T)
      nECLrain[i,j,4]=mean(dayrain_future[-I,i-3],na.rm=T)
    }
   
  }







boxplot((ECLrain[,,1]/(ECLrain[,,1]+nECLrain[,,1]))*100,ylab="% of ESB total rainfall")
abline(h=21,col="red")

boxplot(ECLrain[,,1]/nECLrain[,,1])

ECLchange=((ECLdays[,,4]/ECLdays[,,3])-1)*100
nECLchange=((nECLrain[,,2]/nECLrain[,,1])-1)*100

boxplot(ECLchange,ylim=c(-50,50))
abline(h=0,col="red")
boxplot(nECLchange,ylim=c(-25,25))
abline(h=0,col="red")

now_m=array(0,c(12,15,5,2))
dayspa=22

for(i in 1:15)
  for(j in 1:5)
  {
    b=order(now[,i,j],decreasing=T)
    thresh=now[b[20*dayspa],i,j]
    if(is.na(thresh)) thresh=min(now[,i,j],na.rm=T)
    
    for(m in 1:12)
    {
      I=unique(c(which(now[,i,j]>=thresh & month_9009==m),which(now[,i,j]>=thresh & month_9009==m)-1))
      now_m[m,i,j,1]=length(I)
      if(length(I)>0) now_m[m,i,j,2]=sum(dayrain_now[I,i],na.rm=T)/sum(dayrain_now[month_9009==m,i],na.rm=T)
    }
    
  }

now_m=now_m[,,c(1,2,4),]

now_mN<-apply(now_m[,1:3,,2],c(1,3),mean)
now_mC<-apply(now_m[,4:15,,2],c(1,3),mean)

clist=c("black","blue","red","orange")

plot(1:12,now_mN[,1]*100,col="blue",lwd=4,type="l",xlab="Month",ylab="% ESB rain due to ECLs",ylim=range(0,monthly[,,2]*100))
for(i in 2:3) lines(1:12,now_mN[,i]*100,col=clist[i+1],lwd=4)
for(i in 1:3) lines(1:12,now_mC[,i]*100,col=clist[i+1],lwd=4,lty=2)
legend("bottomleft",c("GV","LAP","PG"),col=clist[2:4],lwd=4,ncol=3,bty="n")
legend("topright",c("NCEP","CMIP"),col=1,lwd=4,lty=c(1,2),bty="n")

##### Okay, want change in total ECL rain split by season - May-Oct & Nov-Apr

Rain<-array(NaN,c(12,5,2))
dayspa=22
dimnames(Rain)[[1]]=dimnames(now)[[2]][4:15]
dimnames(Rain)[[2]]=dimnames(now)[[3]]
dimnames(Rain)[[3]]=c("May-October","Nov-April")

Ecount<-Rain

for(i in 1:12)
  for(j in 1:5)
  {
    b=order(now[,i+3,j],decreasing=T)
    thresh=now[b[20*dayspa],i+3,j]
    if(is.na(thresh)) thresh=min(now[,i+3,j],na.rm=T)
    
    ##May-October
    I=which(month_9009>=5 & month_9009<=10)
    d1=now[I,i+3,j]
    J=which(month_6079>=5 & month_6079<=10)
    d2=future[J,i,j]    
    I1=unique(c(which(d1>=thresh),which(d1>=thresh)-1))
    J1=unique(c(which(d2>=thresh),which(d2>=thresh)-1))    
    Rain[i,j,1]=100*((sum(dayrain_future[J[J1],i],na.rm=T)/sum(dayrain_now[I[I1],i+3],na.rm=T))-1)
    Ecount[i,j,1]=100*((length(which(d2>=thresh))/length(which(d1>=thresh)))-1)

    ##November-April
    I=which(month_9009>=11 | month_9009<=4)
    d1=now[I,i+3,j]
    J=which(month_6079>=11 | month_6079<=4)
    d2=future[J,i,j]    
    I1=unique(c(which(d1>=thresh),which(d1>=thresh)-1))
    J1=unique(c(which(d2>=thresh),which(d2>=thresh)-1))    
    Rain[i,j,2]=100*((sum(dayrain_future[J[J1],i],na.rm=T)/sum(dayrain_now[I[I1],i+3],na.rm=T))-1)
    Ecount[i,j,2]=100*((length(which(d2>=thresh))/length(which(d1>=thresh)))-1)
  }




#What about pdf of daily rain?
#v1 = rain categories
rbreak=c(0,1,5,10,25,Inf)

ECL_dist<-nECL_dist<-array(NaN,c(length(rbreak)-1,15,5,2))
dimnames(ECL_dist)[[1]]<-dimnames(nECL_dist)[[1]]<-rbreak[1:(length(rbreak)-1)]
dimnames(ECL_dist)[[2]]<-dimnames(nECL_dist)[[2]]<-dimnames(now)[[2]]
dimnames(ECL_dist)[[3]]<-dimnames(nECL_dist)[[3]]<-c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
dimnames(ECL_dist)[[4]]<-dimnames(nECL_dist)[[4]]<-c("1990-2009","2060-2079")


for(i in 1:15)
  for(j in 1:5)
  {
    b=order(now[,i,j],decreasing=T)
    thresh=now[b[20*dayspa],i,j]
    I=which(now[,i,j]>=thresh & (month_9009>=11 | month_9009<=4))
    now1<-dayrain_now[I-1,i]
    I=which((now[,i,j]<thresh | is.na(now[,i,j])==TRUE) & (month_9009>=11 | month_9009<=4))
    now2<-dayrain_now[I-1,i]
    
    if(i>3)
    {
      I=which(future[,i-3,j]>=thresh & (month_6079>=11 | month_6079<=4))
      fut1<-dayrain_future[I,i-3]
      I=which((future[,i-3,j]<thresh | is.na(future[,i-3,j])==TRUE) & (month_6079>=11 | month_6079<=4))
      fut2<-dayrain_future[I,i-3]
    }
    
    for(m in 1:(length(rbreak)-1))
    {
      I=which(now1>=rbreak[m] & now1<rbreak[m+1])
      ECL_dist[m,i,j,1]<-length(I)/20
      I=which(now2>=rbreak[m] & now2<rbreak[m+1])
      nECL_dist[m,i,j,1]<-length(I)/20
      
      if(i>3)
      {
        I=which(fut1>=rbreak[m] & fut1<rbreak[m+1])
        ECL_dist[m,i,j,2]<-length(I)/20
        I=which(fut2>=rbreak[m] & fut2<rbreak[m+1])
        nECL_dist[m,i,j,2]<-length(I)/20
      }
    }
    
  }

ECLprop=100*(ECL_dist/(ECL_dist+nECL_dist))
apply(ECLprop[,1:3,,1],c(1,3),mean)


ECLchange=100*((ECL_dist[,,,2]/ECL_dist[,,,1])-1)
ECLchange2=apply(ECLchange[,4:15,],c(1,3),median)
nECLchange=100*((nECL_dist[,,,2]/nECL_dist[,,,1])-1)
nECLchange2=apply(nECLchange[,4:15,],c(1,3),median)

##v2 - actual pdfs

rainobs=density(AWAPrain2[,4],from=0,to=100)
rdensity=array(NaN,c(512,12,5,3,2))
dimnames(rdensity)[[2]]=dimnames(now)[[2]][3:15]
dimnames(rdensity)[[3]]=dimnames(now)[[3]]
dimnames(rdensity)[[4]]=c("All","ECL","non-ECL")
dimnames(rdensity)[[5]]=c("1990-2009","2060-2079")

for(i in 1:12)
  for(j in 1:5)
  {
    b=order(now[,i+3,j],decreasing=T)
    thresh=now[b[20*dayspa],i+3,j]
    if(is.na(thresh)) thresh=min(now[,i+3,j],na.rm=T)
    
    I=unique(c(which(now[,i+3,j]>=thresh),which(now[,i+3,j]>=thresh)-1))
    rdensity[,i,j,1,1]=density(dayrain_now[,i+3],na.rm=T,from=0,to=100)$y
    rdensity[,i,j,2,1]=density(dayrain_now[I,i+3],na.rm=T,from=0,to=100)$y
    rdensity[,i,j,3,1]=density(dayrain_now[-I,i+3],na.rm=T,from=0,to=100)$y
    
    I=unique(c(which(future[,i,j]>=thresh),which(future[,i,j]>=thresh)-1))
    rdensity[,i,j,1,2]=density(dayrain_future[,i],na.rm=T,from=0,to=100)$y
    rdensity[,i,j,2,2]=density(dayrain_future[I,i],na.rm=T,from=0,to=100)$y
    rdensity[,i,j,3,2]=density(dayrain_future[-I,i],na.rm=T,from=0,to=100)$y
  }

plot(rainobs$x,apply(rdensity[,,1,1,1],1,mean),log="x",type="l",lwd=3,col=1,ylim=c(0,0.51),
     xlab="ESB rainfall",ylab="Density",main="Distribution of rainfall in CMIP-WRF models")
for(i in 1:15) lines(rainobs$x,rdensity[,i,1,1,2],type="l",col="darkgray")
lines(rainobs$x,apply(rdensity[,,1,1,1],1,mean),lwd=3)

plot(rainobs$x,apply(rdensity[,,1,1,2]-rdensity[,,1,1,1],1,mean),log="x",type="l",lwd=3,col=1,ylim=c(-0.1,0.1),
     xlab="ESB rainfall",ylab="Change in density",main="Change in distribution of rainfall in CMIP-WRF models")
for(i in 1:15) lines(rainobs$x,rdensity[,i,1,1,2]-rdensity[,i,1,1,1],type="l",col="darkgray")
lines(rainobs$x,apply(rdensity[,,1,1,2]-rdensity[,,1,1,1],1,mean),lwd=3)
abline(h=0,col="red")

a=100*((rdensity[,,1,1,2]/rdensity[,,1,1,1])-1)
plot(rainobs$x,apply(a,1,mean),log="x",type="l",lwd=3,col=1,ylim=c(-100,100),
     xlab="ESB rainfall",ylab="Change in density",main="% change in distribution of rainfall in CMIP-WRF models")
for(i in 1:15) lines(rainobs$x,a[,i],type="l",col="darkgray")
lines(rainobs$x,apply(a,1,mean),lwd=3)
abline(h=0,col="red")

collist=c("black","red","orange","blue","cyan")
plot(rainobs$x,apply(rdensity[,,1,3,2]-rdensity[,,1,3,1],1,mean),log="x",type="l",lwd=3,col=1,
     xlab="ESB rainfall",ylab="Change in density",main="Change in distribution of non-ECL rainfall in CMIP-WRF models",ylim=c(-0.02,0.02))
for(i in 2:5) lines(rainobs$x,apply(rdensity[,,i,3,2]-rdensity[,,i,3,1],1,mean),col=collist[i],lwd=3)
abline(h=0,col="gray",lty=2,lwd=2)
legend("bottomright",c("GV","UM 150km","UM 50km","PG 150km","PG 50km"),ncol=2,lwd=3,col=collist,bty="n")


########### Compare to observations
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->Useful

dates=seq.Date(as.Date("1969-12-30"),as.Date("2010-12-31"),by=1)
#AWAPrain=data.frame(Decade=floor(as.numeric(format.Date(dates,"%Y"))/10),Year=as.numeric(format.Date(dates,"%Y")),Date=as.numeric(format.Date(dates,"%Y%m%d")),Rain=rep(0,length(dates)))
AWAPrain=read.csv("DailyRain/AWAPrain.csv")
AWAPrain=AWAPrain[,2:5]

# for(i in 14613:length(dates))
# {
#   if(i%%100==0) print(AWAPrain[i,3])
#   
#   dir=paste("/media/Seagate Expansion Drive/Data/daily rainfall/rainfall_",AWAPrain[i,1],"0-",AWAPrain[i,1],"9/",sep="")
#   fname<-paste(dir,'/rainfall-',AWAPrain[i,2],'/r',AWAPrain[i,3],'.txt',sep="")
#   read.table(fname, sep="",skip=6,nrows=691)->data
#   as.matrix(data)->data
#   data[data<0]=0
#   data<-data[nrow(data):1,]
#   AWAPrain[i,4]=mean(data*Useful$mask,na.rm=T)
# }
# write.csv(AWAPrain,file="DailyRain/AWAPrain.csv")

### Now, daily ECL from ERAI and from MLDB


filelist=c("~/Documents/ECLs/Algorithm Comparison/MLDB.csv",
           "Fei/vor_NCEP_1990-2010_Andrew.csv",
           "outputUM/proj100/outputUM_ncep_rad2cv06/ECLfixes_umelb_ncep_proj100_rad2cv06_9009.csv",
           "Alejandro/ECLfixes_Alejandro_ncep_res150_9009_pg0.5.csv")

intcol=c(8,5,10,9)
nowA<-array(NaN,c(length(dates_9009),4))
for(ff in 1:4)
{
  data=read.csv(filelist[ff])
  if(ff==3) data$Time=(as.numeric(data$Time)-1)*6

   data$Location2=0
   I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
   data$Location2[I]<-1
   I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
   data$Location2[I]<-1
   I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
   data$Location2[I]<-1
  
  
  for(m in 1:length(dates_9009))
  {
    if(ff==2) I=which(data$Date==dates_9009[m]) else I=which(data$Date==dates_9009[m] & data$Location2==1)
    if(length(I)>0) nowA[m,ff]=mean(data[I,intcol[ff]])
  }
}

now2=array(0,dim(nowA))
I=unique(c(which(!is.na(nowA[,1])),which(!is.na(nowA[,1]))+1))
now2[I,1]=1

for(i in 2:4)
{
b=order(nowA[,i],decreasing=T)
thresh=nowA[b[20*22],i]
if(is.na(thresh)) thresh=min(nowA[,i],na.rm=T)
I=unique(c(which(nowA[,i]>=thresh),which(nowA[,i]>=thresh)+1)) ###ECL on day or ECL on prev day, but +1 because AWAP is rain for prev day
now2[I,i]=1
}

AWAPrain2=AWAPrain[AWAPrain[,2]>=1990 & AWAPrain[,2]<=2009,]
mm=floor(AWAPrain2[,3]/100)%%100
I=which(AWAPrain2[,2]<=2006)
apply(now2[I,],2,sum)/17

for(i in 1:4) print(sum(AWAPrain2[,4]*now2[,i])/sum(AWAPrain2[,4]))

monthly=array(0,c(12,4,2))
for(i in 1:12)
  for(j in 1:4)
  {
    if(j==1) I=which(AWAPrain2[,2]<=2006 & mm==i) else I=which(mm==i)
    monthly[i,j,1]=sum(now2[I,j])
    if(length(I)>0) monthly[i,j,2]=sum(AWAPrain2[I,4]*now2[I,j])/sum(AWAPrain2[I,4])
  }

clist=c("black","blue","red","orange")
plot(1:12,monthly[,1,2]*100,col=1,lwd=4,type="l",xlab="Month",ylab="% ESB rain due to ECLs",ylim=range(0,monthly[,,2]*100))
for(i in 2:4) lines(1:12,monthly[,i,2]*100,col=clist[i],lwd=4)
legend("bottomleft",c("MLDB","GV","LAP","PG"),col=clist,lwd=4,ncol=4,bty="n")




