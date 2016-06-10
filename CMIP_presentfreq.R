rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

freq<-matrix(0,15,5)
rownames(freq)<-rep("aaa",15)
colnames(freq)<-c("UM 150km","UM 50km","PG 150km","PG 50km")

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    rownames(freq)[n]<-paste(cmip[i],wrf[j])
    dir=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""))
    freq[n,1]=length(data[,1])
    dir=paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""))
    freq[n,2]=length(data[,1])
    
    data=read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""))  
    freq[n,3]=length(data[,1])
    data=read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))  
    freq[n,4]=length(data[,1])
    n=n+1
  }

freq=freq/20

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

freq<-matrix(0,15,3)
rownames(freq)<-rep("aaa",15)
colnames(freq)<-c("GV","LAP","PG")

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    rownames(freq)[n]<-paste(cmip[i],wrf[j])
    data=read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""))
    data=data[data$MaxGV>=5.1,]
    freq[n,1]=length(data[,1])
    dir=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv1_9009.csv",sep=""))
    freq[n,2]=length(data[,1])
    data=read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.95.csv",sep=""))  
    freq[n,3]=length(data[,1])
    n=n+1
  }

freq=freq/20

boxplot(freq,xlab="",ylab="ECLs p.a.")
abline(h=22,col="gray",lwd=2,lty=2)
write.csv(freq,"tmp.csv")

######### Okay, re-do monthly, for NCEP 22 events
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
monthly<-array(0,dim=c(15,12,5,5))
dimnames(monthly)[[1]]<-rep("aaa",15)
dimnames(monthly)[[2]]<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(monthly)[[3]]<-c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
thresh=c(22,15,10,5,2)
dimnames(monthly)[[4]]<-paste(thresh,"p.a.")

for(t in 1:5)
{
  n=1
for(i in 1:5)
  for(j in 1:3)
  {
    if(t==1) dimnames(monthly)[[1]][n]<-paste(cmip[i],wrf[j])
    
    data=read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""))
    b=order(data$MaxGV,decreasing=T)
    thresh2=data$MaxGV[b[20*thresh[t]]]
    data=data[data$MaxGV>=thresh2,]
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      monthly[n,m,1,t]=length(I)
    }    
    
    dir=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""))
    b=order(data$CVmax,decreasing=T)
    thresh2=data$CVmax[b[20*thresh[t]]]
    data=data[data$CVmax>=thresh2,]
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      monthly[n,m,2,t]=length(I)
    }    
    
    dir=paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""))
    b=order(data$CVmax,decreasing=T)
    thresh2=data$CVmax[b[20*thresh[t]]]
    data=data[data$CVmax>=thresh2,]
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      monthly[n,m,3,t]=length(I)
    }
    
    data=read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""))  
    b=order(data$PG2,decreasing=T)
    thresh2=data$PG2[b[20*thresh[t]]]
    data=data[data$PG2>=thresh2,]
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      monthly[n,m,4,t]=length(I)
    }
    
    data=read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))  
    b=order(data$PG2,decreasing=T)
    thresh2=data$PG2[b[20*thresh[t]]]
    data=data[data$PG2>=thresh2,]
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      monthly[n,m,5,t]=length(I)
    }
    n=n+1
  }
}

bycmip=array(0,c(5,12,5))
dimnames(bycmip)[[1]]=c("ncep","echam5","csiromk3","miroc","cccma")
dimnames(bycmip)[[2]]<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(bycmip)[[3]]=paste(thresh,"p.a.")
for(i in 1:5) bycmip[i,,]<-apply(monthly[(3*(i-1)+1):(3*i),,2:5,],c(2,4),mean)

erai=list(read.csv('~/Documents/ECLs/Algorithm Comparison/events_Fei_Dowdy.csv'),read.csv('~/Documents/ECLs/Algorithm Comparison/events_UM_rad2_p100.csv'),read.csv('~/Documents/ECLs/Algorithm Comparison/events_Ale_v15CTL.csv'))
monthly2=matrix(0,12,3)
for(i in 1:3)
{
  data=erai[[i]]
  I=which(data$Date1>19900000 & data$Date1<=20100000)
  mm=floor((data$Date1[I]%%10000)/100)
  for(m in 1:12)
  {
    I=which(mm==m)
    monthly2[m,i]=length(I)
  }
}


yL=c(4,3,2,1,0.5)
for(t in 1:5)
{
  pdf(file=paste("Figures/ECLeventsALL_thresh",thresh[t],"_monthlydist_v3.pdf",sep=""),height=4,width=6)
plot(1:12,rep(NaN,12),col=1,type="l",ylim=c(0,yL[t]),xlab="Month",ylab="ECLs p.a.",main=paste("Seasonal distribution of top",thresh[t],"ECLs p.a., 1990-2009"))
n=1
for(j in c(1,2,4))
  {
    lines(1:12,apply(monthly[1:3,,j,t],2,mean)/20,col=j,lwd=2,lty=1)
    lines(1:12,apply(monthly[4:15,,j,t],2,mean)/20,col=j,lwd=2,lty=2)
    lines(1:12,monthly2[,n]/20,col=j,lwd=2,lty=3)
    n=n+1
}
legend("topright",legend=c("GV","LAP","PG"),col=c(1,2,4),lwd=2,bty="n",cex=0.75)
legend("topleft",legend=c("NCEP-WRF","CMIP-WRF","ERAI"),col="black",lwd=2,lty=c(1,2,3),bty="n",cex=0.75)
dev.off()
}
for(t in 1:5)
{
  pdf(file=paste("Figures/ECLeventsALL_thresh",thresh[t],"_monthlydist_bycmip.pdf",sep=""),height=4,width=6)
  plot(1:12,rep(NaN,12),col=1,type="l",ylim=c(0,yL[t]),xlab="Month",ylab="ECLs p.a.",main=paste("Seasonal distribution of top",thresh[t],"ECLs p.a., 1990-2009"))
  for(j in 1:5) lines(1:12,bycmip[j,,t]/20,col=j,lwd=2)
  legend("topright",legend=c("NNRP","ECHAM5","CSIROMk3.0","MIROC3.2","CCCMA3.1"),col=1:5,lwd=2,bty="n",cex=0.75,ncol=2)
  dev.off()
}


coolprop=apply(monthly[,5:10,],c(1,3),sum)/apply(monthly,c(1,3),sum)

############
##
##  Okay, we want to try to do some more rigorous comparisons between the ncep & cmip data - 
##  in terms of intensity distributions, duration, central pressures, seasonality, etc etc. 


rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

freq<-matrix(0,15,5)
Rnames<-rep("aaa",15)

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

eventsGV<-eventsUM150<-eventsUM50<-eventsPG150<-eventsPG50<-list()

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    eventsGV[[n]]<-read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""))
    eventsUM150[[n]]<-read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""))
    eventsUM50[[n]]<-read.csv(paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""))
    eventsPG150[[n]]<-read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""))
    eventsPG50[[n]]<-read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))    
    Rnames[n]<-paste(cmip[i],wrf[j])
    n=n+1
  }
all=list(eventsGV,eventsUM150,eventsUM50,eventsPG150,eventsPG50)

length=c(seq(1,17,4),100)
thresh=c(22,15,10,5,2)
lengthD<-array(0,c(5,5,15,5))
dimnames(lengthD)[[1]]=c("<1","1+","2+","3+","4+")
dimnames(lengthD)[[2]]=c("22 p.a.","15 p.a.","10 p.a.","5 p.a.","2 p.a.")
dimnames(lengthD)[[3]]=Rnames
dimnames(lengthD)[[4]]=c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
intcol=c(8,10,10,10,10)
lencol=c(3,8,8,8,8)

for(i in 1:5)
  for(j in 1:15)
  {
    data=all[[i]][[j]]
    for(k in 1:5)
    {
      b=order(data[,intcol[i]],decreasing=T)
      thresh2=data[b[20*thresh[k]],intcol[i]]
      data2=data[(data[,intcol[i]]>=thresh2),]
      for(l in 1:5)
      {
        I=which(data2[,lencol[i]]>=length[l] & data2[,lencol[i]]<length[l+1] )
        lengthD[l,k,j,i]=length(I)
      }
    }   
  }

meanL<-array(0,c(5,15,5))
dimnames(meanL)=dimnames(apply(lengthD,c(2,3,4),sum))
for(i in 1:5)
  for(j in 1:15)
  {
    data=all[[i]][[j]]
    for(k in 1:5)
    {
      b=order(data[,intcol[i]],decreasing=T)
      thresh2=data[b[20*thresh[k]],intcol[i]]
      data2=data[(data[,intcol[i]]>=thresh2),]
      meanL[k,j,i]=mean(data2[,lencol[i]])
    }   
  }

Cchange<-Cchangeabs<-meanL[,4:15,]
for(i in 1:3) 
  {
  inds=seq(i,12,3)
  for(j in 1:4){
    Cchange[,inds[j],]=100*(meanL[,3+inds[j],]-meanL[,i,])/meanL[,i,]
    Cchangeabs[,inds[j],]=(meanL[,3+inds[j],]-meanL[,i,])
  }  
}

boxplot(Cchangeabs[4,,]/4,xlab="Method",ylab="Difference in average ECL length (days)")
abline(h=0,col="red")
boxplot(Cchange[4,,],xlab="Method",ylab="% change in average ECL length")
abline(h=0,col="red")

count1=apply(lengthD[3:5,,,],c(2,3,4),sum)/20
count1a=apply(lengthD[3:5,,,],c(2,3,4),sum)/apply(lengthD,c(2,3,4),sum)

C1change<-C1changeabs<-count1[,4:15,]
for(i in 1:3) 
{
  inds=seq(i,12,3)
  for(j in 1:4){
    C1change[,inds[j],]=100*(count1[,3+inds[j],]-count1[,i,])/count1[,i,]
    C1changeabs[,inds[j],]=(count1[,3+inds[j],]-count1[,i,])
  }  
}
boxplot(C1change[1,,],xlab="Method",ylab="Change in events of at least 2 day")
abline(h=0,col="red")

##### What about playing with distributions?

a=density(eventsGV[[1]][,3]/4,from=0,to=10)
LenDens<-array(0,c(40,15,5))
dimnames(LenDens)[[2]]<-Rnames
dimnames(LenDens)[[3]]<-c("GV","UM 150km","UM 50km","PG 150km","PG 50km")

for(i in 1:15)
  for(j in 1:5)
  {
    data=all[[j]][[i]]
    b=order(data[,intcol[j]],decreasing=T)
    thresh2=data[b[20*22],intcol[j]]
    data2=data[(data[,intcol[j]]>=thresh2),]
    a=density(data2[,lencol[j]]/4,from=0,to=10,n=40)
    LenDens[,i,j]<-a$y
    }   

for(j in 1:5)
{
plot(1,NaN,xlim=c(0,4),ylim=c(0,max(LenDens[,,j])),xlab="Length (days)",ylab="Density")
for(i in 1:3) lines(a$x,LenDens[,i,j],col="blue")
for(i in 4:15) lines(a$x,LenDens[,i,j],col="red")
lines(a$x,apply(LenDens[,1:3,j],1,mean),col="blue",lwd=3)
lines(a$x,apply(LenDens[,4:15,j],1,mean),col="red",lwd=3)
}

LenDens3<-array(0,c(20,15,5))
dimnames(LenDens3)[[2]]<-Rnames
dimnames(LenDens3)[[3]]<-c("GV","UM 150km","UM 50km","PG 150km","PG 50km")

for(k in 1:20)
for(i in 1:15)
  for(j in 1:5)
  {
    data=all[[j]][[i]]
    b=order(data[,intcol[j]],decreasing=T)
    thresh2=data[b[20*22],intcol[j]]
    data2=data[(data[,intcol[j]]>=thresh2),]
    
    I=which(data2[,lencol[j]]==k)
    LenDens3[k,i,j]<-length(I)
  }  

plot(1,NaN,xlim=c(0,3),ylim=c(0,2.5),xlab="Length (days)",ylab="Density")
for(j in 2:5)
{
  lines(a$x,apply(LenDens[,1:3,j],1,mean),col=j,lwd=3)
  lines(a$x,apply(LenDens[,4:15,j],1,mean),col=j,lwd=3,lty=2)
}
legend("topright",c("UM 150km","UM 50km","PG 150km","PG 50km"),col=2:5,lwd=3,bty="n")
legend("topleft",c("NCEP-WRF","CMIP3-WRF"),col="black",lwd=3,lty=1:2,bty="n")


###### Ooh, testing radii with my method?

RDens<-array(0,c(100,15,2))
dimnames(RDens)[[2]]<-Rnames
dimnames(RDens)[[3]]<-c("UM 150km","UM 50km")

for(i in 1:15)
  for(j in 2:3)
  {
    data=all[[j]][[i]]
    b=order(data$CVmax,decreasing=T)
    thresh2=data$CVmax[b[20*22]]
    data2=data[(data$CVmax>=thresh2),]
    a=density(data2$Rmean,from=0,to=10,n=100)
    RDens[,i,j-1]<-a$y
  }   

plot(1,NaN,xlim=c(0,8),ylim=c(0,1.5),xlab="Length (days)",ylab="Density")
for(j in 1:2)
{
  lines(a$x,apply(RDens[,1:3,j],1,mean),col=j,lwd=3)
  lines(a$x,apply(RDens[,4:15,j],1,mean),col=j,lwd=3,lty=2)
}
legend("topright",c("UM 150km","UM 50km"),col=1:2,lwd=3,bty="n")
legend("topleft",c("NCEP-WRF","CMIP3-WRF"),col="black",lwd=3,lty=1:2,bty="n")

###### Ooh, just for fun, do density for warm vs cold
RDens<-array(0,c(100,15,4))
dimnames(RDens)[[2]]<-Rnames
dimnames(RDens)[[3]]<-c("150km Warm","50km Warm","150km Cool","50km Cool")

for(i in 1:15)
  for(j in 2:3)
  {
    data=all[[j]][[i]]
    b=order(data$CVmax,decreasing=T)
    thresh2=data$CVmax[b[20*22]]
    data2=data[(data$CVmax>=thresh2),]
    
    mm=floor((data2$Date1%%10000)/100)
    I=which(mm>=5 & mm<=10)
    a=density(data2$Rmean[-I],from=0,to=10,n=100)
    RDens[,i,j-1]<-a$y
    a=density(data2$Rmean[I],from=0,to=10,n=100)
    RDens[,i,j+1]<-a$y
  }
plot(1,NaN,xlim=c(0,8),ylim=c(0,1.5),xlab="Length (days)",ylab="Density")
for(j in 1:2)
{
  lines(a$x,apply(RDens[,,j],1,mean),col="red",lwd=3,lty=j)
  lines(a$x,apply(RDens[,,j+2],1,mean),col="blue",lwd=3,lty=j)
}
legend("topright",c("Warm","Cool"),col=c("red","blue"),lwd=3,bty="n")
legend("topleft",c("UM 150km","UM 50km"),col="black",lwd=3,lty=1:2,bty="n")


### Using my method, want to look at the distribution of various characteristics for NCEP v the four GCMs

RDens<-array(0,c(100,5,4))
dimnames(RDens)[[2]]<-cmip
dimnames(RDens)[[3]]<-c("UM 150km","UM 50km","UM 150km shared","UM 50km shared")
CVDens1<-CVDens2<-PDens<-LenDens<-RDens
xlist<-matrix(0,100,5)
colnames(xlist)=c("CV max","CV mean","MSLP min","Length (days)","Radius (deg.lat)")

for(i in 1:5)
  for(j in 2:3)
  {
    a=(i-1)*3+1
    data=all[[j]][[a]]
    b=order(data$CVmax,decreasing=T)
    thresh2=data$CVmax[b[20*22]]
    data2=data[(data$CVmax>=thresh2),]
    data3=data    
    for(k in 1:2)
    {
      data=all[[j]][[a+k]]
      b=order(data$CVmax,decreasing=T)
      thresh2=data$CVmax[b[20*22]]
      data2=rbind(data2,data[(data$CVmax>=thresh2),])
      data3=rbind(data3,data)
    }    
    b=order(data3$CVmax,decreasing=T)
    thresh2=data3$CVmax[b[60*22]]
    data3=data3[(data$CVmax>=thresh2),]
    
    a=density(data2$CVmax,from=0,to=5,n=100)
    CVDens1[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,1]<-a$x
    a=density(data2$CVmean,from=0,to=5,n=100)
    CVDens2[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,2]<-a$x
    a=density(data2$MSLP2,from=950,to=1050,n=100)
    PDens[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,3]<-a$x
    a=density(data2$Length2/4,from=0,to=10,n=100)
    LenDens[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,4]<-a$x
    a=density(data2$Rmean,from=0,to=10,n=100)
    RDens[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,5]<-a$x
    
    a=density(data3$CVmax,from=0,to=5,n=100)
    CVDens1[,i,j+1]<-a$y
    a=density(data3$CVmean,from=0,to=5,n=100)
    CVDens2[,i,j+1]<-a$y
    a=density(data3$MSLP2,from=950,to=1050,n=100)
    PDens[,i,j+1]<-a$y
    a=density(data3$Length2/4,from=0,to=10,n=100)
    LenDens[,i,j+1]<-a$y
    a=density(data3$Rmean,from=0,to=10,n=100)
    RDens[,i,j+1]<-a$y
  }   

plot(1,NaN,xlim=c(0,8),ylim=c(0,1.5),xlab="Mean radius",ylab="Density")
for(j in 1:5) {
  lines(xlist[,5],RDens[,j,1],col=j,lwd=3)
  lines(xlist[,5],RDens[,j,2],col=j,lwd=3,lty=2)
}
legend("topright",cmip,col=1:5,lwd=3,bty="n")
legend("topleft",c("UM 150km","UM 50km"),col="black",lwd=3,lty=1:2,bty="n")

######### What about comparing the three different WRF versions, for both CMIP & NCEP

RDens<-array(0,c(100,6,2))
dimnames(RDens)[[2]]<-c("NCEP R1","NCEP R2","NCEP R2","CMIP R1","CMIP R2","CMIP R3")
dimnames(RDens)[[3]]<-c("UM 150km","UM 50km")
CVDens1<-CVDens2<-PDens<-LenDens<-RDens
xlist<-matrix(0,100,5)
colnames(xlist)=c("CV max","CV mean","MSLP min","Length (days)","Radius (deg.lat)")
LenDens2<-array(0,c(40,6,2))

for(i in 1:3)
  for(j in 2:3)
  {
    data=all[[j]][[i]]
    b=order(data$CVmax,decreasing=T)
    thresh2=data$CVmax[b[20*22]]
    data2=data[(data$CVmax>=thresh2),]
    
    a=density(data2$CVmax,from=0,to=5,n=100)
    CVDens1[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,1]<-a$x
    a=density(data2$CVmean,from=0,to=5,n=100)
    CVDens2[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,2]<-a$x
    a=density(data2$MSLP2,from=950,to=1050,n=100)
    PDens[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,3]<-a$x
    a=density(data2$Length2/4,from=0,to=10,n=40)
    LenDens2[,i,j-1]<-a$y
    if(i==1 & j==2) xx<-a$x
    a=density(data2$Rmean,from=0,to=10,n=100)
    RDens[,i,j-1]<-a$y
    if(i==1 & j==2) xlist[,5]<-a$x
    
    data=all[[j]][[i+3]]
    b=order(data$CVmax,decreasing=T)
    thresh2=data$CVmax[b[20*22]]
    data2=data[(data$CVmax>=thresh2),] 
    for(k in c(6,9,12))
    {
      data=all[[j]][[i+k]]
      b=order(data$CVmax,decreasing=T)
      thresh2=data$CVmax[b[20*22]]
      data2=rbind(data2,data[(data$CVmax>=thresh2),])
    }    
    
    a=density(data2$CVmax,from=0,to=5,n=100)
    CVDens1[,i+3,j-1]<-a$y
    a=density(data2$CVmean,from=0,to=5,n=100)
    CVDens2[,i+3,j-1]<-a$y
    a=density(data2$MSLP2,from=950,to=1050,n=100)
    PDens[,i+3,j-1]<-a$y
    a=density(data2$Length2/4,from=0,to=10,n=40)
    LenDens2[,i+3,j-1]<-a$y
    a=density(data2$Rmean,from=0,to=10,n=100)
    RDens[,i+3,j-1]<-a$y
    
  }   

plot(1,NaN,xlim=c(0,3),ylim=c(0,0.8),xlab="Length in ESB (days)",ylab="Density")
for(j in 1:3) {
  lines(xx,LenDens2[,j,1],col=j,lwd=3)
  lines(xx,LenDens2[,j,2],col=j,lwd=3,lty=2)
}
legend("topright",wrf,col=1:3,lwd=3,bty="n")
legend("topleft",c("UM 150km","UM 50km"),col="black",lwd=3,lty=1:2,bty="n")

### Duration of 24 hours (5 fixes)

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

freq<-matrix(0,15,5)
Rnames<-rep("aaa",15)

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

eventsGV<-eventsUM150<-eventsUM50<-eventsPG150<-eventsPG50<-list()

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    eventsGV[[n]]<-read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""))
    eventsUM150[[n]]<-read.csv(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_9009.csv",sep=""))
    eventsUM50[[n]]<-read.csv(paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_9009.csv",sep=""))
    eventsPG150[[n]]<-read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.5.csv",sep=""))
    eventsPG50[[n]]<-read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg0.6.csv",sep=""))    
    Rnames[n]<-paste(cmip[i],wrf[j])
    n=n+1
  }
all=list(eventsGV,eventsUM150,eventsUM50,eventsPG150,eventsPG50)

freq<-matrix(0,15,5)
rownames(freq)<-rep("aaa",15)
colnames(freq)<-c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
freq4<-freq

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
lencol=c(3,8,8,8,8)

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    rownames(freq)[n]<-rownames(freq4)[n]<-paste(cmip[i],wrf[j])
    for(k in 1:5)
    {
      data=all[[k]][[n]]
      freq[n,k]=length(data[,1])/20
      I=which(data[,lencol[k]]>4)
      freq4[n,k]=length(I)/20
    }
    n=n+1
  }

########## For CCCMA and ECHAM5, compare CMIP to the 3 RCMs
########## All is for proj 60, rad2cv1

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

freq<-array(0,c(4,3,12))
cmip=c("ncep","echam5","cccma")
wrf=c(NA,"R1","R2","R3")
dimnames(freq)[[1]]=cmip
dimnames(freq)[[2]]=c("GCM","R1","R2","R3")

n=1
thresh=10
for(i in 1:3)
{
  dir=paste("outputUM/proj100/outputUM_",cmip[i],"_rad2cv06/",sep="")
  data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_proj100_rad2cv06_9009.csv",sep=""))
  b=order(data$CV2,decreasing=T)
  thresh2=data$CV2[b[20*thresh]]
  data=data[data$CV2>=thresh2,]
  mm=floor((data$Date1%%10000)/100)
  for(m in 1:12)
  {
    I=which(mm==m)
    freq[1,i,m]=length(I)
  } 
                  
  for(j in 2:4)
  {
        dir=paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/",sep="")
    data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_9009.csv",sep=""))
    b=order(data$CV2,decreasing=T)
    thresh2=data$CV2[b[20*thresh]]
    data=data[data$CV2>=thresh2,]
    mm=floor((data$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      freq[j,i,m]=length(I)
    } 
  }
}
freq=freq/20

cool=apply(freq[,,5:10],c(1,2),sum)/apply(freq,c(1,2),sum)
dimnames(cool)[[2]]=cmip
dimnames(cool)[[1]]=c("GCM","R1","R2","R3")

### Testing NCEP different projs

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
monthly<-array(0,dim=c(15,12,4,5))
dimnames(monthly)[[1]]<-rep("aaa",15)
dimnames(monthly)[[2]]<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(monthly)[[3]]<-c("150km, r200km","50km, r200km","150km, r500km","50km, r500km")
cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
thresh=c(22,15,10,5,2)
dimnames(monthly)[[4]]<-paste(thresh,"p.a.")
future<-monthly

plist=c(100,240,100,240)
ilist=c("rad2cv06","rad2cv1","rad5cv015","rad5cv04")

for(t in 1:5)
  for(k in 1:4)
{
  n=1
  for(i in 1:5)
    for(j in 1:3)
    {
      if(t==1) dimnames(monthly)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j])
      
      dir=paste("outputUM/proj",plist[k],"/outputUM_",cmip[i],"_WRF",wrf[j],"_50_",ilist[k],"/",sep="")
      data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj",plist[k],"_",ilist[k],"_9009.csv",sep=""))
      b=order(data$CV2,decreasing=T)
      thresh2=data$CV2[b[20*thresh[t]]]
      data=data[data$CV2>=thresh2,]
      mm=floor((data$Date1%%10000)/100)
      for(m in 1:12)
      {
        I=which(mm==m)
        monthly[n,m,k,t]=length(I)/20
      }    
      
      if(i>1)
      {
        data=read.csv(paste(dir,"ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj",plist[k],"_",ilist[k],"_6079.csv",sep=""))
        data=data[data$CV2>=thresh2,]
        mm=floor((data$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          future[n,m,k,t]=length(I)/20
        }     
      }
      n=n+1
    }
}

annualC<-(apply(future,c(1,3,4),sum)/apply(monthly,c(1,3,4),sum))-1
coolC<-(apply(future[,5:10,,],c(1,3,4),sum)/apply(monthly[,5:10,,],c(1,3,4),sum))-1
warmC<-(apply(future[,c(1:4,11:12),,],c(1,3,4),sum)/apply(monthly[,c(1:4,11:12),,],c(1,3,4),sum))-1
coolC[1:3,,]<-warmC[1:3,,]<-NaN

boxplot(coolC[,,1]*100,ylab="Cool season change (%)",ylim=c(-50,50),main="22 events p.a.")
abline(h=0,col="red")
boxplot(warmC[,,1]*100,ylab="Warm season change (%)",ylim=c(-50,50),main="22 events p.a.")
abline(h=0,col="red")

mmean<-apply(monthly[4:15,,,],c(2,3,4),mean)
mmean2<-apply(monthly[1:3,,,],c(2,4),mean)
yrange1=ceiling(apply(mmean,3,max)*2)/2
for(t in 1:5)
{
  pdf(file=paste("ECLevents_thresh",thresh[t],"_monthlyLAP_radius.pdf",sep=""),height=4,width=6)
  plot(1:12,mmean2[,t],col=1,type="l",lwd=2,xlim=c(0,12),ylim=c(0,yrange1[t]),xlab="Month",ylab="Events p.a.",main=paste(thresh[t],"ECLs p.a."))
  for(n in 1:4) lines(1:12,mmean[,n,t],col=n+1,lwd=2)
  legend("bottomleft",c("150km, r200km","50km, r200km","150km, r500km","50km, r500km","NCEP"),lwd=2,col=c(2:5,1),ncol=3,bty="n")
  dev.off()
}

yrange=ceiling(apply(monthly,3,max)*2)/2
clist=c(rep(0,3),rep(1,3),rep(2,3),rep(3,3))


apply(monthly[10:12,5:10,],3,sum)/apply(monthly[10:12,,],3,sum)

####### Similar for Ale's method
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
monthly<-array(0,dim=c(15,12,4,5))
dimnames(monthly)[[1]]<-rep("aaa",15)
dimnames(monthly)[[2]]<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(monthly)[[3]]<-c("150km, r200km","50km, r200km","150km, r500km","50km, r500km")
cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
thresh=c(22,15,10,5,2)
dimnames(monthly)[[4]]<-paste(thresh,"p.a.")
future<-monthly

n=1
    for(i in 1:5)
      for(j in 1:3)
      {
        if(i==1 & j==1) dimnames(monthly)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j])
        eventsPG50<-read.csv(paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_PG2.csv",sep=""))    
        eventsPG150<-read.csv(paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_PG2.csv",sep="")) 
        for(t in 1:5)
        {
          data=list()
          b=order(eventsPG150$PG200,decreasing=T)
          thresh2=eventsPG150$PG200[b[20*thresh[t]]]
          data[[1]]=eventsPG150[eventsPG150$PG200>=thresh2,]
          b=order(eventsPG50$PG200,decreasing=T)
          thresh2=eventsPG50$PG200[b[20*thresh[t]]]
          data[[2]]=eventsPG50[eventsPG50$PG200>=thresh2,]
          b=order(eventsPG150$PG500,decreasing=T)
          thresh2=eventsPG150$PG500[b[20*thresh[t]]]
          data[[3]]=eventsPG150[eventsPG150$PG500>=thresh2,]
          b=order(eventsPG50$PG500,decreasing=T)
          thresh2=eventsPG50$PG500[b[20*thresh[t]]]
          data[[4]]=eventsPG50[eventsPG50$PG500>=thresh2,]
          
          for(k in 1:4)
          {
          mm=floor((data[[k]]$Date1%%10000)/100)
          for(m in 1:12)
          {
            I=which(mm==m)
            monthly[n,m,k,t]=length(I)/20
          }    
          }
        }
        n=n+1
      }
  }

mmean<-apply(monthly[4:15,,,],c(2,3,4),mean)
mmean2<-apply(monthly[1:3,,,],c(2,4),mean)
yrange1=ceiling(apply(mmean,3,max)*2)/2
for(t in 1:5)
{
  pdf(file=paste("ECLevents_thresh",thresh[t],"_monthlyPG_radius.pdf",sep=""),height=4,width=6)
  plot(1:12,mmean2[,t],col=1,type="l",lwd=2,xlim=c(0,12),ylim=c(0,yrange1[t]),xlab="Month",ylab="Events p.a.",main=paste(thresh[t],"ECLs p.a."))
  for(n in 1:4) lines(1:12,mmean[,n,t],col=n+1,lwd=2)
  legend("bottomleft",c("150km, r200km","50km, r200km","150km, r500km","50km, r500km","NCEP"),lwd=2,col=c(2:5,1),ncol=3,bty="n")
  dev.off()
}

yrange=ceiling(apply(monthly,3,max)*2)/2
clist=c(rep(0,3),rep(1,3),rep(2,3),rep(3,3))
apply(monthly[10:12,5:10,],3,sum)/apply(monthly[10:12,,],3,sum)


################################################
#### Assessing only the 50km ones, cool season proportion as a function of intensity threshold

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

freq<-matrix(0,15,5)
Rnames<-rep("aaa",15)

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

eventsUM50<-eventsPG50<-list()

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    eventsUM50[[n]]<-read.csv(paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_9009.csv",sep=""))
    eventsPG50[[n]]<-read.csv(paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg0.95.csv",sep=""))    
    Rnames[n]<-paste(cmip[i],wrf[j])
    n=n+1
  }

freq<-matrix(0,15,2)
for(i in 1:15)
{
  freq[i,1]=length(eventsUM50[[i]][,1])
  freq[i,2]=length(eventsPG50[[i]][,1])
}

thresh=c(NaN,50,40,30,20,10,NaN)
count<-array(NaN,c(15,6,2))
dimnames(count)[[1]]<-Rnames
dimnames(count)[[2]]<-thresh[1,6]
dimnames(count)[[3]]<-c("UM","PG")
cool<-count

for(i in 1:15)
{
  a=eventsUM50[[i]]
  a=a[order(a[,7],decreasing=T),]
  m1=(floor(a$Date1/100)%%100)
  
  b=eventsPG50[[i]]
  b=b[order(b[,7],decreasing=T),]
  m2=(floor(b$Date1/100)%%100)
  
  for(t in 1:6)
  {
    if(is.na(thresh[t])) thresh1=a[length(a[,1]),7] else thresh1=a[20*thresh[t],7]
    if(is.na(thresh[t+1])) thresh2=a[1,7] else thresh2=a[20*thresh[t+1],7]
    I=which(a[,7]>=thresh1 & a[,7]<thresh2)
    count[i,t,1]=length(I)
    I=which(a[,7]>=thresh1 & a[,7]<thresh2 & m1>=5 & m1<=10)
    cool[i,t,1]=length(I)
    
    if(is.na(thresh[t])) thresh1=b[length(a[,1]),7] else thresh1=b[20*thresh[t],7]
    if(is.na(thresh[t+1])) thresh2=b[1,7] else thresh2=b[20*thresh[t+1],7]
    I=which(b[,7]>=thresh1 & b[,7]<thresh2)
    count[i,t,2]=length(I)
    I=which(b[,7]>=thresh1 & b[,7]<thresh2 & m2>=5 & m2<=10)
    cool[i,t,2]=length(I)
  }
}

prop=(cool/count)*100
dimnames(prop)[[1]]<-Rnames
dimnames(prop)[[2]]<-c(">50","40-50","30-40","20-30","10-20","<10")
dimnames(prop)[[3]]<-c("UM","PG")

boxplot(prop[4:15,,1],xlab="Threshold - events p.a.",ylab="% during May-October",main="LAP method (50km) - seasonality",ylim=c(25,75))
abline(h=50,col="red")
points(1:6,apply(prop[1:3,,1],2,mean),pch=4,cex=2,lwd=3)
legend("topleft",legend=c("CMIP-WRF","NCL-WRF"),pch=c(0,4),pt.cex=2,pt.lwd=3)
boxplot(prop[,,2],xlab="Threshold - events p.a.",ylab="% during May-October",main="PG method (50km) - seasonality",ylim=c(25,75))
abline(h=50,col="red")
points(1:6,apply(prop[1:3,,2],2,mean),pch=4,cex=2,lwd=3)
legend("topleft",legend=c("CMIP-WRF","NCL-WRF"),pch=c(0,4),pt.cex=2,pt.lwd=3)

###
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
monthly<-array(0,dim=c(15,12))
dimnames(monthly)[[1]]<-rep("aaa",15)
dimnames(monthly)[[2]]<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")


n=1
  for(i in 1:5)
    for(j in 1:3)
    {
      data=read.csv(paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_PG2.csv",sep=""))  
      mm=floor((data$Date1%%10000)/100)
      for(m in 1:12)
      {
        I=which(mm==m)
        monthly[n,m]=length(I)
      }
      n=n+1
    }
monthly2=apply(monthly,2,range)/20

names1=c("ERAI","MERRA","JRA55")
monthly3=matrix(0,3,12)
for(j in 1:3)
{
  data=read.csv(paste("Alejandro/ECLevents_Alejandro_",names1[j],"_res50_9009_PG2.csv",sep=""))
  mm=floor((data$Date1%%10000)/100)
  for(m in 1:12)
  {
    I=which(mm==m & data$Date1>=19900000 & data$Date1<=20100000)
    monthly3[j,m]=length(I)/20
  }
}
clist=c("blue","red","black")


plot(NA,xlim=c(1,12),ylim=c(0,20),xlab="Month",ylab="Average number of ECLs")
polygon(c(1:12,12:1), c(monthly2[1,],rev(monthly2[2,])), col = 'grey80', border = NA)
for(j in 1:3) lines(1:12,monthly3[j,],col=clist[j],lwd=2)
legend("topleft",legend=names1,lwd=2,col=clist,bty="n")




