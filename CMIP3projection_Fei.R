##############33

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Fei/")

cmip=c("ECHAM5","MK30","MIROC","CCCMA")
cmipUM=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

now<-future<-nowE<-futureE<-list()
now2<-future2<-matrix(0,12,12)
rownames(now2)<-rownames(future2)<-rep("aaa",12)

n=1
for(i in 1:4)
  for(j in 1:3)
  {
    rownames(now2)[n]<-rownames(future2)[n]<-paste(cmip[i],wrf[j])
    now[[n]]=read.csv(paste("vor_",cmip[i],"_",wrf[j],"_1990-2010_Andrew.csv",sep=""))
    future[[n]]=read.csv(paste("vor_",cmip[i],"_",wrf[j],"_2060-2080_Andrew.csv",sep=""))
    n=n+1
  }

n=1
for(k in 1:4)
  for(j in 1:3)
{
  a=now[[n]]
  x<-rle(a$ID)
  events<-data.frame(ID=x$values,Length=x$lengths,
                     Date1=rep(0,length(x$values)),Year1=rep(0,length(x$values)),Month1=rep(0,length(x$values)),
                     Date2=rep(0,length(x$values)),MaxGV=rep(0,length(x$values)))
  
  for(i in 1:length(events$ID)) 
  {
    I<-which(a$ID==events$ID[i])
    events$Date1[i]=min(a$Date[I]) 
    events$Year1[i]=floor(events$Date1[i]/10000)
    events$Month1[i]=floor((events$Date1[i] %% 10000)/100)
    events$Date2[i]=max(a$Date[I])
    events$MaxGV[i]=max(a$vor[I])
  }
  write.csv(events,file=paste("GVevents_",cmipUM[k],"_",wrf[j],"_9009.csv",sep=""))
  nowE[[n]]=events
  
  a=future[[n]]
  x<-rle(a$ID)
  events<-data.frame(ID=x$values,Length=x$lengths,
                     Date1=rep(0,length(x$values)),Year1=rep(0,length(x$values)),Month1=rep(0,length(x$values)),
                     Date2=rep(0,length(x$values)),MaxGV=rep(0,length(x$values)))
  
  for(i in 1:length(events$ID)) 
  {
    I<-which(a$ID==events$ID[i])
    events$Date1[i]=min(a$Date[I]) 
    events$Year1[i]=floor(events$Date1[i]/10000)
    events$Month1[i]=floor((events$Date1[i] %% 10000)/100)
    events$Date2[i]=max(a$Date[I])
    events$MaxGV[i]=max(a$vor[I])
  }
  futureE[[n]]=events
  write.csv(events,file=paste("GVevents_",cmipUM[k],"_",wrf[j],"_6079.csv",sep=""))
  
  n=n+1
}

for(n in 1:12)
{
  
  a=nowE[[n]]$Month1
  b=futureE[[n]]$Month1
  
  for(m in 1:12)
  {
    I=which(a==m)
    now2[n,m]=length(I)
    I=which(b==m)
    future2[n,m]=length(I)
  }
}

now3<-future3<-now2
for(n in 1:12)
{
  a=nowE[[n]]$Month1[nowE[[n]]$MaxGV>=5.1]
  b=futureE[[n]]$Month1[nowE[[n]]$MaxGV>=5.1]
  
  for(m in 1:12)
  {
    I=which(a==m)
    now3[n,m]=length(I)
    I=which(b==m)
    future3[n,m]=length(I)
  }
}


plot(seq(1:12),now2[1,],ylim=c(0,90),xlab="Month",ylab="Number of ECLs",type="l")
for(i in 2:12) lines(seq(1:12),now2[i,])
colnames(now2)=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
boxplot(now2/20,main="Number of ECLs p.a.")

change=((future2/now2)-1)*100
colnames(change)=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
boxplot(change,main="% Change in ECL numbers between 1990-2009 and 2060-2079")
abline(h=0,col="red")

bymodel=matrix(0,4,12)
rownames(bymodel)=cmip
for(i in 1:4) bymodel[i,]=apply(change[((i-1)*3+1):(i*3),],2,mean)

bywrf=matrix(0,3,12)
rownames(bywrf)=wrf
for(i in 1:3) bywrf[i,]=apply(change[seq(i,12,3),],2,mean)

colnames(bywrf)<-colnames(bymodel)<-colnames(change)

plot(seq(1:12),bymodel[1,],col=1,type="p",pch=19,ylim=c(-60,110),
     xlab="Month",ylab="",main="% Change in ECL numbers between 1990-2009 and 2060-2079 - by model")
for(i in 2:4) points(seq(1:12),bymodel[i,],col=i,pch=19)
abline(h=0,col="red")
legend("topright",pch=19,col=1:4,legend=cmip)

plot(seq(1:12),bywrf[1,],col=1,type="p",pch=19,ylim=c(-60,110),
     xlab="Month",ylab="",main="% Change in ECL numbers between 1990-2009 and 2060-2079 - by WRF version")
for(i in 2:3) points(seq(1:12),bywrf[i,],col=i,pch=19)
abline(h=0,col="red")
legend("topright",pch=19,col=1:3,legend=wrf)


##########Now, restrict to the event CV threshold that gives ~ 5 or ~ 2 p.a. in historical

now4<-future4<-array(0,c(12,12,5))
dimnames(now4)[[1]]<-dimnames(future4)[[1]]<-rownames(now2)

thresh=c(20,10,5,2,1)

for(t in 1:5)
 for(n in 1:12)
  {
    data=nowE[[n]]
    b=order(data$MaxGV,decreasing=T)
    thresh2=data$MaxGV[b[20*thresh[t]]]
    data2=data[data$MaxGV>=thresh2,]    
    for(m in 1:12)
    {
      I=which(data2$Month1==m)
      now4[n,m,t]=length(I)
    }
    
    data=futureE[[n]]
    data2=data[data$MaxGV>=thresh2,]    
    for(m in 1:12)
    {
      I=which(data2$Month1==m)
      future4[n,m,t]=length(I)
    }
}

change=((future4/now4)-1)*100
dimnames(change)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

for(i in 1:5)
{
boxplot(change[,,i],main=paste("% Change in ECLs for thresh of",thresh[i],"events p.a. between 1990-2009 and 2060-2079"),ylim=c(-200,200))
abline(h=0,col="red")
}

nowY=apply(now4,c(1,3),sum)
futureY=apply(future4,c(1,3),sum)
changeY=((futureY/nowY)-1)*100
dimnames(changeY)[[2]]=thresh
boxplot(changeY,main="% Change in ECL events p.a. between 1990-2009 and 2060-2079",cex.main=1,ylim=c(-50,50),xlab="Threshold: ECLs p.a.")
abline(h=0,col="red")

###Oooh, try to make fun density plots?

for(thresh in c(10,5,2))
{
n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      data=nowE[[n]]
      b=order(data$MaxGV,decreasing=T)
      thresh2=data$MaxGV[b[20*thresh]]  
      data1=data[data$MaxGV>=thresh2,]
      if(n==1) now5=data1 else now=rbind(now5,data1)   
      data=futureE[[n]]
      b=order(data$MaxGV,decreasing=T)
      thresh2=data$MaxGV[b[20*thresh]]  
      data2=data[data$MaxGV>=thresh2,]
      if(n==1) future5=data2 else now=rbind(future5,data2)         
      
      a=density(data1$MaxGV,na.rm=T)
      b=density(data2$MaxGV,na.rm=T)
      
      lims=c(0,10)
      
      pdf(file=paste("~/Documents/ECLs/WRFCMIP/Fei/ECL_GVchange_pdf_",cmip[i],"_",wrf[j],"_thresh",thresh,".pdf",sep=""),width=5,height=4)

      plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
           xlab="Intensity",ylab="Frequency",cex.main=1.2,
           main="ECL intensity distribution")
      polygon(a,col=rgb(0,0,1,1/4),density=-1)
      polygon(b,col=rgb(1,0,0,1/4),density=-1)
      legend("topright",legend=c("1990-2009","2060-2079"),
             col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
      dev.off()
      
      n=n+1
    }

a=density(now5$MaxGV,na.rm=T)
b=density(future5$MaxGV,na.rm=T)

pdf(file=paste("~/Documents/ECLs/WRFCMIP/Fei/ECL_GVchange_pdf_thresh",thresh,".pdf",sep=""),width=5,height=4)
plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
     xlab="Intensity",ylab="Frequency",cex.main=1.2,
     main="ECL intensity distribution")
polygon(a,col=rgb(0,0,1,1/4),density=-1)
polygon(b,col=rgb(1,0,0,1/4),density=-1)
legend("topright",legend=c("1990-2009","2060-2079"),
       col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
dev.off()
}


###### Distribution of GV


rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

now<-array(0,c(12,12,5))
dimnames(now)[[1]]<-rep("aaa",12)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
future<-now
intcol=c(8,10,10,10,10)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
n=1
  for(i in 1:4)
    for(j in 1:3)
    {      
      if(n==1) data1=read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep="")) else data1=rbind(data1,read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep="")))
      if(n==1) data2=read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep="")) else data2=rbind(data2,read.csv(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep="")))
    n=n+1
    }            

range(data1$MaxGV,data2$MaxGV)

bk=quantile(data1$MaxGV,seq(0,1,0.1))
bk[1]=0
bk[11]=100

a=hist(data1$MaxGV,breaks=bk,plot=F)
b=hist(data2$MaxGV,breaks=bk,plot=F)
barplot(((b$counts-a$counts)/12),xlab="Event max GV",ylab="Average (absolute) change in occurrence (all runs)",names.arg=paste("Decile",seq(1,10)))

## Aaand split by season

mm1=floor(data1$Date1/10000)%%100
mm2=floor(data1$Date1/10000)%%100
I=which(mm1>=5 & mm1<=10)
J=which(mm2>=5 & mm2<=10)

bk=quantile(data1$MaxGV[-I],seq(0,1,0.1),na.rm=T)
bk[1]=0
bk[11]=100

a=hist(data1$MaxGV[-I],breaks=bk,plot=F)
b=hist(data2$MaxGV[-J],breaks=bk,plot=F)
barplot(((b$counts-a$counts)/12),xlab="Event max GV",ylab="Average (absolute) change in occurrence (all runs)",names.arg=paste("Decile",seq(1,10)))

