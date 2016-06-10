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
thresh=22
for(thresh in c(2,5,10,22))
{
  n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j])
      
      filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
      filelist2=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
      
      for(ff in 1:5)
      {
        data=read.csv(filelist1[ff])
        b=order(data[,intcol[ff]],decreasing=T)
        thresh2=data[b[20*thresh],intcol[ff]]
        data2=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          now[n,m,ff]=length(I)
        }
        
        data=read.csv(filelist2[ff])
        data2=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          future[n,m,ff]=length(I)
        }        
      }
      
      n=n+1  
    }
  
  change=((future/now)-1)*100
  change2=change[,,1]
  for(i in 2:5) change2=rbind(change2,change[,,i])
  
  tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_bymonth_boxplot.tiff",sep=""),width=600,height=400)
  boxplot(change2,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
  abline(h=0,col="red",lwd=2)
  dev.off()
  now2=now[,,1]
  for(i in 2:5) now2=rbind(now2,now[,,i])
  now3=data.frame(MAM=apply(now2[,3:5],1,sum),JJA=apply(now2[,6:8],1,sum),
                  SON=apply(now2[,9:11],1,sum),DJF=apply(now2[,c(1,2,12)],1,sum))
  
  future2=future[,,1]
  for(i in 2:5) future2=rbind(future2,future[,,i])
  future3=data.frame(MAM=apply(future2[,3:5],1,sum),JJA=apply(future2[,6:8],1,sum),
                     SON=apply(future2[,9:11],1,sum),DJF=apply(future2[,c(1,2,12)],1,sum))
  change3=((future3/now3)-1)*100
  
  tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_byseason_boxplot.tiff",sep=""),width=400,height=400)
  boxplot(change3,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
  abline(h=0,col="red",lwd=2)
  dev.off()
  
  now4=apply(now,c(1,3),sum)
  future4=apply(future,c(1,3),sum)
  change4=((future4/now4)-1)*100
  tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_bymethod_boxplot.tiff",sep=""),width=500,height=400)
  boxplot(change4,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
  abline(h=0,col="red",lwd=2)
  dev.off()
  
  now4=apply(now[,c(1:4,11:12),],c(1,3),sum)
  future4=apply(future[,c(1:4,11:12),],c(1,3),sum)
  change4=((future4/now4)-1)*100
  tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_bymethod_boxplot_warm.tiff",sep=""),width=500,height=400)
  boxplot(change4,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
  abline(h=0,col="red",lwd=2)
  dev.off()
  
  now4=apply(now[,5:10,],c(1,3),sum)
  future4=apply(future[,5:10,],c(1,3),sum)
  change4=((future4/now4)-1)*100
  tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_bymethod_boxplot_cool.tiff",sep=""),width=500,height=400)
  boxplot(change4,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
  abline(h=0,col="red",lwd=2)
  dev.off()
  
  bysource=matrix(NaN,4,3)
  colnames(bysource)=c("GCM","RCM","Method")
  bycmip=matrix(NaN,9,4)
  colnames(bycmip)=cmip
  bywrf<-bymethod<-matrix(NaN,12,3)
  colnames(bywrf)=wrf
  colnames(bymethod)=c("GV","UM 150km","PG 150km")
  
  type=c("","_warm","_cool")
  mlist=list(seq(1,12),c(1:4,11:12),seq(5,10))
  
  for(ss in 1:3)
  {
    now4=apply(now[,mlist[[ss]],],c(1,3),sum)
    future4=apply(future[,mlist[[ss]],],c(1,3),sum)
    change4=((future4/now4)-1)*100
    
    tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_bymethod_boxplot",type[ss],".tiff",sep=""),width=500,height=400)
    boxplot(change4,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
    abline(h=0,col="red",lwd=2)
    dev.off()
    
    change4=change4[,c(1,2,4)]
    
    for(i in 1:4)
      for(j in 1:3)
        for(k in 1:3)
          bycmip[3*(k-1)+j,i]<-bywrf[3*(i-1)+k,j]<-bymethod[3*(i-1)+j,k]<-change4[3*(i-1)+j,k]
    
    for(i in 1:4)
    {
      a=change4[seq(3*(i-1)+1,3*i),]
      bysource[i,1]=mean(a)
    }
    for(j in 1:3)
    {
      a=change4[seq(j,12,3),]
      bysource[j,2]=mean(a)
    }
    for(k in 1:3)
    {
      a=change4[,k]
      bysource[k,3]=mean(a)
    }
    
    tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_boxplot_bysource",type[ss],".tiff",sep=""),width=500,height=400)
    boxplot(bysource,main="Range of trends by source of uncertainty",ylab="% change in ECLs",ylim=c(-30,30))
    abline(h=0,col="red",lwd=2)
    dev.off()
    tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_boxplot_byGCM",type[ss],".tiff",sep=""),width=500,height=400)
    boxplot(bycmip,main="Range of trends by GCM",ylab="% change in ECLs")
    abline(h=0,col="red",lwd=2)
    dev.off()
    tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_boxplot_byRCM",type[ss],".tiff",sep=""),width=500,height=400)
    boxplot(bywrf,main="Range of trends by RCM",ylab="% change in ECLs")
    abline(h=0,col="red",lwd=2)
    dev.off()
    tiff(file=paste("Figures/ECLeventsALL_thresh",thresh,"_change_boxplot_bymethod",type[ss],".tiff",sep=""),width=500,height=400)
    boxplot(bymethod,main="Range of trends by ECL detection method",ylab="% change in ECLs")
    abline(h=0,col="red",lwd=2)
    dev.off()
  }
  
}



####### Statistical significance - need number of events by year!


rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
intcol=c(8,10,10,10,10)
thresh=c(2,5,10,15,22)

changeC<-array(NaN,dim=c(12,5,5))
dimnames(changeC)[[1]]<-rep("aaa",12)
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(changeC)[[1]][n]<-paste(cmip[i],wrf[j])
    n=n+1
  }
dimnames(changeC)[[2]]=c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
dimnames(changeC)[[3]]=c(2,5,10,15,22)
changeW<-changeCsig<-changeWsig<-changeC

sigens<-matrix(NaN,5,2)
colnames(sigens)=c("Cool","Warm")
rownames(sigens)=c(2,5,10,15,22)

for(t in 1:5)
{
  n=1
  nowC<-nowW<-futureC<-futureW<-array(NaN,c(12,20,5))
  for(i in 1:4)
    for(j in 1:3)
    { 
      filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
      filelist2=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
      
      for(ff in 1:5)
      {
        data=read.csv(filelist1[ff])
        b=order(data[,intcol[ff]],decreasing=T)
        thresh2=data[b[20*thresh[t]],intcol[ff]]
        data2=data[data[,intcol[ff]]>=thresh2,]
        
        yy=floor(data2$Date1/10000)
        mm=floor((data2$Date1%%10000)/100)
        yy2=yy
        yy2[mm<=4]=yy2[mm<=4]-1
        years=unique(yy)
        for(y in 1:20)
        {
          I=which(mm>=5 & mm<=10 & yy==years[y])
          nowC[n,y,ff]=length(I)
          I=which((mm>=11 | mm<=4) & yy2==years[y])
          if(y<20) nowW[n,y,ff]=length(I)
        }
        
        data=read.csv(filelist2[ff])
        data2=data[data[,intcol[ff]]>=thresh2,]
        yy=floor(data2$Date1/10000)
        mm=floor((data2$Date1%%10000)/100)
        yy2=yy
        yy2[mm<=4]=yy2[mm<=4]-1
        years=unique(yy)
        for(y in 1:20)
        {
          I=which(mm>=5 & mm<=10 & yy==years[y])
          futureC[n,y,ff]=length(I)
          I=which((mm>=11 | mm<=4) & yy2==years[y])
          if(y<20) futureW[n,y,ff]=length(I)
        }        
      }
      
      n=n+1  
    }
  changeC[,,t]=((apply(futureC,c(1,3),mean)/apply(nowC,c(1,3),mean))-1)*100
  changeW[,,t]=((apply(futureW,c(1,3),mean,na.rm=T)/apply(nowW,c(1,3),mean,na.rm=T))-1)*100
  for(i in 1:12)
    for(j in 1:5)
    {
      a=t.test(nowC[i,,j],futureC[i,,j])
      changeCsig[i,j,t]=a$p.value
      a=t.test(nowW[i,,j],futureW[i,,j])
      changeWsig[i,j,t]=a$p.value
    }
  
  a=t.test(apply(futureC,c(1,3),mean),apply(nowC,c(1,3),mean))
  sigens[t,1]=a$p.value
  a=t.test(apply(futureW,c(1,3),mean,na.rm=T),apply(nowW,c(1,3),mean,na.rm=T))
  sigens[t,2]=a$p.value
}


enschangeC<-matrix(NaN,5,5)
rownames(enschangeC)<-thresh
colnames(enschangeC)<-c("Min","Mean","Max","Std","Decline %")
enschangeW<-enschangeC

enschangeC[,1]=apply(changeC,3,min,na.rm=T)
enschangeC[,2]=apply(changeC,3,mean,na.rm=T)
enschangeC[,3]=apply(changeC,3,max,na.rm=T)
enschangeW[,1]=apply(changeW,3,min,na.rm=T)
enschangeW[,2]=apply(changeW,3,mean,na.rm=T)
enschangeW[,3]=apply(changeW,3,max,na.rm=T)

for(i in 1:5)
{
  enschangeC[i,4]=sd(as.vector(changeC[,,i]))
  I=which(changeC[,,i]<0)
  enschangeC[i,5]=length(I)/60
  enschangeW[i,4]=sd(as.vector(changeW[,,i]))
  I=which(changeW[,,i]<0)
  enschangeW[i,5]=length(I)/60
}

######### Change in ECLs of various intensities - version two, for subcats 20-15,15-10,10-5,5-0

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
intcol=c(8,10,10,10,10)
thresh=c(20,15,10,5,0)

changeC<-array(NaN,dim=c(12,5,4))
dimnames(changeC)[[1]]<-rep("aaa",12)
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(changeC)[[1]][n]<-paste(cmip[i],wrf[j])
    n=n+1
  }
dimnames(changeC)[[2]]=c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
dimnames(changeC)[[3]]=c(20,15,10,5)
changeW<-changeCsig<-changeWsig<-changeC

sigens<-matrix(NaN,4,2)
colnames(sigens)=c("Cool","Warm")
rownames(sigens)=c("20-15","15-10","10-5","<5")

##What about for different duration thresholds?

length=4
  for(t in 1:4)
  {
    n=1
    nowC<-nowW<-futureC<-futureW<-array(NaN,c(12,20,5))
    for(i in 1:4)
      for(j in 1:3)
      { 
        filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                    paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                    paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                    paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                    paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
        filelist2=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""),
                    paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                    paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                    paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                    paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
        
        for(ff in 1:5)
        {
          data=read.csv(filelist1[ff])
#          data=data[data$Length>=length,]
          b=order(data[,intcol[ff]],decreasing=T)
          thresh1=data[b[20*thresh[t]],intcol[ff]]
          if(t<4) thresh2=data[b[20*thresh[t+1]],intcol[ff]] else thresh2=9999999
          data2=data[(data[,intcol[ff]]>=thresh1 & data[,intcol[ff]]<thresh2),]
          
          yy=floor(data2$Date1/10000)
          mm=floor((data2$Date1%%10000)/100)
          yy2=yy
          yy2[mm<=4]=yy2[mm<=4]-1
          years=unique(yy)
          for(y in 1:20)
          {
            I=which(mm>=5 & mm<=10 & yy==years[y])
            nowC[n,y,ff]=length(I)
            I=which((mm>=11 | mm<=4) & yy2==years[y])
            if(y<20) nowW[n,y,ff]=length(I)
          }
          
          data=read.csv(filelist2[ff])
#          data=data[data$Length>=length,]
          data2=data[(data[,intcol[ff]]>=thresh1 & data[,intcol[ff]]<thresh2),]
          yy=floor(data2$Date1/10000)
          mm=floor((data2$Date1%%10000)/100)
          yy2=yy
          yy2[mm<=4]=yy2[mm<=4]-1
          years=unique(yy)
          for(y in 1:20)
          {
            I=which(mm>=5 & mm<=10 & yy==years[y])
            futureC[n,y,ff]=length(I)
            I=which((mm>=11 | mm<=4) & yy2==years[y])
            if(y<20) futureW[n,y,ff]=length(I)
          }        
        }
        
        n=n+1  
      }
    changeC[,,t]=((apply(futureC,c(1,3),mean)/apply(nowC,c(1,3),mean))-1)*100
    changeW[,,t]=((apply(futureW,c(1,3),mean,na.rm=T)/apply(nowW,c(1,3),mean,na.rm=T))-1)*100
    for(i in 1:12)
      for(j in 1:5)
      {
        a=t.test(nowC[i,,j],futureC[i,,j])
        changeCsig[i,j,t]=a$p.value
        a=t.test(nowW[i,,j],futureW[i,,j])
        changeWsig[i,j,t]=a$p.value
      }
    
    a=t.test(apply(futureC,c(1,3),mean),apply(nowC,c(1,3),mean))
    sigens[t,1]=a$p.value
    a=t.test(apply(futureW,c(1,3),mean,na.rm=T),apply(nowW,c(1,3),mean,na.rm=T))
    sigens[t,2]=a$p.value
  }
  
  enschangeC<-matrix(NaN,4,5)
  rownames(enschangeC)<-c("20-15","15-10","10-5","<5")
  colnames(enschangeC)<-c("Min","Mean","Max","Std","Decline %")
  enschangeW<-enschangeC
  
  enschangeC[,1]=apply(changeC,3,min,na.rm=T)
  enschangeC[,2]=apply(changeC,3,mean,na.rm=T)
  enschangeC[,3]=apply(changeC,3,max,na.rm=T)
  enschangeW[,1]=apply(changeW,3,min,na.rm=T)
  enschangeW[,2]=apply(changeW,3,mean,na.rm=T)
  enschangeW[,3]=apply(changeW,3,max,na.rm=T)
  
  for(i in 1:4)
  {
    enschangeC[i,4]=sd(as.vector(changeC[,,i]))
    I=which(changeC[,,i]<0)
    enschangeC[i,5]=length(I)/60
    enschangeW[i,4]=sd(as.vector(changeW[,,i]))
    I=which(changeW[,,i]<0)
    enschangeW[i,5]=length(I)/60
  }
  
  ensmeanW=apply(changeW,c(2,3),mean)
  
  bythreshC<-bythreshW<-matrix(NaN,60,4)
  colnames(bythreshC)<-colnames(bythreshW)<-c("20-15","15-10","10-5","<5")
  
  for(i in 1:5)
    for(j in 1:4)
    {
      bythreshC[((12*(i-1)+1):(12*i)),j]<-changeC[,i,j]
      bythreshW[((12*(i-1)+1):(12*i)),j]<-changeW[,i,j]
    }
  
  tiff(file=paste("ECLeventsALL_change_boxplot_bythresh_cool.tiff",sep=""),width=500,height=400)
  boxplot(bythreshC,main="Range of trends by threshold (ECLs p.a.)",ylab="% change in ECLs",ylim=c(-100,100))
  abline(h=0,col="red",lwd=2)
  dev.off()
  tiff(file=paste("ECLeventsALL_change_boxplot_bythresh_warm.tiff",sep=""),width=500,height=400)
  boxplot(bythreshW,main="Range of trends by threshold (ECLs p.a.)",ylab="% change in ECLs",ylim=c(-100,100))
  abline(h=0,col="red",lwd=2)
  dev.off()
#}

############## Looking at distribution of daily max instead of events
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
daylist=read.csv("daylist.csv")
mm=(floor(daylist/100)%%100)

now=array(0,c(length(daylist[,1]),12,5))
dimnames(now)[[2]]<-rep("aaa",12)
dimnames(now)[[1]]=daylist[,1]
dimnames(now)[[3]]=c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
future<-now
dimnames(future)[[1]]=daylist[,2]

cmip=c("echam5","csiromk3","miroc","cccma")
cmipF=c("ECHAM5","MK30","MIROC","CCCMA")
wrf=c("R1","R2","R3")

intcol=c(5,10,10,9,9)

n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(now)[[2]][n]<-dimnames(future)[[2]][n]<-paste(cmip[i],wrf[j])
    
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
      for(t in 1:length(daylist[,1]))
      {
        if(ff==1) I=which(data$Date==daylist[t,1]) else I=which(data$Date==daylist[t,1] & data$Location==1)
        if(length(I)>0) now[t,n,ff]=max(data[I,intcol[ff]])
      }
      
      data=read.csv(filelist2[ff])
      for(t in 1:length(daylist[,1]))
      {
        if(ff==1) I=which(data$Date==daylist[t,2]) else I=which(data$Date==daylist[t,2] & data$Location==1)
        if(length(I)>0) future[t,n,ff]=max(data[I,intcol[ff]])
      }      
    }
    
    n=n+1  
  }

I=which(mm[,1]>=5 & mm[,1]<=10)
J=which(mm[,2]>=5 & mm[,2]<=10)
threshC=apply(now[I,,],c(2,3),quantile,seq(0.95,1,0.01)) ##Anything lower than 0.95 has some thresholds of 0
threshW=apply(now[-I,,],c(2,3),quantile,seq(0.95,1,0.01))
threshC[6,,]<-threshW[6,,]<-100 #For max

countCn<-array(0,c(5,12,5))
dimnames(countCn)[[1]]=paste(">",dimnames(threshC)[[1]][1:5])
dimnames(countCn)[[2]]<-dimnames(threshC)[[2]]
dimnames(countCn)[[3]]<-dimnames(threshC)[[3]]
countCf<-countWn<-countWf<-countCn

for(i in 1:5)
  for(j in 1:12)
    for(k in 1:5)
    {
      K=which(now[I,j,k]>=threshC[i,j,k] & now[I,j,k]<threshC[i+1,j,k])
      countCn[i,j,k]=length(K)
      K=which(now[-I,j,k]>=threshW[i,j,k] & now[-I,j,k]<threshW[i+1,j,k])
      countWn[i,j,k]=length(K)
      K=which(future[J,j,k]>=threshC[i,j,k] & future[J,j,k]<threshC[i+1,j,k])
      countCf[i,j,k]=length(K)
      K=which(future[-J,j,k]>=threshW[i,j,k] & future[-J,j,k]<threshW[i+1,j,k])
      countWf[i,j,k]=length(K)      
    }

diffC=100*(countCf-countCn)/countCn
diffW=100*(countWf-countWn)/countWn

dC2=diffC[,1,]
for(i in 2:12) dC2=cbind(dC2,diffC[,i,])
dW2=diffW[,1,]
for(i in 2:12) dW2=cbind(dW2,diffW[,i,])

a=boxplot(t(dC2),xlab="Quantile ranges (daily max intensity)",ylab="% change in number of ECL days",boxwex=0.5)
abline(h=0,col="red")
b=boxplot(t(dW2),xlab="Quantile ranges (daily max intensity)",ylab="% change in number of ECL days",boxwex=0.5,ylim=c(-100,200))
abline(h=0,col="red")
I=which(b$out>200)
a=unique(b$group[I])
text(a,rep(205,length(a)),"^",cex=2)


diffC2=100*(apply(countCf,c(2,3),sum)-apply(countCn,c(2,3),sum))/apply(countCn,c(2,3),sum)
diffW2=100*(apply(countWf,c(2,3),sum)-apply(countWn,c(2,3),sum))/apply(countWn,c(2,3),sum)
a=boxplot(diffC2,xlab="ECL method",ylab="% change in number of ECL days",boxwex=0.5,ylim=c(-75,75))
abline(h=0,col="red")
b=boxplot(diffW2,xlab="ECLmethod",ylab="% change in number of ECL days",boxwex=0.5,ylim=c(-100,200))
abline(h=0,col="red")
I=which(b$out>200)
a=unique(b$group[I])
text(a,rep(205,length(a)),"^",cex=2)


## Okay, do for intensity thresh 5 + duration thresh >= 5

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")

now<-array(0,c(12,12,5))
dimnames(now)[[1]]<-rep("aaa",12)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
future<-now
intcol=c(8,10,10,10,10)
lencol=c(3,8,8,8,8)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
ithresh=5
dthresh=5

n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j])
      
      filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_9009.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.5.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg0.6.csv",sep=""))
      filelist2=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_6079.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.5.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg0.6.csv",sep=""))
      
      for(ff in 1:5)
      {
        data=read.csv(filelist1[ff])
        I=which(data[,lencol[ff]]>=dthresh)
        data=data[I,]
        b=order(data[,intcol[ff]],decreasing=T)
        thresh2=data[b[20*ithresh],intcol[ff]]
        data2=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          now[n,m,ff]=length(I)
        }
        
        data=read.csv(filelist2[ff])
        I=which(data[,lencol[ff]]>=dthresh)
        data=data[I,]
        data2=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          future[n,m,ff]=length(I)
        }        
      }
      
      n=n+1  
    }
  
  change=((future/now)-1)*100
  change2=change[,,1]
  for(i in 2:5) change2=rbind(change2,change[,,i])

tiff(file=paste("Figures/ECLevents_D5_thresh5_change_bymonth_boxplot.tiff",sep=""),width=600,height=400)
boxplot(change2,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
abline(h=0,col="red",lwd=2)
dev.off()
now2=now[,,1]
for(i in 2:5) now2=rbind(now2,now[,,i])
now3=data.frame(MAM=apply(now2[,3:5],1,sum),JJA=apply(now2[,6:8],1,sum),
                SON=apply(now2[,9:11],1,sum),DJF=apply(now2[,c(1,2,12)],1,sum))

future2=future[,,1]
for(i in 2:5) future2=rbind(future2,future[,,i])
future3=data.frame(MAM=apply(future2[,3:5],1,sum),JJA=apply(future2[,6:8],1,sum),
                   SON=apply(future2[,9:11],1,sum),DJF=apply(future2[,c(1,2,12)],1,sum))
change3=((future3/now3)-1)*100

tiff(file=paste("Figures/ECLevents_D5_thresh5_change_byseason_boxplot.tiff",sep=""),width=400,height=400)
boxplot(change3,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
abline(h=0,col="red",lwd=2)
dev.off()

now4=apply(now,c(1,3),sum)
future4=apply(future,c(1,3),sum)
change4=((future4/now4)-1)*100
tiff(file=paste("Figures/ECLevents_D5_thresh5_change_bymethod_boxplot.tiff",sep=""),width=500,height=400)
boxplot(change4,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
abline(h=0,col="red",lwd=2)
dev.off()

now4=apply(now[,c(1:4,11:12),],c(1,3),sum)
future4=apply(future[,c(1:4,11:12),],c(1,3),sum)
change4=((future4/now4)-1)*100
tiff(file=paste("Figures/ECLevents_D5_thresh5_change_bymethod_boxplot_warm.tiff",sep=""),width=500,height=400)
boxplot(change4,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
abline(h=0,col="red",lwd=2)
dev.off()

now4=apply(now[,5:10,],c(1,3),sum)
future4=apply(future[,5:10,],c(1,3),sum)
change4=((future4/now4)-1)*100
tiff(file=paste("Figures/ECLevents_D5_thresh5_change_bymethod_boxplot_cool.tiff",sep=""),width=500,height=400)
boxplot(change4,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
abline(h=0,col="red",lwd=2)
dev.off()

bysource=matrix(NaN,4,3)
colnames(bysource)=c("GCM","RCM","Method")
bycmip=matrix(NaN,9,4)
colnames(bycmip)=cmip
bywrf<-bymethod<-matrix(NaN,12,3)
colnames(bywrf)=wrf
colnames(bymethod)=c("GV","UM 150km","PG 150km")

type=c("","_warm","_cool")
mlist=list(seq(1,12),c(1:4,11:12),seq(5,10))

for(ss in 1:3)
{
  now4=apply(now[,mlist[[ss]],],c(1,3),sum)
  future4=apply(future[,mlist[[ss]],],c(1,3),sum)
  change4=((future4/now4)-1)*100
  
  tiff(file=paste("Figures/ECLevents_D5_thresh5_change_bymethod_boxplot",type[ss],".tiff",sep=""),width=500,height=400)
  boxplot(change4,main=paste("% change in ECLs from 1990-2009 to 2060-2079"))
  abline(h=0,col="red",lwd=2)
  dev.off()
  
  change4=change4[,c(1,2,4)]
  
  for(i in 1:4)
    for(j in 1:3)
      for(k in 1:3)
        bycmip[3*(k-1)+j,i]<-bywrf[3*(i-1)+k,j]<-bymethod[3*(i-1)+j,k]<-change4[3*(i-1)+j,k]
  
  for(i in 1:4)
  {
    a=change4[seq(3*(i-1)+1,3*i),]
    bysource[i,1]=mean(a)
  }
  for(j in 1:3)
  {
    a=change4[seq(j,12,3),]
    bysource[j,2]=mean(a)
  }
  for(k in 1:3)
  {
    a=change4[,k]
    bysource[k,3]=mean(a)
  }
  
  tiff(file=paste("Figures/ECLevents_D5_thresh5_change_boxplot_bysource",type[ss],".tiff",sep=""),width=500,height=400)
  boxplot(bysource,main="Range of trends by source of uncertainty",ylab="% change in ECLs",ylim=c(-30,30))
  abline(h=0,col="red",lwd=2)
  dev.off()
  tiff(file=paste("Figures/ECLevents_D5_thresh5_change_boxplot_byGCM",type[ss],".tiff",sep=""),width=500,height=400)
  boxplot(bycmip,main="Range of trends by GCM",ylab="% change in ECLs")
  abline(h=0,col="red",lwd=2)
  dev.off()
  tiff(file=paste("Figures/ECLevents_D5_thresh5_change_boxplot_byRCM",type[ss],".tiff",sep=""),width=500,height=400)
  boxplot(bywrf,main="Range of trends by RCM",ylab="% change in ECLs")
  abline(h=0,col="red",lwd=2)
  dev.off()
  tiff(file=paste("Figures/ECLevents_D5_thresh5_change_boxplot_bymethod",type[ss],".tiff",sep=""),width=500,height=400)
  boxplot(bymethod,main="Range of trends by ECL detection method",ylab="% change in ECLs")
  abline(h=0,col="red",lwd=2)
  dev.off()
}
