## Using ONLY the 50 km resolution, comparing 500 & 200km radius for projections
## For both UM and Ale
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
intcol=c(10,10,17,18)
thresh=c(50,40,30,20,10,NaN)
nowC<-nowW<-futureC<-futureW<-array(NaN,c(12,20,4,5))

 n=1
  
  for(i in 1:4)
    for(j in 1:3)
    { 
      filelist1=c(paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad5cv04/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad5cv04_9009.csv",sep=""),
                  paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_PG2.csv",sep=""),
                  paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_PG2.csv",sep=""))
      
      filelist2=c(paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad5cv04/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad5cv04_6079.csv",sep=""),
                  paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_PG2.csv",sep=""),
                  paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_PG2.csv",sep=""))
      

      for(ff in 1:4) 
      {
        data=read.csv(filelist1[ff])
        b=order(data[,intcol[ff]],decreasing=T)
        dataF=read.csv(filelist2[ff])
        #          data=data[data$Length>=length,]
        
        for(t in 1:5)        
        {
        
        thresh1=data[b[20*thresh[t]],intcol[ff]]
        if(is.na(thresh1)) thresh1=data[length(data[,1]),intcol[ff]]
        if(t<5) thresh2=data[b[20*thresh[t+1]],intcol[ff]] else thresh2=9999999
        data2=data[(data[,intcol[ff]]>=thresh1 & data[,intcol[ff]]<thresh2),]
        
        yy=floor(data2$Date1/10000)
        mm=floor((data2$Date1%%10000)/100)
        yy2=yy
        yy2[mm<=4]=yy2[mm<=4]-1
        years=unique(yy)
        for(y in 1:20)
        {
          I=which(mm>=5 & mm<=10 & yy==years[y])
          nowC[n,y,ff,t]=length(I)
          I=which((mm>=11 | mm<=4) & yy2==years[y])
          if(y<20) nowW[n,y,ff,t]=length(I)
        }
        
        data2=dataF[(dataF[,intcol[ff]]>=thresh1 & dataF[,intcol[ff]]<thresh2),]
        yy=floor(data2$Date1/10000)
        mm=floor((data2$Date1%%10000)/100)
        yy2=yy
        yy2[mm<=4]=yy2[mm<=4]-1
        years=unique(yy)
        for(y in 1:20)
        {
          I=which(mm>=5 & mm<=10 & yy==years[y])
          futureC[n,y,ff,t]=length(I)
          I=which((mm>=11 | mm<=4) & yy2==years[y])
          if(y<20) futureW[n,y,ff,t]=length(I)
        }   
        }
      }
      
      n=n+1  
    }

changeC=((apply(futureC,c(1,3,4),mean,na.rm=T)/apply(nowC,c(1,3,4),mean,na.rm=T))-1)*100
changeW=((apply(futureW,c(1,3,4),mean,na.rm=T)/apply(nowW,c(1,3,4),mean,na.rm=T))-1)*100

dimnames(changeC)[[2]]<-dimnames(changeW)[[2]]<-c("UM 200km","UM 500km","PG 200km","PG 500km")
dimnames(changeC)[[1]]<-rep("aaa",12)
boxplot(changeC[,,5])
boxplot(changeW[,,5])

########## Four different intensity methods:
## PG200, PG500, MSLP, LAP, Depth

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
intcol=c(17,18,16,20,21)
thresh=c(50,22,10,NaN)
nowC<-nowW<-futureC<-futureW<-array(NaN,c(12,20,5,3))

n=1

for(i in 1:4)
  for(j in 1:3)
  { 
    filelist1=paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_PG2.csv",sep="")    
    filelist2=paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_PG2.csv",sep="")
    data=read.csv(filelist1)
    dataF=read.csv(filelist2)

    for(ff in 1:5) 
    {
      if(ff==3) b=order(data[,intcol[ff]],decreasing=F) else b=order(data[,intcol[ff]],decreasing=T)
    
      for(t in 1:3)        
      {
        
        thresh1=data[b[20*thresh[t]],intcol[ff]]
        if(is.na(thresh1)) thresh1=data[length(data[,1]),intcol[ff]]
        data2=data[(data[,intcol[ff]]>=thresh1),]
        
        yy=floor(data2$Date1/10000)
        mm=floor((data2$Date1%%10000)/100)
        yy2=yy
        yy2[mm<=4]=yy2[mm<=4]-1
        years=unique(yy)
        for(y in 1:20)
        {
          I=which(mm>=5 & mm<=10 & yy==years[y])
          nowC[n,y,ff,t]=length(I)
          I=which((mm>=11 | mm<=4) & yy2==years[y])
          if(y<20) nowW[n,y,ff,t]=length(I)
        }
        
        data2=dataF[(dataF[,intcol[ff]]>=thresh1),]
        yy=floor(data2$Date1/10000)
        mm=floor((data2$Date1%%10000)/100)
        yy2=yy
        yy2[mm<=4]=yy2[mm<=4]-1
        years=unique(yy)
        for(y in 1:20)
        {
          I=which(mm>=5 & mm<=10 & yy==years[y])
          futureC[n,y,ff,t]=length(I)
          I=which((mm>=11 | mm<=4) & yy2==years[y])
          if(y<20) futureW[n,y,ff,t]=length(I)
        }   
      }
    }
    
    n=n+1  
  }

changeC=((apply(futureC,c(1,3,4),mean,na.rm=T)/apply(nowC,c(1,3,4),mean,na.rm=T))-1)*100
changeW=((apply(futureW,c(1,3,4),mean,na.rm=T)/apply(nowW,c(1,3,4),mean,na.rm=T))-1)*100

dimnames(changeC)[[2]]<-dimnames(changeW)[[2]]<-c("PG200","PG500","MSLP","Laplacian","Depth")
boxplot(changeC[,,3],ylim=c(-50,100),main="Trend in cool season PG ECLs, 50 events p.a.",ylab="% change")
abline(h=0,col="red")
boxplot(changeW[,,3],ylim=c(-50,100),main="Trend in warm season PG ECLs, 50 events p.a.",ylab="% change")
abline(h=0,col="red")


########
## Okay, now trends as a function of radius, duration, etc.
## Or maybe just changes in the distribution across ALL models

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

n=1

for(i in 1:4)
  for(j in 1:3)
  { 
    filelist1=paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_PG2.csv",sep="")    
    filelist2=paste("Alejandro/v2/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_PG2.csv",sep="")
    if(n==1)
      {
      dataN=read.csv(filelist1)
      dataF=read.csv(filelist2)
    } else
    {
      dataN=rbind(dataN,read.csv(filelist1))
      dataF=rbind(dataF,read.csv(filelist2))
    }
    n=n+1
  }

makePDF = function(data1,data2,xlabel,tit="",labloc="topleft") {
  a=density(data1,na.rm=T)
  b=density(data2,na.rm=T)
  
  lims=range(data1,data2,na.rm=T)
  if((lims[2]-lims[1])<10)
  {
    lims[1]=floor(lims[1])
    lims[2]=ceiling(lims[2])
  } else {
    lims[1]=floor(lims[1]/5)*5
    lims[2]=ceiling(lims[2]/5)*5
  }
  
  plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
       xlab=xlabel,ylab="Frequency",cex.main=1.2,
       main=tit)
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=c("1990-2009","2060-2079"),
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
}

mm1=floor(dataN$Date1/100)%%100
mm2=floor(dataF$Date1/100)%%100
makePDF(dataN$MSLP2[mm1>=5 & mm1<=10],dataF$MSLP2[mm2>=5 & mm2<=10],"Maximum central pressure","Distribution of ECL intensity(cool)","topright")
makePDF(dataN$MSLP2[mm1>=11 | mm1<=4],dataF$MSLP2[mm2>=11 | mm2<=4],"Maximum central pressure","Distribution of ECL intensity(warm)","topright")
