rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(16,12,5))
dimnames(now)[[1]]<-rep("ncep",16)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now[1:12,,]
intcol=c(8,10,10,13,13)

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")

thresh=22
for(thresh in c(22,10,5,2,1))
{
  n=1
  for(i in 1:5)
    for(j in 1:3)
    {
      if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
      
      filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                  paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_final.csv",sep=""),
                  paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_final.csv",sep=""))
      filelist2=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                  paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.7_final.csv",sep=""),
                  paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.2_final.csv",sep=""))
      
      for(ff in 1:5)
      {
        data=read.csv(filelist1[ff])
        b=order(data[,intcol[ff]],decreasing=T)
        thresh2=data[b[20*thresh],intcol[ff]]
        if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
        data2=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          now[n,m,ff]=length(I)
        }
        
        if(i<5)
        {
          data=read.csv(filelist2[ff])
          data2=data[data[,intcol[ff]]>=thresh2,]
          mm=floor((data2$Date1%%10000)/100)
          for(m in 1:12)
          {
            I=which(mm==m)
            future[n,m,ff]=length(I)
          } 
        }
      }
      
      n=n+1  
    }
  
  
  ##### Add NCEP!
  
  filelist=c("Fei/GVevents_ncep_9009.csv",
             "outputUM/proj100/outputUM_ncep_rad2cv06/ECLevents_umelb_ncep_proj100_rad2cv06_9009.csv","",
             "Alejandro/ECLevents_Alejandro_ncep_res150_9009_pg0.5.csv")
  
  for(ff in c(1,2,4))
  {
    data=read.csv(filelist[ff])
    b=order(data[,intcol[ff]],decreasing=T)
    thresh2=data[b[20*thresh],intcol[ff]]
    if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
    data2=data[data[,intcol[ff]]>=thresh2,]
    mm=floor((data2$Date1%%10000)/100)
    for(m in 1:12)
    {
      I=which(mm==m)
      now[n,m,ff]=length(I)
    }
  }
  
  prop=100*apply(now[,5:10,],c(1,3),sum)/apply(now,c(1,3),sum)
  prop=as.data.frame(prop)
  prop$Type=c(rep("CMIP",12),rep("NCEP",3),NaN)
  
  data2=melt(prop[1:15,])  
  pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_coolprop_bytype.pdf",sep=""),height=5,width=8)
  ggplot(data2, aes(x = variable, y = value, fill = Type)) +
    geom_boxplot() +
    scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
    theme_bw() + ylab("Proportion of ECLs") + xlab("") + geom_hline(yintercept = 57.5)
  dev.off()
  
  now2=abind(apply(now[1:12,5:10,],c(1,3),sum),apply(now[1:12,c(1:4,11:12),],c(1,3),sum),along=3)
  future2=abind(apply(future[,5:10,],c(1,3),sum),apply(future[,c(1:4,11:12),],c(1,3),sum),along=3)
  dimnames(now2)[[3]]<-dimnames(future2)[[3]]<-c("Cool","Warm")
  names(dimnames(now2))<-names(dimnames(future2))<-c("Source","Method","Season")
  
  data2=melt(100*((future2/now2)-1))
  pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_trend_seasonmethod.pdf",sep=""),height=5,width=8)
  ggplot(data2, aes(x = Method, y = value, fill = Season)) +
    geom_boxplot() +
    scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
    scale_y_continuous(breaks=seq(-50, 100, 25)) +
    theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0)
  dev.off()
  
  pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_trend_seasonmethod_limited.pdf",sep=""),height=5,width=8)
  ggplot(data2, aes(x = Method, y = value, fill = Season)) +
    geom_boxplot() +
    scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
    scale_y_continuous(breaks=seq(-50, 100, 25),limits = c(-50, 100)) +
    theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0)
  dev.off()
  
  change=100*((future2/now2)-1)
  change2=apply(change,c(1,3),mean)
  clist=rep(c("black","red","blue"),4)
  plist=c(1,1,1,2,2,2,3,3,3,4,4,4)
  
  pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_trend_vGCMRCM.pdf",sep=""),height=5,width=6)
  plot(change2[,1],change2[,2],type="p",pch=plist,col=clist,cex=2,lwd=2,
       xlim=c(-30,10),ylim=c(-30,80),xlab="% change in cool season",ylab="% change in warm season")
  abline(h=0,col="gray")
  abline(v=0,col="gray")
  legend("topleft",c("RCM1","RCM2","RCM3","","ECHAM5","CSIROMK3.0","MIROC3.2","CGCM3.1"),
         pch=c(rep(4,4),1:4),col=c("black","red","blue","white",rep("black",4)),
         pt.lwd=2,pt.cex=1.5,bty="n",ncol=2)
  dev.off()
  
}

#####ECL days by location


intcol=c(10,10,9,9)
intcolE=rep(10,4)

ECLdaysN=array(NaN,c(12,12,4,3))
dimnames(ECLdaysN)[[1]]<-rep("ncep",12)
dimnames(ECLdaysN)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(ECLdaysN)[[3]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
dimnames(ECLdaysN)[[4]]=c("Whole domain","Coastal region 1","Coastal region 2")
names(dimnames(ECLdaysN))<-c("Source","Month","Method","Region")
ECLdaysF<-ECLdaysN

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
names=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
thresh=22
n=1

for(i in 1:4)
  for(j in 1:3)
  {    
    rownames(ECLdaysN)[n]<-rownames(ECLdaysF)[n]<-paste(cmip[i],wrf[j])
    
    filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
    
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.8.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.3.csv",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelistE[ff])
      b=order(data[,intcolE[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcolE[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcolE[ff]])
      
      data=read.csv(filelist1[ff])
      data$Location2=0
      I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
      data$Location2[I]<-1
      I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
      data$Location2[I]<-1
      I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
      data$Location2[I]<-1
      data$Location3=0
      I<-which(data$Lon>=149 & data$Lon<=154 & data$Lat<(-37) & data$Lat>=-40)
      data$Location3[I]<-1
      I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(154+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
      data$Location3[I]<-1
      I<-which(data$Lon>=152 & data$Lon<=157 & data$Lat<=(-25) & data$Lat>=-31)
      data$Location3[I]<-1
      
      dataN=data[data[,intcol[ff]]>=thresh2,]     
      mmN=floor(dataN$Date/100)%%100
      locN=cbind(dataN$Location,dataN$Location2,dataN$Location3)
      
      data=read.csv(filelist2[ff])
      data$Location2=0
      I<-which(data$Lon>=149 & data$Lon<=156 & data$Lat<(-37) & data$Lat>=-40)
      data$Location2[I]<-1
      I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(156+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
      data$Location2[I]<-1
      I<-which(data$Lon>=152 & data$Lon<=159 & data$Lat<=(-25) & data$Lat>=-31)
      data$Location2[I]<-1
      data$Location3=0
      I<-which(data$Lon>=149 & data$Lon<=154 & data$Lat<(-37) & data$Lat>=-40)
      data$Location3[I]<-1
      I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(154+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
      data$Location3[I]<-1
      I<-which(data$Lon>=152 & data$Lon<=157 & data$Lat<=(-25) & data$Lat>=-31)
      data$Location3[I]<-1
      
      dataF=data[data[,intcol[ff]]>=thresh2,]
      mmF=floor(dataF$Date/100)%%100
      locF=cbind(dataF$Location,dataF$Location2,dataF$Location3)
      
      for(m in 1:12)
        for(l in 1:3)
        {
          I=which(mmN==m & locN[,l]==1)
          ECLdaysN[n,m,ff,l]=length(unique(dataN$Date[I]))
          I=which(mmF==m & locF[,l]==1)
          ECLdaysF[n,m,ff,l]=length(unique(dataF$Date[I]))
        }
    }      
    n=n+1  
  }

changeC=100*((apply(ECLdaysF[,5:10,,],c(1,3,4),sum)/apply(ECLdaysN[,5:10,,],c(1,3,4),sum))-1)
changeW=100*((apply(ECLdaysF[,c(11:12,1:4),,],c(1,3,4),sum)/apply(ECLdaysN[,c(11:12,1:4),,],c(1,3,4),sum))-1)

names(dimnames(changeC))<-c("Source","Method","Region")
data=melt(changeC[,,c(1,3)])

######## ARI by season - basic GEV
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(extRemes)
rlev=20
nowC.ari<-array(0,c(12,5))
dimnames(nowC.ari)[[1]]<-rep("aaa",12)
dimnames(nowC.ari)[[2]]<-c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
futureW.ari<-nowW.ari<-futureC.ari<-nowC.ari

daylist=read.csv("daylist.csv")
mm=(floor(daylist/100)%%100)
cmip=c("echam5","csiromk3","miroc","cccma")
cmipF=c("ECHAM5","MK30","MIROC","CCCMA")
wrf=c("R1","R2","R3")

intcol=c(5,10,10,9,9)

for(rlev in c(2,10,20,50))
{
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
        yy=floor(data$Date/10000)
        mm=floor((data$Date%%10000)/100)
        years=cbind(unique(yy),matrix(0,length(unique(yy)),2))
        for(t in 1:length(years[,1]))
        {
          if(ff==1) I=which(yy==years[t,1] & mm>=5 & mm <=10) else I=which(yy==years[t,1] & mm>=5 & mm <=10 & data$Location==1)
          years[t,2]=max(data[I,intcol[ff]])
          if(ff==1) I=which(yy==years[t,1] & (mm>=11 | mm <=4)) else I=which(yy==years[t,1] & (mm>=11 | mm <=4) & data$Location==1)
          years[t,3]=max(data[I,intcol[ff]])
        }
        a=fevd(years[,2])
        nowC.ari[n,ff]=return.level(a,rlev)
        a=fevd(years[,3])
        nowW.ari[n,ff]=return.level(a,rlev)
        
        data=read.csv(filelist2[ff])
        yy=floor(data$Date/10000)
        mm=floor((data$Date%%10000)/100)
        years=cbind(unique(yy),matrix(0,length(unique(yy)),2))
        for(t in 1:length(years[,1]))
        {
          if(ff==1) I=which(yy==years[t,1] & mm>=5 & mm <=10) else I=which(yy==years[t,1] & mm>=5 & mm <=10 & data$Location==1)
          years[t,2]=max(data[I,intcol[ff]])
          if(ff==1) I=which(yy==years[t,1] & (mm>=11 | mm <=4)) else I=which(yy==years[t,1] & (mm>=11 | mm <=4) & data$Location==1)
          years[t,3]=max(data[I,intcol[ff]])
        }
        a=fevd(years[,2])
        futureC.ari[n,ff]=return.level(a,rlev)
        a=fevd(years[,3])
        futureW.ari[n,ff]=return.level(a,rlev)      
      }    
      n=n+1  
    }
  
  changeC=100*((futureC.ari/nowC.ari)-1)
  changeW=100*((futureW.ari/nowW.ari)-1)
  
  change2=as.data.frame(rbind(changeC,changeW))
  change2$Season=c(rep("Cool",12),rep("Warm",12))
  write.csv(change2,"tmp.csv")
  
  data2=melt(change2)
  pdf(paste("ECLeventsALL_ARI",rlev,"gev_ggplot2_trend_seasonmethod.pdf",sep=""),height=5,width=8)
  print(ggplot(data2, aes(x = variable, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 150, 25)) +
          theme_bw() + ylab(paste("Percentage change in intensity of 1-in-",rlev," year ECL",sep="")) + xlab("") +  geom_hline(yintercept = 0))
  dev.off()
}


####### ECL mean duration in location

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(15,5))
dimnames(now)[[1]]<-rep("ncep",15)
dimnames(now)[[2]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
future<-now[1:12,]
intcol=c(8,10,10,10,10)

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")

thresh=22
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
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
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      data2=data[data[,intcol[ff]]>=thresh2,]
      if(ff==1) now[n,ff]=mean(data2$Length*6) else now[n,ff]=mean(data2$Length2*6)
      
      if(i<5)
      {
        data=read.csv(filelist2[ff])
        data2=data[data[,intcol[ff]]>=thresh2,]
        if(ff==1) future[n,ff]=mean(data2$Length*6) else future[n,ff]=mean(data2$Length2*6)  
      }
    }
    
    n=n+1  
  }  

now2=as.data.frame(now)
now2$Type=c(rep("CMIP",12),rep("NCEP",3))
data2=melt(now2)  
pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_duration_bytype.pdf",sep=""),height=5,width=8)
ggplot(data2, aes(x = variable, y = value, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  theme_bw() + ylab("Mean duration (hours)") + xlab("")
dev.off()



######## Sensitivity - different radii & durations

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(12,12,4))
dimnames(now)[[1]]<-rep("aaa",12)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("200km, 6hr","500km, 6hr","200km, 24hr","500km, 24hr")
names(dimnames(now))<-c("Source","Month","Method")
future<-now

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

for(thresh in c(2,5,10,22))
{
  n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
      
      filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_9009.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad5cv015/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad5cv015_9009.csv",sep=""))
      filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_6079.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad5cv015/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad5cv015_6079.csv",sep=""))
      
      for(ff in 1:2)
      {
        data=read.csv(filelist1[ff])
        b=order(data$CV2,decreasing=T)
        thresh2=data$CV2[b[20*thresh]]
        if(is.na(thresh2)) thresh2=min(data$CV2,na.rm=T)
        data2=data[data$CV2>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          now[n,m,ff]=length(I)
        }
        
        data=data[data$Length>=5,]
        b=order(data$CV2,decreasing=T)
        thresh3=data$CV2[b[20*thresh]]
        if(is.na(thresh3)) thresh3=min(data$CV2,na.rm=T)
        data2=data[data$CV2>=thresh3,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          now[n,m,ff+2]=length(I)
        }
        
        data=read.csv(filelist2[ff])
        data2=data[data$CV2>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          future[n,m,ff]=length(I)
        } 
        
        data=data[data$Length>=5,]
        data2=data[data$CV2>=thresh3,]
        mm=floor((data2$Date1%%10000)/100)
        for(m in 1:12)
        {
          I=which(mm==m)
          future[n,m,ff+2]=length(I)
        } 
      }
      n=n+1
    }
  
  
  now2=abind(apply(now[,5:10,],c(1,3),sum),apply(now[,c(1:4,11:12),],c(1,3),sum),along=3)
  future2=abind(apply(future[,5:10,],c(1,3),sum),apply(future[,c(1:4,11:12),],c(1,3),sum),along=3)
  dimnames(now2)[[3]]<-dimnames(future2)[[3]]<-c("Cool","Warm")
  names(dimnames(now2))<-names(dimnames(future2))<-c("Source","Method","Season")
  
  data2=melt(100*((future2/now2)-1))
  pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_trend_season_sensitivity_limited.pdf",sep=""),height=5,width=8)
  print(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 100, 25),limits=c(-50,100)) +
          theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))
  dev.off()
}

######## By year, vs. var

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(15,20,12,5))
ylist<-array(0,c(20,5))
dimnames(now)[[1]]<-rep("ncep",15)
dimnames(now)[[3]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[4]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now[1:12,,,]
intcol=c(8,10,10,13,13)

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")

thresh=5
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
    year=seq(1990,2009)
    filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_final.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_final.csv",sep=""))
    filelist2=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.7_final.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.2_final.csv",sep=""))
    
    for(ff in 1:5)
    {
      data=read.csv(filelist1[ff])
      year=unique(floor(data$Date1/10000))
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      data2=data[data[,intcol[ff]]>=thresh2,]
      mm=floor((data2$Date1%%10000)/100)
      yy=floor(data2$Date1/10000)
      if(i==1 & j==1) ylist[,ff]=year
      for(y in 1:20)
        for(m in 1:12)
        {
          I=which(mm==m & yy==year[y])
          now[n,y,m,ff]=length(I)
        }
      
      if(i<5)
      {
        data=read.csv(filelist2[ff])
        year=unique(floor(data$Date1/10000))
        data2=data[data[,intcol[ff]]>=thresh2,]
        mm=floor((data2$Date1%%10000)/100)
        yy=floor(data2$Date1/10000)
        if(i==1 & j==1) ylist[,ff]=year
        for(y in 1:20)
          for(m in 1:12)
          {
            I=which(mm==m & yy==year[y])
            future[n,y,m,ff]=length(I)
          } 
      }
    }
    
    n=n+1  
  }

library(abind)
now2=abind(apply(now[1:12,,5:10,],c(1,2,4),sum),apply(now[1:12,,c(1:4,11:12),],c(1,2,4),sum),along=4)
future2=abind(apply(future[,,5:10,],c(1,2,4),sum),apply(future[,,c(1:4,11:12),],c(1,2,4),sum),along=4)

std=apply(now2,c(1,3,4),sd,na.rm=T)
change=100*((apply(future2,c(1,3,4),mean)/apply(now2,c(1,3,4),mean))-1)
change2=(apply(future2,c(1,3,4),mean)-apply(now2,c(1,3,4),mean))/std
dimnames(change2)[[3]]<-c("Cool","Warm")
names(dimnames(change2))<-c("Source","Method","Season")
data2=melt(change2)

pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_trend_seasonmethod_stdev.pdf",sep=""),height=5,width=10)
ggplot(data2, aes(x = Method, y = value, fill = Season)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  scale_y_continuous(breaks=seq(-1.5,4,0.5)) +
  theme_bw() + ylab("Change in mean ECLs (standard deviations)") + xlab("") +  geom_hline(yintercept = 0)
dev.off()

pval=array(0,c(12,5,2))
for(i in 1:12)
  for(j in 1:5)
    for(k in 1:2)
    {
      a=t.test(now2[i,,j,k],y=future2[i,,j,k])
      pval[i,j,k]=a$p.value
    }

pval2<-matrix(0,5,2)
for(i in 1:5)
  for(j in 1:2)
    pval2[i,j]=mean(pval[,i,j]<0.05)
#Stat sig difference in ~35% of individual members

###### Most intense ECL per year

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(15,20,5))
dimnames(now)[[1]]<-rep("ncep",15)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now[1:12,,]
intcol=c(8,10,10,13,13)

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")

n=1
for(i in 1:5)
  for(j in 1:3)
  {
    if(i<5) dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) else dimnames(now)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
    filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_final.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_final.csv",sep=""))
    filelist2=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.7_final.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.2_final.csv",sep=""))
    
    for(ff in 1:5)
    {
      data=read.csv(filelist1[ff])
      b=order(data[,intcol[ff]],decreasing=T)
      now[n,,ff]=data[b[1:20],intcol[ff]]
      #         year=floor(data$Date1/10000)
      #         ylist=unique(year)
      #         for(m in 1:20)
      #         {
      #           I=which(year==ylist[m])
      #           now[n,m,ff]=max(data[I,intcol[ff]])
      #         }
      
      if(i<5)
      {
        data=read.csv(filelist2[ff])
        b=order(data[,intcol[ff]],decreasing=T)
        future[n,,ff]=data[b[1:20],intcol[ff]]
        #           year=floor(data$Date1/10000)
        #           ylist=unique(year)
        #           for(m in 1:20)
        #           {
        #             I=which(year==ylist[m])
        #             future[n,m,ff]=max(data[I,intcol[ff]])
        #           }
      }
    }
    
    n=n+1  
  }


tt=matrix(0,12,5)
for(i in 1:12)
  for(j in 1:5)
  {
    a=t.test(future[i,,j],now[i,,j])
    tt[i,j]=a$p.value
  }


### What are my thresholds?

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(15,5))
dimnames(now)[[1]]<-rep("ncep",15)
dimnames(now)[[2]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
intcol=c(8,10,10,13,13)

cmip=c("ncep","echam5","csiromk3","miroc","cccma")
cmip2=c("NCEP","ECHAM5","CCCMA3.1","MIROC3.2","CSIROMk3.0")
wrf=c("R1","R2","R3")

thresh=22
n=1
for(i in 1:5)
  for(j in 1:3)
  {
    dimnames(now)[[1]][n]<-paste(cmip2[i],"-",wrf[j],sep="")
    year=seq(1990,2009)
    filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_final.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_final.csv",sep=""))
    
    for(ff in 1:5)
    {
      data=read.csv(filelist1[ff])
      year=unique(floor(data$Date1/10000))
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      now[n,ff]=thresh2
    }
    
    n=n+1  
  }

############## Change for Location 2

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)
cmip=c("echam5","csiromk3","miroc","cccma","ncep")
wrf=c("R1","R2","R3")
namelist=rep("aaa",12)

now<-array(NaN,c(12,2,4))
dimnames(now)[[1]]<-rep("ncep",12)
dimnames(now)[[2]]=c("Cool","Warm")
dimnames(now)[[3]]=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now

thresh=22
intcol=c(10,10,13,13) # CVmax/PG200max
intcolF=c(10,10,10,10)

n=1
for(i in 1:4)
  for(j in 1:3)
  {
    namelist[n]<-paste(cmip[i],wrf[j])
    filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.2_impacts_v2.csv",sep=""))
    filelist2=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_6079_impacts_v2.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_6079_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.7_impacts_v2.csv",sep=""),
                paste("Alejandro/final/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg1.2_impacts_v2.csv",sep=""))
    
    for(ff in 1:4)
    {
      data=read.csv(filelistE[ff])
      b=order(data[,intcol[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
      
      data=read.csv(filelist1[ff])
      data$Location2=0
      I<-which(data$Lon>=149 & data$Lon<=155 & data$Lat<(-37) & data$Lat>=-40)
      data$Location2[I]<-1
      I<-which(data$Lon>=(149+(37+data$Lat)/2) & data$Lon<=(155+(37+data$Lat)/2) & data$Lat<(-31) & data$Lat>=-37)
      data$Location2[I]<-1
      I<-which(data$Lon>=152 & data$Lon<=158 & data$Lat<=(-25) & data$Lat>=-31)
      data$Location2[I]<-1
      data=data[(data[,intcolF[ff]]>=thresh2 & data$Location2==1),]
      mm=floor((data$Date%%10000)/100)
      
      data2=read.csv(filelist2[ff])
      data2$Location2=0
      I<-which(data2$Lon>=149 & data2$Lon<=155 & data2$Lat<(-37) & data2$Lat>=-40)
      data2$Location2[I]<-1
      I<-which(data2$Lon>=(149+(37+data2$Lat)/2) & data2$Lon<=(155+(37+data2$Lat)/2) & data2$Lat<(-31) & data2$Lat>=-37)
      data2$Location2[I]<-1
      I<-which(data2$Lon>=152 & data2$Lon<=158 & data2$Lat<=(-25) & data2$Lat>=-31)
      data2$Location2[I]<-1
      data2=data2[(data2[,intcolF[ff]]>=thresh2 & data2$Location2==1),]
      mm2=floor((data2$Date%%10000)/100)
      data[data==-Inf | data==Inf]=NaN
      data2[data2==-Inf | data2==Inf]=NaN
      
      ### Density
      I=which(mm>=5 & mm<=10)
      now[n,1,ff]=length(unique(data$ID[I]))
      now[n,2,ff]=length(unique(data$ID[-I]))     
      
      I=which(mm2>=5 & mm2<=10)
      future[n,1,ff]=length(unique(data2$ID[I]))
      future[n,2,ff]=length(unique(data2$ID[-I]))       
    }
    n=n+1  
  }


change=100*((future/now)-1)
names(dimnames(change))=c("Source","Season","Method")

data2=melt(change)
pdf(paste("ECLeventsCLOSE_thresh",thresh,"_ggplot2_trend_seasonmethod.pdf",sep=""),height=5,width=8)
ggplot(data2, aes(x = Method, y = value, fill = Season)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  scale_y_continuous(breaks=seq(-50, 100, 25)) +
  theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0)
dev.off()

change2=change
