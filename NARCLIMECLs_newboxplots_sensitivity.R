#################################
#####
#####  Sensitivity testing - different ways to do box plots of projections.
#####
#################################

## Version 1 - standard, threshold depends on method

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

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
threshlist=matrix(0,12,5)

thresh=22
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
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
      threshlist[n,ff]=thresh2
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

now2<-array(0,c(12,5,2,4))
dimnames(now2)[[1]]<-dimnames(now)[[1]][1:12]
dimnames(now2)[[2]]<-dimnames(now)[[3]]
dimnames(now2)[[3]]<-c("Cool","Warm")
dimnames(now2)[[4]]<-c("Individual","Mean","Median","Pooled")
names(dimnames(now2))<-c("Source","Method","Season","Threshold")
future2=now2

now2[,,,1]=abind(apply(now[1:12,5:10,],c(1,3),sum),apply(now[1:12,c(1:4,11:12),],c(1,3),sum),along=3)
future2[,,,1]=abind(apply(future[,5:10,],c(1,3),sum),apply(future[,c(1:4,11:12),],c(1,3),sum),along=3)

## Version 2 - one threshold for each method, which is the median/mean of the different model thresholds. 

nowA=now
futureA=future

for(ff in 1:5)
{
  data1<-data2<-list()
  threshes=rep(0,12)
  n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) 
      
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
      
      
      data1[[n]]=read.csv(filelist1[ff])
      data2[[n]]=read.csv(filelist2[ff])
      
      b=order(data1[[n]][,intcol[ff]],decreasing=T)
      a=data1[[n]][b[20*thresh],intcol[ff]]
      if(is.na(a)) threshes[n]=min(data[,intcol[ff]],na.rm=T) else threshes[n]=a
      n=n+1
    }
  
  thresh1=mean(threshes)
  thresh2=median(threshes)
  
  for(n in 1:12)
  {
    tmp=data1[[n]]
    tmp$Month=floor((tmp$Date1%%10000)/100)
    tmp1n=tmp[tmp[,intcol[ff]]>=thresh1,]
    tmp2n=tmp[tmp[,intcol[ff]]>=thresh2,]
    
    tmp=data2[[n]]
    tmp$Month=floor((tmp$Date1%%10000)/100)
    tmp1f=tmp[tmp[,intcol[ff]]>=thresh1,]
    tmp2f=tmp[tmp[,intcol[ff]]>=thresh2,]
    
    
    for(m in 1:12)
    {
      I=which(tmp1n$Month==m)
      now[n,m,ff]=length(I)
      I=which(tmp2n$Month==m)
      nowA[n,m,ff]=length(I)
      I=which(tmp1f$Month==m)
      future[n,m,ff]=length(I)
      I=which(tmp2f$Month==m)
      futureA[n,m,ff]=length(I)
    }
    
  }
}


now2[,,,2]=abind(apply(now[1:12,5:10,],c(1,3),sum),apply(now[1:12,c(1:4,11:12),],c(1,3),sum),along=3)
future2[,,,2]=abind(apply(future[,5:10,],c(1,3),sum),apply(future[,c(1:4,11:12),],c(1,3),sum),along=3)

now2[,,,3]=abind(apply(nowA[1:12,5:10,],c(1,3),sum),apply(nowA[1:12,c(1:4,11:12),],c(1,3),sum),along=3)
future2[,,,3]=abind(apply(futureA[,5:10,],c(1,3),sum),apply(futureA[,c(1:4,11:12),],c(1,3),sum),along=3)

########### Version 3 - Lumping them all together for each dataset to calculate threshold

for(ff in 1:5)
{
  data1<-data2<-list()
  n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) 
      
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
      
      
      data1[[n]]=read.csv(filelist1[ff])
      data2[[n]]=read.csv(filelist2[ff])
      
      if(n==1) data=data1[[n]] else data=rbind(data,data1[[n]])
      n=n+1
    }
  
  b=order(data[,intcol[ff]],decreasing=T)
  thresh1=data[b[12*20*thresh],intcol[ff]]
  
  for(n in 1:12)
  {
    tmp=data1[[n]]
    tmp$Month=floor((tmp$Date1%%10000)/100)
    tmp1n=tmp[tmp[,intcol[ff]]>=thresh1,]
    
    tmp=data2[[n]]
    tmp$Month=floor((tmp$Date1%%10000)/100)
    tmp1f=tmp[tmp[,intcol[ff]]>=thresh1,]
    
    for(m in 1:12)
    {
      I=which(tmp1n$Month==m)
      now[n,m,ff]=length(I)
      I=which(tmp1f$Month==m)
      future[n,m,ff]=length(I)
    }
    
  }
}

now2[,,,4]=abind(apply(now[1:12,5:10,],c(1,3),sum),apply(now[1:12,c(1:4,11:12),],c(1,3),sum),along=3)
future2[,,,4]=abind(apply(future[,5:10,],c(1,3),sum),apply(future[,c(1:4,11:12),],c(1,3),sum),along=3)

##Individual plots
surnames=c("indivthresh","meanthresh","medianthresh","pooledthresh")

for(n in 1:4)
{
data2=melt(100*((future2[,,,n]/now2[,,,n])-1))
pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_trend_seasonmethod_",surnames[n],".pdf",sep=""),width=8,height=4,pointsize=10)
print(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  scale_y_continuous(breaks=seq(-50, 100, 25)) +
  theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))
dev.off()

pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_trend_seasonmethod_",surnames[n],"_limits.pdf",sep=""),width=8,height=4,pointsize=10)
print(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  scale_y_continuous(breaks=seq(-50, 125, 25),limits = c(-50, 125)) +
  theme_bw() + ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))
dev.off()
}

### Panel plot

data2=melt(100*((future2[,,,1]/now2[,,,1])-1))
plot1<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("a) Individual thresholds") +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 125, 25),limits = c(-50, 125)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))

data2=melt(100*((future2[,,,2]/now2[,,,2])-1))
plot2<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("b) Mean threshold") +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 125, 25),limits = c(-50, 125)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))

data2=melt(100*((future2[,,,3]/now2[,,,3])-1))
plot3<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("c) Median threshold") +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 125, 25),limits = c(-50, 125)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))

data2=melt(100*((future2[,,,4]/now2[,,,4])-1))
plot4<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("d) Pooled threshold") +
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-50, 125, 25),limits = c(-50, 125)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))


library(grid)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf(file=paste("ECLevents_ALL_thresh",thresh,"_ggplot2_trend_seasonmethod_threshmethod_panel2.pdf",sep=""),width=11,height=6,pointsize=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(2, 1))
print(plot3, vp = vplayout(1, 2))
print(plot4, vp = vplayout(2, 2))
dev.off()


### Analyses
change=100*((future2/now2)-1)
apply(change,c(3,4),median)



########## Also, panel figure for diff thresholds

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(12,12,5))
dimnames(now)[[1]]<-rep("ncep",12)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now[1:12,,]
intcol=c(8,10,10,13,13)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
thresh=c(25,20,15,10,5)

now2<-array(0,c(12,5,2,5))
dimnames(now2)[[1]]<-dimnames(now)[[1]]
dimnames(now2)[[2]]<-dimnames(now)[[3]]
dimnames(now2)[[3]]<-c("Cool","Warm")
dimnames(now2)[[4]]<-thresh
names(dimnames(now2))<-c("Source","Method","Season","Threshold")
future2=now2

for(t in 1:5)
{
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(now)[[1]][n]<-dimnames(now2)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
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
      thresh2=data[b[20*thresh[t]],intcol[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcol[ff]],na.rm=T)
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
now2[,,,t]=abind(apply(now[1:12,5:10,],c(1,3),sum),apply(now[1:12,c(1:4,11:12),],c(1,3),sum),along=3)
future2[,,,t]=abind(apply(future[,5:10,],c(1,3),sum),apply(future[,c(1:4,11:12),],c(1,3),sum),along=3)
}

### Panel plot

data2=melt(100*((future2[,,,1]/now2[,,,1])-1))
plot1<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("a) 25 p.a.") + 
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-100, 350, 50)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))

data2=melt(100*((future2[,,,2]/now2[,,,2])-1))
plot2<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("b) 20 p.a.") + 
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-100, 350, 50)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))

data2=melt(100*((future2[,,,3]/now2[,,,3])-1))
plot3<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("c) 15 p.a.") + 
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-100, 350, 50)) + 
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))

data2=melt(100*((future2[,,,4]/now2[,,,4])-1))
plot4<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("d) 10 p.a.") + 
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-100, 350, 50)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))

data2=melt(100*((future2[,,,5]/now2[,,,5])-1))
plot5<-(ggplot(data2, aes(x = Method, y = value, fill = Season)) +
          geom_boxplot() + ggtitle("e) 5 p.a.") + 
          scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
          scale_y_continuous(breaks=seq(-100, 350, 50)) +
          theme_bw() + theme(plot.title = element_text(hjust=0)) +
          ylab("Percentage change in ECLs") + xlab("") +  geom_hline(yintercept = 0))

library(grid)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf(file=paste("ECLevents_ALL_threshALL_ggplot2_trend_seasonmethod_panel4.pdf",sep=""),width=11,height=9,pointsize=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(2, 1))
print(plot3, vp = vplayout(3, 1))
print(plot4, vp = vplayout(1, 2))
print(plot5, vp = vplayout(2, 2))
dev.off()

change=100*((future2/now2)-1)
apply(change,c(3,4),median)

################
## Threshold sensitivity - comparing a thresh of 0.86 & 1.12 for all in UM 150

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(16,12,2))
dimnames(now)[[1]]<-rep("ncep",16)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now[1:12,,]
intcol=c(8,10,10,13,13)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
threshlist=matrix(0,12,5)

thresh=c(0.86,1.12)
ff=2
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
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
    
    for(t in 1:2)
    {
      data=read.csv(filelist1[ff])
      data2=data[data[,intcol[ff]]>=thresh[t],]
      mm=floor((data2$Date1%%10000)/100)
      for(m in 1:12)
      {
        I=which(mm==m)
        now[n,m,t]=length(I)
      }
      
      data=read.csv(filelist2[ff])
      data2=data[data[,intcol[ff]]>=thresh[t],]
      mm=floor((data2$Date1%%10000)/100)
      for(m in 1:12)
      {
        I=which(mm==m)
        future[n,m,t]=length(I)
      } 
      
    }
    n=n+1  
  }

now2=abind(apply(now[1:12,5:10,],c(1,3),sum),apply(now[1:12,c(1:4,11:12),],c(1,3),sum),along=3)
future2=abind(apply(future[,5:10,],c(1,3),sum),apply(future[,c(1:4,11:12),],c(1,3),sum),along=3)

change=100*((future2/now2)-1)


########################
## Interann Var

rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)

now<-array(0,c(12,20,5))
dimnames(now)[[1]]<-rep("ncep",12)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","LAP 150km","LAP 50km","PG 150km","PG 50km")
names(dimnames(now))<-c("Source","Month","Method")
future<-now[1:12,,]
intcol=c(8,10,10,13,13)

cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

thresh=22
n=1
for(i in 1:4)
  for(j in 1:3)
  {
    dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-paste(cmip[i],wrf[j]) 
    
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
      yy=floor(data2$Date1/10000)
      years=unique(yy)
      for(y in 1:20)
      {
        I=which(yy==years[y])
        now[n,y,ff]=length(I)
      }
      
      data=read.csv(filelist2[ff])
      data2=data[data[,intcol[ff]]>=thresh2,]
      yy=floor(data2$Date1/10000)
      years=unique(yy)
      for(y in 1:20)
      {
        I=which(yy==years[y])
        future[n,y,ff]=length(I)
      }
      
    }
    n=n+1  
  }