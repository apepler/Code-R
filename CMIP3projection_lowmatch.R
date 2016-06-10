rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(sp)
library(abind)
library(ggplot2)
library(reshape)

intcol=c(10,10,9,9)
intcolE=rep(10,4)
cmip=c("ncep","echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")
names=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
res=c(10,50,100,150,250,500,1000)
thresh=22

for(thresh in c(22,10,5,2))
{
n=1

matchL<-array(NaN,c(12,4,7))
dimnames(matchL)[[1]]<-rep("aaaa",12)
dimnames(matchL)[[2]]<-c("LAP 150km","LAP 50km","PG 150km","PG 50km")
dimnames(matchL)[[3]]<-paste(res,"km")
matchE<-matchT<-matchL

for(i in 1:4)
  for(j in 1:3)
  {    
    dimnames(matchL)[[1]][n]<-dimnames(matchE)[[1]][n]<-paste(cmip[i],wrf[j])
    
    filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
    
    filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_9009.csv",sep=""),
                paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_9009.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.8.csv",sep=""),
                paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg1.3.csv",sep=""))
    
    for(ff in c(1,3))
    {
      data=read.csv(filelistE[ff])
      b=order(data[,intcolE[ff]],decreasing=T)
      thresh2=data[b[20*thresh],intcolE[ff]]
      if(is.na(thresh2)) thresh2=min(data[,intcolE[ff]])
      IDs=data$ID[data[,intcolE[ff]]>=thresh2]    
      
      data1=read.csv(filelist1[ff])
      I=which(data1$ID%in%IDs  & data1[,intcol[ff]]>=thresh2 )
      data1=data1[I,]
 
      data=read.csv(filelistE[ff+1])
      b=order(data[,intcolE[ff+1]],decreasing=T)
      thresh2=data[b[20*thresh],intcolE[ff+1]]
      if(is.na(thresh2)) thresh2=min(data[,intcolE[ff+1]])
      IDs=data$ID[data[,intcolE[ff+1]]>=thresh2]    
      
      data2=read.csv(filelist1[ff+1])
      I=which(data2$ID%in%IDs  & data2[,intcol[ff]]>=thresh2 )
      data2=data2[I,]
      
      a=matchlows(data1,data2,res)
      matchL[n,ff,]=a[,2]
      matchT[n,ff,]=a[,3]
      matchE[n,ff,]=a[,4]
      a=matchlows(data2,data1,res)
      matchL[n,ff+1,]=a[,2]
      matchT[n,ff+1,]=a[,3]
      matchE[n,ff+1,]=a[,4]
      
    }      
    n=n+1  
  }


matchE2=abind(apply(matchE[,1:2,],c(1,3),mean),apply(matchE[,3:4,],c(1,3),mean),along=3)*100
matchL2=abind(apply(matchL[,1:2,],c(1,3),mean),apply(matchL[,3:4,],c(1,3),mean),along=3)*100

data2=abind(matchL2[,4,],matchE2[,4,],along=3)
dimnames(data2)[[2]]=c("LAP","PG")
dimnames(data2)[[3]]=c("Lows","Events")
names(dimnames(data2))=c("Source","Method","By")
data2=melt(data2)

pdf(paste("ECLeventsALL_thresh",thresh,"_ggplot2_eventmatch_150km.pdf",sep=""),height=5,width=5)
ggplot(data2, aes(x = By, y = value, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))) + 
  scale_y_continuous(breaks=seq(0, 100, 20),limits = c(0, 100)) +
  theme_bw() + ylab("Percentage matched") + xlab("") +  geom_hline(yintercept = 0)
dev.off()
}


matchlows<-function(d1,d2,res)
{
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  d1=d1[y1>=ymin & y1<=ymax & d1$Location==1,]
  d2=d2[y2>=ymin & y2<=ymax,]
  ll=length(d1[,1])
  match=matrix(NaN,ll)
  
  for(i in 1:length(d1[,1])) 
  {
    I=which(d2$Date==d1$Date[i] & d2$Time==d1$Time[i])
    if(length(I)>0) match[i]=min(spDistsN1(as.matrix(cbind(d2$Lon[I],d2$Lat[I])),as.numeric(c(d1$Lon[i],d1$Lat[i])),longlat=TRUE))
  } 
  
  matchev=cbind(unique(d1$ID),rep(NaN,length(unique(d1$ID))))
  for(i in 1:length(matchev[,1]))
  {
    I=which(d1$ID==matchev[i,1])
    matchev[i,2]=min(match[I],na.rm=T)
  }
                
  
  res=cbind(res,matrix(0,length(res),3))
  for(i in 1:length(res[,1])){
    res[i,2]=length(which(match<=res[i,1]))/ll
    res[i,3]=length(which(match<=res[i,1]))/sum(1-is.na(match))
    res[i,4]=length(which(matchev[,2]<=res[i,1]))/length(matchev[,2])
  } 
  
  return(res)
}

