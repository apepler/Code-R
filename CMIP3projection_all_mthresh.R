rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
mldb<-read.csv("~/Documents/ECLs/mldb_events.csv")
mthresh=rep(0,12)
for(i in 1:12) mthresh[i]=length(which(mldb$Month==i))*20/37

mthresh2=mthresh
for(i in 2:11) mthresh2[i]=length(which(mldb$Month>=(i-1) & mldb$Month<=(i+1)))*20/(37*3)
mthresh2[1]=length(which(mldb$Month<=2 | mldb$Month>=12))*20/(37*3)
mthresh2[12]=length(which(mldb$Month<=1 | mldb$Month>=11))*20/(37*3)
plot(1:12,mthresh,type="l")
lines(1:12,mthresh2,col=2) ## Smoothed out

##Alternately, bootstrap the years? Comes out the same!!

## This, of course, assumes 22 events p.a. So maybe better to only take the "significant" ones? No, then diff season
## So just use a scaled-threshold for ease?

mthresh=round(mthresh*0.7)
mthresh2=round(mthresh2*0.7)

now<-array(0,c(12,12,5))
dimnames(now)[[1]]<-rep("aaa",12)
dimnames(now)[[2]]=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dimnames(now)[[3]]=c("GV","UM 150km","UM 50km","PG 150km","PG 50km")
future2<-now2<-future<-now
intcol=c(8,10,10,10,10)
cmip=c("echam5","csiromk3","miroc","cccma")
wrf=c("R1","R2","R3")

n=1
  for(i in 1:4)
    for(j in 1:3)
    {
      dimnames(now)[[1]][n]<-dimnames(future)[[1]][n]<-dimnames(now2)[[1]][n]<-dimnames(future2)[[1]][n]<-paste(cmip[i],wrf[j])
      
      filelist1=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_9009.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_9009.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_9009.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_9009_pg0.6.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_9009_pg0.95.csv",sep=""))
      filelist2=c(paste("Fei/GVevents_",cmip[i],"_",wrf[j],"_6079.csv",sep=""),
                  paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv06_6079.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1_6079.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_6079_pg0.6.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_6079_pg0.95.csv",sep=""))
      
      for(ff in 1:5)
      {
        data=read.csv(filelist1[ff])
        mm=floor((data$Date1%%10000)/100)
        thresh<-thresh2<-rep(0,12)
        for(m in 1:12)
        {
          data2=data[mm==m,]
          b=order(data2[,intcol[ff]],decreasing=T)
          thresh[m]=data2[b[mthresh[m]],intcol[ff]]
          if(is.na(thresh[m])) thresh[m]=data2[b[length(b)],intcol[ff]] 
          now[n,m,ff]=length(which(data2[,intcol[ff]]>=thresh[m]))
          
          thresh2[m]=data2[b[mthresh2[m]],intcol[ff]]
          if(is.na(thresh2[m])) thresh2[m]=data2[b[length(b)],intcol[ff]] 
          now2[n,m,ff]=length(which(data2[,intcol[ff]]>=thresh2[m]))
        }
        
        data=read.csv(filelist2[ff])
        mm=floor((data$Date1%%10000)/100)
        for(m in 1:12)
        {
          data2=data[mm==m,]
          future[n,m,ff]=length(which(data2[,intcol[ff]]>=thresh[m]))
          future2[n,m,ff]=length(which(data2[,intcol[ff]]>=thresh2[m]))
        }        
      }
      n=n+1  
    }

change=100*(apply(future,c(1,3),sum)-apply(now,c(1,3),sum))/apply(now,c(1,3),sum)
boxplot(change)
abline(h=0,col="red")
change2=100*(apply(future2,c(1,3),sum)-apply(now2,c(1,3),sum))/apply(now2,c(1,3),sum) #Basically the same

##Seasonal using the varying thresh?
change=100*(apply(future[,c(1:4,11:12),],c(1,3),sum)-apply(now[,c(1:4,11:12),],c(1,3),sum))/apply(now[,c(1:4,11:12),],c(1,3),sum)
boxplot(change)
abline(h=0,col="red")
