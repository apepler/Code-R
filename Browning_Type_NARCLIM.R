rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library("R.matlab")
library(RNetCDF)
library(abind)
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

cmip=c("echam5","csiromk3","miroc","cccma","ncep")
cmip2=c("ECHAM5","CSIROMK3","MIROC","CCCMA")
wrf=c("R1","R2","R3")
names=c("LAP 150km","LAP 50km","PG 150km","PG 50km")
time=c(9009,6079)
thresh=22
intcol=c(10,10,9,9)
n=1

a=open.nc("/srv/ccrc/data34/z3478332/WRF/MIROC/WRF_slp_R1_199001_regrid.nc")
lat=var.get.nc(a,"lat0")
lon=var.get.nc(a,"lon0")
lat2=rev(seq(-1.5,-48.5,-1.5))
lon2=rev(seq(178.5,106,-1.5))

DS=8
DL=4

threshholds=matrix(0,5,3)

for(t in 1:2)
  for(i in 1:5)
    for(j in 1:3)
      
    {    
      print(paste(cmip[i],wrf[j],time[t]))
      filelistE=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_",time[t],".csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_",time[t],".csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_",time[t],"_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_",time[t],"_pg1.3.csv",sep=""))
      
      filelist1=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_",time[t],".csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLfixes_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_",time[t],".csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_",time[t],"_pg0.8.csv",sep=""),
                  paste("Alejandro/ECLfixes_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_",time[t],"_pg1.3.csv",sep=""))
      filelistT=c(paste("outputUM/proj100/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv06/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj100_rad2cv0.8_",time[t],"_type.csv",sep=""),
                  paste("outputUM/proj240/outputUM_",cmip[i],"_WRF",wrf[j],"_50_rad2cv1/ECLevents_umelb_",cmip[i],"_wrf",wrf[j],"_proj240_rad2cv1.35_",time[t],"_type.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res150_",time[t],"_pg0.8_type.csv",sep=""),
                  paste("Alejandro/ECLevents_Alejandro_",cmip[i],"_wrf",wrf[j],"_res50_",time[t],"_pg1.3_type.csv",sep=""))
      
      
      for(ff in 1:4)
      {
        data=read.csv(filelist1[ff])      
        events=read.csv(filelistE[ff])     
        
        if(t==1) 
        {  
          b=order(events[,intcol[ff]],decreasing=T)
          threshholds[i,j]=events[b[20*thresh],intcol[ff]]
        }
        IDs=events$ID[which(events[,intcol[ff]]>=threshholds[i,j])]
        
        x=1
        I=which(data$ID==IDs[x] & data$Location==1)
        firstfixes=data[I[1],]
        for(x in 2:length(IDs))
        {
          I=which(data$ID==IDs[x] & data$Location==1)
          firstfixes=rbind(firstfixes,firstfixes=data[I[1],])
        }
        
        if(t==1)
        {
          I=which(firstfixes$Date>=20100000)
          if(length(I)>0) firstfixes$Date[I]=firstfixes$Date[I]-200000
        }
        
        
        Date2=data.frame(Year=floor(firstfixes$Date/10000),Month=floor(firstfixes$Date/100)%%100,Day=firstfixes$Date%%100)
        if(ff<3) Date2$Hour=(as.numeric(firstfixes$Time)-1)*6 else Date2$Hour=firstfixes$Time
        Date2$Date=as.POSIXct(as.character(firstfixes$Date*100+Date2$Hour),tz="GMT",format="%Y%m%d%H")
        Types=data.frame(ID=firstfixes$ID,Fix=firstfixes$Fix,Date=firstfixes$Date,Time=firstfixes$Time,Type=rep("aaa",length(IDs)),stringsAsFactors=F)
        
        backtrack=array(0,c(length(firstfixes[,1]),DS+1,2))
        dimnames(backtrack)[[3]]=c("Lon","Lat")
        backtrack[,1,]=as.matrix(round(c(firstfixes$Lon,firstfixes$Lat)/1.5)*1.5)
        
        for(x in 1:length(IDs))
        {
          if(Date2$Day[x]>=3)
          {
            if(i<5) fname=paste("/srv/ccrc/data34/z3478332/WRF/",cmip2[i],"/WRF_slp_R",j,"_",Date2$Year[x],sprintf("%2.2d",Date2$Month[x]),"_regrid.nc",sep="") else
              fname=paste("/srv/ccrc/data34/z3478332/WRF/R",j,"/WRF_slp_R",j,"_",Date2$Year[x],sprintf("%2.2d",Date2$Month[x]),"_regrid.nc",sep="")
            
            if(Date2$Month[x]<12) wtime=seq.POSIXt(ISOdate(Date2$Year[x],Date2$Month[x],1,0),ISOdate(Date2$Year[x],(Date2$Month[x]+1),1,0),by="6 hours") else
              wtime=seq.POSIXt(ISOdate(Date2$Year[x],Date2$Month[x],1,0),ISOdate((Date2$Year[x]+1),1,1,0),by="6 hours")
            
            a=open.nc(fname)
            I=which(wtime==Date2$Date[x])
            slp=var.get.nc(a,"slp0",start=c(1,1,I-8),count=c(length(lon),length(lat),8))
            close.nc(a)
            
          } else {
            if(i<5) {
              fname2=paste("/srv/ccrc/data34/z3478332/WRF/",cmip2[i],"/WRF_slp_R",j,"_",Date2$Year[x],sprintf("%2.2d",Date2$Month[x]),"_regrid.nc",sep="") 
              if(Date2$Month[x]>1) fname1=paste("/srv/ccrc/data34/z3478332/WRF/",cmip2[i],"/WRF_slp_R",j,"_",Date2$Year[x],sprintf("%2.2d",Date2$Month[x]-1),"_regrid.nc",sep="") else
                fname1=paste("/srv/ccrc/data34/z3478332/WRF/",cmip2[i],"/WRF_slp_R",j,"_",Date2$Year[x]-1,"12_regrid.nc",sep="")
            } else {
              fname2=paste("/srv/ccrc/data34/z3478332/WRF/R",j,"/WRF_slp_R",j,"_",Date2$Year[x],sprintf("%2.2d",Date2$Month[x]),"_regrid.nc",sep="")
              if(Date2$Month[x]>1) fname1=paste("/srv/ccrc/data34/z3478332/WRF/R",j,"/WRF_slp_R",j,"_",Date2$Year[x],sprintf("%2.2d",Date2$Month[x]-1),"_regrid.nc",sep="") else
                fname1=paste("/srv/ccrc/data34/z3478332/WRF/R",j,"/WRF_slp_R",j,"_",Date2$Year[x]-1,"12_regrid.nc",sep="")
            }
            
            if(Date2$Month[x]<12) wtime2=seq.POSIXt(ISOdate(Date2$Year[x],Date2$Month[x],1,0),ISOdate(Date2$Year[x],(Date2$Month[x]+1),1,0),by="6 hours") else
              wtime2=seq.POSIXt(ISOdate(Date2$Year[x],Date2$Month[x],1,0),ISOdate((Date2$Year[x]+1),1,1,0),by="6 hours")
            
            if(Date2$Month[x]>1) wtime1=seq.POSIXt(ISOdate(Date2$Year[x],Date2$Month[x]-1,1,0),ISOdate(Date2$Year[x],(Date2$Month[x]),1,0),by="6 hours") else
              wtime1=seq.POSIXt(ISOdate(Date2$Year[x]-1,12,1,0),ISOdate((Date2$Year[x]),1,1,0),by="6 hours")
            
            I=which(wtime2==Date2$Date[x])
            a=open.nc(fname1)
            b=open.nc(fname2)
            slp=abind(var.get.nc(a,"slp0",start=c(1,1,length(wtime1)+I-9),count=c(length(lon),length(lat),8-I)),
                      var.get.nc(b,"slp0",start=c(1,1,1),count=c(length(lon),length(lat),I)),along=3)
            close.nc(a)
            close.nc(b)
            
          }
          
          ## Regrid as lazy average
          slp2=array(NaN,c(length(lon2),length(lat2),8))
          for(xx in 1:length(lon2))
            for(yy in 1:length(lat2))
            {
              I=which(lon==lon2[xx])
              J=which(lat==lat2[yy])
              slp2[xx,yy,]=apply(slp[(I-1):(I+1),(J-1):(J+1),],3,mean,na.rm=T)
            }
          
          ### Do actual backtracking
          
          for(dd in 1:DS)
          {
            J=which(lon2==backtrack[x,dd,1])
            K=which(lat2==backtrack[x,dd,2])
            
            n=1
            tmp=array(0,c((DL*2+1)^2,5))
            for(p in seq((J-DL),(J+DL)))
              for(q in seq((K-DL),(K+DL)))
              {    
                if(p>1 & q>1 & p<length(lon2) & q<length(lat2))
                {
                  CELL = slp2[p,q,9-dd];
                  
                  # MEASURE PRESSURE GRADIENTS AROUND TEST CELL
                  test = c(slp2[p-1,q,9-dd]-CELL, slp2[p-1,q+1,9-dd]-CELL, slp2[p-1,q-1,9-dd]-CELL, 
                           slp2[p+1,q,9-dd]-CELL, slp2[p+1,q+1,9-dd]-CELL, slp2[p+1,q-1,9-dd]-CELL, 
                           slp2[p,q-1,9-dd]-CELL, slp2[p,q+1,9-dd]-CELL)
                  
                  tmp[n,] = c(length(which(test>0)),max(test),lon2[p],lat2[q],CELL)
                }
                n=n+1
              }  
            a=order(-tmp[,1],tmp[,5])
            backtrack[x,dd+1,]=tmp[a[1],3:4]
          }
        }
        
        motion=array(0,c(length(IDs),8,3))
        for(x in 1:8) 
        {
          motion[,x,1:2]=backtrack[,x,]-backtrack[,x+1,]
          
          ##I'll use my lazy GDR boundary for ECLs from before!
          I<-which(backtrack[,x+1,1]<(151-(25+backtrack[,x+1,2])) & backtrack[,x+1,2]>(-25) & backtrack[,x+1,2]<=-10) ## Over more land somehow?
          if(length(I)>0) motion[I,x,3]<-1
          I<-which(backtrack[,x+1,1]<151 & backtrack[,x+1,2]<=(-25) & backtrack[,x+1,2]>=-31)
          if(length(I)>0) motion[I,x,3]<-1
          I<-which(backtrack[,x+1,2]<(-31) & backtrack[,x+1,2]>=-39 & backtrack[,x+1,1]<(148+(37+backtrack[,x+1,2])/2))
          if(length(I)>0) motion[I,x,3]<-1
        }
        
        typing=cbind(apply(motion[,,1:2]*(motion[,,1:2]>1),c(1,3),sum),apply(motion[,,1:2]*(motion[,,1:2]<(-1)),c(1,3),sum)*-1,apply(motion[,,3],1,sum)/DS,matrix(0,length(IDs),2))
        colnames(typing)=c("E","N","W","S","Land prop","N points","S points")
        for(x in 1:length(IDs)){
          typing[x,6]=length(which(backtrack[x,2:9,2]>=-27))
          typing[x,7]=length(which(backtrack[x,2:9,2]<=-39))
        } 
        
        I=which(typing[,5]<=0.5 & (typing[,4]>typing[,2] | typing[,7]<2) )
        Types$Type[I]="ET"
        I=which(typing[,5]<=0.5 & typing[,2]>=typing[,4] & typing[,7]>=2)
        Types$Type[I]="SSL"
        I=which(typing[,5]>0.5 & typing[,6]>=2)
        Types$Type[I]="IT"
        I=which(typing[,5]>0.5 & typing[,6]<2)
        Types$Type[I]="CL"
        
        write.csv(Types,file=filelistT[ff])
      }
    }

