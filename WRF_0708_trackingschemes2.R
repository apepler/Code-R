####### Okay, taking bits from existing code for ECL typing.
### First, we need to load in wrf and restrict to just one year where a range of MLDB ECLs
rm(list=ls())
library(R.matlab)
library(RNetCDF)
library(abind)
library(sp)
library(geosphere)
library(abind)

setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/")

wrfdirs=c("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007/out/slp/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007/out/slp/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007/out/slp/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007_notopo/out/slp/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007_notopo/out/slp/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007_notopo/out/slp/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN/out/slp/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN/out/slp/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN/out/slp/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN_noeac/out/slp/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_noeac/out/slp/",
          "/srv/ccrc/data37/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN_noeac/out/slp/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007_BRAN_2eac/out/slp/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default_2007_BRAN_2eac/out/slp/",
          "/srv/ccrc/data36/z3478332/WRF/output/ERAI_R3_nudging_default_2007_BRAN_2eac/out/slp/")

dom="d01"
cat="rad2_p100"

wnames=c("R1","R2","R3",
         "R1_notopo","R2_notopo","R3_notopo",
         "R1_BRAN","R2_BRAN","R3_BRAN",
         "R1_BRAN_noeac","R2_BRAN_noeac","R3_BRAN_noeac",
         "R2_BRAN_2eac","R2_BRAN_2eac","R3_BRAN_2eac")

### This is my explicit mask to get rid of the places where topo was doing weird things.
##  It's messy, but it works for now

w=1
file=open.nc(paste(wrfdirs[w],"WRF_d01_gph_0708_regrid.nc",sep=""))
lat=var.get.nc(file,"lat0")
lon=var.get.nc(file,"lon0")
lev<-var.get.nc(file,'lev_2')
time1=seq.POSIXt(as.POSIXct("2007010100",format="%Y%m%d%H",tz="GMT"),as.POSIXct("2008123118",format="%Y%m%d%H",tz="GMT"),21600)
hgt=var.get.nc(file,'Z',start=c(1,1,9,1),count=c(length(lon),length(lat),5,30))
diff=hgt[,,1,]-hgt[,,5,]

mask=matrix(1,length(lon),length(lat))
mask[(apply(diff,c(1,2),max)>2220 | is.na(apply(diff,c(1,2),max)))]=NaN

mask[83:85,17]=NaN
mask[84:85,26:27]=NaN
mask[86,26:28]=NaN
mask[89,27]=NaN
mask[90,28:30]=NaN
mask[90:91,33:35]=NaN
mask[93:94,38:40]=NaN
mask[86,82]=NaN
mask[70,94]=NaN
mask[128,10]=NaN
mask[126,10:11]=NaN
mask[126:132,12]=NaN
mask[129,13]=NaN
mask[135,14]=NaN
mask[133:136,15]=NaN
mask[134:138,16]=NaN
mask[135:138,17]=NaN
mask[136:138,18]=NaN
mask[142:143,22]=NaN
mask[143,22]=NaN

rm(hgt)
rm(diff)
close.nc(file)

for(w in 1:length(wnames))
{
  fname=c(paste("ECLfixes_",dom,"_2007_",wnames[w],"_",cat,".csv",sep=""),
          paste("ECLfixes_",dom,"_2008_",wnames[w],"_",cat,".csv",sep=""))
  wrf=rbind(read.csv(fname[1],stringsAsFactors=F),read.csv(fname[2],stringsAsFactors=F))
  wrf$ID=wrf$ID+1000*floor(wrf$Date/10000)
  wrf=wrf[,2:14]
  
  
  fname=c(paste("ECLevents_",dom,"_2007_",wnames[w],"_",cat,".csv",sep=""),
          paste("ECLevents_",dom,"_2008_",wnames[w],"_",cat,".csv",sep=""))
  wrf_E=rbind(read.csv(fname[1],stringsAsFactors=F),read.csv(fname[2],stringsAsFactors=F))
  wrf_E$ID=wrf_E$ID+1000*floor(wrf_E$Date1/10000)
  wrf_E=wrf_E[,2:11]
  
  I=which(!(wrf$ID %in% wrf_E$ID))
  if(length(I)>0) print("IDs are strange - overlapping years")
  
  # Step 1 - Calculate NDR/bomb for each event
  
  wrf$NDR=NaN
  I=which(wrf$ID[2:length(wrf[,1])]==wrf$ID[1:length(wrf[,1])-1])+1
  wrf$NDR[I]= (wrf$MSLP[I]-wrf$MSLP[I-1])*sin(60*pi/180)/(6*sin(wrf$Lat[I]*pi/180))
  
  wrf_E$Bomb=0
  for(i in 1:length(wrf_E))
  {
    I=which(wrf$ID==wrf_E$ID[i] & wrf$Location==1 & !is.na(wrf$NDR))
    if(length(I)>0) if(max(wrf$NDR[I],na.rm=T)>=1) wrf_E$Bomb[i]=1
  }
  
  # Step 2: Identify the event as entered or formed
  #1: Formed in region 
  #2: Entered & didn't intensify
  #3: Entered and intensified - look at NDR & cv?
  #4: Formed in region, but may need to tie to a different event
  
  wrf$Date2=as.POSIXct(paste(as.character(wrf$Date),substr(wrf$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  
  wrf_E$EnteredFormed=0
  
  for(i in 1:length(wrf_E[,1]))
  {
    tmp=wrf[which(wrf$ID==wrf_E$ID[i]),]
    
    if(tmp$Location[1]==1)
    {
      I=which(wrf$ID!=wrf_E$ID[i] & (wrf$Date2>=tmp$Date2[1]-43200 & wrf$Date2<=tmp$Date2[1]+43200))
      if(length(I)>1) {
        dist=spDistsN1(as.matrix(wrf[I,6:7]),as.numeric(tmp[1,6:7]),longlat=T)
        if(min(dist)<500) 
        {
          k=which(wrf_E$ID==min(wrf$ID[I]))
          wrf_E$EnteredFormed[i]=wrf_E$EnteredFormed[k] } else wrf_E$EnteredFormed[i]=1
      }  else wrf_E$EnteredFormed[i]=1
    } else {
      J=which(tmp$Location==1)
      c2=max(tmp$CV[J])
      n2=max(tmp$NDR[J])
      
      K=max(1,J[1]-3):(J[1]-1)
      c1=max(tmp$CV[K])
      
      if(n2>=1 | c2>=c1+0.1)  wrf_E$EnteredFormed[i]=3 else wrf_E$EnteredFormed[i]=2 
    }
  }
  
  ### Next set: Browning classification
  
  #Set days to backtrack = 2 (i.e. 8 timesteps)
  DS=8
  #Set region for tracking - ~ surrounding 6 degrees, so +- 4 1.5 deg cells
  DL=4
  
  IDs=unique(wrf$ID)
  
  i=1
  I=which(wrf$ID==IDs[i] & wrf$Location==1)
  firstfixes=wrf[I[1],]
  for(i in 2:length(IDs))
  {
    I=which(wrf$ID==IDs[i] & wrf$Location==1)
    firstfixes=rbind(firstfixes,wrf[I[1],])
  }
  backtrack=array(NaN,c(length(firstfixes[,1]),DS+1,2))
  dimnames(backtrack)[[3]]=c("Lon","Lat")
  backtrack[,1,]=as.matrix(round(firstfixes[,6:7]/1.5)*1.5)
  
  ### Okay, now set up to loop through ALL
  ##  Easiest to just load the whole slp dataset?
  
  file=open.nc(paste(wrfdirs[w],"WRF_d01_0708_regrid.nc",sep=""))
  lat=var.get.nc(file,"lat0")
  lon=var.get.nc(file,"lon0")
  dates=seq.POSIXt(as.POSIXct("2007010100",format="%Y%m%d%H",tz="GMT"),as.POSIXct("2008123118",format="%Y%m%d%H",tz="GMT"),21600)
  
  ####Do stuart's lazy backtracking
  
  for(i in 1:length(firstfixes[,1]))
  {
    I=which(dates==firstfixes$Date2[i])
    if(I>8) 
    {
      slp=var.get.nc(file,"slp0",start=c(1,1,I-8),count=c(length(lon),length(lat),8),unpack=T) 
      
      for(j in 1:DS)
      {
        J=which(lon==backtrack[i,j,1])
        K=which(lat==backtrack[i,j,2])
        
        if(!is.na(backtrack[i,j,1]))
          if((J-DL)>=1 & (J+DL)<=length(lon) & (K-DL)>=1 & (K+DL)<=length(lat))
          {
            n=1
            tmp=array(0,c((DL*2+1)^2,5))
            for(p in (J-DL):(J+DL))
              for(q in (K-DL):(K+DL))
              {    
                if(p>2 & q>2 & p<length(lon)-1 & q<length(lat)-1)
                {
                  CELL = slp[p,q,9-j]
                  
                  # MEASURE PRESSURE GRADIENTS AROUND TEST CELL
                  test = c(slp[p-2,q,9-j]-CELL, slp[p-2,q+2,9-j]-CELL, slp[p-2,q-2,9-j]-CELL, 
                           slp[p+2,q,9-j]-CELL, slp[p+2,q+2,9-j]-CELL, slp[p+2,q-2,9-j]-CELL, 
                           slp[p,q-2,9-j]-CELL, slp[p,q+2,9-j]-CELL)
                  
                  tmp[n,] = c(length(which(test>0)),max(test),lon[p],lat[q],CELL)
                }
                n=n+1
              }
            a=order(-tmp[,1],tmp[,5])
            backtrack[i,j+1,]=tmp[a[1],3:4]
          }
      }
    }  else {
      slp=var.get.nc(file,"slp0",start=c(1,1,1),count=c(length(lon),length(lat),(I-1)),unpack=T)
      for(j in 1:(I-1))
      {
        J=which(lon==backtrack[i,j,1])
        K=which(lat==backtrack[i,j,2])
        
        if(!is.na(backtrack[i,j,1]))
          if((J-DL)>=1 & (J+DL)<=length(lon) & (K-DL)>=1 & (K+DL)<=length(lat))
          {
            n=1
            tmp=array(0,c((DL*2+1)^2,5))
            for(p in (J-DL):(J+DL))
              for(q in (K-DL):(K+DL))
              {    
                if(p>2 & q>2 & p<length(lon)-1 & q<length(lat)-1)
                {
                  CELL = slp[p,q,9-j];
                  
                  # MEASURE PRESSURE GRADIENTS AROUND TEST CELL
                  test = c(slp[p-2,q,I-j]-CELL, slp[p-2,q+2,I-j]-CELL, slp[p-2,q-2,I-j]-CELL, 
                           slp[p+2,q,I-j]-CELL, slp[p+2,q+2,I-j]-CELL, slp[p+2,q-2,I-j]-CELL, 
                           slp[p,q-2,I-j]-CELL, slp[p,q+2,I-j]-CELL)
                  
                  tmp[n,] = c(length(which(test>0)),max(test),lon[p],lat[q],CELL)
                }
                n=n+1
              }  
            a=order(-tmp[,1],tmp[,5])
            backtrack[i,j+1,]=tmp[a[1],3:4]
          }
      }
    }
    
  }
  close.nc(file)
  ##### Results are still fairly messy, but let's try to apply Stuart's cats anyway
  ## Need to know:
  ## Whether primary movement is from N, S or W
  ## If (mostly) east or west of GDR
  
  motion=array(0,c(length(IDs),8,3))
  for(i in 1:8) 
  {
    motion[,i,1:2]=backtrack[,i,]-backtrack[,i+1,]
    
    ##I'll use my lazy GDR boundary for ECLs from before!
    I<-which(backtrack[,i+1,1]<(151-(25+backtrack[,i+1,2])) & backtrack[,i+1,2]>(-25) & backtrack[,i+1,2]<=-10) ## Over more land somehow?
    if(length(I)>0) motion[I,i,3]<-1
    I<-which(backtrack[,i+1,1]<151 & backtrack[,i+1,2]<=(-25) & backtrack[,i+1,2]>=-31)
    if(length(I)>0) motion[I,i,3]<-1
    I<-which(backtrack[,i+1,2]<(-31) & backtrack[,i+1,2]>=-39 & backtrack[,i+1,1]<(148+(37+backtrack[,i+1,2])/2))
    if(length(I)>0) motion[I,i,3]<-1
  }
  
  typing=cbind(apply(motion[,,1:2]*(motion[,,1:2]>0),c(1,3),sum,na.rm=T),
               apply(motion[,,1:2]*(motion[,,1:2]<0),c(1,3),sum,na.rm=T)*-1,
               apply(motion[,,3],1,sum,na.rm=T)/DS,matrix(0,length(IDs),2))
  colnames(typing)=c("E","N","W","S","Land prop","N points","S points")
  for(i in 1:length(IDs)){
    typing[i,6]=length(which(backtrack[i,2:9,2]>=-27))
    typing[i,7]=length(which(backtrack[i,2:9,2]<=-39))
  } 
  
  Type=rep("NA",length(IDs))
  I=which(typing[,5]<=0.5 & (typing[,4]>typing[,2] | typing[,7]<2) )
  Type[I]="ET"
  I=which(typing[,5]<=0.5 & typing[,2]>=typing[,4] & typing[,7]>=2)
  Type[I]="SSL"
  I=which(typing[,5]>0.5 & typing[,6]>=2)
  Type[I]="IT"
  I=which(typing[,5]>0.5 & typing[,6]<2)
  Type[I]="CL"
  
  wrf_E$TypeSB=Type
  
  ### Next: Hart phase space
  
  file=open.nc(paste(wrfdirs[w],"WRF_d01_gph_0708_regrid.nc",sep=""))
  lat=var.get.nc(file,"lat0")
  lon=var.get.nc(file,"lon0")
  lev<-var.get.nc(file,'lev_2')
  time1=seq.POSIXt(as.POSIXct("2007010100",format="%Y%m%d%H",tz="GMT"),as.POSIXct("2008123118",format="%Y%m%d%H",tz="GMT"),21600)
  
  wrf$VL<-wrf$VU<-wrf$B<-wrf$Bearing<-NaN
  logp=log(lev)
  
  wrf$Lon[wrf$Lon>180]=wrf$Lon[wrf$Lon>180]-360
  lon[lon>180]=lon[lon>180]-360
  
  for(k in 1:length(wrf[,1]))
  {
    dist=matrix(NaN,length(lon),length(lat))
    for(i in 1:length(lon))
      for(j in 1:length(lat))
        dist[i,j]=spDistsN1(cbind(lon[i],lat[j]),cbind(wrf$Lon[k],wrf$Lat[k]))
      
      loc=which(dist==min(dist),arr.ind=T) ##I.e. the lon/lat location
      I=which(time1==wrf$Date2[k])             ##I.e. the time
      
      if(is.matrix(loc)) loc=loc[1,] ##Simplifies in cases where multiples equidistant
      
      if(loc[1]>7 & loc[2]>7 & loc[1]<=(length(lon)-7) & loc[2]<=(length(lat)-7) & min(dist)<=100)
      {
        lon2=lon[(loc[1]-7):(loc[1]+7)]
        lat2=lat[(loc[2]-7):(loc[2]+7)]
        hgt<-var.get.nc(file,'Z',start=c(loc[1]-7,loc[2]-7,1,I),count=c(15,15,13,1),unpack=T)
        for(z in 10:13) hgt[,,z]=hgt[,,z]*mask[(loc[1]-7):(loc[1]+7),(loc[2]-7):(loc[2]+7)] ##Remove dodgy points
        #image(hgt[,,13])
        
        #delZ=(apply(hgt,3,max,na.rm=T)-apply(hgt,3,min,na.rm=T))
        #dd=rep(0,12)
        #for(i in 1:12) dd[i]=(delZ[i+1]-delZ[i])/(logp[i+1]-logp[i]) 
        #wrf$VU[k]=sum(dd[1:6])
        #wrf$VL[k]=sum(dd[7:12])
        
        #delZ=(apply(hgt,3,max,na.rm=T)-apply(hgt,3,min,na.rm=T))/logp
        #wrf$VU[k]=((300-600)/6)*(delZ[1]/2 + delZ[6]/2 + sum(delZ[2:5]))
        #wrf$VL[k]=((600-900)/6)*(delZ[6]/2 + delZ[12]/2 + sum(delZ[7:11]))
        
        delZ=(apply(hgt,3,max,na.rm=T)-apply(hgt,3,min,na.rm=T))
        a=lm(delZ[1:7]~logp[1:7])
        wrf$VU[k]=a$coefficients[2]
        a=lm(delZ[7:13]~logp[7:13])
        wrf$VL[k]=a$coefficients[2]
        
        ##Bearing of low - in past 6 hours if possible, otherwise next six hours
        if(wrf$Fix[k]==1 | k==1)
          wrf$Bearing[k]=bearingRhumb(wrf[k,6:7],wrf[k+1,6:7]) else
            wrf$Bearing[k]=bearingRhumb(wrf[k-1,6:7],wrf[k,6:7])
        
        bearing2=matrix(0,15,15)
        for(i in 1:15)
          for(j in 1:15)
            bearing2[i,j]=bearingRhumb(wrf[k,6:7],c(lon2[i],lat2[j]))
        
        L=which((bearing2<wrf$Bearing[k] & bearing2>=wrf$Bearing[k]-180) | (bearing2<wrf$Bearing[k]+360 & bearing2>=wrf$Bearing[k]+180))
        R=which((bearing2>=wrf$Bearing[k]-360 & bearing2<wrf$Bearing[k]-180) | (bearing2<wrf$Bearing[k]+180 & bearing2>=wrf$Bearing[k]))
        
        Thick=hgt[,,which(lev==600)]-hgt[,,which(lev==900)]
        wrf$B[k]=mean(Thick[L],na.rm=T)-mean(Thick[R],na.rm=T)
      }
  }
  close.nc(file)
  wrf$Type="Mixed"
  wrf$Type[wrf$VL<(-50)]="EC"
  wrf$Type[wrf$VU>-10 & wrf$VL>0]="TC"
  wrf$Type[wrf$VU<(-10) & wrf$VL>0]="SC"
  
  wrf_E$HartType="Unsure"
  for(k in 1:length(wrf_E[,1]))
  {
    I=which(wrf$ID==wrf_E$ID[k] & wrf$Location==1)
    b=unique(wrf$Type[I])
    
    if(length(b)==1) wrf_E$HartType[k]=b else {
      if(length(b)==2 & "Mixed" %in% b) wrf_E$HartType[k]=b[which(b!="Mixed")] else
      {
        count=rep(0,length(b))
        for(j in 1:length(b)) count[j]=length(which(wrf$Type[I]==b[j]))
        J=which(count==max(count))
        if(length(J)==1) wrf_E$HartType[k]=b[J] else wrf_E$HartType[k]="Mixed"
      }
    }
  }
  
  write.csv(wrf,paste("typing/ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_typing.csv",sep=""))
  write.csv(wrf_E,paste("typing/ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing.csv",sep=""))
  
}
