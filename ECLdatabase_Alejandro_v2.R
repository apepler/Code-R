##This is code to take the output from the low tracking software 
##and turn it into a .csv file of ECL track data.

##Note that this script performs only the initial filtering 
##to restrict ECLs based on location, track length and at least one closed fix
##To obtain a better database of ECLs with impacts you can filter for 
##a fix curvature of >= 0.25 & the fix being in the ECL region

##This version of the code requires both the annual version of the data
##as well as that run for the previous summer, to collate the tracks for events
##over January 1.

##Note also that R has issues with large datasets - 
##I tend to separate into sets of 20 years or fewer. Pretty quick. 

##yearS and yearE are start & end years. Must have DJF data for yearS-1
##output is a string to help label the fix & event files (so use "name")

##You will need to set the location to where the track files are stored.
#setwd('~/Charles_bom/output')

cmip=c("ECHAM5","CSIROMk3.0","MIROC3.2","CCCMA3.1","nnrp")
cmipUM=c("echam5","csiromk3","miroc","cccma","ncep")
wrfv=c("R1","R2","R3")
period=c("present","farfuture")
periodUM=c("9009","6079")
res=c("50","150")

indir="/srv/ccrc/data25/z3444417/ECLs/Browning/4_PHASE_SPACE/v14_in-v14_out/1990-2010/PG4/WRF-N50/txt/"
outdir="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Alejandro/v2/"

for(j in 1:3)
  for(l in 1)
  {
    infile=paste(indir,"WRF-50-3-nnrp-",wrfv[j],"-present-",res[l],"-6_v14-D2-PG4.csv",sep="")        
    outfile=paste("Alejandro_ncep_wrf",wrfv[j],"_res",res[l],"_9009_PG4",sep="")
    ECLdatabase(infile,outfile)
  }


for(i in 1:4)
  for(j in 1:3)
    for(k in 1)
      for(l in 1)
      {
        infile=paste(indir,"WRF-50-3-",cmip[i],"-",wrfv[j],"-",period[k],"-",res[l],"-6_v14-D2-PG4.csv",sep="")        
        outfile=paste("Alejandro_",cmipUM[i],"_wrf",wrfv[j],"_res",res[l],"_",periodUM[k],"_PG4",sep="")
        ECLdatabase(infile,outfile)
      }

indir="/srv/ccrc/data25/z3444417/ECLs/Browning/4_PHASE_SPACE/v14_in-v14_out/2060-2080/PG4/WRF-N50/txt/"
cmip=c("ECHAM5","CSIRO-MK3.0","MIROC3.2","CCCMA3.1","nnrp")
for(i in 1:4) 
  for(j in 1:3)
    for(k in 2)
      for(l in 1)
      {
        infile=paste(indir,"WRF-50-3-",cmip[i],"-",wrfv[j],"-",period[k],"-",res[l],"-6_v14-D2-PG4.csv",sep="")        
        outfile=paste("Alejandro_",cmipUM[i],"_wrf",wrfv[j],"_res",res[l],"_",periodUM[k],"_PG4",sep="")
        ECLdatabase(infile,outfile)
      }


ECLdatabase<-function(infile,outfile,thresh=NA)
{
  print(outfile)
  data=read.csv(infile,stringsAsFactors=FALSE)
  
  ##Make similar in structure to my ECL fixes
  data2=data.frame(ID=data[,1],Fix=rep(0,length(data[,1])),Date=data$YEAR*10000+data$MONTH*100+data$DAY,Time=data$HOUR,
                   Lon=data[,8],Lat=data[,7],MSLP=data[,9],PG200=data$PG200,PG500=data$PG500,Location=rep(0,length(data[,1])),
                   Rad=as.numeric(data$SIZE_MEAN_RUD),Lap=as.numeric(data$LAP_MEAN),Depth=as.numeric(data$DEPTH_MEAN),
                   B=data$B500,VL=data$VL500,VU=data$VU500)
  
  data2[1,2]=1
  for(i in 2:length(data2$ID))
  if(data2$ID[i]==data2$ID[i-1]) data2$Fix[i]=data2$Fix[i-1]+1 else data2$Fix[i]=1
  
  if(!is.na(thresh))
  {
    data3=data2[data2$PG200>=thresh/100,]
    data3$Fix2<-data3$ID2<-rep(0,length(data3$ID))
    data3$ID2[1]<-data3$Fix2[1]<-1
    
    for(i in 2:length(data3$ID))
      if(data3$ID[i]==data3$ID[i-1] & data3$Fix[i]==data3$Fix[i-1]+1) {
        data3$ID2[i]=data3$ID2[i-1]
        data3$Fix2[i]=data3$Fix2[i-1]+1
      } else {
        data3$ID2[i]=data3$ID2[i-1]+1
        data3$Fix2[i]=1
      }
    
    data3[,1:2]=data3[,10:11]
    data2=data3[,1:9]    
  }

  I<-which(data2$Lon>=149 & data2$Lon<=161 & data2$Lat<(-37) & data2$Lat>=-41)
  data2$Location[I]<-1
  I<-which(data2$Lon>=(149+(37+data2$Lat)/2) & data2$Lon<=161 & data2$Lat<(-31) & data2$Lat>=-37)
  data2$Location[I]<-1
  I<-which(data2$Lon>=152 & data2$Lon<=161 & data2$Lat<=(-24) & data2$Lat>=-31)
  data2$Location[I]<-1
  
  
  ##Collates information for events
  x<-rle(data2[,1])
  events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Date2=rep(0,length(x$values)),
                     MSLP=rep(0,length(x$values)),Loc=rep(0,length(x$values)),PG200=rep(0,length(x$values)),PG500=rep(0,length(x$values)),
                     Rad=rep(0,length(x$values)),Lap=rep(0,length(x$values)),Depth=rep(0,length(x$values)),
                     B=rep(0,length(x$values)),VL=rep(0,length(x$values)),VU=rep(0,length(x$values)))
  
  ##Single-fix events
  
  I<-which(events[,2]==1)
  if(length(I)>0)
  {
    y<-match(events[I,1],data2[,1])
    events=events[-I,]
    fixes=data2[-y,]
  } else fixes=data2
  
  for(i in 1:length(events[,1])) 
  {
    I<-which(fixes[,1]==events[i,1])
    events[i,3]=min(fixes[I,3]) ##Date1
    events[i,4]=max(fixes[I,3]) ##Date2
    events[i,5]=min(fixes[I,7]) ##Min central pressure
    if(sum(fixes[I,10])>0) events[i,6]=1 ##Ever in ECL region
    events[i,7]=max(fixes[I,8]) ##Max PG200
    events[i,8]=max(fixes[I,9]) ##Max PG500
    events[i,9]=mean(fixes[I,11],na.rm=T)
    events[i,10]=max(fixes[I,12],na.rm=T)
    events[i,11]=max(fixes[I,13],na.rm=T)
    events[i,12:14]=apply(fixes[I,14:16],2,mean,na.rm=T)
  }
  
  ## Take only the events where at least one fix is closed & in the ECL region
  ##And duration is at least 2 consecutive fixes
  
  I<-which(events[,6]>0)
  events2<-events[I,-6]
  include<-match(fixes[,1],events2[,1])
  J<-which(is.na(include)==0)
  fixes2=fixes[J,]

  ##Fix numbering and save final data
  
  eventnum=events2[,1]  
  for(i in 1:length(eventnum))
  {
    I<-which(fixes2[,1]==eventnum[i])
    fixes2[I,1]=i
    events2[i,1]=i
  }

  ##Add extra columns - count,mslp,cv,rad in loc
  events2$Depth2<-events2$Lap2<-events2$Rad2<-events2$PG500a<-events2$PG200a<-events2$MSLP2<-events2$Length2<-rep(0,length(events2$ID))
  for(i in 1:length(events2[,1]))
  {
    I<-which(fixes2$ID==events2[i,1] & fixes2$Location==1)
    events2$Length2[i]=length(I)
    events2[i,15:17]=apply(fixes2[I,7:9],2,max)
    events2[i,18]=mean(fixes2[I,11],na.rm=T)
    events2[i,19:20]=apply(fixes2[I,12:13],2,max,na.rm=T)
  }

  outF=paste(outdir,'ECLfixes_',outfile,'.csv',sep="")
  outE=paste(outdir,'ECLevents_',outfile,'.csv',sep="")
  
  write.csv(fixes2,outF)
  write.csv(events2,outE)
}