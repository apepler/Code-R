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

cmip=c("ECHAM5","CSIROMk3.0","MIROC3.2","CCCMA3.1")
cmipUM=c("echam5","csiromk3","miroc","cccma")
wrfv=c("R1","R2","R3")
period=c("present","farfuture")
periodUM=c("9009","6079")
res=c("50","150")

indir="/srv/ccrc/data25/z3444417/ECLs/Browning/4_FilterEvents/v14_in-v14_out/1990-2010/PG1/WRF-50/txt/"
outdir="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Alejandro/"

for(j in 1:3)
  for(l in 1:2)
  {
    infile=paste(indir,"WRF-50-3-nnrp-",wrfv[j],"-present-",res[l],"-6_v14-D1-PG1.csv",sep="")        
    outfile=paste("Alejandro_ncep_wrf",wrfv[j],"_res",res[l],"_9009_pg0.5",sep="")
    ECLdatabase(infile,outfile)
  }


for(i in 1:4)
  for(j in 1:3)
    for(k in 1:2)
      for(l in 1:2)
      {
        infile=paste(indir,"WRF-50-3-",cmip[i],"-",wrfv[j],"-",period[k],"-",res[l],"-6_v14-D1-PG1.csv",sep="")        
        outfile=paste("Alejandro_",cmipUM[i],"_wrf",wrfv[j],"_res",res[l],"_",periodUM[k],"_pg0.5",sep="")
        ECLdatabase(infile,outfile)
      }

i=3
j=2
k=1
l=2
for(t in seq(0.6,1.4,0.1))
{
  infile=paste(indir,"WRF-50-3-",cmip[i],"-",wrfv[j],"-",period[k],"-",res[l],"-6_v14-D1-PG2.csv",sep="")        
  outfile=paste("Alejandro_",cmipUM[i],"_wrf",wrfv[j],"_res",res[l],"_",periodUM[k],"_PG2_",t,sep="")
  ECLdatabase(infile,outfile,t)
}

count=data.frame(Thresh=seq(0.6,1.4,0.1),Count2=rep(0,9),Count3=rep(0,9))
for(t in 1:9)
{
  a=read.csv(paste(outdir,"ECLevents_Alejandro_",cmipUM[i],"_wrf",wrfv[j],"_res",res[l],"_",periodUM[k],"_PG2_",count[t,1],'.csv',sep=""))
  count[t,2]=length(a[,1])/20
#  a=read.csv(paste(outdir,"ECLevents_Alejandro_",cmipUM[i],"_wrf",wrfv[j],"_res",res[l],"_",periodUM[k],"_PG3_",count[t,1],'.csv',sep=""))
#  count[t,3]=length(a[,1])/20
}

indir="/srv/ccrc/data25/z3444417/ECLs/Browning/3_FILTERED_EVENTS/v14_in-v14_out/1980-2010/PG2/WRF-N50/txt/"
outdir="/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/Alejandro/"

names1=c("CFSR","JRA55","MERRA","ERAI")
names2=c("CFSR-50-1","JRA55-50-6","MERRA-50-1","ERAI-75-6")
for(j in 4)
  for(l in 1:2)
  {
    infile=paste(indir,names2[j],"-global-",res[l],"-6_v14-D2-PG2.csv",sep="")        
    outfile=paste("Alejandro_",names1[j],"_res",res[l],"_9009_PG2",sep="")
    ECLdatabase(infile,outfile)
  }



ECLdatabase<-function(infile,outfile,thresh=NA)
{
  print(outfile)
  data=read.csv(infile)
  
  ##Make similar in structure to my ECL fixes
  data2=data.frame(ID=data[,1],Fix=rep(0,length(data[,1])),Date=data$YEAR*10000+data$MONTH*100+data$DAY,Time=data$HOUR,
                   Lon=data$LON_CEN,Lat=data$LAT_CEN,MSLP=data$P_CEN,PG200=data$PG200,Location=rep(0,length(data[,1])))
  
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
                     MSLP=rep(0,length(x$values)),PG=rep(0,length(x$values)),Loc=rep(0,length(x$values)))
  
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
    events[i,6]=max(fixes[I,8]) ##Max PG
    if(sum(fixes[I,9])>0) events[i,7]=1 ##Ever in ECL region
  }
  
  ## Take only the events where at least one fix is closed & in the ECL region
  ##And duration is at least 2 consecutive fixes
  
  I<-which(events[,7]>0)
  events2<-events[I,1:6]
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
  events2$PG2<-events2$MSLP2<-events2$Length2<-rep(0,length(events2$ID))
  for(i in 1:length(events2[,1]))
  {
    I<-which(fixes2$ID==events2[i,1] & fixes2$Location==1)
    events2$Length2[i]=length(I)
    events2$MSLP2[i]=min(fixes2$MSLP[I])
    events2$PG2[i]=max(fixes2$PG200[I])
  }

  outF=paste(outdir,'ECLfixes_',outfile,'.csv',sep="")
  outE=paste(outdir,'ECLevents_',outfile,'.csv',sep="")
  
  write.csv(fixes2,outF)
  write.csv(events2,outE)
}