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

ECLdatabase<-function(yearS,yearE,output)
{

  ##Files with ECL track data for required years & preceding summers.  
  year<-as.character(seq(yearS,yearE))
  fname1=paste('tracks_',year,'.dat',sep="")
  year2<-as.character(seq(yearS-1,yearE-1))
  fname2=paste('tracks_',year2,'_djf.dat',sep="")
    
  ##collates all the annual tracks into one file
  
  read.table(fname1[1], sep="",skip=1)->fixes
  I<-which(fixes[,7]<0)
  fixes<-fixes[I,]
  
  for(i in 2:length(year))
  {
    read.table(fname1[i], sep="",skip=1)->data
    I<-which(data[,7]<0)
    fixes<-rbind(fixes,data[I,])
  }
  rm(data)
  
  #Unique event numbers (starting with year)
  fixes[,1]<-fixes[,1]+floor(fixes[,3]/10000)*1000000
  
  ##Identify if within the ECL region
  
  fixes[,11]<-0
  I<-which(fixes[,6]>=149 & fixes[,6]<=161 & fixes[,7]<(-37) & fixes[,7]>=-41)
  fixes[I,11]<-1
  I<-which(fixes[,6]>=(149+(37+fixes[,7])/2) & fixes[,6]<=161 & fixes[,7]<(-31) & fixes[,7]>=-37)
  fixes[I,11]<-1
  I<-which(fixes[,6]>=152 & fixes[,6]<=161 & fixes[,7]<=(-24) & fixes[,7]>=-31)
  fixes[I,11]<-1
    
  ##Repeats for DJF
  
  read.table(fname2[1], sep="",skip=1)->fixesDJF
  I<-which(fixesDJF[,7]<0)
  fixesDJF<-fixesDJF[I,]
  
  for(i in 2:length(year2))
  {
    read.table(fname2[i], sep="",skip=1)->data
    I<-which(data[,7]<0)
    fixesDJF<-rbind(fixesDJF,data[I,])
  }
  
  fixesDJF[,1]<-fixesDJF[,1]+floor(fixesDJF[,3]/10000)*1000000
  fixesDJF[,11]<-0
  I<-which(fixesDJF[,6]>=149 & fixesDJF[,6]<=161 & fixesDJF[,7]<(-37) & fixesDJF[,7]>=-41)
  fixesDJF[I,11]<-1
  I<-which(fixesDJF[,6]>=(149+(37+fixesDJF[,7])/2) & fixesDJF[,6]<=161 & fixesDJF[,7]<(-31) & fixesDJF[,7]>=-37)
  fixesDJF[I,11]<-1
  I<-which(fixesDJF[,6]>=152 & fixesDJF[,6]<=161 & fixesDJF[,7]<=(-24) & fixesDJF[,7]>=-31)
  fixesDJF[I,11]<-1
  
  ##Collates information for events
  x<-rle(fixes[,1])
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=8))
  
  ##Single-fix events
  
  I<-which(events[,2]==1)
  y<-match(events[I,1],fixes[,1])
  events[I,3]<-events[I,4]<-fixes[y,3]
  events[I,5]<-fixes[y,8]
  events[I,6]<-fixes[y,9]
  events[I,7]<-events[I,8]<-fixes[y,2]
  events[I,9]<-fixes[y,11] 
  events[I,10]<-fixes[y,5]
  
  ##Multi-fix events
  
  J<-which(events[,2]>1)
  
  for(i in 1:length(J)) 
  {
    I<-which(fixes[,1]==events[J[i],1])
    events[J[i],3]=min(fixes[I,3]) ##Date1
    events[J[i],4]=max(fixes[I,3]) ##Date2
    events[J[i],5]=min(fixes[I,8]) ##Min central pressure
    events[J[i],6]=max(fixes[I,9]) ##Max curvature
    events[J[i],7]=min(fixes[I,2]) ##First fix # - for cross-year events
    events[J[i],8]=max(fixes[I,2]) ##Last fix #
    if(sum(fixes[I,11])>0) events[J[i],9]=1 ##Ever in ECL region
    events[J[i],10]<-min(fixes[I,5]) ##Ever closed  
  }
  
  ## Try to identify cross-year events
  ## For each event that starts higher than 1, find the first row
  ## Then identify the row in the DJF file that has the same values for columns 2-6
  ## Next find a row in the fixes file that matches the row one higher in DJF
  ## Finally, apply the eventID of this event to all of the values for the second event
  
  events2=events;
  
  I<-which(events[,7]>1 & floor(events[,3]/10000)>1950)
  ##Any case where the first fix is not fix #1
  I2=matrix(0,nrow=length(I),ncol=1)
  for(i in 1:length(I))
  {
    J<-which(fixes[,1]==events[I[i],1])
    K<-which(fixesDJF[,3]==fixes[J[1],3] & fixesDJF[,4]==fixes[J[1],4] & fixesDJF[,7]==fixes[J[1],7] & fixesDJF[,6]==fixes[J[1],6])
    L<-which(fixes[,3]==fixesDJF[K-1,3] & fixes[,4]==fixesDJF[K-1,4] & fixes[,7]==fixesDJF[K-1,7] & fixes[,6]==fixesDJF[K-1,6])
    if(length(L)==1)
    {
      ID<-fixes[L,1]
      fixes[J,1]<-ID
      I2[i]=1
      
      ##Make event1 represent the full event now
      M<-which(events[,1]==ID)
      events2[M,2]=events[M,2]+events[I[i],2]
      events2[M,4]=events[I[i],4]
      events2[M,5]=min(events[I[i],5],events[M,5])
      events2[M,6]=max(events[I[i],6],events[M,6])
      events2[M,8]=events[I[i],8]
      events2[M,9]=max(events[I[i],9],events[M,9])
      events2[M,10]=min(events[I[i],10],events[M,10])
    }
  }
  
  ##And remove the end-year events
  I3<-which(I2==1)
  events2<-events2[-I[I3],]
  rm(fixesDJF)
  
  ## Take only the events where at least one fix is closed & in the ECL region
  ##And duration is at least 2 consecutive fixes
  
  I<-which(events2[,9]>0 & events2[,2]>1 & events2[,10]==0)
  events3<-events2[I,]
  include<-match(fixes[,1],events3[,1])
  J<-which(is.na(include)==0)
  fixes2=fixes[J,]
  
  ##Fix numbering and save final data
  
  eventnum=events3[,1]
  
  for(i in 1:length(eventnum))
  {
    I<-which(fixes2[,1]==eventnum[i])
    fixes2[I,1]=i
    events3[i,1]=i
  }
  
  ##Convert to dataframe
  
  Events<-data.frame(events3[,1:6])
  names(Events)=c('ID','Length','Date1','Date2','MSLP','CV')
  
  Fixes<-data.frame(fixes2[,1:9],fixes2[,11])
  names(Fixes)=c('ID','Fix','Date','Time','Open','Lon','Lat','MSLP','CV','Location')
  
  outF=paste('ECLfixes_',output,'.csv',sep="")
  outE=paste('ECLevents_',output,'.csv',sep="")
  
  write.csv(Fixes,outF)
  write.csv(Events,outE)
}

##Gets events & fixes out of a single trackfile. Writes files to my ECLs folder. 
#Assumes file is in current working directory
ECLyear<-function(fin,output)
{
  read.table(fin, sep="",skip=1)->fixes
  
  ##Identify if within the ECL region
  
  fixes[,11]<-0
  I<-which(fixes[,6]>=149 & fixes[,6]<=161 & fixes[,7]<(-37) & fixes[,7]>=-41)
  fixes[I,11]<-1
  I<-which(fixes[,6]>=(149+(37+fixes[,7])/2) & fixes[,6]<=161 & fixes[,7]<(-31) & fixes[,7]>=-37)
  fixes[I,11]<-1
  I<-which(fixes[,6]>=152 & fixes[,6]<=161 & fixes[,7]<=(-24) & fixes[,7]>=-31)
  fixes[I,11]<-1 
  
  ##Collates information for events
  x<-rle(fixes[,1])
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=8))
  
  ##Single-fix events
  
  I<-which(events[,2]==1)
  y<-match(events[I,1],fixes[,1])
  events[I,3]<-events[I,4]<-fixes[y,3]
  events[I,5]<-fixes[y,8]
  events[I,6]<-fixes[y,9]
  events[I,7]<-events[I,8]<-fixes[y,2]
  events[I,9]<-fixes[y,11] 
  events[I,10]<-fixes[y,5]
  
  ##Multi-fix events
  
  J<-which(events[,2]>1)
  
  for(i in 1:length(J)) 
  {
    I<-which(fixes[,1]==events[J[i],1])
    events[J[i],3]=min(fixes[I,3]) ##Date1
    events[J[i],4]=max(fixes[I,3]) ##Date2
    events[J[i],5]=min(fixes[I,8]) ##Min central pressure
    events[J[i],6]=max(fixes[I,9]) ##Max curvature
    events[J[i],7]=min(fixes[I,2]) ##First fix # - for cross-year events
    events[J[i],8]=max(fixes[I,2]) ##Last fix #
    if(sum(fixes[I,11])>0) events[J[i],9]=1 ##Ever in ECL region
    events[J[i],10]<-min(fixes[I,5]) ##Ever closed  
  }
  
  ## Take only the events where at least one fix is closed & in the ECL region
  I<-which(events[,9]>0 & events[,10]==0)
  events2<-events[I,]
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
  
  ##Convert to dataframe
  
  Events<-data.frame(events2[,1:6])
  names(Events)=c('ID','Length','Date1','Date2','MSLP','CV')
  
  Fixes<-data.frame(fixes2[,1:9],fixes2[,11])
  names(Fixes)=c('ID','Fix','Date','Time','Open','Lon','Lat','MSLP','CV','Location')
  
  outF=paste('~/Documents/ECLs/ECLfixes_',output,'.csv',sep="")
  outE=paste('~/Documents/ECLs/ECLevents_',output,'.csv',sep="")
  
  write.csv(Fixes,outF)
  write.csv(Events,outE)
}
