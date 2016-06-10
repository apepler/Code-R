##This is code to take the output from the low tracking software 
##and turn it into a .csv file of ECL track data.

##This version is about extracting any fix that it is within the ECL region & 
##has curvature >= 0.25, even if it's a single-fix event
##Outputs an alternative list of ECL days

#setwd('~/Charles_bom/output')

ECLdatabase<-function(yearS,yearE,output)
{

  ##Files with ECL track data for required years 
  year<-as.character(seq(yearS,yearE))
  fname1=paste('tracks_',year,'.dat',sep="")
    
  ##collates all the annual tracks into one file
  
  read.table(fname1[1], sep="",skip=1)->fixes
  I<-which(fixes[,7]<0)
  fixes<-fixes[I,]
  
  for(i in 2:length(year)) ##Only takes broader E Aus box to make faster
  {
    read.table(fname1[i], sep="",skip=1)->data
    I<-which(data[,7]<(-20) & data[,7]>(-50) & data[,6]>140 & data[,7]<170)
    fixes<-rbind(fixes,data[I,])
  }
  rm(data)
  
  ##Identify if within the ECL region
  
  fixes[,11]<-0
  I<-which(fixes[,6]>=149 & fixes[,6]<=161 & fixes[,7]<(-37) & fixes[,7]>=-41)
  fixes[I,11]<-1
  I<-which(fixes[,6]>=(149+(37+fixes[,7])/2) & fixes[,6]<=161 & fixes[,7]<(-31) & fixes[,7]>=-37)
  fixes[I,11]<-1
  I<-which(fixes[,6]>=152 & fixes[,6]<=161 & fixes[,7]<=(-24) & fixes[,7]>=-31)
  fixes[I,11]<-1
  
  I<-which(fixes[,11]==1 & fixes[,9]>=0.25)
  fixes2=fixes[I,]
  
  Fixes<-data.frame(fixes2[,3:9],fixes2[,11])
  names(Fixes)=c('Date','Time','Open','Lon','Lat','MSLP','CV','Location') 
  outF=paste('ECLfixes_extra_',output,'.csv',sep="")
  write.csv(Fixes,outF)
}
