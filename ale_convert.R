library("RNetCDF")
setwd("~/Documents/ECLs/Algorithm Comparison/Alejandro/")

fnames=c("ERAI-150-6-global-150-6_v15-1-1.csv","ERAI-150-6-global-150-6_v15L-1-1.csv","ERAI-150-6-global-150-6_v15CTL-1-1.csv",
         "ERAI-150-6-global-150-6_v15CTLNEW-1-1.csv","ERAI-150-6-global-150-6_v15LHR-1-1.csv")

fout=c("Ale_ERAI150_v15.csv","Ale_ERAI150_v15L.csv","Ale_ERAI150_v15CTL.csv","Ale_ERAI150_v15CTLNEW.csv","Ale_ERAI150_v15LHR.csv")

fnames=c("ERAI-75-6-global-150-6_v12-1-1.csv","ERAI-75-6-global-150-6_v18-1-1.csv")
fout=c("Ale_ERAI75_v12.csv","Ale_ERAI75_v118.csv")
for(i in 1:2) 
  
readA(fnames[4],fout[4])

readA<-function(fname,fout)
{
  data<-read.csv(fname)
  date=data$YEAR*10000+data$MONTH*100+data$DAY
  data2=data.frame(ID=data[,1],Length=data$DURATION,Date=date,Time=data$HOUR,Lon=data$LON_CEN,Lat=data$LAT_CEN,MSLP=data$P_CEN,Grad=data$PG200)
  data2$Location=0
  I<-which(data2[,5]>=149 & data2[,5]<=161 & data2[,6]<(-37) & data2[,6]>=-41)
  data2[I,9]<-1
  I<-which(data2[,5]>=(149+(37+data2[,6])/2) & data2[,5]<=161 & data2[,6]<(-31) & data2[,6]>=-37)
  data2[I,9]<-1
  I<-which(data2[,5]>=152 & data2[,5]<=161 & data2[,6]<=(-24) & data2[,6]>=-31)
  data2[I,9]<-1
  
  write.csv(data2,file=fout)
}

data<-read.csv(fname)
sort(unique(data$DURATION))

