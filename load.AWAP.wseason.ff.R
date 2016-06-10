##Load data from text files
setwd('/media/Seagate Expansion Drive/monthly rainfall')

##First need months<-c("06","07","08")

getRainWff<-function(m1a,m1b)
{
  rain<-ff(0,dim=c(691,886,111))
  year<-seq(1900,2011,1)
  
  for(i in 1:(length(year)-1))
  {
    for(j in 1:length(m1a)) 
    {
      fname=paste(year[i],m1a[j],'.grid',sep="")
      read.table(fname, sep="",skip=6,nrows=691)->data
      as.matrix(data)->data
      data[data<0]=NaN
      data<-data[nrow(data):1,]
      rain[,,i]=(rain[,,i]+data)
    }
    for(j in 1:length(m1b)) 
    {
      fname=paste(year[i+1],m1b[j],'.grid',sep="")
      read.table(fname, sep="",skip=6,nrows=691)->data
      as.matrix(data)->data
      data[data<0]=NaN
      data<-data[nrow(data):1,]
      rain[,,i]=(rain[,,i]+data)
    }
  }
  
  return(rain)
}
