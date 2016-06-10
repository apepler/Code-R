##Load data from text files

setwd('/media/Seagate Expansion Drive/monthly rainfall')
##First need months<-c("06","07","08")

getRainff<-function(months)
{
  fname<-character(length(months))
  
  rain<-ff(0,dim=c(691,886,112))
  
  year<-seq(1900,2011,1)
  
  for(i in 1:length(year))
    for(j in 1:length(months)) 
    {
      fname[j]=paste(year[i],months[j],'.grid',sep="")
      read.table(fname[j], sep="",skip=6,nrows=691)->data
      as.matrix(data)->data
      data[data<0]=NaN
      data<-data[nrow(data):1,]
      rain[,,i]=(rain[,,i]+data)
    }
  
  return(rain)
}
