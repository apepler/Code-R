##A series of R functions to calculate useful stats for month-to-date rain/temperature
## Written by Acacia Pepler, July 2014.
## Updated by RJS, August 2014.
## Updated by ASP, August 2014 - Added new functions (rain, tables of ranks)
## Updated by RJS, September 2015 - Added function for hadcrut
## Updated by ASP, September 2015 - Added table output
##                                  Added a version that calculates scenarios based on individual
##                                  Changed the "global" function to take either NOAA or HADCRUT as input
##                                  Changed axis scale to reflect data
##                                  Changed upper/lower bounds to exclude current year
## run by opening the file source("<PATH>/YTD_v3.R") then calling function of interest
## IMPORTANT - National state is "aus" for temp and "aust" for rain

##Functions included:
## YTD_temp("mean","aus")/YTD_rain("aust")
##  (return a figure of YTD temp/rain & potential ranks for the chosen state & variable)
##  Optional parameters:
##    outfile="figname.ext"
##    Accepts tiff/tif/jpeg/jpg/png/bmp currently
##    If unspecified, the figure is saved in current working directory as e.g. YTDTmean_aus.png
##    e.g. YTD_temp("mean","aus",outfile="dir/fig.tiff") 
##    TAB=T -> also outputs a table using the same filename but .csv

##
## YTD_tempRanks("mean","aus",c(1,2,3,5,10))/YTD_rainRanks("aust",c(1,2,3,5,10))  
##  (returns a table of anomalies required for certain annual ranks)
##  Optional parameters:
##    outfile="file.csv"/"file.txt" (also saves to a file)
##    RankFrom="high"/"low" (default high) - says if rank is relative to highest or lowest
##    e.g. YTD_tempRanks("mean","aus",c(1,2,3,5,10),outfile="dir/filename.csv")

##
##YTD_globaltemp<-function(VAR,STATE)
# only works for VAR=globalt and ns, hs, and global, e.g. YTD_globaltemp("globalt","global")
YTD_temp<-function(VAR,STATE,outfile=NA,TAB=F,lim=T,internal=F)
{
 if(internal==T) file=paste("http://cmap/climatechange/timeseries/acorn/sat/t",VAR,"/allmonths/",STATE,"/latest.txt",sep="") else
   file=paste("http://www.bom.gov.au/web01/ncc/www/cli_chg/timeseries/t",VAR,"/allmonths/",STATE,"/latest.txt",sep="")
  data=read.table(file)
  
  ##Dates in various formats
  yyyymm=floor(data[,1]/1000000)
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=mean(data2[i,1:j])
  
  ##So, want to plot the historical range, the highest on record, 2014
  
  aves=matrix(0,4,12)
  aves[1,]=apply(data3[1:(length(years)-1),],2,min,na.rm=T) #Min
  aves[2,]=apply(data3[1:(length(years)-1),],2,max,na.rm=T) #Max
  I=which(data3[,12]==max(data3[,12],na.rm=T))
  Mname=years[I]
  aves[3,]=data3[I,] #Year
  I=which(data3[,12]==min(data3[,12],na.rm=T))
  Mname2=years[I]
  aves[4,]=data3[I,] #Year
  
  ###Now a more complex version with scenarios
  curM=mm[length(mm)]
  curY=yy[length(yy)]
  if(curM<11) EOY=apply(data2[,(curM+1):12],1,mean,na.rm=T) else EOY=data2[,12]
 ##Remaining x months
  
  ##Scenarios
  scen=matrix(NaN,6,3)
  for(i in 1:6) scen[i,1]=data3[length(years),curM]
  
  scen[2,2]<-scen[2,3]<-scen[1,1] #Same anomaly 
  scen[1,2]=(curM*scen[2,1]+(12-curM)*max(EOY,na.rm=T))/12 #HOR
  scen[1,3]=max(EOY,na.rm=T) 
  I=which(years>=2000 & years<curY)
  scen[3,2]=(curM*scen[2,1]+(12-curM)*mean(EOY[I]))/12 #2004-2013s
  scen[3,3]=mean(EOY[I])
  I=which(years>=1981 & years<=2010)
  scen[4,2]=(curM*scen[2,1]+(12-curM)*mean(EOY[I]))/12 #1981-2010
  scen[4,3]=mean(EOY[I])
  scen[5,2]=(curM*scen[2,1]+(12-curM)*0)/12 #1961-1990 average
  scen[5,3]=0
  scen[6,2]=(curM*scen[2,1]+(12-curM)*min(EOY,na.rm=T))/12 #LOR
  scen[6,3]=min(EOY,na.rm=T)
 
  names=c("Highest-on-record","Persistence","average since 2000","1981-2010 average","1961-1990 average","Lowest-on-record")
  clist=c("orange","green","yellow","purple","blue","darkblue")
    
  ###Still, want to say how fits into list
  names2=names
  for(i in 1:6)
  {
    I=which(data3[,12]>scen[i,2]) ##How many years HIGHER
    names2[i]=paste(names[i], ": rank ",length(I)+1, sep="")
  }
  
  if(is.na(outfile)==T)
    png(file=paste("YTDT",VAR,"_",STATE,".png",sep=""), height=500, width=1000) else
  {
    a=strsplit(outfile,"\\.")
    ext=a[[1]][2]
    if(ext=="png") png(outfile, height=500, width=1000) else
      if(ext=="tiff" | ext=="tif") tiff(outfile, height=500, width=1000) else
        if(ext=="jpg" | ext=="jpeg") jpeg(outfile, height=500, width=1000) else
          if(ext=="bmp") bmp(outfile, height=500, width=1000) else
        {
          print('Only accepts TIFF/JPEG/PNG/BMP; saving as a png file')
          png(file=paste(a[[1]][1],".png",sep=""), height=500, width=1000)
        }

  }
 
 if(lim==T) lim=ceiling(2*max(abs(range(scen,aves,na.rm=T))))/2 
  mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  plot(seq(1,12),data3[length(years),],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-lim,lim),axes=F)
  axis(1, at = seq(1:12), labels=mnames,cex.axis=1.5) 
  axis(2, at = seq(-3,3),cex.axis=1.5)
  title(main=paste("Monthly year-to-date ",VAR," temperature: ",STATE,sep=""),cex.main=2)
  polygon(c(rev(seq(1:12)),seq(1:12)),c(rev(aves[1,]),aves[2,]),col='grey80',border=NA)
  lines(seq(1,12),aves[3,],col="red",type="b",lwd=4)

  lines(seq(1,12),data3[length(years),],col="black",type="b",lwd=4)
  for(i in 1:6) lines(c(curM,12),scen[i,1:2],col=clist[i],lwd=4)  
  legend("bottomright",legend=c(paste("warmest year:",Mname),paste("year-to-date:",curY),names2),col=c("red","black",clist),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()
  print(paste("ran t", VAR, " for ", STATE, sep=""))
 
 
  if(TAB!=F)
  {
    if(is.na(outfile)==T)
      tabname=paste("YTDT",VAR,"_",STATE,".csv",sep="") else
      {
        a=strsplit(outfile,"\\.")
        tabname=paste(a[[1]][1],".csv",sep="")
      }
    
    if(TAB==1 | TAB==T)
    {
    table=array(NaN,c(11,12))
    colnames(table)=mnames
    rownames(table)=c("Lower Bound","Upper Bound",paste("Warmest year:",Mname),paste("Coldest year:",Mname2),paste("Year-to-date:",curY),names2)
    table[1:4,]=aves[1:4,]
    table[5,]=data3[length(years),]
    table[6:11,12]=scen[,2]
    write.csv(table,file=tabname)
    } else if(TAB==2)  
      {
      table=array(NaN,c(6,3))
      table[,1]=scen[,3]
      table[,2]=scen[,2]
      
      for(i in 1:6)
      {
        I=which(data3[,12]>table[i,2]) ##How many years HIGHER
        J=which(data3[,12]==table[i,2]) ##How many years EQUAL        
        if(length(J)==0) table[i,3]=length(I)+1 else table[i,3]=length(I)+1.5
      }
      
      I=which(EOY==max(EOY,na.rm=T))
      Mname=years[I]
      I=which(EOY==min(EOY,na.rm=T))
      Mname2=years[I]
      
      rownames(table)=c(paste("Warmest on record (",Mname,")",sep=""),"Persistence of anomaly", 
                        "Average since 2000","1981-2010 average",
                        "1961-1990 average",paste("Coldest on record (",Mname2,")",sep=""))    
      colnames(table)=c("EOY anomaly","Annual anomaly","Historical rank")
      print(table,digits=2)
      write.csv(table,file=tabname)
    }
  }

}

YTD_temp_monthly<-function(VAR,STATE,outfile=NA,TAB=F,lim=T,internal=F)
{
  if(internal==T) file=paste("http://cmap/climatechange/timeseries/acorn/sat/t",VAR,"/allmonths/",STATE,"/latest.txt",sep="") else
    file=paste("http://www.bom.gov.au/web01/ncc/www/cli_chg/timeseries/t",VAR,"/allmonths/",STATE,"/latest.txt",sep="")
  data=read.table(file)
  
  ##Dates in various formats
  yyyymm=floor(data[,1]/1000000)
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=mean(data2[i,1:j])
  
  ##So, want to plot the historical range, the highest on record, 2014
  
  aves=matrix(0,4,12)
  aves[1,]=apply(data3[1:(length(years)-1),],2,min,na.rm=T) #Min
  aves[2,]=apply(data3[1:(length(years)-1),],2,max,na.rm=T) #Max
  I=which(data3[,12]==max(data3[,12],na.rm=T))
  Mname=years[I]
  aves[3,]=data3[I,] #Year
  I=which(data3[,12]==min(data3[,12],na.rm=T))
  Mname2=years[I]
  aves[4,]=data3[I,] #Year
  
  ###Now a more complex version with scenarios
  curM=mm[length(mm)]
  curY=yy[length(yy)]
  if(curM<11) EOY=apply(data2[,(curM+1):12],1,mean,na.rm=T) else EOY=data2[,12]
  ##Remaining x months
  
  ##Scenarios
  scen=matrix(NaN,6,12)
  for(i in 1:6) scen[i,curM]=data3[length(years),curM]
  for(i in (curM+1):12)
  {
    scen[2,i]=scen[2,i-1] #Same anomaly
    scen[1,i]=((i-1)*scen[1,i-1]+max(data2[,i],na.rm=T))/i #HOR
    I=which(years>=2000 & years<curY)
    scen[3,i]=((i-1)*scen[3,i-1]+mean(data2[I,i]))/i #2004-2013s
    I=which(years>=1981 & years<=2010)
    scen[4,i]=((i-1)*scen[4,i-1]+mean(data2[I,i]))/i #1981-2010
    scen[5,i]=((i-1)*scen[5,i-1])/i #1961-1990 average
    scen[6,i]=((i-1)*scen[6,i-1]+min(data2[,i],na.rm=T))/i #LOR
  }
    
  names=c("Highest-on-record","Persistence","average since 2000","1981-2010 average","1961-1990 average","Lowest-on-record")
  clist=c("orange","green","yellow","purple","blue","darkblue")
  
  ###Still, want to say how fits into list
  names2=names
  for(i in 1:6)
  {
    I=which(data3[,12]>scen[i,12]) ##How many years HIGHER
    names2[i]=paste(names[i], ": rank ",length(I)+1, sep="")
  }
  
  if(is.na(outfile)==T)
    png(file=paste("YTDT",VAR,"_",STATE,"_v2.png",sep=""), height=500, width=1000) else
    {
      a=strsplit(outfile,"\\.")
      ext=a[[1]][2]
      if(ext=="png") png(outfile, height=500, width=1000) else
        if(ext=="tiff" | ext=="tif") tiff(outfile, height=500, width=1000) else
          if(ext=="jpg" | ext=="jpeg") jpeg(outfile, height=500, width=1000) else
            if(ext=="bmp") bmp(outfile, height=500, width=1000) else
            {
              print('Only accepts TIFF/JPEG/PNG/BMP; saving as a png file')
              png(file=paste(a[[1]][1],".png",sep=""), height=500, width=1000)
            }
      
    }
  
  if(lim==T) lim=ceiling(2*max(abs(range(scen,aves,na.rm=T))))/2 
  mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  plot(seq(1,12),data3[length(years),],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-lim,lim),axes=F)
  axis(1, at = seq(1:12), labels=mnames,cex.axis=1.5) 
  axis(2, at = seq(-3,3),cex.axis=1.5)
  title(main=paste("Monthly year-to-date ",VAR," temperature: ",STATE,sep=""),cex.main=2)
  polygon(c(rev(seq(1:12)),seq(1:12)),c(rev(aves[1,]),aves[2,]),col='grey80',border=NA)
  lines(seq(1,12),aves[3,],col="red",type="b",lwd=4)
  
  lines(seq(1,12),data3[length(years),],col="black",type="b",lwd=4)
  for(i in 1:6) lines(seq(1:12),scen[i,],col=clist[i],lwd=4)  
  legend("bottomright",legend=c(paste("warmest year:",Mname),paste("year-to-date:",curY),names2),col=c("red","black",clist),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()
  print(paste("ran t", VAR, " for ", STATE, sep=""))
  
  
  if(TAB==T)
  {
    if(is.na(outfile)==T)
      tabname=paste("YTDT",VAR,"_",STATE,"_v2.csv",sep="") else
      {
        a=strsplit(outfile,"\\.")
        tabname=paste(a[[1]][1],".csv",sep="")
      }
    
    table=array(NaN,c(11,12))
    colnames(table)=mnames
    rownames(table)=c("Lower Bound","Upper Bound",paste("Warmest year:",Mname),paste("Coldest year:",Mname2),paste("Year-to-date:",curY),names2)
    table[1:4,]=aves[1:4,]
    table[5,]=data3[length(years),]
    table[6:11,]=scen
    write.csv(table,file=tabname)
  }
  
}

YTD_globaltemp<-function(SOURCE,STATE="global",outfile=NA,TAB=F,lim=T)
# SOURCE is NOAA or HACRUT
# STATE only works for HADCRUT - nh, sh, and global, e.g. YTD_globaltemp(globalt,global)
{
  if(SOURCE=="NOAA" | SOURCE=="noaa")
  {
    file="http://www.ncdc.noaa.gov/cag/time-series/global/globe/land_ocean/p12/12/1880-2015.csv"
    data=read.csv(file,skip=2)
  } else if(SOURCE=="HAD" | SOURCE=="HADCRUT" | SOURCE=="had" | SOURCE=="hadcrut")
  {
    file=paste("http://cmap/climatechange/timeseries/global_t/allmonths/",STATE,"/latest.txt",sep="")
    data=read.table(file)
  }
  
  ##Dates in various formats
  yyyymm=data[,1]
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=mean(data2[i,1:j])
  
  ##So, want to plot the historical range, the highest on record, 2014
  
  aves=matrix(0,4,12)
  aves[1,]=apply(data3[1:(length(years)-1),],2,min,na.rm=T) #Min
  aves[2,]=apply(data3[1:(length(years)-1),],2,max,na.rm=T) #Max
  I=which(data3[,12]==max(data3[,12],na.rm=T))
  Mname=years[I]
  aves[3,]=data3[I,] #Year
  I=which(data3[,12]==min(data3[,12],na.rm=T))
  Mname2=years[I]
  aves[4,]=data3[I,] #Year
  
  ###Now a more complex version with scenarios
  curM=mm[length(mm)]
  curY=yy[length(yy)]
  if(curM<11) EOY=apply(data2[,(curM+1):12],1,mean,na.rm=T) else EOY=data2[,12]
  ##Remaining x months
  
  ##Scenarios
  scen=matrix(NaN,6,3)
  for(i in 1:6) scen[i,1]=data3[length(years),curM]
  
  scen[2,2]<-scen[2,3]<-scen[1,1] #Same anomaly 
  scen[1,2]=(curM*scen[2,1]+(12-curM)*max(EOY,na.rm=T))/12 #HOR
  scen[1,3]=max(EOY,na.rm=T) 
  I=which(years>=2000 & years<curY)
  scen[3,2]=(curM*scen[2,1]+(12-curM)*mean(EOY[I]))/12 #2004-2013s
  scen[3,3]=mean(EOY[I])
  I=which(years>=1981 & years<=2010)
  scen[4,2]=(curM*scen[2,1]+(12-curM)*mean(EOY[I]))/12 #1981-2010
  scen[4,3]=mean(EOY[I])
  I=which(years>=1961 & years<=1990)
  scen[5,2]=(curM*scen[2,1]+(12-curM)*mean(EOY[I]))/12 #1981-2010
  scen[5,3]=mean(EOY[I])
  scen[6,2]=(curM*scen[2,1]+(12-curM)*min(EOY,na.rm=T))/12 #LOR
  scen[6,3]=min(EOY,na.rm=T)
  
  names=c("Highest-on-record","Persistence","average since 2000","1981-2010 average","1961-1990 average","Lowest-on-record")
  clist=c("orange","green","yellow","purple","blue","darkblue")
    
  ###Still, want to say how fits into list
  names2=names
  for(i in 1:6)
  {
    I=which(data3[,12]>scen[i,2]) ##How many years HIGHER
    names2[i]=paste(names[i], ": rank ",length(I)+1, sep="")
  }
  
  if(is.na(outfile)==T)
    png(file=paste("YTD_",SOURCE,"_",STATE,".png",sep=""), height=500, width=1000) else
  {
    a=strsplit(outfile,"\\.")
    ext=a[[1]][2]
    if(ext=="png") png(outfile, height=500, width=1000) else
      if(ext=="tiff" | ext=="tif") tiff(outfile, height=500, width=1000) else
        if(ext=="jpg" | ext=="jpeg") jpeg(outfile, height=500, width=1000) else
          if(ext=="bmp") bmp(outfile, height=500, width=1000) else
        {
          print('Only accepts TIFF/JPEG/PNG/BMP; saving as a png file')
          png(file=paste(a[[1]][1],".png",sep=""), height=500, width=1000)
        }

  }
 
  if(lim==T) lim=ceiling(2*max(abs(range(scen,aves,na.rm=T))))/2 
  mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  plot(seq(1,12),data3[length(years),],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-lim,lim),axes=F)
  axis(1, at = seq(1:12), labels=mnames,cex.axis=1.5) 
  axis(2, at = seq(-3,3),cex.axis=1.5)
  title(main=paste("Monthly year-to-date ",SOURCE," temperature: ",STATE,sep=""),cex.main=2)
  polygon(c(rev(seq(1:12)),seq(1:12)),c(rev(aves[1,]),aves[2,]),col='grey80',border=NA)
  lines(seq(1,12),aves[3,],col="red",type="b",lwd=4)

  lines(seq(1,12),data3[length(years),],col="black",type="b",lwd=4)
  for(i in 1:6) lines(c(curM,12),scen[i,1:2],col=clist[i],lwd=4)  
  legend("bottomright",legend=c(paste("warmest year: ",Mname,sep=""),"year-to-date: 2015",names2),col=c("red","black",clist),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()
  print(paste("ran ", SOURCE, " for ", STATE, sep=""))
  
  
  if(TAB!=F)
  {
    if(is.na(outfile)==T)
      tabname=paste("YTDT",VAR,"_",STATE,".csv",sep="") else
      {
        a=strsplit(outfile,"\\.")
        tabname=paste(a[[1]][1],".csv",sep="")
      }
    
    if(TAB==1 | TAB==T)
    {
      table=array(NaN,c(11,12))
      colnames(table)=mnames
      rownames(table)=c("Lower Bound","Upper Bound",paste("Warmest year:",Mname),paste("Coldest year:",Mname2),paste("Year-to-date:",curY),names2)
      table[1:4,]=aves[1:4,]
      table[5,]=data3[length(years),]
      table[6:11,12]=scen[,2]
      write.csv(table,file=tabname)
    } else if(TAB==2)  
    {
      table=array(NaN,c(6,3))
      table[,1]=scen[,3]
      table[,2]=scen[,2]
      
      for(i in 1:6)
      {
        I=which(data3[,12]>table[i,2]) ##How many years HIGHER
        J=which(data3[,12]==table[i,2]) ##How many years EQUAL        
        if(length(J)==0) table[i,3]=length(I)+1 else table[i,3]=length(I)+1.5
      }
      
      I=which(EOY==max(EOY,na.rm=T))
      Mname=years[I]
      I=which(EOY==min(EOY,na.rm=T))
      Mname2=years[I]
      
      rownames(table)=c(paste("Warmest on record (",Mname,")",sep=""),"Persistence of anomaly", 
                        "Average since 2000","1981-2010 average",
                        "1961-1990 average",paste("Coldest on record (",Mname2,")",sep=""))    
      colnames(table)=c("EOY anomaly","Annual anomaly","Historical rank")
      print(table,digits=2)
      write.csv(table,file=tabname)
    }
  }

}

YTD_globaltemp_monthly<-function(SOURCE,STATE="global",outfile=NA,TAB=F,lim=T)
  # SOURCE is NOAA or HACRUT
  # STATE only works for HADCRUT - nh, sh, and global, e.g. YTD_globaltemp(globalt,global)
{
  if(SOURCE=="NOAA" | SOURCE=="noaa")
  {
    file="https://www.ncdc.noaa.gov/cag/time-series/global/globe/land_ocean/p12/12/1880-2015.csv"
    data=read.csv(file,skip=2)
  } else if(SOURCE=="HAD" | SOURCE=="HADCRUT" | SOURCE=="had" | SOURCE=="hadcrut")
  {
    file=paste("http://cmap/climatechange/timeseries/global_t/allmonths/",STATE,"/latest.txt",sep="")
    data=read.table(file)
  }
  
  ##Dates in various formats
  yyyymm=data[,1]
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=mean(data2[i,1:j])
  
  ##So, want to plot the historical range, the highest on record, 2014
  
  aves=matrix(0,4,12)
  aves[1,]=apply(data3[1:(length(years)-1),],2,min,na.rm=T) #Min
  aves[2,]=apply(data3[1:(length(years)-1),],2,max,na.rm=T) #Max
  I=which(data3[,12]==max(data3[,12],na.rm=T))
  Mname=years[I]
  aves[3,]=data3[I,] #Year
  I=which(data3[,12]==min(data3[,12],na.rm=T))
  Mname2=years[I]
  aves[4,]=data3[I,] #Year
  
  ###Now a more complex version with scenarios
  curM=mm[length(mm)]
  curY=yy[length(yy)]
  if(curM<11) EOY=apply(data2[,(curM+1):12],1,mean,na.rm=T) else EOY=data2[,12] ##Remaining x months
  
  ##Scenarios
  scen=matrix(NaN,6,12)
  for(i in 1:6) scen[i,curM]=data3[length(years),curM]
  for(i in (curM+1):12)
  {
    scen[2,i]=scen[2,i-1] #Same anomaly
    scen[1,i]=((i-1)*scen[1,i-1]+max(data2[,i],na.rm=T))/i #HOR
    I=which(years>=2000 & years<curY)
    scen[3,i]=((i-1)*scen[3,i-1]+mean(data2[I,i]))/i #2004-2013s
    I=which(years>=1981 & years<=2010)
    scen[4,i]=((i-1)*scen[4,i-1]+mean(data2[I,i]))/i #1981-2010
    I=which(years>=1961 & years<=1990)
    scen[5,i]=((i-1)*scen[5,i-1]+mean(data2[I,i]))/i  #1961-1990 average
    scen[6,i]=((i-1)*scen[6,i-1]+min(data2[,i],na.rm=T))/i #LOR
  }
  

  names=c("Highest-on-record","Persistence","average since 2000","1981-2010 average","1961-1990 average","Lowest-on-record")
  clist=c("orange","green","yellow","purple","blue","darkblue")
  
  ###Still, want to say how fits into list
  names2=names
  for(i in 1:6)
  {
    I=which(data3[,12]>scen[i,12]) ##How many years HIGHER
    names2[i]=paste(names[i], ": rank ",length(I)+1, sep="")
  }
  
  if(is.na(outfile)==T)
    png(file=paste("YTD_",SOURCE,"_",STATE,"_v2.png",sep=""), height=500, width=1000) else
    {
      a=strsplit(outfile,"\\.")
      ext=a[[1]][2]
      if(ext=="png") png(outfile, height=500, width=1000) else
        if(ext=="tiff" | ext=="tif") tiff(outfile, height=500, width=1000) else
          if(ext=="jpg" | ext=="jpeg") jpeg(outfile, height=500, width=1000) else
            if(ext=="bmp") bmp(outfile, height=500, width=1000) else
            {
              print('Only accepts TIFF/JPEG/PNG/BMP; saving as a png file')
              png(file=paste(a[[1]][1],".png",sep=""), height=500, width=1000)
            }
      
    }
  
  if(lim==T) lim=ceiling(2*max(abs(range(scen,aves,na.rm=T))))/2 
  mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  plot(seq(1,12),data3[length(years),],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-lim,lim),axes=F)
  axis(1, at = seq(1:12), labels=mnames,cex.axis=1.5) 
  axis(2, at = seq(-3,3),cex.axis=1.5)
  title(main=paste("Monthly year-to-date ",SOURCE," temperature: ",STATE,sep=""),cex.main=2)
  polygon(c(rev(seq(1:12)),seq(1:12)),c(rev(aves[1,]),aves[2,]),col='grey80',border=NA)
  lines(seq(1,12),aves[3,],col="red",type="b",lwd=4)
  
  lines(seq(1,12),data3[length(years),],col="black",type="b",lwd=4)
  for(i in 1:6) lines(seq(1,12),scen[i,],col=clist[i],lwd=4)  
  legend("bottomright",legend=c(paste("warmest year: ",Mname,sep=""),"year-to-date: 2015",names2),col=c("red","black",clist),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()
  print(paste("ran ", SOURCE, " for ", STATE, sep=""))
  
  if(TAB==T)
  {
    if(is.na(outfile)==T)
      tabname=paste("YTD_",SOURCE,"_",STATE,"_v2.csv",sep="") else
      {
        a=strsplit(outfile,"\\.")
        tabname=paste(a[[1]][1],".csv",sep="")
      }
    
    table=array(NaN,c(11,12))
    colnames(table)=mnames
    rownames(table)=c("Lower Bound","Upper Bound",paste("Warmest year:",Mname),paste("Coldest year:",Mname2),paste("Year-to-date:",curY),names2)
    table[1:4,]=aves[1:4,]
    table[5,]=data3[length(years),]
    table[6:11,]=scen
    write.csv(table,file=tabname)
  }
  
}


### In this case, the user inputs a variable, state, and list of desired ranks
### The results are printed to the screen
### If a file is specified, they are also printed to a .csv file
### Returns an error if the filename is not .csv or .txt (I hope)
### Usage: YTD_tempRanks("mean","aus",c(1,2,3,5,10))
### Optional variables: 
###   outfile="blah.csv" - returns a csv table as well as printing
###   RankFrom="high"/"low" (default "high") - whether you want ranks relative to HOR (true) or LOR (false)

YTD_tempRanks<-function(VAR,STATE,RANKS,outfile=NA,RankFrom="high",internal=F)
{
  if(internal==T) file=paste("http://cmap/climatechange/timeseries/acorn/sat/t",VAR,"/allmonths/",STATE,"/latest.txt",sep="") else
    file=paste("http://www.bom.gov.au/web01/ncc/www/cli_chg/timeseries/t",VAR,"/allmonths/",STATE,"/latest.txt",sep="")
  data=read.table(file)
  
  ##Dates in various formats
  yyyymm=floor(data[,1]/1000000)
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=mean(data2[i,1:j])
   
  curM=mm[length(mm)]
  curY=yy[length(yy)]
  YTD=data3[length(years),curM] ##YTD
  if(curM<11) EOY=apply(data2[,(curM+1):12],1,mean,na.rm=T)  else EOY=data2[,12] ##Remaining x months, all cases
  
  ###Return ranks of what would EQUAL current or EXCEED?
  ## Note that there are some equal-highest
  
  OutTable=data.frame(Annual.rank=RANKS,Annual.anomaly=rep(0,length(RANKS)),EOY.anomaly=rep(0,length(RANKS)),EOY.rank=rep(0,length(RANKS)))
  
  if(RankFrom=="high")
  {
    ann.sorted=sort(data3[,12],decreasing=T) ##Sorted list of all the anomalies, Sorted[1] is HOR
    for(rr in 1:length(RANKS))
    {
      OutTable[rr,2]=ann.sorted[RANKS[rr]] ##Anomaly of current year in position
      OutTable[rr,3]=(12*OutTable[rr,2] - curM*YTD)/(12-curM) ##Just through backtracking math from v1
      I=which(EOY>OutTable[rr,3]) ##How many of existing EOYs are higher
      OutTable[rr,4]=length(I)+1
    }
  } else if(RankFrom=="low") {
    ann.sorted=sort(data3[,12],decreasing=F) ##Sorted list of all the anomalies, Sorted[1] is HOR
    for(rr in 1:length(RANKS))
    {
      OutTable[rr,2]=ann.sorted[RANKS[rr]] ##Anomaly of current year in position
      OutTable[rr,3]=(12*OutTable[rr,2] - curM*YTD)/(12-curM) ##Just through backtracking math from v1
      I=which(EOY<OutTable[rr,3]) ##How many of existing EOYs are lower
      OutTable[rr,4]=length(I)+1
    }
  } else { 
    stop ('RankFrom must be high, low, or left blank (high)')
  }
  print(paste("T",VAR," Anomalies for ",STATE," required to equal years with chosen ranks (",RankFrom,")",sep=""))
  print(OutTable)
  
  ###Writing to file
  if(is.na(outfile)==F)
  {
    outtype=substr(outfile,nchar(outfile)-3,nchar(outfile)) ##Last 4 characters 
    if(outtype==".txt") write.table(OutTable,file=outfile,row.names=F) else
      if(outtype==".csv") write.csv(OutTable,file=outfile,row.names=F) else
        stop ('outfile must be .csv, .txt, or left blank')
  }
}

##This should find the right input file corresponding to the 
##state (aust,nsw,nt,...) & variable (max,min,mean)
##Find the current year & final month & calc the month-to-date stuff
##Then return a pretty figure

YTD_rain<-function(STATE="aust",outfile=NA,TAB=F,internal=F)
{

  if (internal==T) file=paste("http://cmap/opr/temp_highqual.dir/rainfall.dir/rr.allmonths.",STATE,"av",sep="") else
    file=paste("http://www.bom.gov.au/web01/ncc/www/cli_chg/timeseries/rain/allmonths/",STATE,"/latest.txt",sep="")
  data=read.table(file)
  
  ##Dates in various formats
  yyyymm=floor(data[,1]/1000000)
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=sum(data2[i,1:j])
  
  ##So, want to plot the historical range, the highest on record, 2014
  
  aves=matrix(0,4,12)
  aves[1,]=apply(data3[1:(length(years)-1),],2,min,na.rm=T) #Min
  aves[2,]=apply(data3[1:(length(years)-1),],2,max,na.rm=T) #Max
  I=which(data3[,12]==max(data3[,12],na.rm=T))
  Mname=years[I]
  aves[3,]=data3[I,] #Year
  I=which(data3[1:(length(years)-1),12]==min(data3[1:(length(years)-1),12],na.rm=T))
  Mname2=years[I]
  aves[4,]=data3[I,] #Year
  
  ###Now a more complex version with scenarios
  curM=mm[length(mm)]
  curY=yy[length(yy)]
  if(curM<11) EOY=apply(data2[,(curM+1):12],1,sum,na.rm=T)  else EOY=data2[,12]##Remaining x months
  
  ##Scenarios
  scen=matrix(NaN,5,3)
  for(i in 1:5) scen[i,1]=data3[length(years),curM]
  
  scen[1,2]=scen[2,1]+max(EOY,na.rm=T) #HOR
  scen[1,3]=max(EOY,na.rm=T)
  I=which(years>=2000 & years<curY)
  scen[2,2]=scen[2,1]+(mean(EOY[I])) #2000s
  scen[2,3]=mean(EOY[I])
  I=which(years>=1981 & years<=2010)
  scen[3,2]=scen[2,1]+(mean(EOY[I])) #1981-2010 average
  scen[3,3]=mean(EOY[I])
  I=which(years>=1961 & years<=1990)
  scen[4,2]=scen[2,1]+(mean(EOY[I])) #1961-1990 average
  scen[4,3]=mean(EOY[I])
  scen[5,2]=scen[2,1]+min(EOY[1:(length(years)-1)],na.rm=T) #LOR
  scen[5,3]=min(EOY[1:(length(years)-1)],na.rm=T)
  
  names=c("Highest-on-record","average since 2000","1981-2010 average","1961-1990 average","Lowest-on-record")
  clist=c("orange","yellow","purple","blue","darkblue")
  
  ###Still, want to say how fits into list
  names2=names
  for(i in 1:5)
  {
    I=which(data3[,12]>scen[i,2])
    names2[i]=paste(names[i], ": rank ",length(I)+1, sep="")
  }
  
  ymax=(floor(max(data3[,12],na.rm=T)/100)+1)*100
  
  if(is.na(outfile)==T)
    png(file=paste("YTDRain_",STATE,".png",sep=""), height=500, width=1000) else
    {
      a=strsplit(outfile,"\\.")
      ext=a[[1]][2]
      if(ext=="png") png(outfile, height=500, width=1000) else
        if(ext=="tiff" | ext=="tif") tiff(outfile, height=500, width=1000) else
          if(ext=="jpg" | ext=="jpeg") jpeg(outfile, height=500, width=1000) else
            if(ext=="bmp") bmp(outfile, height=500, width=1000) else
            {
              print('Only accepts TIFF/JPEG/PNG/BMP; saving as a png file')
              png(file=paste(a[[1]][1],".png",sep=""), height=500, width=1000)
            }
      
    }
  
  mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  plot(seq(1,12),data3[length(years),],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(0,ymax),axes=F)
  axis(1, at = seq(1:12), labels=mnames,cex.axis=1.5) 
  axis(2, at = seq(0,ymax,100),cex.axis=1.5)
  title(main=paste("YTD rain for ",STATE,sep=""),cex.main=2)
  polygon(c(rev(seq(1:12)),seq(1:12)),c(rev(aves[1,]),aves[2,]),col='grey80',border=NA)
  lines(seq(1,12),aves[3,],col="red",lwd=4,type="b")
  lines(seq(1,12),data3[length(years),],col="black",lwd=4,type="b")
  for(i in 1:5) lines(c(curM,12),scen[i,1:2],col=clist[i],lwd=4)  
  legend("topleft",legend=c(paste("Max: ",Mname,sep=""),"2014",names2),col=c("red","black",clist),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()
  
  print(paste("ran rain for ", STATE, sep=""))
  
  if(TAB!=F)
  {
    if(is.na(outfile)==T)
      tabname=paste("YTDRain_",STATE,".csv",sep="") else
      {
        a=strsplit(outfile,"\\.")
        tabname=paste(a[[1]][1],".csv",sep="")
      }
    
    if(TAB==1 | TAB==T)
    {
      table=array(NaN,c(10,12))
      colnames(table)=mnames
      rownames(table)=c("Lower Bound","Upper Bound",paste("Wettest year:",Mname),paste("Driest year:",Mname2),paste("Year-to-date:",curY),names2)
      table[1:4,]=aves
      table[5,]=data3[length(years),]
      table[6:10,12]=scen[,2]
      write.csv(table,file=tabname)
    } else if(TAB==2)  
    {
      table=array(NaN,c(5,3))
      I=which(years>=1961 & years<=1990)
      table[,1]=100*((scen[,3]/mean(EOY[I]))-1)
      table[,2]=100*((scen[,2]/mean(data3[I,12]))-1)
      
      for(i in 1:5)
      {
        I=which(data3[,12]>scen[i,2]) ##How many years HIGHER
        J=which(data3[,12]==scen[i,2]) ##How many years EQUAL        
        if(length(J)==0) table[i,3]=length(I)+1 else table[i,3]=length(I)+1.5
      }
      
      I=which(EOY==max(EOY,na.rm=T))
      Mname=years[I]
      I=which(EOY==min(EOY[1:(length(years)-1)],na.rm=T))
      Mname2=years[I]
      
      rownames(table)=c(paste("Wettest on record (",Mname,")",sep=""),
                        "Average since 2000","1981-2010 average",
                        "1961-1990 average",paste("Driest on record (",Mname2,")",sep=""))    
      colnames(table)=c("EOY anomaly (%)","Annual anomaly (%)","Historical rank")
      print(table,digits=2)
      write.csv(table,file=tabname)
    }
  }
  
}

YTD_rain_monthly<-function(STATE="aust",outfile=NA,TAB=F,internal=F)
{
  
  if (internal==T) file=paste("http://cmap/opr/temp_highqual.dir/rainfall.dir/rr.allmonths.",STATE,"av",sep="") else
    file=paste("http://www.bom.gov.au/web01/ncc/www/cli_chg/timeseries/rain/allmonths/",STATE,"/latest.txt",sep="")
  data=read.table(file)
  
  ##Dates in various formats
  yyyymm=floor(data[,1]/1000000)
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=sum(data2[i,1:j])
  
  ##So, want to plot the historical range, the highest on record, 2014
  
  aves=matrix(0,4,12)
  aves[1,]=apply(data3[1:(length(years)-1),],2,min,na.rm=T) #Min
  aves[2,]=apply(data3[1:(length(years)-1),],2,max,na.rm=T) #Max
  I=which(data3[,12]==max(data3[,12],na.rm=T))
  Mname=years[I]
  aves[3,]=data3[I,] #Year
  I=which(data3[,12]==min(data3[1:(length(years)-1),12],na.rm=T))
  Mname2=years[I]
  aves[4,]=data3[I,] #Year
    
  ###Now a more complex version with scenarios
  curM=mm[length(mm)]
  curY=yy[length(yy)]
  if(curM<11) EOY=apply(data2[,(curM+1):12],1,sum,na.rm=T) else EOY=data2[,12]##Remaining x months
  
  ##Scenarios
  scen=matrix(NaN,5,12)
  for(i in 1:5) scen[i,curM]=data3[length(years),curM]
  for(i in (curM+1):12)
  {
    scen[1,i]=scen[1,i-1]+max(data2[,i],na.rm=T) #HOR
    I=which(years>=2000 & years<curY)
    scen[2,i]=scen[2,i-1]+mean(data2[I,i]) #2004-2013s
    I=which(years>=1981 & years<=2010)
    scen[3,i]=scen[3,i-1]+mean(data2[I,i]) #1981-2010
    I=which(years>=1961 & years<=1990)
    scen[4,i]=scen[4,i-1]+mean(data2[I,i]) #1961-1990 average
    scen[5,i]=scen[5,i-1]+min(data2[,i],na.rm=T) #LOR
  }
  
  
  names=c("Highest-on-record","average since 2000","1981-2010 average","1961-1990 average","Lowest-on-record")
  clist=c("orange","yellow","purple","blue","darkblue")
  
  ###Still, want to say how fits into list
  names2=names
  for(i in 1:5)
  {
    I=which(data3[,12]>scen[i,12])
    names2[i]=paste(names[i], ": rank ",length(I)+1, sep="")
  }
  
  ymax=(floor(max(data3[,12],na.rm=T)/100)+1)*100
  
  if(is.na(outfile)==T)
    png(file=paste("YTDRain_",STATE,"_v2.png",sep=""), height=500, width=1000) else
    {
      a=strsplit(outfile,"\\.")
      ext=a[[1]][2]
      if(ext=="png") png(outfile, height=500, width=1000) else
        if(ext=="tiff" | ext=="tif") tiff(outfile, height=500, width=1000) else
          if(ext=="jpg" | ext=="jpeg") jpeg(outfile, height=500, width=1000) else
            if(ext=="bmp") bmp(outfile, height=500, width=1000) else
            {
              print('Only accepts TIFF/JPEG/PNG/BMP; saving as a png file')
              png(file=paste(a[[1]][1],".png",sep=""), height=500, width=1000)
            }
      
    }
  
  mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  plot(seq(1,12),data3[length(years),],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(0,ymax),axes=F)
  axis(1, at = seq(1:12), labels=mnames,cex.axis=1.5) 
  axis(2, at = seq(0,ymax,100),cex.axis=1.5)
  title(main=paste("YTD rain for ",STATE,sep=""),cex.main=2)
  polygon(c(rev(seq(1:12)),seq(1:12)),c(rev(aves[1,]),aves[2,]),col='grey80',border=NA)
  lines(seq(1,12),aves[3,],col="red",lwd=4,type="b")
  lines(seq(1,12),data3[length(years),],col="black",lwd=4,type="b")
  for(i in 1:5) lines(seq(1,12),scen[i,],col=clist[i],lwd=4)  
  legend("topleft",legend=c(paste("Max: ",Mname,sep=""),"2014",names2),col=c("red","black",clist),lwd=4,bty="n",cex=1.5,ncol=2)
  dev.off()
  
  print(paste("ran rain for ", STATE, sep=""))
  
  if(TAB==T)
  {
    if(is.na(outfile)==T)
      tabname=paste("YTDRain_",STATE,"_v2.csv",sep="") else
      {
        a=strsplit(outfile,"\\.")
        tabname=paste(a[[1]][1],".csv",sep="")
      }
    
    table=array(NaN,c(10,12))
    colnames(table)=mnames
    rownames(table)=c("Lower Bound","Upper Bound",paste("Wettest year:",Mname),paste("Driest year:",Mname2),paste("Year-to-date:",curY),names2)
    table[1:4,]=aves[1:4,]
    table[5,]=data3[length(years),]
    table[6:10,]=scen
    write.csv(table,file=tabname)
  }
}


### In this case, the user inputs a variable, state, and list of desired ranks
### The results are printed to the screen
### If a file is specified, they are also printed to a .csv file
### Returns an error if the filename is not .csv or .txt (I hope)
### Usage: YTD_rainRanks("aus",c(1,2,3,5,10))
### Optional variables: 
###   outfile="blah.csv" - returns a csv table as well as printing
###   RankFrom="high"/"low" (default "high") - whether you want ranks relative to HOR (true) or LOR (false)

YTD_rainRanks<-function(STATE,RANKS,outfile=NA,RankFrom="high",internal=F)
{
  
  if (internal==T) file=paste("http://cmap/opr/temp_highqual.dir/rainfall.dir/rr.allmonths.",STATE,"av",sep="") else
    file=paste("http://www.bom.gov.au/web01/ncc/www/cli_chg/timeseries/rain/allmonths/",STATE,"/latest.txt",sep="")
  data=read.table(file)
  
  ##Dates in various formats
  yyyymm=floor(data[,1]/1000000)
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=sum(data2[i,1:j])
  
  curM=mm[length(mm)]
  curY=yy[length(yy)]
  YTD=data3[length(years),curM] ##YTD
  if(curM<11) EOY=apply(data2[,(curM+1):12],1,sum,na.rm=T)  else EOY=data2[,12] ##Remaining x months
  
  ###Return ranks of what would EQUAL current or EXCEED?
  ## Note that there are some equal-highest
  
  OutTable=data.frame(Annual.rank=RANKS,Annual.anomaly=rep(0,length(RANKS)),EOY.anomaly=rep(0,length(RANKS)),EOY.rank=rep(0,length(RANKS)))
  
  if(RankFrom=="high")
  {
    ann.sorted=sort(data3[,12],decreasing=T) ##Sorted list of all the anomalies, Sorted[1] is HOR
    for(rr in 1:length(RANKS))
    {
      OutTable[rr,2]=ann.sorted[RANKS[rr]] ##Total of current year in position
      OutTable[rr,3]=OutTable[rr,2] - YTD  ##Difference
      I=which(EOY>OutTable[rr,3]) ##How many of existing EOYs are higher
      OutTable[rr,4]=length(I)+1
    }
  } else if(RankFrom=="low") {
    ann.sorted=sort(data3[,12],decreasing=F) ##Sorted list of all the anomalies, Sorted[1] is HOR
    for(rr in 1:length(RANKS))
    {
      OutTable[rr,2]=ann.sorted[RANKS[rr]] ##Total of current year in position
      OutTable[rr,3]=OutTable[rr,2] - YTD  ##Difference
      I=which(EOY<OutTable[rr,3]) ##How many of existing EOYs are lower
      OutTable[rr,4]=length(I)+1
    }
  } else { 
    stop ('RankFrom must be high, low, or left blank (high)')
  }
  print(paste("Rainfall totals for ",STATE," required to equal years with chosen ranks (",RankFrom,")",sep=""))
  print(OutTable)
  
  ###Writing to file
  if(is.na(outfile)==F)
  {
    outtype=substr(outfile,nchar(outfile)-3,nchar(outfile)) ##Last 4 characters 
    if(outtype==".txt") write.table(OutTable,file=outfile,row.names=F) else
      if(outtype==".csv") write.csv(OutTable,file=outfile,row.names=F) else
        stop ('outfile must be .csv, .txt, or left blank')
  }
}

### This version creates a histogram of output

YTD_temp_hist<-function(VAR,STATE,outfile=NA,int=seq(-2.05,2.05,0.1),type="h",STATE2=NA,internal=F)
{
  if(internal==T) file=paste("http://cmap/climatechange/timeseries/acorn/sat/t",VAR,"/allmonths/",STATE,"/latest.txt",sep="") else
    file=paste("http://www.bom.gov.au/web01/ncc/www/cli_chg/timeseries/t",VAR,"/allmonths/",STATE,"/latest.txt",sep="")
  data=read.table(file)
  if(is.na(STATE2)) STATE2=STATE
  ##Dates in various formats
  yyyymm=floor(data[,1]/1000000)
  yy=floor(yyyymm/100) 
  mm=yyyymm%%100
  
  ##Create a new table of yearsx12
  years=unique(yy)
  data2=matrix(NaN,length(years),12)
  n=1
  for(i in 1:length(years))
    for(j in 1:12)
    {
      data2[i,j]=data[n,2]
      n=n+1
    }
  
  ##Create a third table of YTD
  data3=data2
  for(i in 1:length(years))
    for(j in 1:12)
      data3[i,j]=mean(data2[i,1:j])
  
  ###Now a more complex version with scenarios
  curM=mm[length(mm)]
  curY=yy[length(yy)]
 if(curM==11) EOY=data2[,12] else EOY=apply(data2[,(curM+1):12],1,mean,na.rm=T) ##Remaining x months
  
  FINAL=((12-curM)*EOY + curM*data3[length(years),curM])/12
  ANNUAL=apply(data2,1,mean)


    
  if(is.na(outfile)==T)
    png(file=paste("YTDT",VAR,"_",STATE,"_dens.png",sep=""), height=400, width=600, pointsize=16) else
    {
      oo=strsplit(outfile,"\\.")
      ext=oo[[1]][2]
      if(ext=="png") png(outfile, height=400, width=600, pointsize=16) else
        if(ext=="tiff" | ext=="tif") tiff(outfile, height=400, width=600, pointsize=16) else
          if(ext=="jpg" | ext=="jpeg") jpeg(outfile, height=400, width=600, pointsize=16) else
            if(ext=="bmp") bmp(outfile, height=400, width=600, pointsize=16) else
            {
              print('Only accepts TIFF/JPEG/PNG/BMP; saving as a png file')
              png(file=paste(oo[[1]][1],".png",sep=""), height=400, width=600, pointsize=16)
            }
      
    }
  
  if(type=="d")
  {  
    b=density(ANNUAL,na.rm=T)
    a=density(FINAL,na.rm=T)
    plot(b,col=rgb(0,0,1,1/4),xlim=c(-2,2),ylim=range(0,a$y,b$y),
         xlab=paste("Annual ",VAR," temperature anomaly",sep=""),ylab="Frequency",cex.main=1,
         main=paste("Distribution of annual ",VAR," temperature anomalies in ",STATE2,sep=""))    
    polygon(b,col=rgb(0,0,1,1/4),density=-1)
    polygon(a,col=rgb(1,0,0,1/4),density=-1)
    mod=c(a$x[which(a$y==max(a$y))],a$y[which(a$y==max(a$y))])
    points(mod[1],mod[2],pch=4,lwd=2)
    if(mod[1]<0)
    {
      I=which(ANNUAL<mod[1])
      name=paste(round(mod[1],1),", ",length(I)+1,"th lowest",sep="")
    } else {
      I=which(ANNUAL>mod[1])
      name=paste("+",round(mod[1],1),", ",length(I)+1,"th highest",sep="")
    }
    text(mod[1],mod[2],pos=4,offset=1,labels=name,cex=0.8)
  } else if(type=="h")
  {
    b=hist(apply(data2,1,mean),breaks=int,plot=F)
    a=hist(FINAL,breaks=int,plot=F)
    mx=max(max(a$counts),max(b$counts))
    mx2=(floor(mx/5)+1)*5
    plot(b,col=rgb(0,0,1,1/4),xlim=c(-2,2),ylim=range(0,mx2),
         xlab=paste("Annual ",VAR," temperature anomaly",sep=""),ylab="Frequency",cex.main=1,
         main=paste("Distribution of annual ",VAR," temperature anomalies in ",STATE2,sep=""))
    plot(a,col=rgb(1,0,0,1/4),add=T)
  }
  
  legend("topleft",legend=c("Historical",paste("Outcomes for",curY)),
         col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=4,cex=1,bty="n")         
  dev.off()
  print(paste("ran t", VAR, " for ", STATE, " at month ",curM,sep=""))
  
  
}

