####### This is a collection of functions that are used in analysing ECLs

## Code to make a colourbar for figures
ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))

### Compares the PDF of two datasets
makePDF = function(data1,data2,xlabel="Intensity",labloc="topright",leg=c("Data1","Data2"),tit="") {
  a=density(data1,na.rm=T)
  b=density(data2,na.rm=T)
  
  lims=range(data1,data2,na.rm=T)
  if((lims[2]-lims[1])<10)
  {
    lims[1]=floor(lims[1])
    lims[2]=ceiling(lims[2])
  } else {
    lims[1]=floor(lims[1]/5)*5
    lims[2]=ceiling(lims[2]/5)*5
  }
  
  plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
       xlab=xlabel,ylab="Frequency",main=tit)
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=leg,
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
  
  return(ks.test(data1,data2))
}


### Identifies what proportion of events in database 1 are matched by database 2
eventmatch <- function(events1,fixes1,events2,fixes2,GV=F)
{
  if(GV==F)
  {
  match<-array(NaN,c(length(events1$ID),8))
  dimnames(match)[[2]]=c("CV","MSLP","Length","MatchEvents","MatchHours","CV2","MSLP2","Length2")
  match[,1]=events1$CV2
  match[,2]=events1$MSLP2
  match[,3]=events1$Length2
  
  fixes1$Date2=as.POSIXct(paste(as.character(fixes1$Date),substr(fixes1$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  fixes2$Date2=as.POSIXct(paste(as.character(fixes2$Date),substr(fixes2$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  
  for(i in 1:length(events1$ID))
  {
    tmp=fixes1[(fixes1$ID==events1$ID[i] & fixes1$Location==1),]
    rn=range(tmp$Date2)
    
    I=which(fixes2$Date2<=rn[2]+(60*60*6) & fixes2$Date2>=rn[1]-(60*60*6) & fixes2$Location==1)
    if(length(I)>0)
    {
      J=unique(fixes2$ID[I])
      match[i,4]=length(J) #All events that match
      match[i,5]=length(which(fixes2$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
      
      K=which(events2$ID %in% J)
      match[i,6]=max(events2$CV2[K])
      match[i,7]=min(events2$MSLP2[K])
      match[i,8]=min(events2$Length2[K])
    } else match[i,4]=0
  }
  } else {
    {
      match<-array(NaN,c(length(events1$ID),10))
      dimnames(match)[[2]]=c("CV","MSLP","Length","GV","MatchEvents","MatchHours","CV2","MSLP2","Length2","GV2")
      match[,1]=events1$CV2
      match[,2]=events1$MSLP2
      match[,3]=events1$Length2
      match[,4]=events1$GV
      
      fixes1$Date2=as.POSIXct(paste(as.character(fixes1$Date),substr(fixes1$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
      fixes2$Date2=as.POSIXct(paste(as.character(fixes2$Date),substr(fixes2$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
      
      for(i in 1:length(events1$ID))
      {
        tmp=fixes1[(fixes1$ID==events1$ID[i] & fixes1$Location==1),]
        rn=range(tmp$Date2)
        
        I=which(fixes2$Date2<=rn[2]+(60*60*6) & fixes2$Date2>=rn[1]-(60*60*6) & fixes2$Location==1)
        if(length(I)>0)
        {
          J=unique(fixes2$ID[I])
          match[i,5]=length(J) #All events that match
          match[i,6]=length(which(fixes2$Date2[I] %in% tmp$Date2)) ## Hours where both have an ECL at some time
          
          K=which(events2$ID %in% J)
          match[i,7]=max(events2$CV2[K])
          match[i,8]=min(events2$MSLP2[K])
          match[i,9]=min(events2$Length2[K])
          match[i,10]=max(events2$GV[K])
        } else match[i,4]=0
      }
    }
  }
  
  return(match)
}

### Looking at HR for days

CSI_days <- function(fixes1,fixes2,day1=NaN,day2=NaN)
{
  if(is.nan(day1)) day1=min(c(fixes1$Date,fixes2$Date))
  if(is.nan(day2)) day2=max(c(fixes1$Date,fixes2$Date))
  
  dates=seq(as.Date(as.character(day1),format="%Y%m%d"),
            as.Date(as.character(day2),format="%Y%m%d"),
            by="day",format="%Y%m%d")
  dates2=cbind(as.numeric(format.Date(dates,"%Y%m%d")),matrix(0,length(dates),2))
  
  tmp=sort(unique(fixes1$Date[fixes1$Location==1]))
  I=which(dates2[,1] %in% tmp)
  dates2[I,2]=1
  tmp=sort(unique(fixes2$Date[fixes2$Location==1]))
  I=which(dates2[,1] %in% tmp)
  dates2[I,3]=1
  
  CSI=matrix(0,1,3)
  colnames(CSI)=c("HR","FAR","CSI")
  CSI[1]=length(which(dates2[,2]==1 & dates2[,3]==1))/length(which(dates2[,2]==1)) ## HR - Proportion of data1 matched by data2
  CSI[2]=length(which(dates2[,3]==1 & dates2[,2]==0))/length(which(dates2[,3]==1)) ## FAR - Proportion of data2 not in data1
  CSI[3]=length(which(dates2[,2]==1 & dates2[,3]==1))/length(which(dates2[,3]==1 | dates2[,2]==1)) ## CSI
  
  return(CSI)
}

