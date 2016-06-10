rm(list=ls())
setwd("~/Documents")
data=read.csv("snowdata.csv")

data2=data.frame(Year=data$YR,Depth=data$MAX,ENSO=rep("aaa",length(data$YR)),IOD=rep("aaa",length(data$YR)),stringsAsFactors=F)

I=which(data$ENSO=="EN")
data2$ENSO[I]="El Nino"
I=which(data$ENSO=="LN")
data2$ENSO[I]="La Nina"
I=which(data$ENSO=="")
data2$ENSO[I]="Neutral"

I=which(data$IOD..1958..=="P")
data2$IOD[I]="Positive"
I=which(data$IOD..1958..=="N")
data2$IOD[I]="Negative"
I=which(data$IOD..1958..=="")
data2$IOD[I]="Neutral"

data2$ENSO2=factor(data2$ENSO,levels=c("El Nino","Neutral","La Nina"))
data2$IOD2=factor(data2$IOD,levels=c("Positive","Neutral","Negative"))

library(reshape)
library(ggplot2)

pdf(file="snow_enso_boxplot.pdf",height=6,width=4,pointsize=12)
ggplot(data2, aes(x=ENSO2, y=Depth)) + geom_boxplot(fill=c(rgb(1,0,0,1/4),rgb(1,1,0,1/4),rgb(0,0,1,1/4))) +
  theme_bw() + theme(text=element_text(size=20)) + ylab("Maximum snow depth (cm)") + xlab("")
dev.off()
pdf(file="snow_iod_boxplot.pdf",height=6,width=4,pointsize=12)
ggplot(data2[5:61,], aes(x=IOD2, y=Depth)) + geom_boxplot(fill=c(rgb(1,0,0,1/4),rgb(1,1,0,1/4),rgb(0,0,1,1/4))) +
  theme_bw() + theme(text=element_text(size=20)) + ylab("Maximum snow depth (cm)") + xlab("")
dev.off()

### V2
rm(list=ls())
setwd("~/Documents")
data=read.csv("snowdata.csv")

data2=data.frame(Year=data$YR,Depth=data$MAX,ENSO=rep("aaa",length(data$YR)),IOD=rep("aaa",length(data$YR)),SAM=rep("aaa",length(data$YR)),stringsAsFactors=F)

cols=rep("aaa",length(data2$ENSO))
I=which(data$ENSO=="EN")
data2$ENSO[I]=1
cols[I]="red"
I=which(data$ENSO=="LN")
data2$ENSO[I]=3
cols[I]="blue"
I=which(data$ENSO=="")
data2$ENSO[I]=2
cols[I]="orange"
plot(data2$ENSO,data2$Depth,type="p",pch=19,axes=F,xlim=c(0.5,3.5),ylim=c(0,400),col=cols,
     xlab="",ylab="Maximum snow depth (cm)",cex=1,lwd=2)
axis(1,at=c(1,2,3),labels=c("El Nino","Neutral","La Nina"))
axis(2,at=seq(50,350,100))
plot(data2$Year,data2$Depth,type="l",col="black",xlab="Year",ylab="Maximum snow depth (cm)")
points(data2$Year,data2$Depth,col=cols,pch=19)
legend("topright",c("El Nino","Neutral","La Nina"),pch=19,col=c("red","orange","blue"),bty="n",ncol=3)

I=which(data$IOD..1958..=="P")
data2$IOD[I]=1
cols[I]="red"
I=which(data$IOD..1958..=="N")
data2$IOD[I]=3
cols[I]="blue"
I=which(data$IOD..1958..=="")
data2$IOD[I]=2
cols[I]="orange"
plot(data2$IOD[5:61],data2$Depth[5:61],type="p",pch=4,axes=F,xlim=c(0.5,3.5),ylim=c(0,400),col=cols[5:61],
     xlab="",ylab="Maximum snow depth (cm)",cex=1.5,lwd=2)
axis(1,at=c(1,2,3),labels=c("Positive","Neutral","Negative"))
axis(2,at=seq(50,350,100))


plot(data2$Year[5:61],data2$Depth[5:61],type="l",lwd=2,col="black",xlab="Year",ylab="Maximum snow depth (cm)")
points(data2$Year[5:61],data2$Depth[5:61],col=cols[5:61],pch=19)
legend("topright",c("Positive","Neutral","Negative"),pch=19,col=c("red","orange","blue"),bty="n",ncol=3)

cols=rep("NA",length(data2$ENSO))
I=which(data$JJA.SAM>=0.7)
data2$SAM[I]=1
cols[I]="red"
I=which(data$JJA.SAM<=-0.7)
data2$SAM[I]=3
cols[I]="blue"
I=which(abs(data$JJA.SAM)<0.7)
data2$SAM[I]=2
cols[I]="orange"
plot(data2$SAM[4:61],data2$Depth[4:61],type="p",pch=4,axes=F,xlim=c(0.5,3.5),ylim=c(0,400),col=cols[4:61],
     xlab="",ylab="Maximum snow depth (cm)",cex=1.5,lwd=2)
axis(1,at=c(1,2,3),labels=c("Positive","Neutral","Negative"))
axis(2,at=seq(50,350,100))

plot(data$JJA.SAM[4:61],data$MAX[4:61],type="p",pch=4,cex=1.5,lwd=2,ylim=c(0,400),col=cols[4:61],xlab="JJA mean SAM",ylab="Maximum snow depth (cm)")



