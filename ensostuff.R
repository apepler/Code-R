setwd("~/Documents/ECLs/Algorithm Comparison/")
read.csv('events_ann.csv')->ann
read.csv('events_warm.csv')->warm
read.csv('events_cool.csv')->cool

cA<-cW<-cC<-matrix(0,5,5)
for(i in 1:5)
  for(j in 1:5)
  {
    cA[i,j]=cor(ann[,i+1],ann[,j+1],use="pairwise.complete.obs")
    cW[i,j]=cor(warm[,i+1],warm[,j+1],use="pairwise.complete.obs")
    cC[i,j]=cor(cool[,i+1],cool[,j+1],use="pairwise.complete.obs")
  }

for(i in 1:5) print(summary(lm(warm[,i+1]~warm[,1])))

read.table('/home/nfs/z3478332/Documents/Data/soi.txt',header=T)->soi
read.table('/home/nfs/z3478332/Documents/Data/n34.txt',header=T)->n34
read.csv('/home/nfs/z3478332/Documents/Data/enso.csv',header=T)->enso
enso=enso[(enso[,1]>=1980 & enso[,1]<=2009),]

I=which(soi[,1]>=1980 & soi[,1]<=2009)
J=which(n34[,1]>=1980 & n34[,1]<=2009)

cool$SOI=apply(soi[I,6:11],1,mean)
cool$N34=apply(n34[J,6:11],1,mean)
warm$SOI=apply(cbind(soi[min(I):(max(I)-1),12:13],soi[(min(I)+1):(max(I)),2:5]),1,mean)
warm$N34=apply(cbind(n34[min(J):(max(J)-1),12:13],n34[(min(J)+1):(max(J)),2:5]),1,mean)

corels=matrix(0,5,7)
for(i in 1:5)
  {
  corels[i,1]=cor(cool[,i+1],cool[,7],use="pairwise.complete.obs")
  corels[i,2]=cor(cool[,i+1],cool[,8],use="pairwise.complete.obs")
  corels[i,3]=cor(warm[,i+1],warm[,7],use="pairwise.complete.obs")
  corels[i,4]=cor(warm[,i+1],warm[,8],use="pairwise.complete.obs")
  corels[i,5]=mean(cool[enso[,2]=="EN",i+1],na.rm=T)
  corels[i,6]=mean(cool[enso[,2]=="",i+1],na.rm=T)
  corels[i,7]=mean(cool[enso[,2]=="LN",i+1],na.rm=T)
}

enso2=data.frame(year=enso[,1],state=factor(enso[,2],levels=c("EN","","LN")))
