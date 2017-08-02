## Load libraries


```r
library(dnar)
packageVersion('dnar')
```

```
## [1] '0.1'
```

```r
source('functions.R')
```

## Load data


```r
source("loadData.R")
```

```
## matK
```

```
## Cache work/matK_rarefyOtus.Rdat does exist. Loading data
```

```
## rbcL
```

```
## Cache work/rbcL_rarefyOtus.Rdat does exist. Loading data
```

### Make heatmap of OTU abundances

```r
sortSamples<-samples[order(samples$Species,samples$area2,samples$malaria,samples$Code),]
cols<-c('white',tail(rev(heat.colors(110)),99)) 

minProp<-.02
for(ii in names(swarmData)){
  if(ii=='matK')subLetters<-c('A','B')
  else subLetters<-c('C','D')
  message(ii)
  plotProp<-swarmData[[ii]][['props']][swarmData[[ii]][['isEnough']],]
  plotProp<-plotProp[,apply(plotProp,2,max)>minProp]
  colnames(plotProp)<-swarmData[[ii]][['taxa']][colnames(plotProp),'bestId']
  plotProp<-plotProp[,order(apply(plotProp,2,mean),decreasing=TRUE)]
  plotProp<-plotProp[orderIn(rownames(plotProp),rownames(sortSamples)),]
  maxProp<-apply(plotProp,2,function(x)x/max(x))
  maxTree<-hclust(dist(t(maxProp[,])))
  maxProp<-maxProp[,rev(maxTree$labels[maxTree$order])]
  plotProp<-plotProp[,rev(maxTree$labels[maxTree$order])]
  breaks<-c(-1e-6,seq(min(plotProp[plotProp>0])-1e-10,max(plotProp)+1e-10,length.out=100))
  breaks2<-c(-1e-6,seq(min(maxProp[maxProp>0])-1e-10,max(maxProp)+1e-10,length.out=100))
  par(mfrow=c(2,1),mar=c(12,.1,3.5,15))
  image(1:ncol(plotProp),1:nrow(plotProp),t(plotProp),col=cols,breaks=breaks,xlab='',ylab='',xaxt='n',yaxt='n',main=ii)
  text(grconvertX(.005, "nfc", "user"),grconvertY(.995, "nfc", "user"),subLetters[1],xpd=NA,cex=3,adj=0:1)
  insetScale(round(breaks2,6),cols,c(.97,.58,.98,.83),label='Proportion of OTU within each sample')
  box()
  axis(1,1:ncol(plotProp),colnames(plotProp),cex.axis=.7,las=2,tcl=-.1,mgp=c(0,.3,0))
  metadata<-sortSamples[rownames(plotProp),c('chimpBonobo','area2','plasmoPM','Code')]
  colnames(metadata)<-c('Species','Site','Laverania','Sample')
  addMetaData(metadata,cex=.75)
  abline(h=1:nrow(plotProp)-.5,v=1:ncol(plotProp)+.5,col='#00000011')
  image(1:ncol(maxProp),1:nrow(maxProp),t(maxProp),col=cols,breaks=breaks2,xlab='',ylab='',xaxt='n',yaxt='n',main=ii)
  text(grconvertX(.005, "nfc", "user"),grconvertY(.995, "nfc", "user"),subLetters[2],xpd=NA,cex=3,adj=0:1)
  box()
  insetScale(round(breaks2,6),cols,c(.97,.58,.98,.83),label='Proportion of OTU maximum')
  axis(1,1:ncol(maxProp),colnames(maxProp),cex.axis=.7,las=2,tcl=-.1,mgp=c(0,.3,0))
  abline(h=1:nrow(maxProp)-.5,v=1:ncol(maxProp)+.5,col='#00000011')
  metadata<-sortSamples[rownames(maxProp),c('chimpBonobo','area2','plasmoPM','Code')]
  colnames(metadata)<-c('Species','Site','Laverania','Sample')
  addMetaData(metadata,cex=.75)
}
```

```
## matK
```

```
## rbcL
```

![plot of chunk plantHeatmap](figure/plantHeatmap-1.png)![plot of chunk plantHeatmap](figure/plantHeatmap-2.png)