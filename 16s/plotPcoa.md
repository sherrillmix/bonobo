## Load libraries

```r
library(phyloseq)
packageVersion('phyloseq')
```

```
## [1] '1.20.0'
```

```r
library(Rtsne)
packageVersion('Rtsne')
```

```
## [1] '0.13'
```

```r
library(ape)
packageVersion('ape')
```

```
## [1] '4.1'
```

## Load the data

```r
source('readData.R')
```

```
## ape package version 4.1
```

```
## phyloseq package version 1.20.0
```

```
## Requiring samples to have 15000 reads
```

```r
source('../functions.R')
```



## Calculate distance matrices from OTU table

```r
  selectSamples<-samples[isEnough[samples$Code],]
  subsampledOtus<-apply(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples$Code],2,rarefyCounts,nRequiredReads)
  phyOtu<-otu_table(subsampledOtus,taxa_are_rows=TRUE)
  qiimeData<-phyloseq(otu_table=phyOtu,phy_tree=tree)
  uniDist<-UniFrac(qiimeData,weighted=FALSE)
  save(uniDist,file='work/uniDist.csv')
  phyOtuW<-otu_table(otuProp[rownames(otuTab) %in% tree$tip.label,selectSamples$Code],taxa_are_rows=TRUE)
  qiimeDataW<-phyloseq(otu_table=phyOtuW,phy_tree=tree)
  uniDistW<-UniFrac(qiimeDataW,weighted=TRUE)
  brayDist<-distance(qiimeDataW,'bray')
  brayDistUW<-distance(qiimeData,'bray',binary=TRUE)
  uniPca<-pcoa(uniDist)
  mantels<-list(
    'uniW'=ade4::mantel.rtest(uniDist,uniDistW,nrepet=1e4),
    'brayW'=ade4::mantel.rtest(uniDist,brayDist,nrepet=1e4),
    'brayUW'=ade4::mantel.rtest(uniDist,brayDistUW,nrepet=1e4)
  )
  predictors<-model.matrix(~0+Species+malaria+SIV+area,selectSamples)
  colnames(predictors)<-sub('^Species','',colnames(predictors))
  colnames(predictors)[colnames(predictors)=='malariaTRUE']<-'malariaPos'
```


## Set up plotting parameters

```r
colorBrew<-c('#e41a1cBB','#377eb8BB','#4daf4aBB','#984ea3BB','#ff7f00BB','#ffff33BB','#a65628BB','#f781bfBB','#999999BB','#88ddffBB')
nArea<-length(unique(selectSamples$area2))
if(nArea>length(colorBrew))stop('Need to adjust colors for more areas')
areaCols<-colorBrew[1:nArea]
names(areaCols)<-unique(selectSamples$area2)
areaPch<-sapply(names(areaCols),function(x)mostAbundant(selectSamples$chimpBonobo[selectSamples$area2==x]))
malariaCols<-c('#00000022','#000000CC')
mediumMalariaCol<-'#00000077'
malariaCols3<-c('#00000022','#000000CC')
names(malariaCols3)<-names(malariaCols)<-c('Laverania negative','Laverania positive')
speciesPch<-20+1:length(unique(selectSamples$chimpBonobo))
speciesCols<-rainbow.lab(length(unique(selectSamples$chimpBonobo)),start=-2,end=1,alpha=.8,lightMultiple=.8)
names(speciesCols)<-names(speciesPch)<-sort(unique(selectSamples$chimpBonobo))
importance<-uniPca$values$Relative_eig
colnames(uniPca$vectors)<-sprintf('Principal coordinate %d (%d%% of variance)',1:length(importance),round(importance*100))[1:ncol(targetPca$vectors)]
```

```
## Error in ncol(targetPca$vectors): object 'targetPca' not found
```
## Plot 16s PCoA

```r
  axes<-1:2
  #malaria/species
  pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[selectSamples$chimpBonobo],col=malariaCols3[selectSamples$malaria+1],cex=2.25,lwd=4,arrowsFilter=Inf,las=1,mgp=c(2.75,.75,0),sameAxis=FALSE,bty='l',type='n')
  points(pos[!selectSamples$malaria,],col=malariaCols3[1],cex=2.25,lwd=4,bg=speciesCols[selectSamples[!selectSamples$malaria,'chimpBonobo']],pch=21)
  points(pos[selectSamples$malaria,],col=malariaCols3[2],cex=2.25,lwd=4,bg=speciesCols[selectSamples[selectSamples$malaria,'chimpBonobo']],pch=21)
  title(main=sprintf('16S rRNA',axes[1],axes[2]))
```

![plot of chunk bacteriaPCoA](figure/bacteriaPCoA-1.png)

## Plot 16s t-SNE

```r
tsne<-Rtsne(uniDist,is_distance=TRUE,verbose=TRUE,perplexity=15,max_iter=5000)
```

```
## Read the 94 x 94 data matrix successfully!
## Using no_dims = 2, perplexity = 15.000000, and theta = 0.500000
## Computing input similarities...
## Building tree...
##  - point 0 of 94
## Done in 0.01 seconds (sparsity = 0.616342)!
## Learning embedding...
## Iteration 50: error is 58.289865 (50 iterations in 0.04 seconds)
## Iteration 100: error is 59.650702 (50 iterations in 0.03 seconds)
## Iteration 150: error is 58.595945 (50 iterations in 0.03 seconds)
## Iteration 200: error is 59.967717 (50 iterations in 0.02 seconds)
## Iteration 250: error is 56.462338 (50 iterations in 0.02 seconds)
## Iteration 300: error is 1.886983 (50 iterations in 0.01 seconds)
## Iteration 350: error is 1.326446 (50 iterations in 0.01 seconds)
## Iteration 400: error is 0.967031 (50 iterations in 0.01 seconds)
## Iteration 450: error is 0.666075 (50 iterations in 0.01 seconds)
## Iteration 500: error is 0.443368 (50 iterations in 0.01 seconds)
## Iteration 550: error is 0.392551 (50 iterations in 0.01 seconds)
## Iteration 600: error is 0.382026 (50 iterations in 0.01 seconds)
## Iteration 650: error is 0.378080 (50 iterations in 0.01 seconds)
## Iteration 700: error is 0.332201 (50 iterations in 0.01 seconds)
## Iteration 750: error is 0.321300 (50 iterations in 0.01 seconds)
## Iteration 800: error is 0.315114 (50 iterations in 0.01 seconds)
## Iteration 850: error is 0.281121 (50 iterations in 0.01 seconds)
## Iteration 900: error is 0.280487 (50 iterations in 0.01 seconds)
## Iteration 950: error is 0.278743 (50 iterations in 0.01 seconds)
## Iteration 1000: error is 0.278635 (50 iterations in 0.01 seconds)
## Iteration 1050: error is 0.279910 (50 iterations in 0.01 seconds)
## Iteration 1100: error is 0.279977 (50 iterations in 0.01 seconds)
## Iteration 1150: error is 0.278770 (50 iterations in 0.01 seconds)
## Iteration 1200: error is 0.280171 (50 iterations in 0.01 seconds)
## Iteration 1250: error is 0.280643 (50 iterations in 0.01 seconds)
## Iteration 1300: error is 0.280017 (50 iterations in 0.01 seconds)
## Iteration 1350: error is 0.281352 (50 iterations in 0.01 seconds)
## Iteration 1400: error is 0.280629 (50 iterations in 0.01 seconds)
## Iteration 1450: error is 0.278825 (50 iterations in 0.01 seconds)
## Iteration 1500: error is 0.280534 (50 iterations in 0.01 seconds)
## Iteration 1550: error is 0.281185 (50 iterations in 0.01 seconds)
## Iteration 1600: error is 0.283002 (50 iterations in 0.01 seconds)
## Iteration 1650: error is 0.279119 (50 iterations in 0.01 seconds)
## Iteration 1700: error is 0.282295 (50 iterations in 0.01 seconds)
## Iteration 1750: error is 0.278931 (50 iterations in 0.01 seconds)
## Iteration 1800: error is 0.280701 (50 iterations in 0.01 seconds)
## Iteration 1850: error is 0.279879 (50 iterations in 0.01 seconds)
## Iteration 1900: error is 0.279457 (50 iterations in 0.01 seconds)
## Iteration 1950: error is 0.279548 (50 iterations in 0.01 seconds)
## Iteration 2000: error is 0.282116 (50 iterations in 0.01 seconds)
## Iteration 2050: error is 0.281311 (50 iterations in 0.01 seconds)
## Iteration 2100: error is 0.279824 (50 iterations in 0.01 seconds)
## Iteration 2150: error is 0.278998 (50 iterations in 0.01 seconds)
## Iteration 2200: error is 0.280140 (50 iterations in 0.01 seconds)
## Iteration 2250: error is 0.280266 (50 iterations in 0.01 seconds)
## Iteration 2300: error is 0.282527 (50 iterations in 0.01 seconds)
## Iteration 2350: error is 0.276199 (50 iterations in 0.02 seconds)
## Iteration 2400: error is 0.278470 (50 iterations in 0.02 seconds)
## Iteration 2450: error is 0.278977 (50 iterations in 0.02 seconds)
## Iteration 2500: error is 0.279900 (50 iterations in 0.02 seconds)
## Iteration 2550: error is 0.279521 (50 iterations in 0.02 seconds)
## Iteration 2600: error is 0.280701 (50 iterations in 0.02 seconds)
## Iteration 2650: error is 0.280156 (50 iterations in 0.02 seconds)
## Iteration 2700: error is 0.279635 (50 iterations in 0.02 seconds)
## Iteration 2750: error is 0.281502 (50 iterations in 0.02 seconds)
## Iteration 2800: error is 0.281417 (50 iterations in 0.02 seconds)
## Iteration 2850: error is 0.280907 (50 iterations in 0.02 seconds)
## Iteration 2900: error is 0.280152 (50 iterations in 0.02 seconds)
## Iteration 2950: error is 0.280882 (50 iterations in 0.02 seconds)
## Iteration 3000: error is 0.281527 (50 iterations in 0.02 seconds)
## Iteration 3050: error is 0.282565 (50 iterations in 0.02 seconds)
## Iteration 3100: error is 0.279673 (50 iterations in 0.02 seconds)
## Iteration 3150: error is 0.279873 (50 iterations in 0.02 seconds)
## Iteration 3200: error is 0.279581 (50 iterations in 0.02 seconds)
## Iteration 3250: error is 0.282249 (50 iterations in 0.02 seconds)
## Iteration 3300: error is 0.281236 (50 iterations in 0.02 seconds)
## Iteration 3350: error is 0.279163 (50 iterations in 0.02 seconds)
## Iteration 3400: error is 0.281716 (50 iterations in 0.02 seconds)
## Iteration 3450: error is 0.282428 (50 iterations in 0.02 seconds)
## Iteration 3500: error is 0.280189 (50 iterations in 0.02 seconds)
## Iteration 3550: error is 0.280365 (50 iterations in 0.02 seconds)
## Iteration 3600: error is 0.279144 (50 iterations in 0.02 seconds)
## Iteration 3650: error is 0.280452 (50 iterations in 0.02 seconds)
## Iteration 3700: error is 0.280388 (50 iterations in 0.02 seconds)
## Iteration 3750: error is 0.279483 (50 iterations in 0.02 seconds)
## Iteration 3800: error is 0.281393 (50 iterations in 0.02 seconds)
## Iteration 3850: error is 0.282006 (50 iterations in 0.02 seconds)
## Iteration 3900: error is 0.278511 (50 iterations in 0.02 seconds)
## Iteration 3950: error is 0.282289 (50 iterations in 0.02 seconds)
## Iteration 4000: error is 0.281744 (50 iterations in 0.02 seconds)
## Iteration 4050: error is 0.280475 (50 iterations in 0.02 seconds)
## Iteration 4100: error is 0.280037 (50 iterations in 0.02 seconds)
## Iteration 4150: error is 0.280059 (50 iterations in 0.02 seconds)
## Iteration 4200: error is 0.280619 (50 iterations in 0.02 seconds)
## Iteration 4250: error is 0.280058 (50 iterations in 0.02 seconds)
## Iteration 4300: error is 0.281159 (50 iterations in 0.02 seconds)
## Iteration 4350: error is 0.280558 (50 iterations in 0.02 seconds)
## Iteration 4400: error is 0.281114 (50 iterations in 0.02 seconds)
## Iteration 4450: error is 0.279915 (50 iterations in 0.02 seconds)
## Iteration 4500: error is 0.280025 (50 iterations in 0.02 seconds)
## Iteration 4550: error is 0.279443 (50 iterations in 0.02 seconds)
## Iteration 4600: error is 0.279553 (50 iterations in 0.02 seconds)
## Iteration 4650: error is 0.279743 (50 iterations in 0.02 seconds)
## Iteration 4700: error is 0.281385 (50 iterations in 0.02 seconds)
## Iteration 4750: error is 0.279032 (50 iterations in 0.02 seconds)
## Iteration 4800: error is 0.278856 (50 iterations in 0.02 seconds)
## Iteration 4850: error is 0.281935 (50 iterations in 0.02 seconds)
## Iteration 4900: error is 0.280185 (50 iterations in 0.02 seconds)
## Iteration 4950: error is 0.280052 (50 iterations in 0.02 seconds)
## Iteration 5000: error is 0.278943 (50 iterations in 0.02 seconds)
## Fitting performed in 1.55 seconds.
```

```r
par(mar=c(4,4,1.5,9.5))
plot(tsne$Y,pch=speciesPch[selectSamples$chimpBonobo],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5,ylab='t-SNE 2',xlab='t-SNE 1',main='',bty='l')
legend(
  par('usr')[2]+.01*diff(par('usr')[1:2]), 
  mean(par('usr')[3:4]),
  c(names(malariaCols),names(areaCols),names(speciesPch)),
  col=c(malariaCols,rep(c(malariaCols[1],mediumMalariaCol),c(length(areaCols),length(speciesPch)))),
  pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch),
  pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch))),
  inset=.01,pt.lwd=3,pt.cex=2.5,
  xjust=0,xpd=NA,bty='n'
)
```

![plot of chunk bacteria_tSNE](figure/bacteria_tSNE-1.png)
