if(!exists('swarmData'))source("loadData.R")
library(phyloseq)
library(ape)
library(vipor)
library(Rtsne)
library(vegan)
source('16s/myBiplot.R')

interactAdonis<-plantAdonis<-chimpAdonis<-bonoboAdonis<-list()
for(ii in names(swarmData)){
  plotProp<-swarmData[[ii]][['props']][swarmData[[ii]][['isEnough']]&rownames(swarmData[[ii]][['props']]) %in% rownames(samples),]
  plotProp2<-swarmData[[ii]][['rare']][swarmData[[ii]][['isEnough']]&rownames(swarmData[[ii]][['rare']]) %in% rownames(samples),]
  phyOtuW<-otu_table(plotProp,taxa_are_rows=FALSE)
  phyOtuU<-otu_table(plotProp2,taxa_are_rows=FALSE)
  qiimeDataW<-phyloseq(otu_table=phyOtuW)
  qiimeDataU<-phyloseq(otu_table=phyOtuW,phy_tree=swarmData[[ii]][['tree']])
  brayDist<-distance(qiimeDataW,'bray',binary=TRUE)
  uniDist<-UniFrac(qiimeDataU,weighted=FALSE)
  brayPca<-pcoa(brayDist)
  uniPca<-pcoa(uniDist)
  tsneBray<-Rtsne(brayDist,is_distance=TRUE,verbose=TRUE,perplexity=15,max_iter=3000)
  tsneUni<-Rtsne(uniDist,is_distance=TRUE,verbose=TRUE,perplexity=15,max_iter=3000)
  selectPca<-uniPca
  selectTsne<-tsneUni

  importance<-selectPca$values$Relative_eig
  colnames(selectPca$vectors)<-sprintf('Principal coordinate %d (%d%% of variance)',1:length(importance),round(importance*100))[1:ncol(selectPca$vectors)]
  selectSamples<-samples[rownames(plotProp),]
  selectSamples<-selectSamples[order(selectSamples$bonobo,selectSamples$area2,selectSamples$malaria),]
  colorBrew<-c('#e41a1cBB','#377eb8BB','#4daf4aBB','#984ea3BB','#ff7f00BB','#ffff33BB','#a65628BB','#f781bfBB','#999999BB','#88ddffBB')
  nArea<-length(unique(selectSamples$area2))
  if(nArea>length(colorBrew))stop('Need to adjust colors for more areas')
  areaCols<-colorBrew[1:nArea]
  names(areaCols)<-unique(selectSamples$area2)
  areaPch<-sapply(names(areaCols),function(x)mostAbundant(selectSamples$chimpBonobo[selectSamples$area2==x]))
  #malariaCols3<-c(NA,'#000000CC')
  malariaCols3<-c('#00000022','#000000CC')
  mediumMalariaCol<-'#00000077'
  malariaCols<-c('#00000022','#000000CC')
  malariaCols2<-rainbow.lab(2,alpha=.9,lightMultiple=.7)
  names(malariaCols3)<-names(malariaCols2)<-names(malariaCols)<-c('Laverania negative','Laverania positive')
  speciesPch<-20+1:length(unique(selectSamples$Species))
  #speciesCols<-c('#FF0000CC','#0000FFCC')
  #speciesCols<-rainbow.lab(length(unique(selectSamples$Species)),start=-2,end=2,alpha=.9,lightMultiple=.8)
  speciesCols<-rainbow.lab(length(unique(selectSamples$Species)),start=-2,end=1,alpha=.8,lightMultiple=.8)
  names(speciesCols)<-names(speciesPch)<-sort(unique(selectSamples$chimpBonobo))

  predictors<-model.matrix(~0+Species+malaria+SIV+area,selectSamples)
  #pos<-my.biplot.pcoa(selectPca,predictors,plot.axes=1:2,pch=speciesPch[selectSamples$chimpBonobo],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5,mar=c(4,4,1.5,10),arrowsFilter=Inf)
  #legend(
    #par('usr')[2]+.01*diff(par('usr')[1:2]), 
    #mean(par('usr')[3:4]),
    #c(names(malariaCols),names(areaCols),names(speciesPch)),
    #col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch)))),
    #pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch),
    #pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch))),
    #inset=.01,pt.lwd=3,pt.cex=2.5,
    #xjust=0,xpd=NA
  #)

  print(table(selectSamples$malaria,selectSamples$chimpBonobo))
  pdf(sprintf('out/pcoa_%s.pdf',ii),width=6,height=6)
    pos<-my.biplot.pcoa(selectPca,predictors,plot.axes=1:2,pch=21,bg=speciesCols[selectSamples$chimpBonobo],col=malariaCols3[selectSamples$malaria+1],cex=2.25,lwd=4,arrowsFilter=Inf,las=1,mgp=c(2.75,.75,0),sameAxis=FALSE,bty='l',type='n')
    points(pos[!selectSamples$malaria,],col=malariaCols3[1],cex=2.25,lwd=4,bg=speciesCols[selectSamples[!selectSamples$malaria,'chimpBonobo']],pch=21)
    points(pos[selectSamples$malaria,],col=malariaCols3[2],cex=2.25,lwd=4,bg=speciesCols[selectSamples[selectSamples$malaria,'chimpBonobo']],pch=21)
    #legend('bottomright',as.vector(outer(names(speciesCols),names(malariaCols3),paste,sep=' ')),col=as.vector(malariaCols3[outer(names(speciesCols),names(malariaCols3),function(x,y)y)]),pch=21,pt.bg=as.vector(speciesCols[outer(names(speciesCols),names(malariaCols3),function(x,y)x)]),inset=.01,pt.lwd=4,pt.cex=2.5,bty='n')
    title(main=sprintf('%s',ii,1,2))
  dev.off()

  pdf(sprintf('out/tsne_%s.pdf',ii),width=10,height=8)
    par(mar=c(4,4,1.5,10))
    plot(selectTsne$Y,pch=speciesPch[selectSamples$chimpBonobo],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.5,lwd=3,ylab='t-SNE 2',xlab='t-SNE 1',main=sprintf('%s',ii),bty='l',las=1)
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
  dev.off()

  ss<-samples[labels(brayDist),]
  plantAdonis[[ii]]<-cacheOperation(sprintf('work/adonis_%s.Rdat',ii),adonis,brayDist~bonobo+area2+malaria,data=ss,permutations=1e7,parallel=5)
  interactAdonis[[ii]]<-cacheOperation(sprintf('work/interactAdonis_%s.Rdat',ii),adonis,brayDist~bonobo+area2+malaria*bonobo,data=ss,permutations=1e7,parallel=5)
  chimpDist<-as.matrix(brayDist)
  chimpDist<-as.dist(chimpDist[!ss$bonobo,!ss$bonobo])
  chimpAdonis[[ii]]<-cacheOperation(sprintf('work/adonisChimp_%s.Rdat',ii),adonis,chimpDist~area2+malaria,data=ss[!ss$bonobo,],permutations=1e7,parallel=5)
  bonoboDist<-as.matrix(brayDist)
  bonoboDist<-as.dist(bonoboDist[ss$bonobo,ss$bonobo])
  bonoboAdonis[[ii]]<-cacheOperation(sprintf('work/adonisBonobo_%s.Rdat',ii),adonis,bonoboDist~area2+malaria,data=ss[ss$bonobo,],permutations=1e7,parallel=5)
}
