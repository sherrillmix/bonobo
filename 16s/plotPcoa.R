if(!exists('otuTab'))source('runQiime.R')

expandedHull<-function(xys,magnification=1,type=c('convex','ellipse')){
  library(cluster)
  type<-match.arg(type)
  centroid<-apply(xys,2,mean)
  magnified<-t(apply(xys,1,function(xy)(xy-centroid)*magnification+centroid))
  if(type=='convex')return(magnified[chull(magnified),])
  else if(type=='ellipse')return(predict(ellipsoidhull(magnified)))
  else stop("Unknown type")
}

library(phyloseq)
library(ape)
library(vipor)

if(!exists('uniPca')){
  selectSamples<-samples[samples$isEnough,]
  tree<-multi2di(read_tree('work/qiime/rep_set.tre'))
  subsampledOtus<-cacheOperation('work/rarefyOtus.Rdat',apply,otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples$name],2,rarefyCounts,nRequiredReads)
  #phyOtu<-otu_table(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples$name],taxa_are_rows=TRUE)
  phyOtu<-otu_table(subsampledOtus,taxa_are_rows=TRUE)
  qiimeData<-phyloseq(otu_table=phyOtu,phy_tree=tree)
  #make sure tree is bifurcating or breaks UniFrac without error
  uniDist<-UniFrac(qiimeData,weighted=FALSE)
  phyOtuW<-otu_table(otuProp[rownames(otuTab) %in% tree$tip.label,selectSamples$name],taxa_are_rows=TRUE)
  qiimeDataW<-phyloseq(otu_table=phyOtuW,phy_tree=tree)
  uniDistW<-UniFrac(qiimeDataW,weighted=TRUE)
  brayDist<-distance(qiimeDataW,'bray')
  brayDistUW<-distance(qiimeDataW,'bray',binary=TRUE)
  jacDist<-distance(qiimeData,'jaccard')
  uniPca<-pcoa(uniDist)
  uniPcaW<-pcoa(uniDistW)
  brayPca<-pcoa(brayDist)
  brayPcaUW<-pcoa(brayDistUW)
  jacPca<-pcoa(jacDist)
  predictors<-model.matrix(~0+Species+malaria+SIV+area,selectSamples)
  colnames(predictors)<-sub('^Species','',colnames(predictors))
  colnames(predictors)[colnames(predictors)=='malariaTRUE']<-'malariaPos'
}


#areaCols<-rainbow.lab(length(unique(samples$area)),alpha=.8)
source('myBiplot.R')
#last color made up
colorBrew<-c('#e41a1cBB','#377eb8BB','#4daf4aBB','#984ea3BB','#ff7f00BB','#ffff33BB','#a65628BB','#f781bfBB','#999999BB','#88ddffBB')
nArea<-length(unique(selectSamples$area2))
if(nArea>length(colorBrew))stop('Need to adjust colors for more areas')
areaCols<-colorBrew[1:nArea]
names(areaCols)<-unique(selectSamples$area2)
areaPch<-sapply(names(areaCols),function(x)mostAbundant(selectSamples$Species[selectSamples$area2==x]))
malariaCols<-c('#00000022','#000000CC')
malariaCols2<-rainbow.lab(2,alpha=.9,lightMultiple=.7)
#malariaCols3<-c(NA,'#000000CC')
malariaCols3<-c('#00000022','#000000CC')
names(malariaCols3)<-names(malariaCols2)<-names(malariaCols)<-c('Laverania negative','Laverania positive')
speciesPch<-20+1:length(unique(selectSamples$Species))
#speciesCols<-c('#FF0000CC','#0000FFCC')
speciesCols<-rainbow.lab(length(unique(selectSamples$Species)),start=-2,end=1,alpha=.8,lightMultiple=.8)
names(speciesCols)<-names(speciesPch)<-sort(unique(selectSamples$chimpBonobo))
#split out TL2-E and -W, bonobos and KR and chimps
splitAreas<-ifelse(selectSamples$area2 %in% c('TL-E','TL-W','KR','TL-NE'),selectSamples$area2,selectSamples$Species)
splitAreaCols<-rainbow.lab(length(unique(splitAreas)))
names(splitAreaCols)<-unique(splitAreas)[order(unique(splitAreas) %in% c('P.t. schweinfurthii'))]
targetPca<-uniPca
importance<-targetPca$values$Relative_eig
colnames(targetPca$vectors)<-sprintf('Principal coordinate %d (%d%% of variance)',1:length(importance),round(importance*100))[1:ncol(targetPca$vectors)]
pdf('out/pcoa.pdf',height=6,width=6)
  sapply(list(1:2,3:4,5:6),function(axes){
    pos<-my.biplot.pcoa(targetPca,predictors,plot.axes=axes,pch=speciesPch[selectSamples$chimpBonobo],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5)
    points(pos[selectSamples$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
    legend(
      'bottomright',
      c(names(malariaCols),names(areaCols),names(speciesPch),'SIVPos'),
      col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch))),'#FF000099'),
      pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch,1),
      pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch)),NA),
      inset=.01,pt.lwd=3,pt.cex=2.5)
    title(main=sprintf('All variables PC %d and %d',axes[1],axes[2]))
    text(pos,selectSamples$Code,cex=.25)
    #species
    pos<-my.biplot.pcoa(targetPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[selectSamples$chimpBonobo],col="#00000077",cex=1.8,lwd=2.5)
    legend('bottomright',names(speciesCols),col='#00000077',pch=21,pt.bg=speciesCols,inset=.01,pt.lwd=3,pt.cex=2,bty='n')
    title(main=sprintf('Species PC %d and %d',axes[1],axes[2]))
    #malaria
    pos<-my.biplot.pcoa(targetPca,predictors,plot.axes=axes,pch=21,bg=malariaCols2[selectSamples$malaria+1],col="#00000077",cex=1.8,lwd=2.5)
    legend('bottomright',names(malariaCols),col='#00000077',pch=21,pt.bg=malariaCols2,inset=.01,pt.lwd=3,pt.cex=2)
    title(main=sprintf('Malaria PC %d and %d',axes[1],axes[2]))
    #malaria/species
    pos<-my.biplot.pcoa(targetPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[selectSamples$chimpBonobo],col=malariaCols3[selectSamples$malaria+1],cex=2.25,lwd=4,arrowsFilter=Inf,las=1,mgp=c(2.75,.75,0),sameAxis=FALSE,bty='l',type='n')
    #legend('bottomleft',as.vector(outer(names(speciesCols),names(malariaCols3),paste,sep=' ')),col=as.vector(malariaCols3[outer(names(speciesCols),names(malariaCols3),function(x,y)y)]),pch=21,pt.bg=as.vector(speciesCols[outer(names(speciesCols),names(malariaCols),function(x,y)x)]),inset=.01,pt.lwd=4,pt.cex=2.5,bty='n')
    points(pos[!selectSamples$malaria,],col=malariaCols3[1],cex=2.25,lwd=4,bg=speciesCols[selectSamples[!selectSamples$malaria,'chimpBonobo']],pch=21)
    points(pos[selectSamples$malaria,],col=malariaCols3[2],cex=2.25,lwd=4,bg=speciesCols[selectSamples[selectSamples$malaria,'chimpBonobo']],pch=21)
    print(table(selectSamples$malaria,selectSamples$chimpBonobo))
    title(main=sprintf('16S rRNA',axes[1],axes[2]))
    for(ii in unique(selectSamples$chimpBonobo)){ 
      hull<-expandedHull(pos[selectSamples$chimpBonobo==ii,],1.02,'ellipse')
      #polygon(hull,border=speciesCols[ii],col=NA,lwd=2.4)
    }
    #labelled
    pos<-my.biplot.pcoa(targetPca,predictors,plot.axes=axes,pch=21,bg=splitAreaCols[splitAreas],col=malariaCols[1],cex=2.25,lwd=4)
    legend('bottomright',names(splitAreaCols),col=malariaCols[1],pch=21,pt.bg=splitAreaCols,inset=.01,pt.lwd=4,pt.cex=2.5)
    title(main=sprintf('TL and KR PC %d and %d',axes[1],axes[2]))
    #show sample name text
    #bak<-par('cex')
    #par('cex'=.65)
    #biplot.pcoa(uniPca,predictors,plot.axes=axes)
    #title(main='Sample names')
    #par('cex'=bak)
  })
dev.off()
system('pdftk out/pcoa.pdf cat 4 output out/Fig.5C.pdf')

pdf('out/pcoa_bray.pdf',height=9,width=9)
  sapply(list(1:2,3:4,5:6),function(axes){
    pos<-my.biplot.pcoa(brayPca,predictors,plot.axes=axes,pch=speciesPch[selectSamples$Species],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5)
    points(pos[selectSamples$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
    legend(
      'bottomright',
      c(names(malariaCols),names(areaCols),names(speciesPch),'SIVPos'),
      col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch))),'#FF000099'),
      pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch,1),
      pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch)),NA),
      inset=.01,pt.lwd=3,pt.cex=2.5)
    title(main=sprintf('All variables PC %d and %d',axes[1],axes[2]))
    text(pos,selectSamples$Code,cex=.25)
    #species
    pos<-my.biplot.pcoa(brayPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[selectSamples$Species],col="#00000077",cex=1.8,lwd=2.5)
    legend('bottomright',names(speciesCols),col='#00000077',pch=21,pt.bg=speciesCols,inset=.01,pt.lwd=3,pt.cex=2)
    title(main=sprintf('Species PC %d and %d',axes[1],axes[2]))
    #malaria
    pos<-my.biplot.pcoa(brayPca,predictors,plot.axes=axes,pch=21,bg=malariaCols2[selectSamples$malaria+1],col="#00000077",cex=1.8,lwd=2.5)
    legend('bottomright',names(malariaCols),col='#00000077',pch=21,pt.bg=malariaCols2,inset=.01,pt.lwd=3,pt.cex=2)
    title(main=sprintf('Malaria PC %d and %d',axes[1],axes[2]))
    #malaria/species
    pos<-my.biplot.pcoa(brayPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[selectSamples$Species],col=malariaCols[selectSamples$malaria+1],cex=2.25,lwd=4,arrowsFilter=Inf)
    legend('bottomright',as.vector(outer(names(speciesCols),names(malariaCols),paste,sep=' ')),col=as.vector(malariaCols[outer(names(speciesCols),names(malariaCols),function(x,y)y)]),pch=21,pt.bg=as.vector(speciesCols[outer(names(speciesCols),names(malariaCols),function(x,y)x)]),inset=.01,pt.lwd=4,pt.cex=2.5)
    title(main=sprintf('Species/malaria PC %d and %d',axes[1],axes[2]))
    #for(ii in unique(selectSamples$Species)){ 
      #hull<-expandedHull(pos[selectSamples$Species==ii,],1.02,'ellipse')
      #polygon(hull,border=speciesCols[ii],col=NA,lwd=2.4)
    #}
    #labelled
    pos<-my.biplot.pcoa(brayPca,predictors,plot.axes=axes,pch=21,bg=splitAreaCols[splitAreas],col=malariaCols[1],cex=2.25,lwd=4)
    legend('bottomright',names(splitAreaCols),col=malariaCols[1],pch=21,pt.bg=splitAreaCols,inset=.01,pt.lwd=4,pt.cex=2.5)
    title(main=sprintf('TL and KR PC %d and %d',axes[1],axes[2]))
  })
dev.off()
system('pdftk out/pcoa_bray.pdf cat 4 output out/pcoa_bray_select.pdf')


#just bonobos
selectSamples2<-samples[samples$bonobo&samples$isEnough,]
phyOtu2<-otu_table(apply(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples2$name],2,rarefyCounts,nRequiredReads),taxa_are_rows=TRUE)
qiimeData2<-phyloseq(otu_table=phyOtu2,phy_tree=tree)
#make sure tree is bifurcating or breaks UniFrac without error
phyOtuW2<-otu_table(otuProp[rownames(otuTab) %in% tree$tip.label,selectSamples2$name],taxa_are_rows=TRUE)
qiimeDataW2<-phyloseq(otu_table=phyOtuW2,phy_tree=tree)
uniDist2<-UniFrac(qiimeData2,weighted=FALSE)
uniDistW2<-UniFrac(qiimeDataW2,weighted=TRUE)
brayDist2<-distance(qiimeDataW2,'bray')
uniPca2<-pcoa(uniDist2)
brayPca2<-pcoa(brayDist2)

selectSamples3<-samples[!samples$bonobo&samples$isEnough,]
phyOtu3<-otu_table(apply(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples3$name],2,rarefyCounts,nRequiredReads),taxa_are_rows=TRUE)
qiimeData3<-phyloseq(otu_table=phyOtu3,phy_tree=tree)
phyOtuW3<-otu_table(otuProp[rownames(otuTab) %in% tree$tip.label,selectSamples3$name],taxa_are_rows=TRUE)
qiimeDataW3<-phyloseq(otu_table=phyOtuW3,phy_tree=tree)
uniDist3<-UniFrac(qiimeData3,weighted=FALSE)
uniDistW3<-UniFrac(qiimeDataW3,weighted=TRUE)
brayDist3<-distance(qiimeDataW3,'bray')
predictors2<-model.matrix(~malaria+area+season,selectSamples2[selectSamples2$bonobo,])
colnames(predictors2)[colnames(predictors2)=='malariaTRUE']<-'malariaPos'
brayPca3<-pcoa(brayDist3)

splitAreas<-withAs(selectSamples2=selectSamples2[selectSamples2$bonobo,],ifelse(selectSamples2$area %in% c('TL-E','TL-W'),selectSamples2$area,selectSamples2$Species))
splitAreaCols<-rainbow.lab(length(unique(splitAreas)))
names(splitAreaCols)<-unique(splitAreas)
pdf('out/bonobo_pcoa.pdf',height=9,width=9)
  sapply(list(1:2,3:4,5:6),function(axes){
    pos<-my.biplot.pcoa(uniPca2,predictors2,plot.axes=axes,pch=21,bg=splitAreaCols[splitAreas],col=malariaCols[selectSamples2[selectSamples2$bonobo,'malaria']+1],cex=2.25,lwd=4)
    legend('bottomright',names(splitAreaCols),col=malariaCols[1],pch=21,pt.bg=splitAreaCols,inset=.01,pt.lwd=4,pt.cex=2.5)
    title(main=sprintf('TL and KR PC %d and %d',axes[1],axes[2]))
  })
dev.off()


sampleNs<-apply(otuTab[,samples$name],2,sum)
#apply(otuTab[,samples$name],1,function(x){
  #glm(cbind(x,sampleNs)
#})

library(Rtsne)
#uw unifrac
tsne<-cacheOperation('work/tsne.Rdat',Rtsne,uniDist,is_distance=TRUE,verbose=TRUE,perplexity=15,max_iter=5000)
#save(tsne,file='work/tnse.Rdat')
tsne2<-Rtsne(uniDist2,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
tsne3<-Rtsne(uniDist3,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
#tsne<-Rtsne(t(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples$name]),verbose=TRUE,perplexity=10)
#w unifrac
tsneW<-Rtsne(uniDistW,is_distance=TRUE,verbose=TRUE,perplexity=10,max_iter=3000)
tsneW2<-Rtsne(uniDistW2,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
tsneW3<-Rtsne(uniDistW3,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
#bray
tsneBray<-Rtsne(brayDist,is_distance=TRUE,verbose=TRUE,perplexity=10,max_iter=3000)
tsneBray2<-Rtsne(brayDist2,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
tsneBray3<-Rtsne(brayDist3,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
pdf('out/tsne.pdf',height=8,width=10)
  inputs<-list('Unweighted'=list(tsne,tsne2,tsne3),'Weighted'=list(tsneW,tsneW2,tsneW3),'Bray'=list(tsneBray,tsneBray2,tsneBray3))
  for(tsneType in names(inputs)){
    tsnes<-inputs[[tsneType]]
    par(mar=c(4,4,1.5,9))
    plot(tsnes[[1]]$Y,pch=speciesPch[selectSamples$Species],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5,ylab='t-SNE 2',xlab='t-SNE 1',main=sprintf('Chimp and bonobo %s',tsneType))
    #points(tsnes[[1]]$Y[selectSamples$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
    legend(
      par('usr')[2]+.01*diff(par('usr')[1:2]), 
      mean(par('usr')[3:4]),
      c(names(malariaCols),names(areaCols),names(speciesPch)),
      col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch)))),
      pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch),
      pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch))),
      inset=.01,pt.lwd=3,pt.cex=2.5,
      xjust=0,xpd=NA
    )
    #text(tsnes[[1]]$Y,selectSamples$Code,cex=.25)
    #bonobo only
    plot(tsnes[[2]]$Y,pch=speciesPch[selectSamples2$Species],bg=areaCols[selectSamples2$area2],col=malariaCols[selectSamples2$malaria+1],cex=2.2,lwd=2.5,ylab='t-SNE 2',xlab='t-SNE 1',main='Bonobo')
    points(tsnes[[2]]$Y[selectSamples2$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
    legend(
      par('usr')[2]+.01*diff(par('usr')[1:2]), 
      mean(par('usr')[3:4]),
      c(names(malariaCols),names(areaCols),names(speciesPch),'SIVPos'),
      col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch))),'#FF000099'),
      pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch,1),
      pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch)),NA),
      inset=.01,pt.lwd=3,pt.cex=2.5,
      xjust=0,xpd=NA
    )
    text(tsnes[[2]]$Y,selectSamples2$Code,cex=.25)
    #bonobo only
    plot(tsnes[[3]]$Y,pch=speciesPch[selectSamples3$Species],bg=areaCols[selectSamples3$area2],col=malariaCols[selectSamples3$malaria+1],cex=2.2,lwd=2.5,ylab='t-SNE 2',xlab='t-SNE 1',main='Chimp')
    points(tsnes[[3]]$Y[selectSamples3$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
    legend(
      par('usr')[2]+.01*diff(par('usr')[1:2]), 
      mean(par('usr')[3:4]),
      c(names(malariaCols),names(areaCols),names(speciesPch),'SIVPos'),
      col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch))),'#FF000099'),
      pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch,1),
      pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch)),NA),
      inset=.01,pt.lwd=3,pt.cex=2.5,
      xjust=0,xpd=NA
    )
    text(tsnes[[3]]$Y,selectSamples3$Code,cex=.25)
  }
dev.off()
system('pdftk out/tsne.pdf cat 1 output out/Fig.S6C.pdf')

