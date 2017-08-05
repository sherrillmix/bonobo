if(!exists('swarmData'))source("loadData.R")
library(phyloseq)
library(ape)
library(vipor)
library(Rtsne)
library(vegan)
library(GUniFrac)
library(ade4)
source('16s/myBiplot.R')

tlAdonis<-interactAdonis<-plantAdonis<-chimpAdonis<-bonoboAdonis<-list()
mantels<-list()
unis<-list()
for(ii in names(swarmData)){
  plotProp<-swarmData[[ii]][['props']][swarmData[[ii]][['isEnough']]&rownames(swarmData[[ii]][['props']]) %in% rownames(samples),]
  plotProp2<-swarmData[[ii]][['rare']][swarmData[[ii]][['isEnough']]&rownames(swarmData[[ii]][['rare']]) %in% rownames(samples),]
  phyOtuW<-otu_table(plotProp,taxa_are_rows=FALSE)
  phyOtuU<-otu_table(plotProp2,taxa_are_rows=FALSE)
  qiimeDataW<-phyloseq(otu_table=phyOtuW,phy_tree=swarmData[[ii]][['tree']])
  qiimeDataU<-phyloseq(otu_table=phyOtuU,phy_tree=swarmData[[ii]][['tree']])
  brayDist<-distance(qiimeDataU,'bray',binary=TRUE)
  brayDistW<-distance(qiimeDataW,'bray',binary=FALSE)
  uniDist<-UniFrac(qiimeDataU,weighted=FALSE)
  unis[[ii]]<-uniDist
  uniDistW<-UniFrac(qiimeDataW,weighted=TRUE)
  #uniDistG<-cacheOperation(sprintf('work/gunifrac_%s.Rdat',ii),GUniFrac,plotProp,swarmData[[ii]][['tree']])$unifracs[,,2]
  mantels[[ii]]<-list(
    'uniW'=ade4::mantel.rtest(uniDist,uniDistW,nrepet=1e4),
    'brayW'=ade4::mantel.rtest(uniDist,brayDistW,nrepet=1e4),
    'brayUW'=ade4::mantel.rtest(uniDist,brayDist,nrepet=1e4)
  )
  brayPca<-pcoa(brayDist)
  uniPca<-pcoa(uniDist)
  uniPcaW<-pcoa(uniDistW)
  #uniPcaG<-pcoa(uniDistG)
  tsneBray<-Rtsne(brayDist,is_distance=TRUE,verbose=TRUE,perplexity=10,max_iter=3000)
  tsneUni<-Rtsne(uniDist,is_distance=TRUE,verbose=TRUE,perplexity=15,max_iter=3000)
  tsneUniW<-Rtsne(uniDist,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
  selectDist<-uniDist
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
  names(areaCols)<-unique(selectSamples$area2[order(selectSamples$chimpBonobo)])
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

  pdf(sprintf('work/check_%s.pdf',ii));
    tmp<-0+(plotProp2>0)
    tmp<-tmp[,apply(tmp,2,sum)>1]
    uniTree<-hclust(uniDist)
    heatmap(0+(tmp>0),col=rev(heat.colors(100)),mar=c(5,5),scale='none',Rowv=as.dendrogram(uniTree))
    pos<-my.biplot.pcoa(selectPca,plotProp2>0,plot.axes=1:2,pch=speciesPch[selectSamples$chimpBonobo],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5,mar=c(4,4,1.5,10),arrowsFilter=1)
    text(pos,rownames(plotProp2),cex=.25)
    #plot(pos[,1],plotProp2[,'236'],xlab='PCoA 1',ylab='Morus 235 OTU counts')
    plot(swarmData[[ii]][['tree']])
    #prob236<-swarmData[['matK']][['taxa']]['236','seq']
    #matkSample<-swarmData[['matK']][['taxa']][1:20,'seq']
    #rbclSample<-swarmData[['rbcL']][['taxa']][1:20,'seq'] 
    #matkAlign<-levenAlign(matkSample,prob236)
    #rbclAlign<-levenAlign(rbclSample,prob236)
    #plotDNA(unlist(matkAlign),groups=c('Morus 236',rep('matk',20)))
    #plotDNA(unlist(rbclAlign),groups=c('Morus 236',rep('rbcl',20)))
 dev.off();

  #dists<-leven(swarmData[[1]][['taxa']][,'seq'],swarmData[[1]][['taxa']]['236','seq'],nThreads=50)
  #dists1<-leven(substring(swarmData[[1]][['taxa']][,'seq'],1,220),substring(swarmData[[1]][['taxa']]['236','seq'],1,220),nThreads=50)
  #dists2<-leven(substring(swarmData[[1]][['taxa']][,'seq'],240),substring(swarmData[[1]][['taxa']]['236','seq'],240),nThreads=50)
  #dists3<-leven(substring(swarmData[[1]][['taxa']][,'seq'],1,100),substring(swarmData[[1]][['taxa']]['236','seq'],1,100),nThreads=50)
  #distsRev<-leven(swarmData[[1]][['taxa']][,'seq'],revComp(swarmData[[1]][['taxa']]['236','seq']),nThreads=50)
  #rbcl<-leven(substring(swarmData[[2]][['taxa']][,'seq'],1),substring(swarmData[[ii]][['taxa']]['236','seq'],1),nThreads=50)
  #pdf('test.pdf',height=90,width=30);plot(phyloTree);hist(dists);dev.off()
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

  predictors<-model.matrix(~0+Species+malaria+SIV+area,selectSamples)
  print(table(selectSamples$malaria,selectSamples$chimpBonobo))
  pdf(sprintf('out/pcoa_%s.pdf',ii),width=6,height=6)
    pos<-my.biplot.pcoa(selectPca,predictors,plot.axes=1:2,pch=21,bg=speciesCols[selectSamples$chimpBonobo],col=malariaCols3[selectSamples$malaria+1],cex=2.25,lwd=4,arrowsFilter=Inf,las=1,mgp=c(2.75,.75,0),sameAxis=FALSE,bty='l',type='n')
    points(pos[!selectSamples$malaria,],col=malariaCols3[1],cex=2.25,lwd=4,bg=speciesCols[selectSamples[!selectSamples$malaria,'chimpBonobo']],pch=21)
    points(pos[selectSamples$malaria,],col=malariaCols3[2],cex=2.25,lwd=4,bg=speciesCols[selectSamples[selectSamples$malaria,'chimpBonobo']],pch=21)
    #text(posl,1],pos[,2],1:nrow(selectSamples))
    #tmp<-0+(plotProp2>0)
    #tmp<-plotProp2
    #colnames(tmp)<-swarmData[[ii]][['taxa']][colnames(tmp),'bestId']
    #biplot(selectPca,tmp)
    #biplot(selectPca,predictors)
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

  ss<-samples[labels(selectDist),]
  plantAdonis[[ii]]<-cacheOperation(sprintf('work/adonis_%s.Rdat',ii),adonis,selectDist~bonobo+area2+malaria,data=ss,permutations=1e7,parallel=5)
  interactAdonis[[ii]]<-cacheOperation(sprintf('work/interactAdonis_%s.Rdat',ii),adonis,selectDist~bonobo+area2+malaria*bonobo,data=ss,permutations=1e7,parallel=5)
  chimpDist<-as.matrix(selectDist)
  chimpDist<-as.dist(chimpDist[!ss$bonobo,!ss$bonobo])
  chimpAdonis[[ii]]<-cacheOperation(sprintf('work/adonisChimp_%s.Rdat',ii),adonis,chimpDist~area2+malaria,data=ss[!ss$bonobo,],permutations=1e7,parallel=5)
  bonoboDist<-as.matrix(selectDist)
  bonoboDist<-as.dist(bonoboDist[ss$bonobo,ss$bonobo])
  bonoboAdonis[[ii]]<-cacheOperation(sprintf('work/adonisBonobo_%s.Rdat',ii),adonis,bonoboDist~area2+malaria,data=ss[ss$bonobo,],permutations=1e7,parallel=5)
  tlDist<-as.dist(as.matrix(selectDist)[ss$isTL,ss$isTL])
  tlAdonis[[ii]]<-cacheOperation(sprintf('work/adonisTL_%s.Rdat',ii),adonis,tlDist~malaria,data=ss[ss$isTL,],permutations=1e7,parallel=5)

  uniTree<-hclust(uniDist)
  pcoaTree<-hclust(dist(selectPca$vectors[,1:2]))
  phyloTree<-swarmData[[ii]][['tree']]
  pdf(sprintf('out/tree_%s.pdf',ii),height=12,width=10)
    par(mar=c(4,0,0,4))
    plot(as.dendrogram(uniTree),horiz=TRUE,xlab='UniFrac distance')
    plot(as.dendrogram(pcoaTree),horiz=TRUE,xlab='First two PCoA distance')
    tmp<-0+(plotProp2>0)
    tmp<-tmp[,apply(tmp,2,sum)>1]
    if(ncol(tmp)>2000) tmp<-tmp[,apply(tmp,2,sum)>5]
    heatmap(0+(tmp>0),col=rev(heat.colors(100)),mar=c(5,5),scale='none',Rowv=as.dendrogram(uniTree))
    tmp<-tmp[,orderIn(colnames(tmp),phyloTree$tip.label)]
    heatmap(0+(tmp>0),col=rev(heat.colors(100)),mar=c(5,5),scale='none',Rowv=as.dendrogram(uniTree),Colv=NA)
  dev.off()

  pullDists<-function(xx,distMat){
    isIdentical<-length(xx[[1]])==length(xx[[2]]) && all(xx[[1]]==xx[[2]])
    select<-distMat[xx[[1]],xx[[2]]]
    if(isIdentical)dists<-select[upper.tri(select)]
    else dists<-as.vector(select)
    return(dists)
  }
  isEnough<-swarmData[[ii]][['isEnough']][rownames(samples)]
  comparisons<-withAs(s=samples[isEnough,],list(
    list(
      'Within bonobo'=list(s[s$bonobo,'Code'],s[s$bonobo,'Code']),
      'Within chimpanzee'=list(s[!s$bonobo,'Code'],s[!s$bonobo,'Code']),
      'Between bonobo\nand chimpanzee'=list(s[s$bonobo,'Code'],s[!s$bonobo,'Code'])
    ),list(
      'Within non-\nendemic field sites'=list(0,0),
      'Between non-\nendemic field sites'=list(0,0),
      'Between TL2 and\nnon-endemic field sites'=list(s[s$isTL&s$bonobo,'Code'],s[!s$isTL&s$bonobo,'Code'])
    ),list(
      'Within TL2 Laverania negative'=list(s[s$isTL&!s$malaria,'Code'],s[s$isTL&!s$malaria,'Code']),
      'Within TL2 Laverania positive'=list(s[s$isTL&s$malaria,'Code'],s[s$isTL&s$malaria,'Code']),
      'Between TL2 Laverania\npositive and negative'=list(s[s$isTL&s$malaria,'Code'],s[s$isTL&!s$malaria,'Code'])
    #),list(
    #  'Within BI Laverania negative'=list(s[s$area=='BI'&!s$malaria,'Code'],s[s$area=='BI'&!s$malaria,'Code']),
    #  'Within BI Laverania positive'=list(s[s$area=='BI'&s$malaria,'Code'],s[s$area=='BI'&s$malaria,'Code']),
    #  'Between BI Laverania\npositive and negative'=list(s[s$area=='BI'&s$malaria,'Code'],s[s$area=='BI'&!s$malaria,'Code'])
    #),list(
    #  'Within UB Laverania negative'=list(s[s$area=='UB'&!s$malaria,'Code'],s[s$area=='UB'&!s$malaria,'Code']),
    #  'Within UB Laverania positive'=list(s[s$area=='UB'&s$malaria,'Code'],s[s$area=='UB'&s$malaria,'Code']),
    #  'Between UB Laverania\npositive and negative'=list(s[s$area=='UB'&s$malaria,'Code'],s[s$area=='UB'&!s$malaria,'Code'])
    ),list(
      'Within Laverania\nnegative chimpanzees'=list(s[!s$bonobo&!s$malaria,'Code'],s[!s$bonobo&!s$malaria,'Code']),
      'Within Laverania\npositive chimpanzees'=list(s[!s$bonobo&s$malaria,'Code'],s[!s$bonobo&s$malaria,'Code']),
      'Between Laverania negative\n and positive chimpanzees'=list(s[!s$bonobo&s$malaria,'Code'],s[!s$bonobo&!s$malaria,'Code'])
    )
  ))
  nonTL<-unique(samples[samples$bonobo&!samples$isTL,'area'])
  names(nonTL)<-nonTL
  allCombo<-unique(t(apply(expand.grid(nonTL,nonTL),1,sort)))
  allCombo<-allCombo[allCombo[,1]!=allCombo[,2],]
  betweenSites<-apply(allCombo,1,function(xx,distMat)pullDists(list(samples[isEnough&samples$area==xx[1],'Code'],samples[isEnough&samples$area==xx[2],'Code']),distMat),as.matrix(selectDist))
  withinSites<-lapply(nonTL,function(xx,distMat)pullDists(list(samples[isEnough&samples$area==xx,'Code'],samples[isEnough&samples$area==xx,'Code']),distMat),as.matrix(selectDist))

  distList<-lapply(comparisons,function(xx)lapply(xx,pullDists,as.matrix(selectDist)))
  distList<-lapply(distList,function(xx){
    names(xx)<-ifelse(nchar(names(xx))>20,sub(' vs ',' vs\n',names(xx)),names(xx))
    if(any(names(xx)=='Between non-\nendemic field sites'))xx[['Between non-\nendemic field sites']]<-unlist(betweenSites)
    if(any(names(xx)=='Within non-\nendemic field sites'))xx[['Within non-\nendemic field sites']]<-unlist(withinSites)
    return(xx)
  })

  pVals<-lapply(distList,function(dists)outer(dists,dists,function(xx,yy)mapply(function(xxx,yyy){wilcox.test(xxx,yyy)$p.value},xx,yy)))
  pVals<-do.call(rbind,lapply(pVals,function(xx){
    n<-nrow(xx)
    cols<-colnames(xx)[matrix(1:n,nrow=n,ncol=n,byrow=TRUE)[upper.tri(xx)]]
    rows<-rownames(xx)[matrix(1:n,nrow=n,ncol=n,byrow=FALSE)[upper.tri(xx)]]
    ps<-xx[upper.tri(xx)]
    data.frame('x'=rows,'y'=cols,'p'=ps,stringsAsFactors=FALSE)
  }))
  pVals$sig<- symnum(pVals$p, corr = FALSE, na = FALSE, cutpoints = c(0, 1e-6, 1e-4, 1e-2,1), symbols = c("***", "**", "*",''))
  pVals<-pVals[pVals$p<.01,]

  cols<-rainbow.lab(length(distList),start=1,end=-2)
  groupId<-rep(1:length(distList),sapply(distList,length))
  spacer<-.6
  pdf(sprintf('out/dists_%s.pdf',ii),width=9,height=6)
    par(mar=c(9.2,2.75,1.5,4),lheight=.8)
    compareFactor<-factor(rep(unlist(lapply(distList,names)),unlist(lapply(distList,sapply,length))),levels=unlist(lapply(distList,function(xx)rev(names(xx)))))
    stats<-boxplot(unlist(distList)~compareFactor,range=Inf,notch=TRUE,plot=FALSE)
    betaCI<-tapply(unlist(distList),compareFactor,function(xx)medianCI(xx))
    pos<-sum(sapply(distList,length)):1-rep((1:length(distList)-1)*spacer,sapply(distList,length))
    names(pos)<-levels(compareFactor)
    pVals$top<-pos[pVals$x]
    pVals$bottom<-pos[pVals$y]
    pVals$row<-stackRegions(pVals$top,pVals$bottom)
    pVals$middle<-apply(pVals[,c('bottom','top')],1,mean)
    pVals$xPos<-1.03+.03*(pVals$row-1)
    plot(1,1,type='n',xlim=range(pos)+c(-.5-spacer,.5+spacer),ylim=c(min(unlist(distList)),1.08),xaxt='n',xlab='',ylab='UniFrac distance',mgp=c(1.75,.4,0),tcl=-.3,xaxs='i',las=1,bty='l',main=ii)
    for(ii in ncol(stats$stats):1){
      rawDists<-distList[[groupId[ii]]][[stats$names[ii]]]
      #points(rawDists,pos[ii]+offsetX(rawDists),cex=.5,pch=21,col=NA,bg=cols[groupId[ii]])
      thisCI<-betaCI[[stats$names[ii]]]
      #xCoords<-c(stats$stats[2,ii],stats$conf[1,ii],stats$stats[3,ii],stats$conf[2,ii],stats$stats[4,ii])
      xCoords<-c(stats$stats[2,ii],thisCI[1],stats$stats[3,ii],thisCI[2],stats$stats[4,ii])
      yCoords<-c(.4,.4,.1,.4,.4)
      segments(pos[ii],stats$stats[1,ii],pos[ii],stats$stats[5,ii],lwd=3,col=cols[groupId[ii]])
      polygon(c(yCoords,-rev(yCoords))+pos[ii],c(xCoords,rev(xCoords)),col=cols[groupId[ii]])
      segments(pos[ii]+yCoords[3],xCoords[3],pos[ii]-yCoords[3],xCoords[3])
    }
    text(pVals$middle,pVals$xPos+.005,pVals$sig,adj=c(.5,0.5))
    segments(pVals$bottom,pVals$xPos,pVals$top,pVals$xPos)
    breaks<-which(c(FALSE,pos[-1]-pos[-length(pos)]< -1))
    #abline(h=sapply(breaks,function(xx)mean(c(pos[xx],pos[xx-1]))))
    slantAxis(1,pos,names(pos),srt=-40)
  dev.off()

  if(FALSE){
    tmp<-0+(plotProp2>0)
    tmp<-tmp[,apply(tmp,2,sum)>1]
    if(ncol(tmp)>2000) tmp<-tmp[,apply(tmp,2,sum)>5]
    heatmap(tmp,col=rev(heat.colors(100)),mar=c(5,5),scale='none',Rowv=as.dendrogram(uniTree))
    points.stand <- scale(selectPca$vectors[, 1:2])
    S <- cov(0+(plotProp2>0), points.stand)
    rownames(S)<-swarmData[[ii]][['taxa']][rownames(S),'bestId']
    weight<-apply(S,1,function(x)sqrt(sum(x^2)))
    weightSelect<-which(weight>.05)
    weightSelect<-weightSelect[order(weightSelect[weightSelect])]
    tmp<-plotProp2
    colnames(tmp)<-swarmData[[ii]][['taxa']][colnames(tmp),'bestId']
    heatmap(0+(tmp[,weightSelect]>0),col=rev(heat.colors(100)),mar=c(5,5),scale='none',Rowv=as.dendrogram(pcoaTree))
    biplot(prcomp(plotProp2>0))
  }
}
system('cp out/pcoa_matK.pdf out/Fig.5A.pdf')
system('cp out/pcoa_rbcL.pdf out/Fig.5B.pdf')
system('cp out/tsne_matK.pdf out/Fig.S6A.pdf')
system('cp out/tsne_rbcL.pdf out/Fig.S6B.pdf')

table(samples[names(swarmData[['matK']][['isEnough']])[swarmData[['matK']][['isEnough']]],c('chimpBonobo','plasmoPM')])
table(samples[names(swarmData[['rbcL']][['isEnough']])[swarmData[['rbcL']][['isEnough']]],c('chimpBonobo','plasmoPM')])

tmp<-unis[[1]]
tmp2<-unis[[2]]
shared<-intersect(labels(tmp),labels(tmp2))
tmp<-as.dist(as.matrix(tmp)[shared,shared])
tmp2<-as.dist(as.matrix(tmp2)[shared,shared])
mantel.rtest(tmp,tmp2)
