library(taxonomizr)
source('readSamples.R')
library(vipor)

antiMal<-read.csv('AntiMalPlant_List_021817.csv',stringsAsFactors=FALSE)
if(!exists('taxaNodes')){
  getNamesAndNodes('~/db/taxo/')
  taxaNodes<-read.nodes('~/db/taxo/nodes.dmp')
  taxaNamesAmbig<-read.names('~/db/taxo/names.dmp',FALSE)
  taxaNames<-read.names('~/db/taxo/names.dmp')
}
antiMal$clean<-sub(' sp$','',sub('([A-Za-z]+ [a-z]+).*','\\1',antiMal$Plant))
antiMal$clean[antiMal$clean=="AAloe parvibracteata"]<-"Aloe parvibracteata"
antiMal$clean[antiMal$clean=="Curcuma aromatic"]<-"Curcuma aromatica"
antiMal$clean[antiMal$clean=="Uvariopsis congensis"]<-"Uvariopsis"
antiMal$isGenus<-!grepl(' ',antiMal$clean)
antiMal$id<-getId(antiMal$clean,taxaNamesAmbig)
antiMal$species<-getTaxonomy(antiMal$id,taxaNodes,taxaNames,'species')
antiMal$genus<-getTaxonomy(antiMal$id,taxaNodes,taxaNames,'genus')
write.csv(antiMal,'work/antiMal.csv',row.names=FALSE)

hitCut<-.95
if(!exists('swarmData'))source("loadData.R")
antiMals<-lapply(names(swarmData),function(ii){
  species<-swarmData[[ii]][['taxa']][colnames(swarmData[[ii]][['props']]),'species']
  genus<-swarmData[[ii]][['taxa']][colnames(swarmData[[ii]][['props']]),'genus']
  bestHit<-withAs(sw=swarmData[[ii]][['taxa']][colnames(swarmData[[ii]][['props']]),],sw$bestScore/nchar(sw$seq))
  isAntiMal<-bestHit>hitCut &(species %in% antiMal$species[!is.na(antiMal$species)] | genus %in% antiMal$genus[antiMal$isGenus])
  isPresent<-(antiMal$species %in% species[bestHit>hitCut] &!is.na(antiMal$species))| (antiMal$genus %in% genus[bestHit>hitCut] &antiMal$isGenus)
  antiMalProp<-swarmData[[ii]][['props']][,isAntiMal]
  antiMalTaxa<-swarmData[[ii]][['taxa']][colnames(antiMalProp),]
  antiMalProp<-antiMalProp[swarmData[[ii]][['isEnough']][rownames(antiMalProp)],]
  antiMalPlants<-antiMal[isPresent,]
  return(list('prop'=antiMalProp,'taxa'=antiMalTaxa,'plants'=antiMalPlants))
})
names(antiMals)<-names(swarmData)

pdf('out/allAntiMal.pdf',width=10)
for(ii in names(antiMals)){
  propsAntiMal<-antiMals[[ii]][['prop']]
  propsAntiMal<-propsAntiMal[rownames(propsAntiMal) %in% rownames(samples),]
  sampleInfo<-samples[rownames(samples) %in% rownames(propsAntiMal),]
  sampleInfo$speciesMal<-paste(sampleInfo$Species,ifelse(sampleInfo$malaria,'Positive','Negative'))
  propsAntiMal<-propsAntiMal[rownames(sampleInfo),]
  propAntiMal<-apply(propsAntiMal,1,sum)
  #
  cols<-rainbow.lab(length(unique(sampleInfo$speciesMal)),alpha=.7)
  names(cols)<-unique(sampleInfo$speciesMal)
  #
  print(tapply(propAntiMal,sampleInfo$speciesMal,median)*100)
  tapply(propAntiMal,list(sampleInfo$area2,sampleInfo$speciesMal),mean)*100
  #
  pTL<-wilcox.test(
    propAntiMal[grepl('TL',sampleInfo$area)&sampleInfo$bonobo],
    propAntiMal[!grepl('TL',sampleInfo$area)&sampleInfo$bonobo],
  )$p.value
  print(pTL)
  vpPlot(factor(sampleInfo$area2,levels=unique(sampleInfo$area2)),propAntiMal*100,ylab='Percent antimalarial',las=1,main=sprintf('%s p(TL not equal other bonobos)=%0.3f',ii,pTL),bg=cols[sampleInfo$speciesMal],pch=21,cex=1.5)
  legend('topleft',names(cols),pch=21,pt.bg=cols,pt.cex=1.5)
}
dev.off()


pCut<-Inf #show all
for(ii in names(swarmData)){
  message(ii)
  #otuProp<-swarmData[[ii]][['props']][swarmData[[ii]][['isEnough']],]
  otuProp<-antiMals[[ii]][['prop']]
  otuProp<-otuProp[rownames(otuProp) %in% rownames(samples),]
  ss<-samples[rownames(otuProp),]
  ss<-ss[order(!ss$bonobo,ss$area2,ss$malaria),]
  #
  tlGroups<-withAs(s=ss[ss$bonobo&ss$isTL,],tapply(rownames(s),s$area2,c))
  nonTlGroups<-withAs(s=ss[ss$bonobo&!ss$isTL,],tapply(rownames(s),s$area2,c))
  nonTlGroups<-nonTlGroups[sapply(nonTlGroups,length)>5]
  comparisons<-expand.grid('tl'=names(tlGroups),'nontl'=names(nonTlGroups))
  inBonobo<-otuProp[,apply(otuProp[c(unlist(tlGroups),unlist(nonTlGroups)),],2,max)>0]
  message('nOTU: ',ncol(inBonobo),' nTaxa: ',length(unique(antiMals[[ii]][['taxa']][colnames(inBonobo),'species'])))
  ninePs<-apply(inBonobo,2,function(xx)apply(comparisons,1,function(select)suppressWarnings(wilcox.test(xx[tlGroups[[select[1]]]],xx[nonTlGroups[[select[2]]]],alternative='less'))$p.value))
  condenseP<-p.adjust(apply(ninePs,2,fishers,correct=3),'fdr')
  print(summary(condenseP<.05))
  simpleP<-apply(inBonobo,2,function(xx)suppressWarnings(wilcox.test(xx[unlist(tlGroups)],xx[unlist(nonTlGroups)],alternative='less')$p.value))
  swarmData[[ii]][['taxa']][names(condenseP)[condenseP<.05],]
  print(summary(simpleP<.05))
  if(any(condenseP<pCut)){
    selectPropAll<-apply(inBonobo[ss$Code,condenseP<pCut&apply(inBonobo>0,2,sum)>1],2,function(x)x/max(x))
    colnames(selectPropAll)<-sprintf('%s p=%0.3f',sub('^[a-z]_','',swarmData[[ii]][['taxa']][colnames(selectPropAll),'bestId']),condenseP[colnames(selectPropAll)])
    breaks<-c(-1e-6,seq(min(c(selectPropAll[selectPropAll>0],.01))-1e-10,max(selectPropAll)+1e-10,length.out=100))
    print(breaks)
    cols<-c('white',tail(rev(heat.colors(110)),99)) 
    maxTree<-hclust(dist(t(selectPropAll[,])))
    selectPropAll<-selectPropAll[,rev(maxTree$labels[maxTree$order])]
    selectPropAll<-selectPropAll[,order(colnames(selectPropAll))]
    pdf(sprintf('out/nineCompare_%s.pdf',ii),height=13,width=13)
      par(mar=c(10.5,.1,3,14),lheight=.7)
      metadata<-ss[rownames(selectPropAll),c('chimpBonobo','area2','plasmoPM','Code')]
      colnames(metadata)<-c('Species','Area','Laverania','Sample')
      plotHeat(selectPropAll,breaks,cols,yaxt='n',xaxt='n')
      title(main=ii)
      axis(1,1:ncol(selectPropAll),colnames(selectPropAll),las=2,cex.axis=.7)
      addMetaData(metadata,cex=.75)
    dev.off()
  }
}
system('cp out/nineCompare_matK.pdf out/Fig.S7A.pdf')
system('cp out/nineCompare_rbcL.pdf out/Fig.S7B.pdf')

#check chimp
pCut<-.05
for(ii in names(swarmData)){
  message(ii)
  otuProp<-swarmData[[ii]][['props']][swarmData[[ii]][['isEnough']],]
  #otuProp<-antiMals[[ii]][['prop']]
  otuProp<-otuProp[rownames(otuProp) %in% rownames(samples),]
  ss<-samples[rownames(otuProp),]
  ss<-ss[order(!ss$bonobo,ss$area2,ss$malaria),]
  #
  mal<-rownames(ss)[!ss$bonobo&ss$malaria]
  noMal<-rownames(ss)[!ss$bonobo&!ss$malaria]
  inChimp<-otuProp[,apply(otuProp[c(mal,noMal),],2,max)>0]
  ninePs<-apply(inChimp,2,function(xx)suppressWarnings(wilcox.test(xx[mal],xx[noMal],alternative='less'))$p.value)
  condenseP<-p.adjust(ninePs,'fdr')
  print(summary(condenseP<pCut))
}



