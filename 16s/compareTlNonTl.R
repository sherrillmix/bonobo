
library(dnar)
if(!exists('otuTab'))source('runQiime.R')
source('functions.R')
ss<-samples[order(!samples$bonobo,samples$area2,samples$malaria),]

tlGroups<-withAs(s=ss[ss$bonobo&ss$isEnough&ss$isTL,],tapply(s$name,s$area2,c))
nonTlGroups<-withAs(s=ss[ss$bonobo&ss$isEnough&!ss$isTL,],tapply(s$name,s$area2,c))
nonTlGroups<-nonTlGroups[sapply(nonTlGroups,length)>5]

comparisons<-expand.grid('tl'=names(tlGroups),'nontl'=names(nonTlGroups))

inBonobo<-otuProp[apply(otuProp[,c(unlist(tlGroups),unlist(nonTlGroups))],1,max)>.001,]
ninePs<-apply(inBonobo,1,function(xx)apply(comparisons,1,function(select)suppressWarnings(wilcox.test(xx[tlGroups[[select[1]]]],xx[nonTlGroups[[select[2]]]],alternative='less'))$p.value))

ninePsGt<-apply(inBonobo,1,function(xx)apply(comparisons,1,function(select)suppressWarnings(wilcox.test(xx[tlGroups[[select[1]]]],xx[nonTlGroups[[select[2]]]],alternative='greater')$p.value)))

ninePsBoth<-apply(inBonobo,1,function(xx)apply(comparisons,1,function(select)suppressWarnings(wilcox.test(xx[tlGroups[[select[1]]]],xx[nonTlGroups[[select[2]]]]))$p.value))
ninePsBoth[is.na(ninePsBoth)]<-1

condenseP<-p.adjust(apply(ninePs,2,fishers,correct=3),'fdr')
condensePGt<-p.adjust(apply(ninePsGt,2,fishers,correct=3),'fdr')
condensePBoth<-p.adjust(apply(ninePsBoth,2,fishers,correct=3),'fdr')


pCut<-.05
selectPropAll<-apply(inBonobo[condenseP<pCut,ss$name[ss$isEnough]],1,function(x)x/max(x))
#rownames(selectPropAll)<-sprintf('%s%s',ifelse(ss[rownames(selectPropAll),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(selectPropAll)))
selectPropAll<-selectPropAll[,order(condenseP[condenseP<pCut])]
colnames(selectPropAll)<-sprintf('%s q=%0.3f',taxa[colnames(selectPropAll),'bestId'],condenseP[colnames(selectPropAll)])
#
selectPropAllGt<-apply(inBonobo[condensePGt<pCut,ss$name[ss$isEnough]],1,function(x)x/max(x))
selectPropAllGt<-selectPropAllGt[,order(condensePGt[condensePGt<pCut])]
colnames(selectPropAllGt)<-sprintf('%s q=%0.3f',taxa[colnames(selectPropAllGt),'bestId'],condensePGt[colnames(selectPropAllGt)])
#
selectPropAllBoth<-apply(inBonobo[condensePBoth<pCut,ss$name[ss$isEnough]],1,function(x)x/max(x))
selectPropAllBoth<-selectPropAllBoth[,order(condensePBoth[condensePBoth<pCut])]
colnames(selectPropAllBoth)<-sprintf('%s q=%0.3f',taxa[colnames(selectPropAllBoth),'bestId'],condensePBoth[colnames(selectPropAllBoth)])

breaks<-c(-1e-6,seq(min(selectPropAll[selectPropAll>0])-1e-10,max(selectPropAll)+1e-10,length.out=100))
cols<-c('white',tail(rev(heat.colors(110)),99)) 

pdf('out/nineCompare.pdf',height=13,width=12)
  par(mar=c(11.5,.1,3,14.5),lheight=.7)
  metadata<-ss[rownames(selectPropAll),c('chimpBonobo','area2','plasmoPM','Code')]
  colnames(metadata)<-c('Species','Area','Laverania','Sample')
  plotHeat(selectPropAll,breaks,cols,yaxt='n')
  title(main='Depleted')
  addMetaData(metadata,cex=.75)
  plotHeat(selectPropAllGt,breaks,cols,yaxt='n')
  title(main='Enriched')
  addMetaData(metadata,cex=.75)
  #plotHeat(selectPropAllBoth,breaks,cols,yaxt='n')
  #title(main='Both')
  #addMetaData(metadata,cex=.75)
dev.off()

sigTaxa<-names(condensePBoth)[condensePBoth<pCut]
seqDists<-outer(strsplit(taxa[sigTaxa,'seq'],''),strsplit(taxa[sigTaxa,'seq'],''),function(xx,yy)mapply(function(x,y)sum(x!=y),xx,yy))
#seqDists2<-leven(taxa[sigTaxa,'seq'])
rownames(seqDists)<-colnames(seqDists)<-taxa[sigTaxa,'bestId']
breaks<-seq(min(seqDists),max(seqDists),length.out=100)
rownames(seqDists2)<-colnames(seqDists2)<-taxa[sigTaxa,'bestId']
pdf('out/depleteDist.pdf')
  par(mar=c(4,1,1,11))
  plot(as.dendrogram(hclust(as.dist(seqDists))),horiz=TRUE,xlab='Hamming distance between aligned sequences')
  #plot(as.dendrogram(hclust(as.dist(seqDists2))),horiz=TRUE,xlab='Edit distance between aligned sequences')
  #insetScale(round(breaks,6),cols,c(.97,.01,.98,.25),label='Proportion of OTU maximum')
dev.off()
