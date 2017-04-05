
library(dnar)
if(!exists('otuTab'))source('runQiime.R')
source('functions.R')

tlGroups<-withAs(s=samples[samples$bonobo&samples$isEnough&samples$isTL,],tapply(s$name,s$area2,c))
nonTlGroups<-withAs(s=samples[samples$bonobo&samples$isEnough&!samples$isTL,],tapply(s$name,s$area2,c))
nonTlGroups<-nonTlGroups[sapply(nonTlGroups,length)>5]

comparisons<-expand.grid('tl'=names(tlGroups),'nontl'=names(nonTlGroups))

inBonobo<-otuProp[apply(otuProp[,c(unlist(tlGroups),unlist(nonTlGroups))],1,max)>.001,]
ninePs<-apply(inBonobo,1,function(xx)apply(comparisons,1,function(select)wilcox.test(xx[tlGroups[[select[1]]]],xx[nonTlGroups[[select[2]]]],alternative='less')$p.value))

condenseP<-p.adjust(apply(ninePs,2,fishers,correct=3),'bonferroni')

pCut<-.01
selectPropAll<-apply(inBonobo[condenseP<pCut,samples$name[samples$isEnough]],1,function(x)x/max(x))
rownames(selectPropAll)<-sprintf('%s%s',ifelse(samples[rownames(selectPropAll),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(selectPropAll)))
colnames(selectPropAll)<-sprintf('%s (%s)\nq=%0.3f',colnames(selectPropAll),taxa[colnames(selectPropAll),'bestId'],condenseP[colnames(selectPropAll)])
breaks<-c(-1e-6,seq(min(selectPropAll[selectPropAll>0])-1e-10,max(selectPropAll)+1e-10,length.out=100))
cols<-c('white',tail(rev(heat.colors(110)),99)) 

pdf('out/nineCompare.pdf',height=13,width=12)
  par(mar=c(11.5,.1,3,8.5),lheight=.7)
  plotHeat(selectPropAll,breaks,cols)
dev.off()


