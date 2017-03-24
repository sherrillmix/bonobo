if(!exists('otuTab'))source('plotPcoa.R')
splitDist<-function(dists,splits){
  uniqSplits<-unique(splits)
  names(uniqSplits)<-uniqSplits
  splitDists<-lapply(uniqSplits,function(xx){
    lapply(uniqSplits,function(yy){
      out<-dists[splits==xx,splits==yy]
      if(xx==yy)return(out[upper.tri(out)])
      else return(as.vector(out))
    })
  })
  out<-data.frame('x'=NA,'y'=NA,'dist'=NA,stringsAsFactors=FALSE)[0,]
  for(ii in uniqSplits){
    for(jj in uniqSplits){
      out<-rbind(out,data.frame('x'=ii,'y'=jj,'dist'=splitDists[[ii]][[jj]],stringsAsFactors=FALSE))
    }
  }
  return(out)
}

#BI0331 is from separate study
allEnough<-apply(otuTab,2,sum)>15000 & !grepl('Bonobo_BI0331',colnames(otuTab))
speciesOtu<-otu_table(apply(otuTab[rownames(otuTab) %in% tree$tip.label,allEnough],2,rarefyCounts,nRequiredReads),taxa_are_rows=TRUE)
speciesQiimeData<-phyloseq(otu_table=speciesOtu,phy_tree=tree)
speciesProp<-otu_table(otuProp[rownames(otuTab) %in% tree$tip.label,allEnough],taxa_are_rows=TRUE)
speciesQiimeDataW<-phyloseq(otu_table=speciesProp,phy_tree=tree)
#make sure tree is bifurcating or breaks UniFrac without error
speciesUniDist<-UniFrac(speciesQiimeData,weighted=TRUE)
speciesBray<-distance(speciesQiimeDataW,'bray')
#uniDistMat<-as.matrix(speciesUniDist)
uniDistMat<-as.matrix(speciesBray)
species<-sub('(Toddler)?_.*$','',rownames(uniDistMat))
splits<-splitDist(uniDistMat,species)
meanDists<-withAs(splits=splits[splits$x=="Bonobo"&splits$y!='Primate',],tapply(splits$dist,paste(splits$x,splits$y,sep='\n'),mean))
speciesCols<-rainbow.lab(length(unique(species)))
speciesPch<-21+rep(0:2,length.out=length(unique(species)))
names(speciesPch)<-names(speciesCols)<-unique(species)
speciesTsne<-Rtsne(speciesBray,is_distance=TRUE,verbose=TRUE,perplexity=20,max_iter=3000)
speciesPca<-pcoa(speciesBray)
pdf('out/species.pdf',width=10,height=10)
  par(lheight=.7,mar=c(9,4,.1,.1))
  heatmap(t(speciesProp[apply(speciesProp,1,max)>.05,]),col=rev(heat.colors(100)))
  withAs(splits=splits[splits$x=="Bonobo"&splits$y!='Primate',],vpPlot(factor(paste(splits$x,splits$y,sep='\n'),levels=names(sort(meanDists))),splits$dist,las=2,ylab='Pairwise UniFrac distance',las=2))
  heatmap(tapply(splits$dist,list(splits$x,splits$y),mean),scale='none',col=rev(heat.colors(100)))
  plot(speciesTsne$Y,pch=speciesPch[species],bg=speciesCols[species],col='#00000033',cex=2.2,ylab='t-SNE 2',xlab='t-SNE 1')
  text(speciesTsne$Y,rownames(uniDistMat),cex=.2)
  legend('topright',names(speciesCols),pch=speciesPch,pt.bg=speciesCols,col='#00000033',cex=.8)
  pos<-my.biplot.pcoa(speciesPca,matrix(1,nrow=nrow(as.matrix(speciesBray))),plot.axes=1:2,pch=speciesPch[species],bg=speciesCols[species],col='#00000033',cex=2.25,arrowsFilter=Inf)
  text(pos,rownames(uniDistMat),cex=.2)
  pos<-my.biplot.pcoa(speciesPca,matrix(1,nrow=nrow(as.matrix(speciesBray))),plot.axes=3:4,pch=speciesPch[species],bg=speciesCols[species],col='#00000033',cex=2.25,arrowsFilter=Inf)
  text(pos,rownames(uniDistMat),cex=.2)
dev.off()


