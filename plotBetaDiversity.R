if(!exists('swarmData'))source('loadData.R')

for(ii in names(swarmData)){

  nonTL<-unique(samples[samples$bonobo&!samples$isTL,'area'])
  names(nonTL)<-nonTL
  allCombo<-unique(t(apply(expand.grid(nonTL,nonTL),1,sort)))
  allCombo<-allCombo[allCombo[,1]!=allCombo[,2],]

  plotProp<-swarmData[[ii]][['props']][swarmData[[ii]][['isEnough']]&rownames(swarmData[[ii]][['props']]) %in% rownames(samples),]
  plotProp2<-swarmData[[ii]][['rare']][swarmData[[ii]][['isEnough']]&rownames(swarmData[[ii]][['rare']]) %in% rownames(samples),]
  phyOtuW<-otu_table(plotProp,taxa_are_rows=FALSE)
  qiimeDataW<-phyloseq(otu_table=phyOtuW,phy_tree=swarmData[[ii]][['tree']])
  phyOtuU<-otu_table(plotProp2,taxa_are_rows=FALSE)
  qiimeDataU<-phyloseq(otu_table=phyOtuU,phy_tree=swarmData[[ii]][['tree']])
  uniDist<-UniFrac(qiimeDataU,weighted=FALSE)
  brayDist<-distance(qiimeDataU,'bray',binary=TRUE)
  brayDistW<-distance(qiimeDataW,'bray',binary=FALSE)
  uniDistG<-cacheOperation(sprintf('work/gunifrac_%s.Rdat',ii),GUniFrac,plotProp,swarmData[[ii]][['tree']])$unifracs[,,2]

  betweenSites<-apply(allCombo,1,function(xx,distMat)pullDists(list(samples[isEnough&samples$area==xx[1],'Code'],samples[isEnough&samples$area==xx[2],'Code']),distMat),as.matrix(brayDist))
  withinSites<-lapply(nonTL,function(xx,distMat)pullDists(list(samples[isEnough&samples$area==xx,'Code'],samples[isEnough&samples$area==xx,'Code']),distMat),as.matrix(brayDist))


  distList<-lapply(comparisons,function(xx)lapply(xx,pullDists,as.matrix(brayDist)))
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
}
