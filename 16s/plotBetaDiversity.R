library(vipor)
library(dnar)
if(!exists('brayDist'))source('plotPcoa.R')

# TL-E vs TL-W
# TL vs TL
# TL vs other sites
# other sites vs other sites

pullDists<-function(xx,distMat){
  isIdentical<-length(xx[[1]])==length(xx[[2]]) && all(xx[[1]]==xx[[2]])
  select<-distMat[xx[[1]],xx[[2]]]
  if(isIdentical)dists<-select[upper.tri(select)]
  else dists<-as.vector(select)
  return(dists)
}


comparisons<-withAs(s=samples[samples$isEnough,],list(
  list(
    'Within bonobo'=list(s[s$bonobo,'name'],s[s$bonobo,'name']),
    'Within chimpanzee'=list(s[!s$bonobo,'name'],s[!s$bonobo,'name']),
    'Between bonobo\nand chimpanzee'=list(s[s$bonobo,'name'],s[!s$bonobo,'name'])
  ),list(
    'Within TL2'=list(s[s$isTL,'name'],s[s$isTL,'name']),
    'Between TL2-E\nand TL2-W'=list(s[s$area=='TL-E','name'],s[s$area=='TL-W','name']),
    'Between TL2 plasmodium\npositive and negative'=list(s[s$malaria&s$isTL,'name'],s[s$isTL&!s$malaria,'name'])
  ),list(
    'Within non-\nendemic field sites'=list(0,0),
    'Between non-\nendemic field sites'=list(0,0),
    'Between TL2 and\nnon-endemic field sites'=list(s[s$isTL&s$bonobo,'name'],s[!s$isTL&s$bonobo,'name'])
  ),list(
    'Within plasmodium\nnegative chimpanzees'=list(s[!s$bonobo&!s$malaria,'name'],s[!s$bonobo&!s$malaria,'name']),
    'Within plasmodium\npositive chimpanzees'=list(s[!s$bonobo&s$malaria,'name'],s[!s$bonobo&s$malaria,'name']),
    'Between plasmodium\nnegative and positive\nchimpanzees'=list(s[!s$bonobo&s$malaria,'name'],s[!s$bonobo&!s$malaria,'name'])
  )
))
# list(
#   'Within plasmodium\nnegative bonobos'=list(s[s$bonobo&!s$malaria,'name'],s[s$bonobo&!s$malaria,'name']),
#   'Within plasmodium\npositive bonobos'=list(s[s$bonobo&s$malaria,'name'],s[s$bonobo&s$malaria,'name']),
#   'Between plasmodium\nnegative and\npositive bonobos'=list(s[s$bonobo&s$malaria,'name'],s[s$bonobo&!s$malaria,'name'])
# )
nonTL<-unique(samples[samples$bonobo&!samples$isTL,'area'])
names(nonTL)<-nonTL
allCombo<-unique(t(apply(expand.grid(nonTL,nonTL),1,sort)))
allCombo<-allCombo[allCombo[,1]!=allCombo[,2],]
betweenSites<-apply(allCombo,1,function(xx,distMat)pullDists(list(samples[samples$isEnough&samples$area==xx[1],'name'],samples[samples$isEnough&samples$area==xx[2],'name']),distMat),as.matrix(brayDist))
withinSites<-lapply(nonTL,function(xx,distMat)pullDists(list(samples[samples$isEnough&samples$area==xx,'name'],samples[samples$isEnough&samples$area==xx,'name']),distMat),as.matrix(brayDist))
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
spacer<-.5
pdf('out/dists.pdf',width=6,height=6)
  par(mar=c(2.5,10.7,.1,.1),lheight=.9)
  compareFactor<-factor(rep(unlist(lapply(distList,names)),unlist(lapply(distList,sapply,length))),levels=unlist(lapply(distList,names)))
  stats<-boxplot(unlist(distList)~compareFactor,range=Inf,notch=TRUE,plot=FALSE)
  pos<-sum(sapply(distList,length)):1-rep((1:length(distList)-1)*spacer,sapply(distList,length))
  names(pos)<-levels(compareFactor)
  pVals$top<-pos[pVals$x]
  pVals$bottom<-pos[pVals$y]
  pVals$row<-stackRegions(pVals$bottom,pVals$top)
  pVals$middle<-apply(pVals[,c('bottom','top')],1,mean)
  pVals$xPos<-1.02+.025*(pVals$row-1)
  plot(1,1,type='n',ylim=range(pos)+c(-.5,.5),xlim=c(min(unlist(distList)),1.1),yaxt='n',ylab='',xlab='Bray-Curtis dissimilarity',mgp=c(1.5,.6,0),tcl=-.3,yaxs='i')
  for(ii in ncol(stats$stats):1){
    rawDists<-distList[[groupId[ii]]][[stats$names[ii]]]
    #points(rawDists,pos[ii]+offsetX(rawDists),cex=.5,pch=21,col=NA,bg=cols[groupId[ii]])
    xCoords<-c(stats$stats[2,ii],stats$conf[1,ii],stats$stats[3,ii],stats$conf[2,ii],stats$stats[4,ii])
    yCoords<-c(.4,.4,.1,.4,.4)
    segments(stats$stats[1,ii],pos[ii],stats$stats[5,ii],pos[ii],lwd=3,col=cols[groupId[ii]])
    polygon(c(xCoords,rev(xCoords)),c(yCoords,-rev(yCoords))+pos[ii],col=cols[groupId[ii]])
    segments(xCoords[3],pos[ii]+yCoords[3],xCoords[3],pos[ii]-yCoords[3])
  }
  text(pVals$xPos+.003,pVals$middle,pVals$sig,srt=90,adj=c(0.5,1))
  segments(pVals$xPos,pVals$bottom,pVals$xPos,pVals$top)
  breaks<-which(c(FALSE,pos[-1]-pos[-length(pos)]< -1))
  #abline(h=sapply(breaks,function(xx)mean(c(pos[xx],pos[xx-1]))))
  axis(2,pos,names(pos),las=1,mgp=c(3,.7,0))
dev.off()



#looks like mostly bonobo_lg4300
zz<-which(as.matrix(uniDist)>.85,arr.ind=TRUE)
zzz<-cbind(labels(uniDist)[zz[,1]],labels(uniDist)[zz[,2]])
zzz[apply(zzz,1,function(xx)!any(grepl('LG4300',xx))),]

