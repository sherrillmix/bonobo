library(vipor)
library(dnar)
if(!exists('uniDist'))source('plotPcoa.R')

# TL-E vs TL-W
# TL vs TL
# TL vs other sites
# other sites vs other sites

comparisons<-withAs(s=samples[samples$isEnough,],list(
    'TL vs TL'=list(s[s$isTL,'name'],s[s$isTL,'name']),
    'TL-E vs TL-W'=list(s[s$area=='TL-E','name'],s[s$area=='TL-W','name']),
    'TL vs other bonobo'=list(s[s$isTL&s$bonobo,'name'],s[!s$isTL&s$bonobo,'name']),
    # 'Between non-TL bonobo sites'=lapply(unique(s[!s$isTL&s$bonobo,'area']),function(xx)s[s$area==xx,'name']),
    'Malaria vs non-malaria bonobo'=list(s[s$malaria&s$bonobo,'name'],s[s$bonobo&!s$malaria,'name']),
    'Malaria vs non-TL bonobo'=list(s[s$malaria&s$bonobo,'name'],s[!s$isTL&s$bonobo&!s$malaria,'name']),
    'Bonobo vs bonobo'=list(s[s$bonobo,'name'],s[s$bonobo,'name']),
    'Bonobo vs chimp'=list(s[s$bonobo,'name'],s[!s$bonobo,'name']),
    'Non-malaria vs non-malaria chimp'=list(s[!s$bonobo&!s$malaria,'name'],s[!s$bonobo&!s$malaria,'name']),
    'Malaria vs non-malaria chimp'=list(s[!s$bonobo&s$malaria,'name'],s[!s$bonobo&!s$malaria,'name'])
))

pullDists<-function(xx,distMat){
  isIdentical<-length(xx[[1]])==length(xx[[2]]) && all(xx[[1]]==xx[[2]])
  select<-distMat[xx[[1]],xx[[2]]]
  if(isIdentical)dists<-select[upper.tri(select)]
  else dists<-as.vector(select)
  return(dists)
}

nonTL<-unique(samples[samples$bonobo&!samples$isTL,'area'])
allCombo<-unique(t(apply(expand.grid(nonTL,nonTL),1,sort)))
allCombo<-allCombo[allCombo[,1]!=allCombo[,2],]
betweenSites<-apply(allCombo,1,function(xx,distMat)pullDists(list(samples[samples$isEnough&samples$area==xx[1],'name'],samples[samples$isEnough&samples$area==xx[2],'name']),distMat),as.matrix(uniDist))

distList<-lapply(comparisons,pullDists,as.matrix(uniDist))
names(distList)<-ifelse(nchar(names(distList))>20,sub(' vs ',' vs\n',names(distList)),names(distList))
distList<-c('Between non-TL\nbonobo sites'=list(unlist(betweenSites)),distList)

pdf('out/dists.pdf')
  par(mar=c(9,4,.1,.1))
  vpPlot(factor(rep(names(distList),sapply(distList,length)),levels=names(distList)),unlist(distList),las=2,ylab='Unifrac distance')
dev.off()




