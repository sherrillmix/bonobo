naReplace<-function(x,replace){x[is.na(x)]<-replace;return(x)}

plotHeat<-function(selectProp,breaks,cols,xaxt='',yaxt=''){
  image(1:ncol(selectProp),1:nrow(selectProp),t(selectProp),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='')
  box()
  insetScale(round(breaks,6),cols,c(.97,.01,.98,.25),label='Proportion of OTU maximum')
  if(xaxt!='n')slantAxis(1,1:ncol(selectProp),colnames(selectProp))
  if(yaxt!='n')axis(4,1:nrow(selectProp),rownames(selectProp),las=1,tcl=0,mgp=c(0,.2,0),cex.axis=.7)
  abline(h=1:nrow(selectProp)-.5,v=1:ncol(selectProp)+.5,col='#00000011')
}

addMetaData<-function(metadata,cex=1,...){
  widths<-apply(metadata,2,function(x)max(strwidth(x,cex=cex)))
  headerWidths<-strwidth(colnames(metadata),cex=par('cex.axis'))
  widths<-apply(rbind(headerWidths,widths),2,max)
  spacer<-diff(par('usr')[1:2])*.005
  labX<-cumsum(c(0,widths[-length(widths)])+spacer)
  for(ii in 1:ncol(metadata)){
    text(par('usr')[2]+labX[ii],1:nrow(metadata),metadata[,ii],adj=c(0,.5),xpd=NA,cex=cex,...)
  }
  text(par('usr')[2]+labX,convertLineToUser(.3,3),colnames(metadata),adj=c(0,0),xpd=NA,cex=par('cex.axis'))
}


fishers<-function(ps,correct=1){
  chi2<- -2*sum(log(ps))
  out<-1-pchisq(chi2/correct,2*length(ps)/correct)
  return(out)
}

setupHeat<-function(pos,neg,otuProp,taxa,minProp=.001,pCut=.1,treeScale=1/3){
  inSample<-otuProp[apply(otuProp[,c(pos,neg)],1,max)>minProp,]
  pVals<-p.adjust(apply(inSample,1,function(xx)suppressWarnings(wilcox.test(xx[pos],xx[neg]))$p.value),'fdr')
  effect<-apply(inSample,1,function(xx)median(xx[pos])-median(xx[neg]))
  pVals[is.na(pVals)]<-1
  selectProp<-apply(inSample[pVals<pCut,,drop=FALSE],1,function(x)x/max(x))
  hTree<-hclust(dist(t(selectProp[c(pos,neg),]^treeScale)))
  #selectProp<-selectProp[,hTree$labels[hTree$order]]
  effectP<-((effect>0)*2-1)*(1-pVals)
  selectProp<-selectProp[,order(effect[pVals<pCut]>0,order(orderIn(colnames(selectProp),hTree$labels[hTree$order]))),drop=FALSE]
  #selectProp<-selectProp[,order(effectP[pVals<pCut]),drop=FALSE]
  effectSplit<-min(which(effectP[colnames(selectProp)]>0))
  colnames(selectProp)<-sprintf('%s (q=%0.3f)',taxa[colnames(selectProp),'bestId'],pVals[colnames(selectProp)])
  return(list(selectProp,effectSplit))
}

