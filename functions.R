
runQiime<-function(seqs,storeDir=NULL){
  if(any(is.na(seqs)))stop(simpleError('NAs in seqs'))
  seqIds<-1:length(seqs)
  readDir<-tempfile()
  dir.create(readDir)
  readFile<-file.path(readDir,'XXX.fa')
  outDir<-tempfile()
  seqNames<-sprintf('XXX_%d',seqIds)
  write.fa(seqNames,seqs,readFile)
  #miniconda doesn't like sh so need to use bash
  cmd<-sprintf('echo "source activate qiime1; pick_de_novo_otus.py --input %s --output %s --parallel --jobs_to_start 32 --force"|bash',readFile,outDir)
  message(cmd)
  exit<-system(cmd)
  if(exit!=0)stop(simpleError('Problem running qiime'))
  #get otu assignments
  assigns<-strsplit(readLines(file.path(outDir,'uclust_picked_otus/XXX_otus.txt')),'\t')
  names(assigns)<-sapply(assigns,'[[',1)
  assigns<-lapply(assigns,'[',-1)
  otus<-rep(names(assigns),sapply(assigns,length))
  names(otus)<-unlist(assigns)
  out<-otus[seqNames]
  #get taxa assignments
  taxa<-strsplit(readLines(file.path(outDir,'uclust_assigned_taxonomy/XXX_rep_set_tax_assignments.txt')),'\t')
  names(taxa)<-sapply(taxa,'[[',1)
  taxa<-sapply(taxa,'[[',2)
  #get sequences
  seqs<-read.fa(file.path(outDir,'pynast_aligned_seqs/XXX_rep_set_aligned_pfiltered.fasta'))
  seqs<-structure(seqs$seq,.Names=seqs$name)
  if(!is.null(storeDir)){
    #avoid problems with cross filesystem copies
    tmpDir<-tempfile(tmpdir='.')
    dir.create(tmpDir)
    file.copy(outDir,tmpDir)
    file.rename(file.path(tmpDir,basename(outDir)),storeDir)
  }
  return(list('otus'=out,'seqs'=seqs,'taxa'=taxa))
}
parseQiimeTaxa<-function(taxas,desiredTaxa=c('k','p','c','o','f','g','s'),concatLastTwo=TRUE){
  taxa<-strsplit(taxaRaw$taxa,'[;] ?')
  taxa<-do.call(rbind,lapply(taxa,function(xx){
    y<-xx[grep('^[a-z]__',xx)]
    names(y)<-sub('^([a-z])__.*','\\1',y)
    y<-sub('^[a-z]__','',y)
    y[y=='']<-NA
    out<-y[desiredTaxa]
    if(!is.na(out[length(desiredTaxa)]))out[length(desiredTaxa)]<-paste(out[length(desiredTaxa)-1:0],collapse=' ')
    return(out)
  }))
  return(as.data.frame(taxa,stringsAsFactors=FALSE))
}


runSwarm<-function(seqs,swarmBin='swarm',swarmArgs='-f'){
  if(any(is.na(seqs)))stop(simpleError('NAs in seqs'))
  seqIds<-as.numeric(as.factor(seqs))
  seqCounts<-ave(seqIds,seqIds,FUN=length)
  seqNames<-sprintf('%08d_%d',seqIds,seqCounts)
  readFile<-tempfile()
  outFile<-tempfile()
  seqFile<-tempfile()
  uniqSelector<-!duplicated(seqs)
  write.fa(seqNames[uniqSelector],seqs[uniqSelector],readFile)
  cmd<-sprintf('%s %s %s -o %s -w %s',swarmBin,swarmArgs,readFile,outFile,seqFile)
  system(cmd)
  swarm<-readLines(outFile)
  otuSplit<-strsplit(swarm,' ')
  otus<-rep(1:length(otuSplit),sapply(otuSplit,length))
  names(otus)<-unlist(otuSplit)  
  out<-otus[seqNames]
  seedSeqs<-read.fa(seqFile)
  if(nrow(seedSeqs)!=length(otuSplit))stop('Seqs and OTU assigns not same length')
  if(any(!mapply(function(xx,yy)sub('_[0-9]+$','',xx) %in% sub('_[0-9]+$','',yy),seedSeqs$name,otuSplit)))stop('Mismatch between sequence and OTU names')
  if(any(!mapply(function(xx,yy)as.numeric(sub('^[^_]*_','',xx)) == sum(as.numeric(sub('^[^_]*_','',yy))),seedSeqs$name,otuSplit)))stop('Mismatch between sequence and OTU counts')
  seedSeqs$bak<-seedSeqs$name
  seedSeqs$name<-1:nrow(seedSeqs)
  return(list('otus'=out,'seqs'=seedSeqs))
}


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

medianCI<-function(xx,na.rm=TRUE){
  if(na.rm)xx<-xx[!is.na(xx)]
  return(unname(sort(xx)[qbinom(c(.025,.975), length(xx), 0.5)]))
}

