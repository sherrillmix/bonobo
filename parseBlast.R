library(dnar)
library(taxonomizr)
library(parallel)

sqlFile<-'~/db/taxo/accessionTaxa.sql'
if(!file.exists(sqlFile)){
  tmp<-tempdir()
  #this is a big download
  getNamesAndNodes('~/db/taxo/')
  getAccession2taxid(tmp)
  read.nodes2('~/db/taxo/nodes.dmp',sqlFile)
  read.names2('~/db/taxo/names.dmp',sqlFile)
  read.accession2taxid(list.files(tmp,'accession2taxid.gz$',full.names=TRUE),sqlFile)
  file.remove(list.files(tmp,'accession2taxid.gz$',full.names=TRUE))
}

blastFiles<-list.files('work/swarmPair','\\.blast\\.gz$',full.names=TRUE)

taxas<-lapply(blastFiles,function(ii,taxaNodes,taxaNames,sqlFile,...){
  message(ii)
  outFile<-sub('.blast.gz$','_taxa.csv',ii)
  outFile2<-sub('.blast.gz$','_allHits.csv',ii)
  if(file.exists(outFile)&file.exists(outFile2)){
    message(' Reading ',outFile)
    taxaAssigns<-read.csv(outFile,row.names=1,stringsAsFactors=FALSE)
    taxonomy<-read.csv(outFile2,stringsAsFactors=FALSE)
  }else{
    message(' Creating ',outFile)
    message('  Reading blast')
    x<-read.blast(ii)
    x<-x[x$tName!='ref',]
    x$accession<-sapply(strsplit(x$tName,'\\|'),'[[',3)
    message('  Accession to taxonomy')
    x$taxa<-accessionToTaxa(x$accession,sqlFile)
    x$sumScore<-ave(x$bit,paste(x$tName,x$qName,sep='_-_'),FUN=sum)
    x$maxScore<-ave(x$sumScore,x$qName,FUN=max)
    x<-x[x$sumScore>x$maxScore*.98&!is.na(x$taxa),]
    gc()
    message('  Getting upstream taxonomy')
    taxonomy<-getTaxonomy2(x$taxa,sqlFile)
    taxonomy<-as.data.frame(taxonomy,stringsAsFactors=FALSE)
    message('  Condensing taxonomy')
    taxaAssigns<-condenseTaxa2(taxonomy,x$qName)
    taxaAssigns<-as.data.frame(taxaAssigns,stringsAsFactors=FALSE)
    rownames(taxaAssigns)<-taxaAssigns$id
    taxaAssigns<-taxaAssigns[,colnames(taxaAssigns)!='id']
    taxonomy$qName<-x$qName
    write.csv(taxonomy,outFile2)
    taxaAssigns$best<-apply(taxaAssigns,1,lastNotNa)
    bestScore<-x[x$sumScore==x$maxScore,c('qName','alignLength','percID','sumScore')]
    bestScore<-bestScore[!duplicated(bestScore$qName),]
    rownames(bestScore)<-bestScore$qName
    taxaAssigns<-cbind(taxaAssigns,'bestScore'=bestScore[rownames(taxaAssigns),c('sumScore')])
    write.csv(taxaAssigns,outFile)
  }
  return(list('taxa'=taxaAssigns,'taxonomy'=taxonomy))
},taxaNodes,taxaNames,sqlFile,extraCode='library(dnar);library(taxonomizr);library(parallel)',mc.cores=4,nSplits=40)


#taxonomy<-lapply(taxas,'[[','taxonomy')
#taxas<-lapply(taxas,'[[','taxa')
#names(taxas)<-names(taxonomy)<-basename(blastFiles)

rm(taxas)
#rm(taxaNames)
#rm(taxaNodes)


