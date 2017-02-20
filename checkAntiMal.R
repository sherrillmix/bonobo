library(taxonomizr)
antiMal<-read.csv('AntiMalPlant_List_021817.csv',stringsAsFactors=FALSE)
getNamesAndNodes('~/db/taxo/')
taxaNodes<-read.nodes('~/db/taxo/nodes.dmp')
taxaNamesAmbig<-read.names('~/db/taxo/names.dmp',FALSE)
taxaNames<-read.names('~/db/taxo/names.dmp')
antiMal$clean<-sub(' sp$','',sub('([A-Za-z]+ [a-z]+).*','\\1',antiMal$Plant))
antiMal$clean[antiMal$clean=="AAloe parvibracteata"]<-"Aloe parvibracteata"
antiMal$clean[antiMal$clean=="Curcuma aromatic"]<-"Curcuma aromatica"
antiMal$id<-getId(antiMal$clean,taxaNamesAmbig)
antiMal$species<-getTaxonomy(antiMal$id,taxaNodes,taxaNames,'species')

if(!exists('chimpFilter')){
  source('parseBlast.R')
  chimpFilter<-lapply(taxas,function(x)x[x[,'class']!='Mammalia'|is.na(x[,'class']),])
  #dont need currently
  rm(taxas)
  rm(taxonomy)
}

