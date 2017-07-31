library(ape)
message('ape package version ',packageVersion('ape'))
library(phyloseq)
message('phyloseq package version ',packageVersion('phyloseq'))

source('../readSamples.R',chdir=TRUE)
source('../functions.R')

#require 15000 reads
nRequiredReads<-15000
message('Requiring samples to have ',nRequiredReads,' reads')

otus<-read.csv('work/qiimeOtuIds.csv.gz',stringsAsFactors=FALSE)
otuTab<-as.data.frame.matrix(table(otus$otu,otus$file))

# Discard singletons
otuTab<-otuTab[apply(otuTab,1,sum)>1,]

# Read qiime taxonomy
taxaRaw<-read.csv('work/qiimeOtus.taxa',stringsAsFactors=FALSE)
taxa<-parseQiimeTaxa(taxaRaw$taxa)
rownames(taxa)<-taxaRaw$name
taxa$best<-apply(taxa,1,function(x)ifelse(is.na(x),x,sprintf('%s_%s',names(x),x))[max(c(1,which(!is.na(x))))])
taxa$bestId<-ave(taxa$best,naReplace(taxa$best,'__NAFILLER__'),FUN=function(x){sprintf('%s #%d',ifelse(is.na(x),'Unknown',x),1:length(x))})

# Filter chloroplast reads
isChloro<-taxa[rownames(otuTab),'c']=='Chloroplast'&!is.na(taxa[rownames(otuTab),'c'])
otuTab<-otuTab[!isChloro,]

# Filter samples without enough reads
isEnough<-apply(otuTab,2,sum)>nRequiredReads

# Convert to proportions
otuProp<-apply(otuTab,2,function(x)x/sum(x))

# Read tree
#make sure tree is bifurcating or breaks UniFrac without error
tree<-ape::multi2di(phyloseq::read_tree('work/qiime/rep_set.tre'))
