## Load libraries

```r
library(dnar)
packageVersion('dnar')
```

```
## [1] '0.1'
```

```r
library(taxonomizr)
packageVersion('taxonomizr')
```

```
## [1] '0.3.0'
```

```r
library(parallel)
packageVersion('parallel')
```

```
## [1] '3.4.1'
```

## Download NCBI accession numbers and taxonomy

```r
sqlFile<-'work/accessionTaxa.sql'
getNamesAndNodes('work')
```

```
## work/names.dmp, work/nodes.dmp already exist. Delete to redownload
```

```
## [1] "work/names.dmp" "work/nodes.dmp"
```

```r
#this is a big download
getAccession2taxid('work')
```

```
## This can be a big (several gigabytes) download. Please be patient and use a fast connection.
```

```
## [1] "work/nucl_gb.accession2taxid.gz"  "work/nucl_est.accession2taxid.gz"
## [3] "work/nucl_gss.accession2taxid.gz" "work/nucl_wgs.accession2taxid.gz"
```

## Prepare database of accession numbers and taxonomy

```r
read.nodes2('work/nodes.dmp',sqlFile)
read.names2('work/names.dmp',sqlFile)
read.accession2taxid(list.files('work','accession2taxid.gz$',full.names=TRUE),sqlFile,overwrite=TRUE)
```

```
## Reading work/nucl_est.accession2taxid.gz.
```

```
## Reading work/nucl_gb.accession2taxid.gz.
```

```
## Reading work/nucl_gss.accession2taxid.gz.
```

```
## Reading work/nucl_wgs.accession2taxid.gz.
```

```
## Reading in values. This may take a while.
```

```
## Adding index. This may also take a while.
```

## Assign taxonomy to blast hits

```r
blastFiles<-list.files('work/swarmPair','\\.blast\\.gz$',full.names=TRUE)

mclapply(blastFiles,function(ii,sqlFile,...){
  message('Parsing blast hits in ',ii)
  outFile<-sub('.blast.gz$','_taxa.csv',ii)
  outFile2<-sub('.blast.gz$','_allHits.csv',ii)
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
  return(c('taxa'=outFile,'allHits'=outFile2))
},sqlFile,mc.cores=2)
```

```
## [[1]]
##                              taxa                           allHits 
##    "work/swarmPair/matK_taxa.csv" "work/swarmPair/matK_allHits.csv" 
## 
## [[2]]
##                              taxa                           allHits 
##    "work/swarmPair/rbcL_taxa.csv" "work/swarmPair/rbcL_allHits.csv"
```
