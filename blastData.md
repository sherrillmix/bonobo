## Software versions

```r
system('makeblastdb -version',intern=TRUE)
```

```
## [1] "makeblastdb: 2.2.31+"                             
## [2] "Package: blast 2.2.31, build Jan  7 2016 23:17:17"
```

```r
system('blastn -version',intern=TRUE)
```

```
## [1] "blastn: 2.2.31+"                                  
## [2] "Package: blast 2.2.31, build Jan  7 2016 23:17:17"
```

## Prepare blast databases
rbcL and matk data were downloaded manually from the EBI:
* [http://www.ebi.ac.uk/ena/data/search?query=rbcl](http://www.ebi.ac.uk/ena/data/search?query=rbcl)
* [http://www.ebi.ac.uk/ena/data/search?query=matk](http://www.ebi.ac.uk/ena/data/search?query=matk)


```r
if(!file.exists('ebi/rbcl.fa.gz'))stop('Please download rbcl data from http://www.ebi.ac.uk/ena/data/search?query=rbcl to ebi/rbcl.fa.gz')
if(!file.exists('ebi/matk.fa.gz'))stop('Please download matk data from http://www.ebi.ac.uk/ena/data/search?query=matk to ebi/matk.fa.gz')
system('zcat ebi/matk.fa.gz | makeblastdb -in - -title matk -dbtype nucl -out work/matk')
system('zcat ebi/rbcl.fa.gz | makeblastdb -in - -title rbcl -dbtype nucl -out work/rbcl')
```

## Blast swarm OTUs against rbcl/matk databases

```r
contigFa<-list.files('work/swarmPair','.fa.gz$',full.names=TRUE)
contigFa<-contigFa[!grepl('_align',contigFa)]
for(fasta in contigFa){
  outFile<-sub('fa.gz$','blast.gz',fasta)
  if(grepl('rbcL',fasta))gene<-'rbcl'
  else if(grepl('matK',fasta))gene<-'matk'
  else stop('Unknown gene')
  cmd<-sprintf("zcat %s|parallel --block 1m --recstart '>' -L2 -j 14 --pipe blastn -db work/%s -num_threads 5 -culling_limit 30 -outfmt 6 -query -|gzip > %s",fasta,gene,outFile)
  message(cmd)
  system(cmd)
}
```

```
## zcat work/swarmPair/matK.fa.gz|parallel --block 1m --recstart '>' -L2 -j 14 --pipe blastn -db work/matk -num_threads 5 -culling_limit 30 -outfmt 6 -query -|gzip > work/swarmPair/matK.blast.gz
```

```
## zcat work/swarmPair/rbcL.fa.gz|parallel --block 1m --recstart '>' -L2 -j 14 --pipe blastn -db work/rbcl -num_threads 5 -culling_limit 30 -outfmt 6 -query -|gzip > work/swarmPair/rbcL.blast.gz
```
