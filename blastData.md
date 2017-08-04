## Blast OTUs against rbcL and matK databases


```r
#set seed so reproducible
set.seed(12350)
#stop on errors
knitr::opts_chunk$set(error=FALSE,tidy=TRUE)
```

### Software versions

```r
system("makeblastdb -version", intern = TRUE)
```

```
## [1] "makeblastdb: 2.2.31+"                             
## [2] "Package: blast 2.2.31, build Jan  7 2016 23:17:17"
```

```r
system("blastn -version", intern = TRUE)
```

```
## [1] "blastn: 2.2.31+"                                  
## [2] "Package: blast 2.2.31, build Jan  7 2016 23:17:17"
```

```r
system("parallel --version", intern = TRUE)
```

```
##  [1] "GNU parallel 20141022"                                                                             
##  [2] "Copyright (C) 2007,2008,2009,2010,2011,2012,2013,2014 Ole Tange and Free Software Foundation, Inc."
##  [3] "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>"                     
##  [4] "This is free software: you are free to change and redistribute it."                                
##  [5] "GNU parallel comes with no warranty."                                                              
##  [6] ""                                                                                                  
##  [7] "Web site: http://www.gnu.org/software/parallel"                                                    
##  [8] ""                                                                                                  
##  [9] "When using programs that use GNU Parallel to process data for publication please cite:"            
## [10] ""                                                                                                  
## [11] "O. Tange (2011): GNU Parallel - The Command-Line Power Tool, "                                     
## [12] ";login: The USENIX Magazine, February 2011:42-47."                                                 
## [13] ""                                                                                                  
## [14] "Or you can get GNU Parallel without this requirement by paying 10000 EUR."
```

### Prepare rbcL and matK blast databases from downloaded sequences


```r
if (!file.exists("ebi/rbcl.fa.gz")) stop("Please download rbcl data from http://www.ebi.ac.uk/ena/data/search?query=rbcl to ebi/rbcl.fa.gz")
if (!file.exists("ebi/matk.fa.gz")) stop("Please download matk data from http://www.ebi.ac.uk/ena/data/search?query=matk to ebi/matk.fa.gz")
system("zcat ebi/matk.fa.gz | makeblastdb -in - -title matk -dbtype nucl -out work/matk")
system("zcat ebi/rbcl.fa.gz | makeblastdb -in - -title rbcl -dbtype nucl -out work/rbcl")
```

### Blast swarm OTUs against rbcl/matk databases

```r
contigFa <- list.files("work/swarmPair", ".fa.gz$", full.names = TRUE)
contigFa <- contigFa[!grepl("_align", contigFa)]
for (fasta in contigFa) {
    outFile <- sub("fa.gz$", "blast.gz", fasta)
    if (grepl("rbcL", fasta)) 
        gene <- "rbcl" else if (grepl("matK", fasta)) 
        gene <- "matk" else stop("Unknown gene")
    cmd <- sprintf("zcat %s|parallel --block 1m --recstart '>' -L2 -j 14 --pipe blastn -db work/%s -num_threads 5 -culling_limit 30 -outfmt 6 -query -|gzip > %s", 
        fasta, gene, outFile)
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
