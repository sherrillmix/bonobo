## Download rbcl and matk databases from EBI


```r
#set seed so reproducible
set.seed(12352)
#stop on errors
knitr::opts_chunk$set(error=FALSE,tidy=TRUE)
```

### Load libraries

```r
library(parallel)
packageVersion("parallel")
```

```
## [1] '3.4.1'
```

### Download data from ENA
Download the equivalent of searches for "rbcl" and "matk" from EMBL-EBI European Nucleotide Archive:
* [http://www.ebi.ac.uk/ena/data/search?query=rbcl](http://www.ebi.ac.uk/ena/data/search?query=rbcl)
* [http://www.ebi.ac.uk/ena/data/search?query=matk](http://www.ebi.ac.uk/ena/data/search?query=matk)

```r
if (!dir.exists("ebi")) dir.create("ebi")
# note limit=400000
downloadUrls <- c(`ebi/rbcl.fa.gz` = "http://www.ebi.ac.uk/ena/data/search?query=rbcl&result=sequence_release&display=fasta&download=gzip&limit=400000", 
    `ebi/matk.fa.gz` = "http://www.ebi.ac.uk/ena/data/search?query=matk&result=sequence_release&display=fasta&download=gzip&limit=400000")
mclapply(names(downloadUrls), function(xx) download.file(downloadUrls[xx], xx, 
    mode = "wb"), mc.cores = 2)
```

```
## [[1]]
## [1] 0
## 
## [[2]]
## [1] 0
```
