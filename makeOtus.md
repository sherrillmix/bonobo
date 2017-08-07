## Make rbcL and matK OTUs


```r
#set seed so reproducible
set.seed(12354)
#stop on errors
knitr::opts_chunk$set(error=FALSE,tidy=TRUE)
```

### Load libraries

```r
# installed from https://github.com/sherrillmix/dnar
library(dnar)
packageVersion("dnar")
```

```
## [1] '0.1'
```

```r
library(parallel)
packageVersion("parallel")
```

```
## [1] '3.4.1'
```

```r
source("functions.R")
```
### Check software versions

```r
suppressWarnings(system("mafft --help 2>&1", intern = TRUE))
```

```
##  [1] ""                                                                               
##  [2] "------------------------------------------------------------------------------" 
##  [3] "  MAFFT v7.310 (2017/Mar/17)"                                                   
##  [4] "  http://mafft.cbrc.jp/alignment/software/"                                     
##  [5] "  MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)"                               
##  [6] "------------------------------------------------------------------------------" 
##  [7] "High speed:"                                                                    
##  [8] "  % mafft in > out"                                                             
##  [9] "  % mafft --retree 1 in > out (fast)"                                           
## [10] ""                                                                               
## [11] "High accuracy (for <~200 sequences x <~2,000 aa/nt):"                           
## [12] "  % mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok)"
## [13] "  % mafft --maxiterate 1000 --genafpair  in > out (% einsi in > out)"           
## [14] "  % mafft --maxiterate 1000 --globalpair in > out (% ginsi in > out)"           
## [15] ""                                                                               
## [16] "If unsure which option to use:"                                                 
## [17] "  % mafft --auto in > out"                                                      
## [18] ""                                                                               
## [19] "--op # :         Gap opening penalty, default: 1.53"                            
## [20] "--ep # :         Offset (works like gap extension penalty), default: 0.0"       
## [21] "--maxiterate # : Maximum number of iterative refinement, default: 0"            
## [22] "--clustalout :   Output: clustal format, default: fasta"                        
## [23] "--reorder :      Outorder: aligned, default: input order"                       
## [24] "--quiet :        Do not report progress"                                        
## [25] "--thread # :     Number of threads (if unsure, --thread -1)"                    
## attr(,"status")
## [1] 1
```

```r
system("test -x `which dpkg` && dpkg -s fasttree", intern = TRUE)
```

```
##  [1] "Package: fasttree"                                                                            
##  [2] "Status: install ok installed"                                                                 
##  [3] "Priority: optional"                                                                           
##  [4] "Section: science"                                                                             
##  [5] "Installed-Size: 432"                                                                          
##  [6] "Maintainer: Ubuntu Developers <ubuntu-devel-discuss@lists.ubuntu.com>"                        
##  [7] "Architecture: amd64"                                                                          
##  [8] "Version: 2.1.8-2"                                                                             
##  [9] "Depends: libc6 (>= 2.14), libgomp1 (>= 4.9)"                                                  
## [10] "Description: phylogenetic trees from alignments of nucleotide or protein sequences"           
## [11] " FastTree infers approximately-maximum-likelihood phylogenetic trees from"                    
## [12] " alignments of nucleotide or protein sequences. It handles alignments"                        
## [13] " with up to a million of sequences in a reasonable amount of time and"                        
## [14] " memory. For large alignments, FastTree is 100-1,000 times faster than"                       
## [15] " PhyML 3.0 or RAxML 7."                                                                       
## [16] " ."                                                                                           
## [17] " FastTree is more accurate than PhyML 3 with default settings, and much"                      
## [18] " more accurate than the distance-matrix methods that are traditionally"                       
## [19] " used for large alignments. FastTree uses the Jukes-Cantor or generalized"                    
## [20] " time-reversible (GTR) models of nucleotide evolution and the JTT"                            
## [21] " (Jones-Taylor-Thornton 1992) model of amino acid evolution. To account"                      
## [22] " for the varying rates of evolution across sites, FastTree uses a single"                     
## [23] " rate for each site (the \"CAT\" approximation). To quickly estimate the"                     
## [24] " reliability of each split in the tree, FastTree computes local support"                      
## [25] " values with the Shimodaira-Hasegawa test (these are the same as PhyML 3's"                   
## [26] " \"SH-like local supports\")."                                                                
## [27] " ."                                                                                           
## [28] " This package contains a single threaded version (fasttree) and a"                            
## [29] " parallel version which uses OpenMP (fasttreMP)."                                             
## [30] "Original-Maintainer: Debian Med Packaging Team <debian-med-packaging@lists.alioth.debian.org>"
## [31] "Homepage: http://www.microbesonline.org/fasttree/"
```

```r
suppressWarnings(system("swarm --version 2>&1", intern = TRUE))
```

```
## [1] "Swarm 2.1.1 [Mar 31 2015 15:55:38]"                                     
## [2] "Copyright (C) 2012-2015 Torbjorn Rognes and Frederic Mahe"              
## [3] "https://github.com/torognes/swarm"                                      
## [4] ""                                                                       
## [5] "Please cite: Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2014)"
## [6] "Swarm: robust and fast clustering method for amplicon-based studies."   
## [7] "PeerJ 2:e593 https://dx.doi.org/10.7717/peerj.593"                      
## [8] ""                                                                       
## attr(,"status")
## [1] 1
```

### Primers used

```r
primerSeqs <- list(rbcl = c(rbcL634F = "ATGCGTTGGAGAGACCGTTTC", rbcLbR = "TCGGTYAGAGCRGGCATRTGCCA"), 
    matk = c(matK472F = "CCCRTYCATCTGGAAATCTTGGTTC", matK1248R = "GCTRTRATAATGAGAAAGATTTCTGC"))
```

### Find fastqs

```r
fastqs <- list.files("data/", "_R[12]_.*\\.fastq\\.gz$", recursive = TRUE, full.names = TRUE)
fastqs <- fastqs[!grepl("Undetermined", fastqs)]
primers <- sub(".*(matK|rbcL).*_R([0-9]+)_.*", "\\1\\2", basename(fastqs))
primerBases <- sub("[0-9]$", "", primers)
```


### Run swarm, mafft and fasttree

```r
for (primerBase in unique(primerBases)) {
    message("Working on ", primerBase)
    if (!dir.exists("work/swarmPair")) 
        dir.create("work/swarmPair", recursive = TRUE)
    outMat <- sprintf("work/swarmPair/%s.Rdat", primerBase)
    outFa <- sprintf("work/swarmPair/%s.fa.gz", primerBase)
    outAlign <- sprintf("work/swarmPair/%s_align.fa.gz", primerBase)
    outTree <- sprintf("work/swarmPair/%s_align.tre", primerBase)
    message("Trimming primers off the start of reads")
    trimReads <- lapply(sprintf("%s%d", primerBase, 1:2), function(ii) {
        thisPrimer <- primerSeqs[[sub("[12]$", "", tolower(ii))]][as.numeric(substring(ii, 
            nchar(ii)))]
        thisFiles <- fastqs[primers == ii]
        reads <- mclapply(thisFiles, dnar::read.fastq, mc.cores = 20, mc.preschedule = FALSE)
        if (mean(unlist(lapply(reads, function(xx) substring(xx$seq, 1, nchar(thisPrimer)))) %in% 
            expandAmbiguous(thisPrimer)[[1]]) < 0.75) 
            stop(simpleError("Expected primer does not match read start"))
        trimReads <- lapply(reads, function(xx) {
            xx$primerMatch <- substring(xx$seq, 1, nchar(thisPrimer)) %in% expandAmbiguous(thisPrimer)[[1]]
            xx$seq <- substring(xx$seq, nchar(thisPrimer) + 1)
            return(xx)
        })
        names(trimReads) <- thisFiles
        return(trimReads)
    })
    if (any(names(trimReads[[1]]) != sub("_R2_", "_R1_", names(trimReads[[2]])))) 
        stop("Read 1 vs reads 2 file name mismatch")
    message("Discarding reads with >1 expected error in left or right read and concatenating left-right")
    trimReads <- mcmapply(function(left, right, ...) {
        cat(".")
        if (any(sub(" .*$", "", left$name) != sub(" .*$", "", right$name))) 
            stop("Read 1 vs read 2 name mismatch")
        # last base always low qual so ignore
        q1 <- sapply(qualToInts(substring(left$qual, 1, nchar(left$qual) - 1)), 
            function(xx) sum(10^(-xx/10)))
        q2 <- sapply(qualToInts(substring(right$qual, 1, nchar(right$qual) - 
            1)), function(xx) sum(10^(-xx/10)))
        # less than 1 expected error in both reads, no Ns and primer match
        selector <- q1 < 1 & q2 < 1 & !grepl("[^ACTG]", left$seq) & !grepl("[^ACTG]", 
            right$seq) & left$primerMatch & right$primerMatch
        seqs <- paste(left[selector, "seq"], revComp(right[selector, "seq"]), 
            sep = "")
        return(seqs)
    }, trimReads[[1]], trimReads[[2]], mc.cores = 10, SIMPLIFY = FALSE)
    samples <- rep(basename(names(trimReads)), sapply(trimReads, length))
    message("Running swarm")
    otus <- runSwarm(unlist(trimReads), swarmArgs = "-f -t 40")
    swarmOtus <- as.data.frame.matrix(table(samples, otus[["otus"]]))
    write.fa(otus[["seqs"]]$name, otus[["seqs"]]$seq, outFa)
    save(swarmOtus, file = outMat)
    load(outMat)
    tmpFile <- tempfile()
    message("Discarding singleton OTUs")
    swarmOtus <- swarmOtus[, apply(swarmOtus, 2, sum) > 1]
    seqs <- otus[["seqs"]]
    rownames(seqs) <- seqs$name
    write.fa(colnames(swarmOtus), seqs[colnames(swarmOtus), "seq"], tmpFile)
    message("Aligning with mafft")
    cmd <- sprintf("mafft --thread 50 %s|gzip>%s", tmpFile, outAlign)
    message(cmd)
    system(cmd)
    message("Creating tree with fasttree")
    cmd <- sprintf("zcat %s|fasttree -gtr -nt>%s", outAlign, outTree)
    message(cmd)
    system(cmd)
    message("Output files in:")
    message("    ", outMat)
    message("    ", outFa)
    message("    ", outAlign)
    message("    ", outTree)
}
```

```
## Working on matK
```

```
## Trimming primers off the start of reads
```

```
## Discarding reads with >1 expected error in left or right read and concatenating left-right
```

```
## Running swarm
```

```
## swarm -f -t 40 /tmp/RtmppEgiBM/file29f750e0cfaf -o /tmp/RtmppEgiBM/file29f797b4330 -w /tmp/RtmppEgiBM/file29f73e1183fe
```

```
## Discarding singleton OTUs
```

```
## Aligning with mafft
```

```
## mafft --thread 50 /tmp/RtmppEgiBM/file29f739e6d6d|gzip>work/swarmPair/matK_align.fa.gz
```

```
## Creating tree with fasttree
```

```
## zcat work/swarmPair/matK_align.fa.gz|fasttree -gtr -nt>work/swarmPair/matK_align.tre
```

```
## Output files in:
```

```
##     work/swarmPair/matK.Rdat
```

```
##     work/swarmPair/matK.fa.gz
```

```
##     work/swarmPair/matK_align.fa.gz
```

```
##     work/swarmPair/matK_align.tre
```

```
## Working on rbcL
```

```
## Trimming primers off the start of reads
```

```
## Discarding reads with >1 expected error in left or right read and concatenating left-right
```

```
## Running swarm
```

```
## swarm -f -t 40 /tmp/RtmppEgiBM/file29f7263fd323 -o /tmp/RtmppEgiBM/file29f722ba70c8 -w /tmp/RtmppEgiBM/file29f76248ff3b
```

```
## Discarding singleton OTUs
```

```
## Aligning with mafft
```

```
## mafft --thread 50 /tmp/RtmppEgiBM/file29f72c533751|gzip>work/swarmPair/rbcL_align.fa.gz
```

```
## Creating tree with fasttree
```

```
## zcat work/swarmPair/rbcL_align.fa.gz|fasttree -gtr -nt>work/swarmPair/rbcL_align.tre
```

```
## Output files in:
```

```
##     work/swarmPair/rbcL.Rdat
```

```
##     work/swarmPair/rbcL.fa.gz
```

```
##     work/swarmPair/rbcL_align.fa.gz
```

```
##     work/swarmPair/rbcL_align.tre
```

