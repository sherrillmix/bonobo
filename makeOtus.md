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
        reads <- mclapply(thisFiles, function(xx) {
            library(dnar)
            cat(".")
            read.fastq(xx)
        }, mc.cores = 20, mc.preschedule = FALSE)
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
    readCounts <- as.data.frame(sapply(trimReads, function(xx) sapply(xx, nrow)))
    colnames(readCounts) <- c("raw1", "raw2")
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
    readCounts$filter <- sapply(trimReads, length)
    message("Read counts")
    print(readCounts)
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
## Read counts
```

```
##                                 raw1  raw2 filter
## data//BI0054_matK_R1_.fastq.gz 37716 37716  22419
## data//BI0055_matK_R1_.fastq.gz 29870 29870  19976
## data//BI0093_matK_R1_.fastq.gz 23870 23870   8397
## data//BI0097_matK_R1_.fastq.gz 29597 29597  19234
## data//BI0246_matK_R1_.fastq.gz 23907 23907  15538
## data//BI0248_matK_R1_.fastq.gz 52496 52496  24563
## data//BI0257_matK_R1_.fastq.gz 45035 45035  27103
## data//BI0260_matK_R1_.fastq.gz 38342 38342  13744
## data//BI2414_matK_R1_.fastq.gz 42410 42410  24285
## data//BI2415_matK_R1_.fastq.gz 21547 21547  12159
## data//IK3158_matK_R1_.fastq.gz 22764 22764  14387
## data//IK3276_matK_R1_.fastq.gz 51495 51495  35038
## data//IK3358_matK_R1_.fastq.gz 22084 22084  14917
## data//IK3469_matK_R1_.fastq.gz 35564 35564  23276
## data//IK3513_matK_R1_.fastq.gz 43770 43770  25887
## data//IK3650_matK_R1_.fastq.gz 35772 35772  25011
## data//IK3701_matK_R1_.fastq.gz 53473 53473  35065
## data//IK3777_matK_R1_.fastq.gz 24679 24679  13341
## data//IK4184_matK_R1_.fastq.gz 31518 31518  14978
## data//IK4214_matK_R1_.fastq.gz 13382 13382    620
## data//KR02_matK_R1_.fastq.gz   37289 37289  17316
## data//KR05_matK_R1_.fastq.gz   17096 17096   6343
## data//KR07_matK_R1_.fastq.gz   49192 49192  30653
## data//KR10_matK_R1_.fastq.gz   39057 39057  24279
## data//KR12_matK_R1_.fastq.gz   16755 16755  10786
## data//KR21_matK_R1_.fastq.gz   31486 31486  18668
## data//KR33_matK_R1_.fastq.gz   31554 31554  19568
## data//KR35_matK_R1_.fastq.gz   46935 46935  30066
## data//KR52_matK_R1_.fastq.gz   50819 50819  33481
## data//KR57_matK_R1_.fastq.gz   22953 22953   7905
## data//KR67_matK_R1_.fastq.gz   18119 18119  10759
## data//LG4300_matK_R1_.fastq.gz 31703 31703  16584
## data//LG4314_matK_R1_.fastq.gz 24115 24115  15106
## data//LG4322_matK_R1_.fastq.gz 26690 26690  17135
## data//LG4327_matK_R1_.fastq.gz 36202 36202  23733
## data//LK645_matK_R1_.fastq.gz  28947 28947  19557
## data//LK647_matK_R1_.fastq.gz  39875 39875  27466
## data//LK653_matK_R1_.fastq.gz  40197 40197  25629
## data//LK661_matK_R1_.fastq.gz  33390 33390  14521
## data//LK665_matK_R1_.fastq.gz  54621 54621  35101
## data//LK668_matK_R1_.fastq.gz  30992 30992  20217
## data//LK670_matK_R1_.fastq.gz  43269 43269  27602
## data//LK682_matK_R1_.fastq.gz  32752 32752  20517
## data//LK685_matK_R1_.fastq.gz  26060 26060   7601
## data//LK686_matK_R1_.fastq.gz  24219 24219   8620
## data//PA0367_matK_R1_.fastq.gz 26739 26739  10411
## data//PA0368_matK_R1_.fastq.gz 23805 23805   8210
## data//PA0370_matK_R1_.fastq.gz 19569 19569   2627
## data//PA0456_matK_R1_.fastq.gz 19652 19652   6823
## data//PA1038_matK_R1_.fastq.gz 25960 25960  13353
## data//PA1039_matK_R1_.fastq.gz 26785 26785  13824
## data//PA1044_matK_R1_.fastq.gz 54164 54164  32326
## data//PA1049_matK_R1_.fastq.gz 46651 46651  29706
## data//PA1059_matK_R1_.fastq.gz 37381 37381  19092
## data//PA1065_matK_R1_.fastq.gz 30597 30597  18900
## data//TL3793_matK_R1_.fastq.gz 21275 21275  10999
## data//TL3797_matK_R1_.fastq.gz 11636 11636   6806
## data//TL3814_matK_R1_.fastq.gz 50184 50184  24968
## data//TL3816_matK_R1_.fastq.gz 46883 46883  23842
## data//TL3820_matK_R1_.fastq.gz 66064 66064  17098
## data//TL3821_matK_R1_.fastq.gz 24907 24907  11624
## data//TL3824_matK_R1_.fastq.gz 29323 29323  15046
## data//TL3826_matK_R1_.fastq.gz 61039 61039  31482
## data//TL3838_matK_R1_.fastq.gz 35465 35465  18700
## data//TL3842_matK_R1_.fastq.gz 38016 38016  18161
## data//TL3856_matK_R1_.fastq.gz 39587 39587  20218
## data//TL3862_matK_R1_.fastq.gz 32951 32951  19278
## data//TL3882_matK_R1_.fastq.gz 15760 15760   9306
## data//TL3889_matK_R1_.fastq.gz 37674 37674  20620
## data//TL3905_matK_R1_.fastq.gz 20322 20322  11871
## data//TL3910_matK_R1_.fastq.gz 24363 24363  13210
## data//TL3911_matK_R1_.fastq.gz 19392 19392   7831
## data//TL3915_matK_R1_.fastq.gz 34545 34545  15882
## data//TL3916_matK_R1_.fastq.gz    35    35      4
## data//TL3918_matK_R1_.fastq.gz 19370 19370   6740
## data//TL3925_matK_R1_.fastq.gz 33913 33913  20472
## data//TL3926_matK_R1_.fastq.gz 27740 27740  15625
## data//TL3927_matK_R1_.fastq.gz 18085 18085   9124
## data//TL3929_matK_R1_.fastq.gz 17720 17720  10116
## data//TL3932_matK_R1_.fastq.gz 19395 19395  11229
## data//TL3936_matK_R1_.fastq.gz 13435 13435   8603
## data//TL3939_matK_R1_.fastq.gz 16014 16014   9191
## data//TL3940_matK_R1_.fastq.gz 12612 12612   5663
## data//TL3942_matK_R1_.fastq.gz 65864 65864  40859
## data//TL3943_matK_R1_.fastq.gz 62439 62439  34926
## data//TL3944_matK_R1_.fastq.gz 48857 48857  29247
## data//TL3945_matK_R1_.fastq.gz 51533 51533  31465
## data//TL3946_matK_R1_.fastq.gz 31692 31692  20208
## data//TL3948_matK_R1_.fastq.gz 58040 58040  32802
## data//UB0439_matK_R1_.fastq.gz 33727 33727  20532
## data//UB0445_matK_R1_.fastq.gz 21950 21950   8277
## data//UB0599_matK_R1_.fastq.gz 30433 30433  19637
## data//UB1430_matK_R1_.fastq.gz 36057 36057  21126
## data//UB1435_matK_R1_.fastq.gz 44581 44581  29690
## data//UB1446_matK_R1_.fastq.gz 22388 22388  11675
## data//UB1452_matK_R1_.fastq.gz 27983 27983  17216
## data//UB1454_matK_R1_.fastq.gz 17005 17005   5163
## data//UB2037_matK_R1_.fastq.gz 28886 28886  18178
```

```
## Running swarm
```

```
## swarm -f -t 40 /tmp/RtmpX0M3mf/fileb79b7c57cd81 -o /tmp/RtmpX0M3mf/fileb79b185dcc32 -w /tmp/RtmpX0M3mf/fileb79b12db05da
```

```
## Discarding singleton OTUs
```

```
## Aligning with mafft
```

```
## mafft --thread 50 /tmp/RtmpX0M3mf/fileb79b7edc61ae|gzip>work/swarmPair/matK_align.fa.gz
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
## Read counts
```

```
##                                  raw1   raw2 filter
## data//BI0054_rbcL_R1_.fastq.gz 112777 112777  69585
## data//BI0055_rbcL_R1_.fastq.gz  40341  40341  17662
## data//BI0093_rbcL_R1_.fastq.gz  60295  60295  38983
## data//BI0097_rbcL_R1_.fastq.gz  38815  38815  25219
## data//BI0246_rbcL_R1_.fastq.gz  97069  97069  59088
## data//BI0248_rbcL_R1_.fastq.gz  63847  63847  37019
## data//BI0257_rbcL_R1_.fastq.gz  68795  68795  38721
## data//BI0260_rbcL_R1_.fastq.gz  75418  75418  47727
## data//BI2414_rbcL_R1_.fastq.gz  72678  72678  37522
## data//BI2415_rbcL_R1_.fastq.gz  72083  72083  48003
## data//IK3158_rbcL_R1_.fastq.gz  63877  63877  38759
## data//IK3276_rbcL_R1_.fastq.gz  75130  75130  45827
## data//IK3358_rbcL_R1_.fastq.gz  83741  83741  53066
## data//IK3469_rbcL_R1_.fastq.gz  57184  57184  35361
## data//IK3513_rbcL_R1_.fastq.gz  72466  72466  45628
## data//IK3650_rbcL_R1_.fastq.gz  39570  39570  24195
## data//IK3701_rbcL_R1_.fastq.gz  71382  71382  46583
## data//IK3777_rbcL_R1_.fastq.gz  65948  65948  41943
## data//IK4184_rbcL_R1_.fastq.gz  49004  49004  21269
## data//IK4214_rbcL_R1_.fastq.gz  42566  42566  18519
## data//KR02_rbcL_R1_.fastq.gz    46441  46441  28595
## data//KR05_rbcL_R1_.fastq.gz    14003  14003   5182
## data//KR07_rbcL_R1_.fastq.gz    82888  82888  41950
## data//KR10_rbcL_R1_.fastq.gz    71966  71966  31446
## data//KR12_rbcL_R1_.fastq.gz    65322  65322  37758
## data//KR21_rbcL_R1_.fastq.gz    43489  43489  26051
## data//KR33_rbcL_R1_.fastq.gz    50031  50031  29133
## data//KR35_rbcL_R1_.fastq.gz    64584  64584  39959
## data//KR52_rbcL_R1_.fastq.gz    25519  25519  15825
## data//KR57_rbcL_R1_.fastq.gz    45206  45206  19134
## data//KR67_rbcL_R1_.fastq.gz    74544  74544  44921
## data//LG4300_rbcL_R1_.fastq.gz  64265  64265  21340
## data//LG4314_rbcL_R1_.fastq.gz  73212  73212  43803
## data//LG4322_rbcL_R1_.fastq.gz  20771  20771    490
## data//LG4327_rbcL_R1_.fastq.gz  66867  66867  44228
## data//LK645_rbcL_R1_.fastq.gz   87554  87554  41568
## data//LK647_rbcL_R1_.fastq.gz   62652  62652  30685
## data//LK653_rbcL_R1_.fastq.gz   83147  83147  43127
## data//LK661_rbcL_R1_.fastq.gz   98343  98343  55578
## data//LK665_rbcL_R1_.fastq.gz   71165  71165  41582
## data//LK668_rbcL_R1_.fastq.gz   73765  73765  44586
## data//LK670_rbcL_R1_.fastq.gz   87593  87593  42477
## data//LK682_rbcL_R1_.fastq.gz   22211  22211   6347
## data//LK685_rbcL_R1_.fastq.gz   38567  38567  11181
## data//LK686_rbcL_R1_.fastq.gz   38183  38183  11242
## data//PA0367_rbcL_R1_.fastq.gz  41356  41356  15938
## data//PA0368_rbcL_R1_.fastq.gz  73594  73594  41622
## data//PA0370_rbcL_R1_.fastq.gz  72373  72373  32032
## data//PA0456_rbcL_R1_.fastq.gz  33617  33617  12533
## data//PA1038_rbcL_R1_.fastq.gz  66394  66394  41350
## data//PA1039_rbcL_R1_.fastq.gz  75997  75997  46225
## data//PA1044_rbcL_R1_.fastq.gz  35224  35224  19129
## data//PA1049_rbcL_R1_.fastq.gz  49007  49007  24232
## data//PA1059_rbcL_R1_.fastq.gz  70463  70463  43821
## data//PA1065_rbcL_R1_.fastq.gz  81615  81615  46783
## data//TL3793_rbcL_R1_.fastq.gz  63067  63067  28465
## data//TL3797_rbcL_R1_.fastq.gz  63333  63333  35573
## data//TL3814_rbcL_R1_.fastq.gz  51654  51654  29185
## data//TL3816_rbcL_R1_.fastq.gz  50622  50622  27000
## data//TL3820_rbcL_R1_.fastq.gz  31450  31450  17230
## data//TL3821_rbcL_R1_.fastq.gz  44491  44491  25940
## data//TL3824_rbcL_R1_.fastq.gz  63570  63570  38574
## data//TL3826_rbcL_R1_.fastq.gz  51304  51304  27871
## data//TL3838_rbcL_R1_.fastq.gz  60441  60441  30161
## data//TL3842_rbcL_R1_.fastq.gz  82661  82661  37857
## data//TL3856_rbcL_R1_.fastq.gz  66213  66213  42606
## data//TL3862_rbcL_R1_.fastq.gz  81095  81095  48460
## data//TL3882_rbcL_R1_.fastq.gz  82875  82875  45432
## data//TL3889_rbcL_R1_.fastq.gz  51434  51434  28092
## data//TL3905_rbcL_R1_.fastq.gz  57568  57568  34654
## data//TL3910_rbcL_R1_.fastq.gz  44510  44510  27865
## data//TL3911_rbcL_R1_.fastq.gz  62300  62300  30603
## data//TL3915_rbcL_R1_.fastq.gz  56357  56357  31221
## data//TL3916_rbcL_R1_.fastq.gz  46986  46986  27034
## data//TL3918_rbcL_R1_.fastq.gz  25204  25204  13423
## data//TL3925_rbcL_R1_.fastq.gz  38217  38217  21399
## data//TL3926_rbcL_R1_.fastq.gz  41420  41420  24414
## data//TL3927_rbcL_R1_.fastq.gz  18964  18964  10807
## data//TL3929_rbcL_R1_.fastq.gz  35209  35209  19817
## data//TL3932_rbcL_R1_.fastq.gz  31296  31296  20399
## data//TL3936_rbcL_R1_.fastq.gz  19524  19524  10383
## data//TL3939_rbcL_R1_.fastq.gz  47335  47335  30660
## data//TL3940_rbcL_R1_.fastq.gz  28553  28553  19082
## data//TL3942_rbcL_R1_.fastq.gz  36835  36835  23240
## data//TL3943_rbcL_R1_.fastq.gz  35710  35710  21884
## data//TL3944_rbcL_R1_.fastq.gz  15835  15835   4396
## data//TL3945_rbcL_R1_.fastq.gz  77135  77135  26415
## data//TL3946_rbcL_R1_.fastq.gz  68358  68358  20998
## data//TL3948_rbcL_R1_.fastq.gz  96681  96681  20639
## data//UB0439_rbcL_R1_.fastq.gz  80079  80079  43612
## data//UB0445_rbcL_R1_.fastq.gz  45547  45547  17013
## data//UB0599_rbcL_R1_.fastq.gz  24799  24799  10412
## data//UB1430_rbcL_R1_.fastq.gz  90585  90585  58837
## data//UB1435_rbcL_R1_.fastq.gz  56934  56934  21268
## data//UB1446_rbcL_R1_.fastq.gz  93903  93903  60458
## data//UB1452_rbcL_R1_.fastq.gz  25249  25249  16096
## data//UB1454_rbcL_R1_.fastq.gz  49844  49844  17972
## data//UB2037_rbcL_R1_.fastq.gz  91724  91724  46810
```

```
## Running swarm
```

```
## swarm -f -t 40 /tmp/RtmpX0M3mf/fileb79b4ad55825 -o /tmp/RtmpX0M3mf/fileb79b64a365dd -w /tmp/RtmpX0M3mf/fileb79b7e560527
```

```
## Discarding singleton OTUs
```

```
## Aligning with mafft
```

```
## mafft --thread 50 /tmp/RtmpX0M3mf/fileb79b35d4e73a|gzip>work/swarmPair/rbcL_align.fa.gz
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

