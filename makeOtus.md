## Load libraries

```r
# installed from https://github.com/sherrillmix/dnar 
library(dnar)
packageVersion('dnar')
```

```
## [1] '0.1'
```

```r
library(parallel)
packageVersion('parallel')
```

```
## [1] '3.4.1'
```

```r
source('functions.R')
```
## Check software versions

```r
suppressWarnings(system('mafft --help 2>&1',intern=TRUE))
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
system('dpkg -s fasttree 2>&1',intern=TRUE)
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
suppressWarnings(system('swarm --version 2>&1',intern=TRUE))
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

## Primers used

```r
primerSeqs<-list('rbcl'=c('rbcL634F'='ATGCGTTGGAGAGACCGTTTC','rbcLbR'='TCGGTYAGAGCRGGCATRTGCCA'),'matk'=c('matK472F'='CCCRTYCATCTGGAAATCTTGGTTC','matK1248R'='GCTRTRATAATGAGAAAGATTTCTGC'))
```

## Find fastqs

```r
fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',basename(fastqs))
primerBases<-sub('[0-9]$','',primers)
```


## Run swarm, mafft and fasttree

```r
for(primerBase in unique(primerBases)){
  message('Working on ',primerBase)
  if(!dir.exists('work/swarmPair'))dir.create('work/swarmPair',recursive=TRUE)
  outMat<-sprintf('work/swarmPair/%s.Rdat',primerBase)
  outFa<-sprintf('work/swarmPair/%s.fa.gz',primerBase)
  outAlign<-sprintf('work/swarmPair/%s_align.fa.gz',primerBase)
  outTree<-sprintf('work/swarmPair/%s_align.tre',primerBase)
  message('Trimming primers off the start of reads')
  trimReads<-lapply(sprintf('%s%d',primerBase,1:2),function(ii){
    thisPrimer<-primerSeqs[[sub('[12]$','',tolower(ii))]][as.numeric(substring(ii,nchar(ii)))]
    thisFiles<-fastqs[primers==ii]
    reads<-mclapply(thisFiles,function(xx){library(dnar);cat('.');read.fastq(xx)},mc.cores=20,mc.preschedule=FALSE)
    if(mean(unlist(lapply(reads,function(xx)substring(xx$seq,1,nchar(thisPrimer))))%in% expandAmbiguous(thisPrimer)[[1]])<.75)stop(simpleError('Expected primer does not match read start'))
    trimReads<-lapply(reads,function(xx){xx$seq<-substring(xx$seq,nchar(thisPrimer)+1);return(xx)})
    names(trimReads)<-thisFiles
    return(trimReads)
  })
  if(any(names(trimReads[[1]])!=sub('_R2_','_R1_',names(trimReads[[2]]))))stop('Read 1 vs reads 2 file name mismatch')
  readCounts<-as.data.frame(sapply(trimReads,function(xx)sapply(xx,nrow)))
  colnames(readCounts)<-c('raw1','raw2')
  message('Discarding reads with >1 expected error in left or right read and concatenating left-right')
  trimReads<-mcmapply(function(left,right,...){
    cat('.')
    if(any(sub(' .*$','',left$name)!=sub(' .*$','',right$name)))stop('Read 1 vs read 2 name mismatch')
    #last base always low qual so ignore
    q1<-sapply(qualToInts(substring(left$qual,1,nchar(left$qual)-1)),function(xx)sum(10^(-xx/10)))
    q2<-sapply(qualToInts(substring(right$qual,1,nchar(right$qual)-1)),function(xx)sum(10^(-xx/10)))
    #less than 1 expected error in both reads
    selector<-q1<1&q2<1 & !grepl('[^ACTG]',left$seq)&!grepl('[^ACTG]',right$seq)
    seqs<-paste(left[selector,'seq'],revComp(right[selector,'seq']),sep='')
    return(seqs)
  },trimReads[[1]],trimReads[[2]],mc.cores=10,SIMPLIFY=FALSE)
  readCounts$filter<-sapply(trimReads,length)
  message('Read counts')
  print(readCounts)
  samples<-rep(basename(names(trimReads)),sapply(trimReads,length))
  message('Running swarm')
  otus<-runSwarm(unlist(trimReads),'~/installs/swarm/swarm',swarmArgs='-f -t 40')
  swarmOtus<-as.data.frame.matrix(table(samples,otus[['otus']]))
  write.fa(otus[['seqs']]$name,otus[['seqs']]$seq,outFa)
  save(swarmOtus,file=outMat)
  load(outMat)
  tmpFile<-tempfile()
  message('Discarding singleton OTUs')
  swarmOtus<-swarmOtus[,apply(swarmOtus,2,sum)>1]
  seqs<-otus[['seqs']]
  rownames(seqs)<-seqs$name
  write.fa(colnames(swarmOtus),seqs[colnames(swarmOtus),'seq'],tmpFile)
  message('Aligning with mafft')
  cmd<-sprintf('~/installs/mafft/bin/mafft --thread 50 %s|gzip>%s',tmpFile,outAlign)
  message(cmd)
  system(cmd)
  message('Creating tree with fasttree')
  cmd<-sprintf('zcat %s|fasttree -gtr -nt>%s',outAlign,outTree)
  message(cmd)
  system(cmd)
  message('Output files in:')
  message('    ',outMat)
  message('    ',outFa)
  message('    ',outAlign)
  message('    ',outTree)
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
## data//BI0054_matK_R1_.fastq.gz 37716 37716  26073
## data//BI0055_matK_R1_.fastq.gz 29870 29870  23513
## data//BI0093_matK_R1_.fastq.gz 23870 23870  10092
## data//BI0097_matK_R1_.fastq.gz 29597 29597  22221
## data//BI0246_matK_R1_.fastq.gz 23907 23907  18027
## data//BI0248_matK_R1_.fastq.gz 52496 52496  33150
## data//BI0257_matK_R1_.fastq.gz 45035 45035  31673
## data//BI0260_matK_R1_.fastq.gz 38342 38342  18338
## data//BI2414_matK_R1_.fastq.gz 42410 42410  28511
## data//BI2415_matK_R1_.fastq.gz 21547 21547  14033
## data//IK3158_matK_R1_.fastq.gz 22764 22764  16654
## data//IK3276_matK_R1_.fastq.gz 51495 51495  40161
## data//IK3358_matK_R1_.fastq.gz 22084 22084  17191
## data//IK3469_matK_R1_.fastq.gz 35564 35564  26925
## data//IK3513_matK_R1_.fastq.gz 43770 43770  29588
## data//IK3650_matK_R1_.fastq.gz 35772 35772  28257
## data//IK3701_matK_R1_.fastq.gz 53473 53473  41066
## data//IK3777_matK_R1_.fastq.gz 24679 24679  15380
## data//IK4184_matK_R1_.fastq.gz 31518 31518  18712
## data//IK4214_matK_R1_.fastq.gz 13382 13382    801
## data//KR02_matK_R1_.fastq.gz   37289 37289  19915
## data//KR05_matK_R1_.fastq.gz   17096 17096   7503
## data//KR07_matK_R1_.fastq.gz   49192 49192  35617
## data//KR10_matK_R1_.fastq.gz   39057 39057  27940
## data//KR12_matK_R1_.fastq.gz   16755 16755  12276
## data//KR21_matK_R1_.fastq.gz   31486 31486  21340
## data//KR33_matK_R1_.fastq.gz   31554 31554  22179
## data//KR35_matK_R1_.fastq.gz   46935 46935  35174
## data//KR52_matK_R1_.fastq.gz   50819 50819  38528
## data//KR57_matK_R1_.fastq.gz   22953 22953   9926
## data//KR67_matK_R1_.fastq.gz   18119 18119  12177
## data//LG4300_matK_R1_.fastq.gz 31703 31703  19577
## data//LG4314_matK_R1_.fastq.gz 24115 24115  17599
## data//LG4322_matK_R1_.fastq.gz 26690 26690  20210
## data//LG4327_matK_R1_.fastq.gz 36202 36202  27889
## data//LK645_matK_R1_.fastq.gz  28947 28947  22350
## data//LK647_matK_R1_.fastq.gz  39875 39875  32155
## data//LK653_matK_R1_.fastq.gz  40197 40197  29558
## data//LK661_matK_R1_.fastq.gz  33390 33390  16452
## data//LK665_matK_R1_.fastq.gz  54621 54621  41331
## data//LK668_matK_R1_.fastq.gz  30992 30992  23630
## data//LK670_matK_R1_.fastq.gz  43269 43269  31813
## data//LK682_matK_R1_.fastq.gz  32752 32752  24215
## data//LK685_matK_R1_.fastq.gz  26060 26060   9180
## data//LK686_matK_R1_.fastq.gz  24219 24219  10801
## data//PA0367_matK_R1_.fastq.gz 26739 26739  12832
## data//PA0368_matK_R1_.fastq.gz 23805 23805   9937
## data//PA0370_matK_R1_.fastq.gz 19569 19569   3152
## data//PA0456_matK_R1_.fastq.gz 19652 19652   8279
## data//PA1038_matK_R1_.fastq.gz 25960 25960  17685
## data//PA1039_matK_R1_.fastq.gz 26785 26785  18058
## data//PA1044_matK_R1_.fastq.gz 54164 54164  37088
## data//PA1049_matK_R1_.fastq.gz 46651 46651  34806
## data//PA1059_matK_R1_.fastq.gz 37381 37381  25263
## data//PA1065_matK_R1_.fastq.gz 30597 30597  22017
## data//TL3793_matK_R1_.fastq.gz 21275 21275  14047
## data//TL3797_matK_R1_.fastq.gz 11636 11636   8184
## data//TL3814_matK_R1_.fastq.gz 50184 50184  34233
## data//TL3816_matK_R1_.fastq.gz 46883 46883  32265
## data//TL3820_matK_R1_.fastq.gz 66064 66064  22533
## data//TL3821_matK_R1_.fastq.gz 24907 24907  15491
## data//TL3824_matK_R1_.fastq.gz 29323 29323  20044
## data//TL3826_matK_R1_.fastq.gz 61039 61039  37636
## data//TL3838_matK_R1_.fastq.gz 35465 35465  24724
## data//TL3842_matK_R1_.fastq.gz 38016 38016  23573
## data//TL3856_matK_R1_.fastq.gz 39587 39587  23645
## data//TL3862_matK_R1_.fastq.gz 32951 32951  22043
## data//TL3882_matK_R1_.fastq.gz 15760 15760  10761
## data//TL3889_matK_R1_.fastq.gz 37674 37674  24025
## data//TL3905_matK_R1_.fastq.gz 20322 20322  14548
## data//TL3910_matK_R1_.fastq.gz 24363 24363  16564
## data//TL3911_matK_R1_.fastq.gz 19392 19392   9375
## data//TL3915_matK_R1_.fastq.gz 34545 34545  19411
## data//TL3916_matK_R1_.fastq.gz    35    35      6
## data//TL3918_matK_R1_.fastq.gz 19370 19370   8090
## data//TL3925_matK_R1_.fastq.gz 33913 33913  25148
## data//TL3926_matK_R1_.fastq.gz 27740 27740  19163
## data//TL3927_matK_R1_.fastq.gz 18085 18085  12155
## data//TL3929_matK_R1_.fastq.gz 17720 17720  12543
## data//TL3932_matK_R1_.fastq.gz 19395 19395  13956
## data//TL3936_matK_R1_.fastq.gz 13435 13435  10229
## data//TL3939_matK_R1_.fastq.gz 16014 16014  11328
## data//TL3940_matK_R1_.fastq.gz 12612 12612   7220
## data//TL3942_matK_R1_.fastq.gz 65864 65864  49513
## data//TL3943_matK_R1_.fastq.gz 62439 62439  45675
## data//TL3944_matK_R1_.fastq.gz 48857 48857  35878
## data//TL3945_matK_R1_.fastq.gz 51533 51533  37625
## data//TL3946_matK_R1_.fastq.gz 31692 31692  23760
## data//TL3948_matK_R1_.fastq.gz 58040 58040  39163
## data//UB0439_matK_R1_.fastq.gz 33727 33727  23452
## data//UB0445_matK_R1_.fastq.gz 21950 21950  10108
## data//UB0599_matK_R1_.fastq.gz 30433 30433  23163
## data//UB1430_matK_R1_.fastq.gz 36057 36057  25013
## data//UB1435_matK_R1_.fastq.gz 44581 44581  35384
## data//UB1446_matK_R1_.fastq.gz 22388 22388  14080
## data//UB1452_matK_R1_.fastq.gz 27983 27983  19799
## data//UB1454_matK_R1_.fastq.gz 17005 17005   6383
## data//UB2037_matK_R1_.fastq.gz 28886 28886  21143
```

```
## Running swarm
```

```
## Discarding singleton OTUs
```

```
## Aligning with mafft
```

```
## ~/installs/mafft/bin/mafft --thread 50 /tmp/RtmpLe6yQ0/filea1a6374c999f|gzip>work/swarmPair/matK_align.fa.gz
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
## data//BI0054_rbcL_R1_.fastq.gz 112777 112777  73958
## data//BI0055_rbcL_R1_.fastq.gz  40341  40341  18774
## data//BI0093_rbcL_R1_.fastq.gz  60295  60295  41305
## data//BI0097_rbcL_R1_.fastq.gz  38815  38815  26566
## data//BI0246_rbcL_R1_.fastq.gz  97069  97069  62450
## data//BI0248_rbcL_R1_.fastq.gz  63847  63847  39235
## data//BI0257_rbcL_R1_.fastq.gz  68795  68795  40798
## data//BI0260_rbcL_R1_.fastq.gz  75418  75418  50323
## data//BI2414_rbcL_R1_.fastq.gz  72678  72678  40729
## data//BI2415_rbcL_R1_.fastq.gz  72083  72083  50399
## data//IK3158_rbcL_R1_.fastq.gz  63877  63877  40937
## data//IK3276_rbcL_R1_.fastq.gz  75130  75130  49304
## data//IK3358_rbcL_R1_.fastq.gz  83741  83741  55911
## data//IK3469_rbcL_R1_.fastq.gz  57184  57184  37481
## data//IK3513_rbcL_R1_.fastq.gz  72466  72466  47834
## data//IK3650_rbcL_R1_.fastq.gz  39570  39570  25611
## data//IK3701_rbcL_R1_.fastq.gz  71382  71382  49404
## data//IK3777_rbcL_R1_.fastq.gz  65948  65948  43926
## data//IK4184_rbcL_R1_.fastq.gz  49004  49004  22209
## data//IK4214_rbcL_R1_.fastq.gz  42566  42566  19568
## data//KR02_rbcL_R1_.fastq.gz    46441  46441  30165
## data//KR05_rbcL_R1_.fastq.gz    14003  14003   5546
## data//KR07_rbcL_R1_.fastq.gz    82888  82888  44832
## data//KR10_rbcL_R1_.fastq.gz    71966  71966  33087
## data//KR12_rbcL_R1_.fastq.gz    65322  65322  39892
## data//KR21_rbcL_R1_.fastq.gz    43489  43489  27475
## data//KR33_rbcL_R1_.fastq.gz    50031  50031  30648
## data//KR35_rbcL_R1_.fastq.gz    64584  64584  42490
## data//KR52_rbcL_R1_.fastq.gz    25519  25519  16818
## data//KR57_rbcL_R1_.fastq.gz    45206  45206  20272
## data//KR67_rbcL_R1_.fastq.gz    74544  74544  47472
## data//LG4300_rbcL_R1_.fastq.gz  64265  64265  22515
## data//LG4314_rbcL_R1_.fastq.gz  73212  73212  46293
## data//LG4322_rbcL_R1_.fastq.gz  20771  20771    519
## data//LG4327_rbcL_R1_.fastq.gz  66867  66867  46913
## data//LK645_rbcL_R1_.fastq.gz   87554  87554  43858
## data//LK647_rbcL_R1_.fastq.gz   62652  62652  32151
## data//LK653_rbcL_R1_.fastq.gz   83147  83147  45611
## data//LK661_rbcL_R1_.fastq.gz   98343  98343  59599
## data//LK665_rbcL_R1_.fastq.gz   71165  71165  43637
## data//LK668_rbcL_R1_.fastq.gz   73765  73765  47183
## data//LK670_rbcL_R1_.fastq.gz   87593  87593  44699
## data//LK682_rbcL_R1_.fastq.gz   22211  22211   6683
## data//LK685_rbcL_R1_.fastq.gz   38567  38567  12138
## data//LK686_rbcL_R1_.fastq.gz   38183  38183  12012
## data//PA0367_rbcL_R1_.fastq.gz  41356  41356  17128
## data//PA0368_rbcL_R1_.fastq.gz  73594  73594  43716
## data//PA0370_rbcL_R1_.fastq.gz  72373  72373  33608
## data//PA0456_rbcL_R1_.fastq.gz  33617  33617  13157
## data//PA1038_rbcL_R1_.fastq.gz  66394  66394  43390
## data//PA1039_rbcL_R1_.fastq.gz  75997  75997  48673
## data//PA1044_rbcL_R1_.fastq.gz  35224  35224  20087
## data//PA1049_rbcL_R1_.fastq.gz  49007  49007  25464
## data//PA1059_rbcL_R1_.fastq.gz  70463  70463  46012
## data//PA1065_rbcL_R1_.fastq.gz  81615  81615  49437
## data//TL3793_rbcL_R1_.fastq.gz  63067  63067  30701
## data//TL3797_rbcL_R1_.fastq.gz  63333  63333  37742
## data//TL3814_rbcL_R1_.fastq.gz  51654  51654  31105
## data//TL3816_rbcL_R1_.fastq.gz  50622  50622  28237
## data//TL3820_rbcL_R1_.fastq.gz  31450  31450  18192
## data//TL3821_rbcL_R1_.fastq.gz  44491  44491  27286
## data//TL3824_rbcL_R1_.fastq.gz  63570  63570  40537
## data//TL3826_rbcL_R1_.fastq.gz  51304  51304  29295
## data//TL3838_rbcL_R1_.fastq.gz  60441  60441  31902
## data//TL3842_rbcL_R1_.fastq.gz  82661  82661  39880
## data//TL3856_rbcL_R1_.fastq.gz  66213  66213  45127
## data//TL3862_rbcL_R1_.fastq.gz  81095  81095  51790
## data//TL3882_rbcL_R1_.fastq.gz  82875  82875  47941
## data//TL3889_rbcL_R1_.fastq.gz  51434  51434  29859
## data//TL3905_rbcL_R1_.fastq.gz  57568  57568  36660
## data//TL3910_rbcL_R1_.fastq.gz  44510  44510  29473
## data//TL3911_rbcL_R1_.fastq.gz  62300  62300  32203
## data//TL3915_rbcL_R1_.fastq.gz  56357  56357  32800
## data//TL3916_rbcL_R1_.fastq.gz  46986  46986  28441
## data//TL3918_rbcL_R1_.fastq.gz  25204  25204  14408
## data//TL3925_rbcL_R1_.fastq.gz  38217  38217  23089
## data//TL3926_rbcL_R1_.fastq.gz  41420  41420  26002
## data//TL3927_rbcL_R1_.fastq.gz  18964  18964  11599
## data//TL3929_rbcL_R1_.fastq.gz  35209  35209  21349
## data//TL3932_rbcL_R1_.fastq.gz  31296  31296  22154
## data//TL3936_rbcL_R1_.fastq.gz  19524  19524  11240
## data//TL3939_rbcL_R1_.fastq.gz  47335  47335  32687
## data//TL3940_rbcL_R1_.fastq.gz  28553  28553  20397
## data//TL3942_rbcL_R1_.fastq.gz  36835  36835  24937
## data//TL3943_rbcL_R1_.fastq.gz  35710  35710  23403
## data//TL3944_rbcL_R1_.fastq.gz  15835  15835   4738
## data//TL3945_rbcL_R1_.fastq.gz  77135  77135  29743
## data//TL3946_rbcL_R1_.fastq.gz  68358  68358  23644
## data//TL3948_rbcL_R1_.fastq.gz  96681  96681  23264
## data//UB0439_rbcL_R1_.fastq.gz  80079  80079  45974
## data//UB0445_rbcL_R1_.fastq.gz  45547  45547  18294
## data//UB0599_rbcL_R1_.fastq.gz  24799  24799  11017
## data//UB1430_rbcL_R1_.fastq.gz  90585  90585  62323
## data//UB1435_rbcL_R1_.fastq.gz  56934  56934  22837
## data//UB1446_rbcL_R1_.fastq.gz  93903  93903  65113
## data//UB1452_rbcL_R1_.fastq.gz  25249  25249  16972
## data//UB1454_rbcL_R1_.fastq.gz  49844  49844  19228
## data//UB2037_rbcL_R1_.fastq.gz  91724  91724  49379
```

```
## Running swarm
```

```
## Discarding singleton OTUs
```

```
## Aligning with mafft
```

```
## ~/installs/mafft/bin/mafft --thread 50 /tmp/RtmpLe6yQ0/filea1a64bba67ea|gzip>work/swarmPair/rbcL_align.fa.gz
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

