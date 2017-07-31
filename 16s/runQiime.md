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
source('../functions.R')
```

## Software versions

```r
system('echo "source activate qiime1;print_qiime_config.py -t"|bash 2>&1',intern=TRUE)
```

```
##  [1] ""                                                                                                                                                                            
##  [2] "System information"                                                                                                                                                          
##  [3] "=================="                                                                                                                                                          
##  [4] "         Platform:\tlinux2"                                                                                                                                                  
##  [5] "   Python version:\t2.7.12 |Continuum Analytics, Inc.| (default, Jul  2 2016, 17:42:40)  [GCC 4.4.7 20120313 (Red Hat 4.4.7-1)]"                                             
##  [6] "Python executable:\t/home/scott/miniconda3/envs/qiime1/bin/python"                                                                                                           
##  [7] ""                                                                                                                                                                            
##  [8] "QIIME default reference information"                                                                                                                                         
##  [9] "==================================="                                                                                                                                         
## [10] "For details on what files are used as QIIME's default references, see here:"                                                                                                 
## [11] " https://github.com/biocore/qiime-default-reference/releases/tag/0.1.3"                                                                                                      
## [12] ""                                                                                                                                                                            
## [13] "Dependency versions"                                                                                                                                                         
## [14] "==================="                                                                                                                                                         
## [15] "          QIIME library version:\t1.9.1"                                                                                                                                     
## [16] "           QIIME script version:\t1.9.1"                                                                                                                                     
## [17] "qiime-default-reference version:\t0.1.3"                                                                                                                                     
## [18] "                  NumPy version:\t1.10.4"                                                                                                                                    
## [19] "                  SciPy version:\t0.17.1"                                                                                                                                    
## [20] "                 pandas version:\t0.18.1"                                                                                                                                    
## [21] "             matplotlib version:\t1.4.3"                                                                                                                                     
## [22] "            biom-format version:\t2.1.5"                                                                                                                                     
## [23] "                   h5py version:\t2.6.0 (HDF5 version: 1.8.16)"                                                                                                              
## [24] "                   qcli version:\t0.1.1"                                                                                                                                     
## [25] "                   pyqi version:\t0.3.2"                                                                                                                                     
## [26] "             scikit-bio version:\t0.2.3"                                                                                                                                     
## [27] "                 PyNAST version:\t1.2.2"                                                                                                                                     
## [28] "                Emperor version:\t0.9.51"                                                                                                                                    
## [29] "                burrito version:\t0.9.1"                                                                                                                                     
## [30] "       burrito-fillings version:\t0.1.1"                                                                                                                                     
## [31] "              sortmerna version:\tSortMeRNA version 2.0, 29/11/2014"                                                                                                         
## [32] "              sumaclust version:\tSUMACLUST Version 1.0.00"                                                                                                                  
## [33] "                  swarm version:\tSwarm 1.2.19 [Mar  1 2016 23:41:10]"                                                                                                       
## [34] "                          gdata:\tInstalled."                                                                                                                                
## [35] ""                                                                                                                                                                            
## [36] "QIIME config values"                                                                                                                                                         
## [37] "==================="                                                                                                                                                         
## [38] "For definitions of these settings and to learn how to configure QIIME, see here:"                                                                                            
## [39] " http://qiime.org/install/qiime_config.html"                                                                                                                                 
## [40] " http://qiime.org/tutorials/parallel_qiime.html"                                                                                                                             
## [41] ""                                                                                                                                                                            
## [42] "                     blastmat_dir:\tNone"                                                                                                                                    
## [43] "      pick_otus_reference_seqs_fp:\t/home/scott/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta"               
## [44] "                         sc_queue:\tall.q"                                                                                                                                   
## [45] "      topiaryexplorer_project_dir:\tNone"                                                                                                                                    
## [46] "     pynast_template_alignment_fp:\t/home/scott/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set_aligned/85_otus.pynast.fasta"
## [47] "                  cluster_jobs_fp:\tstart_parallel_jobs.py"                                                                                                                  
## [48] "pynast_template_alignment_blastdb:\tNone"                                                                                                                                    
## [49] "assign_taxonomy_reference_seqs_fp:\t/home/scott/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta"               
## [50] "                     torque_queue:\tfriendlyq"                                                                                                                               
## [51] "                    jobs_to_start:\t1"                                                                                                                                       
## [52] "                       slurm_time:\tNone"                                                                                                                                    
## [53] "            denoiser_min_per_core:\t50"                                                                                                                                      
## [54] "assign_taxonomy_id_to_taxonomy_fp:\t/home/scott/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt"        
## [55] "                         temp_dir:\t/tmp/"                                                                                                                                   
## [56] "                     slurm_memory:\tNone"                                                                                                                                    
## [57] "                      slurm_queue:\tNone"                                                                                                                                    
## [58] "                      blastall_fp:\tblastall"                                                                                                                                
## [59] "                 seconds_to_sleep:\t1"                                                                                                                                       
## [60] ""                                                                                                                                                                            
## [61] "QIIME base install test results"                                                                                                                                             
## [62] "==============================="                                                                                                                                             
## [63] "........."                                                                                                                                                                   
## [64] "----------------------------------------------------------------------"                                                                                                      
## [65] "Ran 9 tests in 0.019s"                                                                                                                                                       
## [66] ""                                                                                                                                                                            
## [67] "OK"
```

```r
system('bbmerge.sh --version 2>&1',intern=TRUE)
```

```
## [1] "java -Djava.library.path=/home/scott/installs/bbmap/jni/ -ea -Xmx1000m -Xmx1000m -cp /home/scott/installs/bbmap/current/ jgi.BBMerge --version"
## [2] "BBMap version 36.77"                                                                                                                           
## [3] "For help, please run the shellscript with no parameters, or look in /docs/."
```

## Require 15000 reads in a sample

```r
nRequiredReads<-15000
```

## Merge reads

```r
  if(!dir.exists('data/joined'))dir.create('data/joined')
  fastq1s<-list.files('data','_R1_.*\\.fastq.gz',full.name=TRUE)
  mclapply(fastq1s,function(ii){
    out<-sprintf('data/joined/%s.fastq',sub('_R1_.*$','',basename(ii)))
    cmd<-sprintf('bbmerge.sh in1=%s in2=%s out=%s t=10 2>%s',ii,sub('_R1_','_R2_',ii),out,sub('fastq$','out',out))
    exit<-system(cmd)
    if(exit!=0)stop("Problem running pairing in",ii)
    system(sprintf('gzip %s',out))
    return(cmd)
  },mc.cores=10)
```

```
## [[1]]
## [1] "bbmerge.sh in1=data/BI0054_16s_R1_.fastq.gz in2=data/BI0054_16s_R2_.fastq.gz out=data/joined/BI0054_16s.fastq t=10 2>data/joined/BI0054_16s.out"
## 
## [[2]]
## [1] "bbmerge.sh in1=data/BI0055_16s_R1_.fastq.gz in2=data/BI0055_16s_R2_.fastq.gz out=data/joined/BI0055_16s.fastq t=10 2>data/joined/BI0055_16s.out"
## 
## [[3]]
## [1] "bbmerge.sh in1=data/BI0093_16s_R1_.fastq.gz in2=data/BI0093_16s_R2_.fastq.gz out=data/joined/BI0093_16s.fastq t=10 2>data/joined/BI0093_16s.out"
## 
## [[4]]
## [1] "bbmerge.sh in1=data/BI0097_16s_R1_.fastq.gz in2=data/BI0097_16s_R2_.fastq.gz out=data/joined/BI0097_16s.fastq t=10 2>data/joined/BI0097_16s.out"
## 
## [[5]]
## [1] "bbmerge.sh in1=data/BI0246_16s_R1_.fastq.gz in2=data/BI0246_16s_R2_.fastq.gz out=data/joined/BI0246_16s.fastq t=10 2>data/joined/BI0246_16s.out"
## 
## [[6]]
## [1] "bbmerge.sh in1=data/BI0248_16s_R1_.fastq.gz in2=data/BI0248_16s_R2_.fastq.gz out=data/joined/BI0248_16s.fastq t=10 2>data/joined/BI0248_16s.out"
## 
## [[7]]
## [1] "bbmerge.sh in1=data/BI0257_16s_R1_.fastq.gz in2=data/BI0257_16s_R2_.fastq.gz out=data/joined/BI0257_16s.fastq t=10 2>data/joined/BI0257_16s.out"
## 
## [[8]]
## [1] "bbmerge.sh in1=data/BI0260_16s_R1_.fastq.gz in2=data/BI0260_16s_R2_.fastq.gz out=data/joined/BI0260_16s.fastq t=10 2>data/joined/BI0260_16s.out"
## 
## [[9]]
## [1] "bbmerge.sh in1=data/BI2414_16s_R1_.fastq.gz in2=data/BI2414_16s_R2_.fastq.gz out=data/joined/BI2414_16s.fastq t=10 2>data/joined/BI2414_16s.out"
## 
## [[10]]
## [1] "bbmerge.sh in1=data/BI2415_16s_R1_.fastq.gz in2=data/BI2415_16s_R2_.fastq.gz out=data/joined/BI2415_16s.fastq t=10 2>data/joined/BI2415_16s.out"
## 
## [[11]]
## [1] "bbmerge.sh in1=data/IK3158_16s_R1_.fastq.gz in2=data/IK3158_16s_R2_.fastq.gz out=data/joined/IK3158_16s.fastq t=10 2>data/joined/IK3158_16s.out"
## 
## [[12]]
## [1] "bbmerge.sh in1=data/IK3276_16s_R1_.fastq.gz in2=data/IK3276_16s_R2_.fastq.gz out=data/joined/IK3276_16s.fastq t=10 2>data/joined/IK3276_16s.out"
## 
## [[13]]
## [1] "bbmerge.sh in1=data/IK3358_16s_R1_.fastq.gz in2=data/IK3358_16s_R2_.fastq.gz out=data/joined/IK3358_16s.fastq t=10 2>data/joined/IK3358_16s.out"
## 
## [[14]]
## [1] "bbmerge.sh in1=data/IK3469_16s_R1_.fastq.gz in2=data/IK3469_16s_R2_.fastq.gz out=data/joined/IK3469_16s.fastq t=10 2>data/joined/IK3469_16s.out"
## 
## [[15]]
## [1] "bbmerge.sh in1=data/IK3513_16s_R1_.fastq.gz in2=data/IK3513_16s_R2_.fastq.gz out=data/joined/IK3513_16s.fastq t=10 2>data/joined/IK3513_16s.out"
## 
## [[16]]
## [1] "bbmerge.sh in1=data/IK3650_16s_R1_.fastq.gz in2=data/IK3650_16s_R2_.fastq.gz out=data/joined/IK3650_16s.fastq t=10 2>data/joined/IK3650_16s.out"
## 
## [[17]]
## [1] "bbmerge.sh in1=data/IK3701_16s_R1_.fastq.gz in2=data/IK3701_16s_R2_.fastq.gz out=data/joined/IK3701_16s.fastq t=10 2>data/joined/IK3701_16s.out"
## 
## [[18]]
## [1] "bbmerge.sh in1=data/IK3777_16s_R1_.fastq.gz in2=data/IK3777_16s_R2_.fastq.gz out=data/joined/IK3777_16s.fastq t=10 2>data/joined/IK3777_16s.out"
## 
## [[19]]
## [1] "bbmerge.sh in1=data/IK4184_16s_R1_.fastq.gz in2=data/IK4184_16s_R2_.fastq.gz out=data/joined/IK4184_16s.fastq t=10 2>data/joined/IK4184_16s.out"
## 
## [[20]]
## [1] "bbmerge.sh in1=data/IK4214_16s_R1_.fastq.gz in2=data/IK4214_16s_R2_.fastq.gz out=data/joined/IK4214_16s.fastq t=10 2>data/joined/IK4214_16s.out"
## 
## [[21]]
## [1] "bbmerge.sh in1=data/KR02_16s_R1_.fastq.gz in2=data/KR02_16s_R2_.fastq.gz out=data/joined/KR02_16s.fastq t=10 2>data/joined/KR02_16s.out"
## 
## [[22]]
## [1] "bbmerge.sh in1=data/KR05_16s_R1_.fastq.gz in2=data/KR05_16s_R2_.fastq.gz out=data/joined/KR05_16s.fastq t=10 2>data/joined/KR05_16s.out"
## 
## [[23]]
## [1] "bbmerge.sh in1=data/KR07_16s_R1_.fastq.gz in2=data/KR07_16s_R2_.fastq.gz out=data/joined/KR07_16s.fastq t=10 2>data/joined/KR07_16s.out"
## 
## [[24]]
## [1] "bbmerge.sh in1=data/KR10_16s_R1_.fastq.gz in2=data/KR10_16s_R2_.fastq.gz out=data/joined/KR10_16s.fastq t=10 2>data/joined/KR10_16s.out"
## 
## [[25]]
## [1] "bbmerge.sh in1=data/KR12_16s_R1_.fastq.gz in2=data/KR12_16s_R2_.fastq.gz out=data/joined/KR12_16s.fastq t=10 2>data/joined/KR12_16s.out"
## 
## [[26]]
## [1] "bbmerge.sh in1=data/KR21_16s_R1_.fastq.gz in2=data/KR21_16s_R2_.fastq.gz out=data/joined/KR21_16s.fastq t=10 2>data/joined/KR21_16s.out"
## 
## [[27]]
## [1] "bbmerge.sh in1=data/KR33_16s_R1_.fastq.gz in2=data/KR33_16s_R2_.fastq.gz out=data/joined/KR33_16s.fastq t=10 2>data/joined/KR33_16s.out"
## 
## [[28]]
## [1] "bbmerge.sh in1=data/KR35_16s_R1_.fastq.gz in2=data/KR35_16s_R2_.fastq.gz out=data/joined/KR35_16s.fastq t=10 2>data/joined/KR35_16s.out"
## 
## [[29]]
## [1] "bbmerge.sh in1=data/KR52_16s_R1_.fastq.gz in2=data/KR52_16s_R2_.fastq.gz out=data/joined/KR52_16s.fastq t=10 2>data/joined/KR52_16s.out"
## 
## [[30]]
## [1] "bbmerge.sh in1=data/KR57_16s_R1_.fastq.gz in2=data/KR57_16s_R2_.fastq.gz out=data/joined/KR57_16s.fastq t=10 2>data/joined/KR57_16s.out"
## 
## [[31]]
## [1] "bbmerge.sh in1=data/KR67_16s_R1_.fastq.gz in2=data/KR67_16s_R2_.fastq.gz out=data/joined/KR67_16s.fastq t=10 2>data/joined/KR67_16s.out"
## 
## [[32]]
## [1] "bbmerge.sh in1=data/LG4300_16s_R1_.fastq.gz in2=data/LG4300_16s_R2_.fastq.gz out=data/joined/LG4300_16s.fastq t=10 2>data/joined/LG4300_16s.out"
## 
## [[33]]
## [1] "bbmerge.sh in1=data/LG4314_16s_R1_.fastq.gz in2=data/LG4314_16s_R2_.fastq.gz out=data/joined/LG4314_16s.fastq t=10 2>data/joined/LG4314_16s.out"
## 
## [[34]]
## [1] "bbmerge.sh in1=data/LG4322_16s_R1_.fastq.gz in2=data/LG4322_16s_R2_.fastq.gz out=data/joined/LG4322_16s.fastq t=10 2>data/joined/LG4322_16s.out"
## 
## [[35]]
## [1] "bbmerge.sh in1=data/LG4327_16s_R1_.fastq.gz in2=data/LG4327_16s_R2_.fastq.gz out=data/joined/LG4327_16s.fastq t=10 2>data/joined/LG4327_16s.out"
## 
## [[36]]
## [1] "bbmerge.sh in1=data/LK645_16s_R1_.fastq.gz in2=data/LK645_16s_R2_.fastq.gz out=data/joined/LK645_16s.fastq t=10 2>data/joined/LK645_16s.out"
## 
## [[37]]
## [1] "bbmerge.sh in1=data/LK647_16s_R1_.fastq.gz in2=data/LK647_16s_R2_.fastq.gz out=data/joined/LK647_16s.fastq t=10 2>data/joined/LK647_16s.out"
## 
## [[38]]
## [1] "bbmerge.sh in1=data/LK653_16s_R1_.fastq.gz in2=data/LK653_16s_R2_.fastq.gz out=data/joined/LK653_16s.fastq t=10 2>data/joined/LK653_16s.out"
## 
## [[39]]
## [1] "bbmerge.sh in1=data/LK661_16s_R1_.fastq.gz in2=data/LK661_16s_R2_.fastq.gz out=data/joined/LK661_16s.fastq t=10 2>data/joined/LK661_16s.out"
## 
## [[40]]
## [1] "bbmerge.sh in1=data/LK665_16s_R1_.fastq.gz in2=data/LK665_16s_R2_.fastq.gz out=data/joined/LK665_16s.fastq t=10 2>data/joined/LK665_16s.out"
## 
## [[41]]
## [1] "bbmerge.sh in1=data/LK668_16s_R1_.fastq.gz in2=data/LK668_16s_R2_.fastq.gz out=data/joined/LK668_16s.fastq t=10 2>data/joined/LK668_16s.out"
## 
## [[42]]
## [1] "bbmerge.sh in1=data/LK670_16s_R1_.fastq.gz in2=data/LK670_16s_R2_.fastq.gz out=data/joined/LK670_16s.fastq t=10 2>data/joined/LK670_16s.out"
## 
## [[43]]
## [1] "bbmerge.sh in1=data/LK682_16s_R1_.fastq.gz in2=data/LK682_16s_R2_.fastq.gz out=data/joined/LK682_16s.fastq t=10 2>data/joined/LK682_16s.out"
## 
## [[44]]
## [1] "bbmerge.sh in1=data/LK685_16s_R1_.fastq.gz in2=data/LK685_16s_R2_.fastq.gz out=data/joined/LK685_16s.fastq t=10 2>data/joined/LK685_16s.out"
## 
## [[45]]
## [1] "bbmerge.sh in1=data/LK686_16s_R1_.fastq.gz in2=data/LK686_16s_R2_.fastq.gz out=data/joined/LK686_16s.fastq t=10 2>data/joined/LK686_16s.out"
## 
## [[46]]
## [1] "bbmerge.sh in1=data/PA0367_16s_R1_.fastq.gz in2=data/PA0367_16s_R2_.fastq.gz out=data/joined/PA0367_16s.fastq t=10 2>data/joined/PA0367_16s.out"
## 
## [[47]]
## [1] "bbmerge.sh in1=data/PA0368_16s_R1_.fastq.gz in2=data/PA0368_16s_R2_.fastq.gz out=data/joined/PA0368_16s.fastq t=10 2>data/joined/PA0368_16s.out"
## 
## [[48]]
## [1] "bbmerge.sh in1=data/PA0370_16s_R1_.fastq.gz in2=data/PA0370_16s_R2_.fastq.gz out=data/joined/PA0370_16s.fastq t=10 2>data/joined/PA0370_16s.out"
## 
## [[49]]
## [1] "bbmerge.sh in1=data/PA0456_16s_R1_.fastq.gz in2=data/PA0456_16s_R2_.fastq.gz out=data/joined/PA0456_16s.fastq t=10 2>data/joined/PA0456_16s.out"
## 
## [[50]]
## [1] "bbmerge.sh in1=data/PA1038_16s_R1_.fastq.gz in2=data/PA1038_16s_R2_.fastq.gz out=data/joined/PA1038_16s.fastq t=10 2>data/joined/PA1038_16s.out"
## 
## [[51]]
## [1] "bbmerge.sh in1=data/PA1039_16s_R1_.fastq.gz in2=data/PA1039_16s_R2_.fastq.gz out=data/joined/PA1039_16s.fastq t=10 2>data/joined/PA1039_16s.out"
## 
## [[52]]
## [1] "bbmerge.sh in1=data/PA1044_16s_R1_.fastq.gz in2=data/PA1044_16s_R2_.fastq.gz out=data/joined/PA1044_16s.fastq t=10 2>data/joined/PA1044_16s.out"
## 
## [[53]]
## [1] "bbmerge.sh in1=data/PA1049_16s_R1_.fastq.gz in2=data/PA1049_16s_R2_.fastq.gz out=data/joined/PA1049_16s.fastq t=10 2>data/joined/PA1049_16s.out"
## 
## [[54]]
## [1] "bbmerge.sh in1=data/PA1059_16s_R1_.fastq.gz in2=data/PA1059_16s_R2_.fastq.gz out=data/joined/PA1059_16s.fastq t=10 2>data/joined/PA1059_16s.out"
## 
## [[55]]
## [1] "bbmerge.sh in1=data/PA1065_16s_R1_.fastq.gz in2=data/PA1065_16s_R2_.fastq.gz out=data/joined/PA1065_16s.fastq t=10 2>data/joined/PA1065_16s.out"
## 
## [[56]]
## [1] "bbmerge.sh in1=data/TL3793_16s_R1_.fastq.gz in2=data/TL3793_16s_R2_.fastq.gz out=data/joined/TL3793_16s.fastq t=10 2>data/joined/TL3793_16s.out"
## 
## [[57]]
## [1] "bbmerge.sh in1=data/TL3797_16s_R1_.fastq.gz in2=data/TL3797_16s_R2_.fastq.gz out=data/joined/TL3797_16s.fastq t=10 2>data/joined/TL3797_16s.out"
## 
## [[58]]
## [1] "bbmerge.sh in1=data/TL3814_16s_R1_.fastq.gz in2=data/TL3814_16s_R2_.fastq.gz out=data/joined/TL3814_16s.fastq t=10 2>data/joined/TL3814_16s.out"
## 
## [[59]]
## [1] "bbmerge.sh in1=data/TL3816_16s_R1_.fastq.gz in2=data/TL3816_16s_R2_.fastq.gz out=data/joined/TL3816_16s.fastq t=10 2>data/joined/TL3816_16s.out"
## 
## [[60]]
## [1] "bbmerge.sh in1=data/TL3820_16s_R1_.fastq.gz in2=data/TL3820_16s_R2_.fastq.gz out=data/joined/TL3820_16s.fastq t=10 2>data/joined/TL3820_16s.out"
## 
## [[61]]
## [1] "bbmerge.sh in1=data/TL3821_16s_R1_.fastq.gz in2=data/TL3821_16s_R2_.fastq.gz out=data/joined/TL3821_16s.fastq t=10 2>data/joined/TL3821_16s.out"
## 
## [[62]]
## [1] "bbmerge.sh in1=data/TL3824_16s_R1_.fastq.gz in2=data/TL3824_16s_R2_.fastq.gz out=data/joined/TL3824_16s.fastq t=10 2>data/joined/TL3824_16s.out"
## 
## [[63]]
## [1] "bbmerge.sh in1=data/TL3826_16s_R1_.fastq.gz in2=data/TL3826_16s_R2_.fastq.gz out=data/joined/TL3826_16s.fastq t=10 2>data/joined/TL3826_16s.out"
## 
## [[64]]
## [1] "bbmerge.sh in1=data/TL3838_16s_R1_.fastq.gz in2=data/TL3838_16s_R2_.fastq.gz out=data/joined/TL3838_16s.fastq t=10 2>data/joined/TL3838_16s.out"
## 
## [[65]]
## [1] "bbmerge.sh in1=data/TL3842_16s_R1_.fastq.gz in2=data/TL3842_16s_R2_.fastq.gz out=data/joined/TL3842_16s.fastq t=10 2>data/joined/TL3842_16s.out"
## 
## [[66]]
## [1] "bbmerge.sh in1=data/TL3856_16s_R1_.fastq.gz in2=data/TL3856_16s_R2_.fastq.gz out=data/joined/TL3856_16s.fastq t=10 2>data/joined/TL3856_16s.out"
## 
## [[67]]
## [1] "bbmerge.sh in1=data/TL3862_16s_R1_.fastq.gz in2=data/TL3862_16s_R2_.fastq.gz out=data/joined/TL3862_16s.fastq t=10 2>data/joined/TL3862_16s.out"
## 
## [[68]]
## [1] "bbmerge.sh in1=data/TL3882_16s_R1_.fastq.gz in2=data/TL3882_16s_R2_.fastq.gz out=data/joined/TL3882_16s.fastq t=10 2>data/joined/TL3882_16s.out"
## 
## [[69]]
## [1] "bbmerge.sh in1=data/TL3889_16s_R1_.fastq.gz in2=data/TL3889_16s_R2_.fastq.gz out=data/joined/TL3889_16s.fastq t=10 2>data/joined/TL3889_16s.out"
## 
## [[70]]
## [1] "bbmerge.sh in1=data/TL3905_16s_R1_.fastq.gz in2=data/TL3905_16s_R2_.fastq.gz out=data/joined/TL3905_16s.fastq t=10 2>data/joined/TL3905_16s.out"
## 
## [[71]]
## [1] "bbmerge.sh in1=data/TL3910_16s_R1_.fastq.gz in2=data/TL3910_16s_R2_.fastq.gz out=data/joined/TL3910_16s.fastq t=10 2>data/joined/TL3910_16s.out"
## 
## [[72]]
## [1] "bbmerge.sh in1=data/TL3911_16s_R1_.fastq.gz in2=data/TL3911_16s_R2_.fastq.gz out=data/joined/TL3911_16s.fastq t=10 2>data/joined/TL3911_16s.out"
## 
## [[73]]
## [1] "bbmerge.sh in1=data/TL3915_16s_R1_.fastq.gz in2=data/TL3915_16s_R2_.fastq.gz out=data/joined/TL3915_16s.fastq t=10 2>data/joined/TL3915_16s.out"
## 
## [[74]]
## [1] "bbmerge.sh in1=data/TL3916_16s_R1_.fastq.gz in2=data/TL3916_16s_R2_.fastq.gz out=data/joined/TL3916_16s.fastq t=10 2>data/joined/TL3916_16s.out"
## 
## [[75]]
## [1] "bbmerge.sh in1=data/TL3918_16s_R1_.fastq.gz in2=data/TL3918_16s_R2_.fastq.gz out=data/joined/TL3918_16s.fastq t=10 2>data/joined/TL3918_16s.out"
## 
## [[76]]
## [1] "bbmerge.sh in1=data/TL3925_16s_R1_.fastq.gz in2=data/TL3925_16s_R2_.fastq.gz out=data/joined/TL3925_16s.fastq t=10 2>data/joined/TL3925_16s.out"
## 
## [[77]]
## [1] "bbmerge.sh in1=data/TL3926_16s_R1_.fastq.gz in2=data/TL3926_16s_R2_.fastq.gz out=data/joined/TL3926_16s.fastq t=10 2>data/joined/TL3926_16s.out"
## 
## [[78]]
## [1] "bbmerge.sh in1=data/TL3927_16s_R1_.fastq.gz in2=data/TL3927_16s_R2_.fastq.gz out=data/joined/TL3927_16s.fastq t=10 2>data/joined/TL3927_16s.out"
## 
## [[79]]
## [1] "bbmerge.sh in1=data/TL3929_16s_R1_.fastq.gz in2=data/TL3929_16s_R2_.fastq.gz out=data/joined/TL3929_16s.fastq t=10 2>data/joined/TL3929_16s.out"
## 
## [[80]]
## [1] "bbmerge.sh in1=data/TL3932_16s_R1_.fastq.gz in2=data/TL3932_16s_R2_.fastq.gz out=data/joined/TL3932_16s.fastq t=10 2>data/joined/TL3932_16s.out"
## 
## [[81]]
## [1] "bbmerge.sh in1=data/TL3936_16s_R1_.fastq.gz in2=data/TL3936_16s_R2_.fastq.gz out=data/joined/TL3936_16s.fastq t=10 2>data/joined/TL3936_16s.out"
## 
## [[82]]
## [1] "bbmerge.sh in1=data/TL3939_16s_R1_.fastq.gz in2=data/TL3939_16s_R2_.fastq.gz out=data/joined/TL3939_16s.fastq t=10 2>data/joined/TL3939_16s.out"
## 
## [[83]]
## [1] "bbmerge.sh in1=data/TL3940_16s_R1_.fastq.gz in2=data/TL3940_16s_R2_.fastq.gz out=data/joined/TL3940_16s.fastq t=10 2>data/joined/TL3940_16s.out"
## 
## [[84]]
## [1] "bbmerge.sh in1=data/TL3942_16s_R1_.fastq.gz in2=data/TL3942_16s_R2_.fastq.gz out=data/joined/TL3942_16s.fastq t=10 2>data/joined/TL3942_16s.out"
## 
## [[85]]
## [1] "bbmerge.sh in1=data/TL3943_16s_R1_.fastq.gz in2=data/TL3943_16s_R2_.fastq.gz out=data/joined/TL3943_16s.fastq t=10 2>data/joined/TL3943_16s.out"
## 
## [[86]]
## [1] "bbmerge.sh in1=data/TL3944_16s_R1_.fastq.gz in2=data/TL3944_16s_R2_.fastq.gz out=data/joined/TL3944_16s.fastq t=10 2>data/joined/TL3944_16s.out"
## 
## [[87]]
## [1] "bbmerge.sh in1=data/TL3945_16s_R1_.fastq.gz in2=data/TL3945_16s_R2_.fastq.gz out=data/joined/TL3945_16s.fastq t=10 2>data/joined/TL3945_16s.out"
## 
## [[88]]
## [1] "bbmerge.sh in1=data/TL3946_16s_R1_.fastq.gz in2=data/TL3946_16s_R2_.fastq.gz out=data/joined/TL3946_16s.fastq t=10 2>data/joined/TL3946_16s.out"
## 
## [[89]]
## [1] "bbmerge.sh in1=data/TL3948_16s_R1_.fastq.gz in2=data/TL3948_16s_R2_.fastq.gz out=data/joined/TL3948_16s.fastq t=10 2>data/joined/TL3948_16s.out"
## 
## [[90]]
## [1] "bbmerge.sh in1=data/UB0439_16s_R1_.fastq.gz in2=data/UB0439_16s_R2_.fastq.gz out=data/joined/UB0439_16s.fastq t=10 2>data/joined/UB0439_16s.out"
## 
## [[91]]
## [1] "bbmerge.sh in1=data/UB0445_16s_R1_.fastq.gz in2=data/UB0445_16s_R2_.fastq.gz out=data/joined/UB0445_16s.fastq t=10 2>data/joined/UB0445_16s.out"
## 
## [[92]]
## [1] "bbmerge.sh in1=data/UB0599_16s_R1_.fastq.gz in2=data/UB0599_16s_R2_.fastq.gz out=data/joined/UB0599_16s.fastq t=10 2>data/joined/UB0599_16s.out"
## 
## [[93]]
## [1] "bbmerge.sh in1=data/UB1430_16s_R1_.fastq.gz in2=data/UB1430_16s_R2_.fastq.gz out=data/joined/UB1430_16s.fastq t=10 2>data/joined/UB1430_16s.out"
## 
## [[94]]
## [1] "bbmerge.sh in1=data/UB1435_16s_R1_.fastq.gz in2=data/UB1435_16s_R2_.fastq.gz out=data/joined/UB1435_16s.fastq t=10 2>data/joined/UB1435_16s.out"
## 
## [[95]]
## [1] "bbmerge.sh in1=data/UB1446_16s_R1_.fastq.gz in2=data/UB1446_16s_R2_.fastq.gz out=data/joined/UB1446_16s.fastq t=10 2>data/joined/UB1446_16s.out"
## 
## [[96]]
## [1] "bbmerge.sh in1=data/UB1452_16s_R1_.fastq.gz in2=data/UB1452_16s_R2_.fastq.gz out=data/joined/UB1452_16s.fastq t=10 2>data/joined/UB1452_16s.out"
## 
## [[97]]
## [1] "bbmerge.sh in1=data/UB1454_16s_R1_.fastq.gz in2=data/UB1454_16s_R2_.fastq.gz out=data/joined/UB1454_16s.fastq t=10 2>data/joined/UB1454_16s.out"
## 
## [[98]]
## [1] "bbmerge.sh in1=data/UB2037_16s_R1_.fastq.gz in2=data/UB2037_16s_R2_.fastq.gz out=data/joined/UB2037_16s.fastq t=10 2>data/joined/UB2037_16s.out"
```

## Read in all sequences

```r
fastqs<-list.files('data/joined','.fastq.gz',full.name=TRUE)
#just get sequences to reduce memory
allSeq<-mclapply(fastqs,function(xx)read.fastq(xx)$seq,mc.cores=12)
message('Read ',length(unlist(allSeq)),' sequences in ',length(allSeq),' fastqs')
```

```
## Read 6574907 sequences in 98 fastqs
```

## Run qiime

```r
if(!dir.exists('work/'))dir.create('work')
out<-runQiime(unlist(allSeq),storeDir='work/qiime')
```

```
## echo "source activate qiime1; pick_de_novo_otus.py --input /tmp/Rtmp0111pK/file450b380b69a6/XXX.fa --output /tmp/Rtmp0111pK/file450b4f21d341 --parallel --jobs_to_start 32 --force"|bash
```

```r
outDf<-data.frame(
  'file'=rep(sub('_16s.fastq.gz','',basename(fastqs)),sapply(allSeq,length)),
  'otu'=out[['otus']],
  stringsAsFactors=FALSE
)
withAs(outFile=gzfile('work/qiimeOtuIds.csv.gz'),write.csv(outDf,outFile,row.names=FALSE))
write.fa(names(out[['seqs']]),out[['seqs']],'work/qiimeOtus.fa.gz')
outTaxa<-data.frame(
  'name'=names(out[['taxa']]),
  'taxa'=out[['taxa']],
  stringsAsFactors=FALSE
)
write.csv(outTaxa,'work/qiimeOtus.taxa',row.names=FALSE)
```
