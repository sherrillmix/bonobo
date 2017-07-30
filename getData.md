## Read in SRA run table
The SRA run table was downloaded from [https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP108776](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP108776).

```r
samples<-read.table('SraRunTable.txt',sep='\t',header=TRUE,stringsAsFactors=FALSE)
rownames(samples)<-samples$Run_s
samples[,c('sample','primer')]<-do.call(rbind,strsplit(samples$Library_Name_s,'_'))
```

## Download data from SRA

```r
system('fastq-dump --version')
```


```r
if(!dir.exists('sra'))dir.create('sra')
for(ii in samples$Run_s){
  outFiles<-sprintf('sra/%s_%d.fastq.gz',ii,1:2)
  cmd<-sprintf('fastq-dump --gzip --split-files --outdir data %s',ii)
  message(cmd)
  if(any(!file.exists(outFiles)))system(cmd)
  if(any(!file.exists(outFiles)))stop('fastq-dump unsucessful')
}
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656515
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656723
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656724
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656721
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656722
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656520
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656517
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656719
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656720
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656717
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656514
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656521
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656581
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656582
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656620
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656621
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656646
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656647
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656699
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656700
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656701
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656702
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656710
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656711
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656712
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656713
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656714
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656715
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656716
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656718
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656432
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656433
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656479
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656604
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656642
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656643
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656644
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656645
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656648
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656649
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656652
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656654
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656655
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656662
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656664
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656665
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656666
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656667
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656674
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656675
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656686
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656687
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656688
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656689
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656690
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656691
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656692
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656693
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656694
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656695
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656696
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656697
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656698
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656703
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656704
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656705
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656706
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656707
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656708
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656709
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656431
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656434
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656440
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656441
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656442
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656451
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656462
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656463
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656464
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656465
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656483
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656488
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656489
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656490
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656491
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656492
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656493
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656497
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656498
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656499
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656500
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656501
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656512
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656513
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656524
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656525
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656528
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656529
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656535
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656544
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656545
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656546
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656548
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656566
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656568
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656579
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656601
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656602
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656603
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656612
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656619
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656629
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656630
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656631
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656632
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656633
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656634
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656635
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656636
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656637
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656638
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656639
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656640
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656641
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656650
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656651
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656653
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656656
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656657
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656658
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656659
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656660
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656661
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656663
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656668
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656669
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656670
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656671
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656672
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656673
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656676
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656677
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656678
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656679
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656680
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656681
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656682
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656683
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656684
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656685
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656435
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656436
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656437
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656438
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656439
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656443
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656444
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656445
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656446
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656447
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656448
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656449
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656450
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656452
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656453
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656454
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656455
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656456
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656457
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656458
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656459
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656460
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656461
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656466
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656467
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656468
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656469
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656470
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656471
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656472
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656473
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656474
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656475
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656476
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656477
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656478
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656480
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656481
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656482
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656484
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656485
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656486
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656487
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656494
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656495
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656496
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656502
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656503
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656504
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656505
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656506
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656507
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656508
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656509
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656510
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656511
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656516
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656518
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656519
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656522
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656523
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656526
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656527
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656530
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656531
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656532
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656533
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656534
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656536
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656537
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656538
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656539
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656540
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656541
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656542
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656543
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656547
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656549
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656550
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656551
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656552
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656553
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656554
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656555
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656556
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656557
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656558
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656559
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656560
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656561
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656562
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656563
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656564
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656565
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656567
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656569
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656570
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656571
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656572
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656573
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656574
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656575
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656576
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656577
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656578
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656580
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656583
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656584
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656585
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656586
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656587
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656588
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656589
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656590
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656591
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656592
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656593
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656594
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656595
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656596
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656597
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656598
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656599
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656600
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656605
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656606
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656607
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656608
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656609
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656610
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656611
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656613
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656614
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656615
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656616
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656617
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656618
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656622
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656623
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656624
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656625
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656626
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656627
```

```
## fastq-dump --gzip --split-files --outdir data SRR5656628
```

## Rename into convenient form

```r
if(!dir.exists('16s/data'))dir.create('16s/data')
if(!dir.exists('data'))dir.create('data')
for(ii in samples$Run_s){
  sraNames<-sprintf('sra/%s_%d.fastq.gz',ii,1:2)
  newNames<-sprintf('%s/%s_%s_R%d_.fastq.gz', ifelse(samples[ii,'primer']=='16s','16s/data','data'),
    samples[ii,'sample'], samples[ii,'primer'], 1:2)
  message('Moving ',paste(sraNames,collapse=' '),' to ',paste(newNames,collapse=' '))
  file.rename(sraNames,newNames)
}
```

```
## Moving sra/SRR5656515_1.fastq.gz sra/SRR5656515_2.fastq.gz to 16s/data/TL3814_16s_R1_.fastq.gz 16s/data/TL3814_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656723_1.fastq.gz sra/SRR5656723_2.fastq.gz to data/TL3814_matK_R1_.fastq.gz data/TL3814_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656724_1.fastq.gz sra/SRR5656724_2.fastq.gz to data/TL3814_rbcL_R1_.fastq.gz data/TL3814_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656721_1.fastq.gz sra/SRR5656721_2.fastq.gz to data/TL3816_matK_R1_.fastq.gz data/TL3816_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656722_1.fastq.gz sra/SRR5656722_2.fastq.gz to data/TL3816_rbcL_R1_.fastq.gz data/TL3816_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656520_1.fastq.gz sra/SRR5656520_2.fastq.gz to 16s/data/TL3816_16s_R1_.fastq.gz 16s/data/TL3816_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656517_1.fastq.gz sra/SRR5656517_2.fastq.gz to 16s/data/TL3793_16s_R1_.fastq.gz 16s/data/TL3793_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656719_1.fastq.gz sra/SRR5656719_2.fastq.gz to data/TL3793_matK_R1_.fastq.gz data/TL3793_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656720_1.fastq.gz sra/SRR5656720_2.fastq.gz to data/TL3793_rbcL_R1_.fastq.gz data/TL3793_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656717_1.fastq.gz sra/SRR5656717_2.fastq.gz to data/TL3797_matK_R1_.fastq.gz data/TL3797_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656514_1.fastq.gz sra/SRR5656514_2.fastq.gz to 16s/data/TL3797_16s_R1_.fastq.gz 16s/data/TL3797_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656521_1.fastq.gz sra/SRR5656521_2.fastq.gz to 16s/data/TL3820_16s_R1_.fastq.gz 16s/data/TL3820_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656581_1.fastq.gz sra/SRR5656581_2.fastq.gz to data/TL3943_rbcL_R1_.fastq.gz data/TL3943_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656582_1.fastq.gz sra/SRR5656582_2.fastq.gz to data/TL3943_matK_R1_.fastq.gz data/TL3943_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656620_1.fastq.gz sra/SRR5656620_2.fastq.gz to data/TL3944_matK_R1_.fastq.gz data/TL3944_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656621_1.fastq.gz sra/SRR5656621_2.fastq.gz to data/TL3944_rbcL_R1_.fastq.gz data/TL3944_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656646_1.fastq.gz sra/SRR5656646_2.fastq.gz to data/IK3358_rbcL_R1_.fastq.gz data/IK3358_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656647_1.fastq.gz sra/SRR5656647_2.fastq.gz to data/IK3358_matK_R1_.fastq.gz data/IK3358_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656699_1.fastq.gz sra/SRR5656699_2.fastq.gz to data/TL3946_matK_R1_.fastq.gz data/TL3946_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656700_1.fastq.gz sra/SRR5656700_2.fastq.gz to data/TL3946_rbcL_R1_.fastq.gz data/TL3946_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656701_1.fastq.gz sra/SRR5656701_2.fastq.gz to data/TL3945_matK_R1_.fastq.gz data/TL3945_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656702_1.fastq.gz sra/SRR5656702_2.fastq.gz to data/TL3945_rbcL_R1_.fastq.gz data/TL3945_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656710_1.fastq.gz sra/SRR5656710_2.fastq.gz to 16s/data/IK3358_16s_R1_.fastq.gz 16s/data/IK3358_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656711_1.fastq.gz sra/SRR5656711_2.fastq.gz to 16s/data/TL3944_16s_R1_.fastq.gz 16s/data/TL3944_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656712_1.fastq.gz sra/SRR5656712_2.fastq.gz to 16s/data/TL3943_16s_R1_.fastq.gz 16s/data/TL3943_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656713_1.fastq.gz sra/SRR5656713_2.fastq.gz to 16s/data/TL3946_16s_R1_.fastq.gz 16s/data/TL3946_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656714_1.fastq.gz sra/SRR5656714_2.fastq.gz to 16s/data/TL3945_16s_R1_.fastq.gz 16s/data/TL3945_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656715_1.fastq.gz sra/SRR5656715_2.fastq.gz to data/TL3820_matK_R1_.fastq.gz data/TL3820_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656716_1.fastq.gz sra/SRR5656716_2.fastq.gz to data/TL3820_rbcL_R1_.fastq.gz data/TL3820_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656718_1.fastq.gz sra/SRR5656718_2.fastq.gz to data/TL3797_rbcL_R1_.fastq.gz data/TL3797_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656432_1.fastq.gz sra/SRR5656432_2.fastq.gz to 16s/data/BI0260_16s_R1_.fastq.gz 16s/data/BI0260_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656433_1.fastq.gz sra/SRR5656433_2.fastq.gz to 16s/data/BI2414_16s_R1_.fastq.gz 16s/data/BI2414_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656479_1.fastq.gz sra/SRR5656479_2.fastq.gz to 16s/data/IK4184_16s_R1_.fastq.gz 16s/data/IK4184_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656604_1.fastq.gz sra/SRR5656604_2.fastq.gz to data/IK4184_matK_R1_.fastq.gz data/IK4184_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656642_1.fastq.gz sra/SRR5656642_2.fastq.gz to data/IK3158_rbcL_R1_.fastq.gz data/IK3158_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656643_1.fastq.gz sra/SRR5656643_2.fastq.gz to data/IK3158_matK_R1_.fastq.gz data/IK3158_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656644_1.fastq.gz sra/SRR5656644_2.fastq.gz to data/IK3276_rbcL_R1_.fastq.gz data/IK3276_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656645_1.fastq.gz sra/SRR5656645_2.fastq.gz to data/IK3276_matK_R1_.fastq.gz data/IK3276_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656648_1.fastq.gz sra/SRR5656648_2.fastq.gz to data/IK3469_rbcL_R1_.fastq.gz data/IK3469_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656649_1.fastq.gz sra/SRR5656649_2.fastq.gz to data/IK3469_matK_R1_.fastq.gz data/IK3469_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656652_1.fastq.gz sra/SRR5656652_2.fastq.gz to data/UB1452_rbcL_R1_.fastq.gz data/UB1452_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656654_1.fastq.gz sra/SRR5656654_2.fastq.gz to data/UB1454_rbcL_R1_.fastq.gz data/UB1454_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656655_1.fastq.gz sra/SRR5656655_2.fastq.gz to data/UB1454_matK_R1_.fastq.gz data/UB1454_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656662_1.fastq.gz sra/SRR5656662_2.fastq.gz to 16s/data/BI2415_16s_R1_.fastq.gz 16s/data/BI2415_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656664_1.fastq.gz sra/SRR5656664_2.fastq.gz to data/UB2037_rbcL_R1_.fastq.gz data/UB2037_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656665_1.fastq.gz sra/SRR5656665_2.fastq.gz to data/UB2037_matK_R1_.fastq.gz data/UB2037_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656666_1.fastq.gz sra/SRR5656666_2.fastq.gz to 16s/data/BI0248_16s_R1_.fastq.gz 16s/data/BI0248_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656667_1.fastq.gz sra/SRR5656667_2.fastq.gz to 16s/data/BI0257_16s_R1_.fastq.gz 16s/data/BI0257_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656674_1.fastq.gz sra/SRR5656674_2.fastq.gz to data/BI0054_matK_R1_.fastq.gz data/BI0054_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656675_1.fastq.gz sra/SRR5656675_2.fastq.gz to data/BI0054_rbcL_R1_.fastq.gz data/BI0054_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656686_1.fastq.gz sra/SRR5656686_2.fastq.gz to 16s/data/UB1454_16s_R1_.fastq.gz 16s/data/UB1454_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656687_1.fastq.gz sra/SRR5656687_2.fastq.gz to 16s/data/UB2037_16s_R1_.fastq.gz 16s/data/UB2037_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656688_1.fastq.gz sra/SRR5656688_2.fastq.gz to data/BI2415_matK_R1_.fastq.gz data/BI2415_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656689_1.fastq.gz sra/SRR5656689_2.fastq.gz to data/BI2415_rbcL_R1_.fastq.gz data/BI2415_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656690_1.fastq.gz sra/SRR5656690_2.fastq.gz to data/BI0260_matK_R1_.fastq.gz data/BI0260_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656691_1.fastq.gz sra/SRR5656691_2.fastq.gz to data/BI0260_rbcL_R1_.fastq.gz data/BI0260_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656692_1.fastq.gz sra/SRR5656692_2.fastq.gz to data/BI2414_matK_R1_.fastq.gz data/BI2414_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656693_1.fastq.gz sra/SRR5656693_2.fastq.gz to data/BI2414_rbcL_R1_.fastq.gz data/BI2414_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656694_1.fastq.gz sra/SRR5656694_2.fastq.gz to data/BI0248_matK_R1_.fastq.gz data/BI0248_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656695_1.fastq.gz sra/SRR5656695_2.fastq.gz to data/BI0248_rbcL_R1_.fastq.gz data/BI0248_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656696_1.fastq.gz sra/SRR5656696_2.fastq.gz to data/BI0257_matK_R1_.fastq.gz data/BI0257_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656697_1.fastq.gz sra/SRR5656697_2.fastq.gz to data/BI0257_rbcL_R1_.fastq.gz data/BI0257_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656698_1.fastq.gz sra/SRR5656698_2.fastq.gz to data/IK4184_rbcL_R1_.fastq.gz data/IK4184_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656703_1.fastq.gz sra/SRR5656703_2.fastq.gz to 16s/data/IK3276_16s_R1_.fastq.gz 16s/data/IK3276_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656704_1.fastq.gz sra/SRR5656704_2.fastq.gz to 16s/data/IK3158_16s_R1_.fastq.gz 16s/data/IK3158_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656705_1.fastq.gz sra/SRR5656705_2.fastq.gz to data/TL3948_matK_R1_.fastq.gz data/TL3948_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656706_1.fastq.gz sra/SRR5656706_2.fastq.gz to data/TL3948_rbcL_R1_.fastq.gz data/TL3948_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656707_1.fastq.gz sra/SRR5656707_2.fastq.gz to 16s/data/BI0054_16s_R1_.fastq.gz 16s/data/BI0054_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656708_1.fastq.gz sra/SRR5656708_2.fastq.gz to 16s/data/TL3948_16s_R1_.fastq.gz 16s/data/TL3948_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656709_1.fastq.gz sra/SRR5656709_2.fastq.gz to 16s/data/IK3469_16s_R1_.fastq.gz 16s/data/IK3469_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656431_1.fastq.gz sra/SRR5656431_2.fastq.gz to 16s/data/PA1038_16s_R1_.fastq.gz 16s/data/PA1038_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656434_1.fastq.gz sra/SRR5656434_2.fastq.gz to 16s/data/PA1059_16s_R1_.fastq.gz 16s/data/PA1059_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656440_1.fastq.gz sra/SRR5656440_2.fastq.gz to data/LK682_rbcL_R1_.fastq.gz data/LK682_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656441_1.fastq.gz sra/SRR5656441_2.fastq.gz to data/LK682_matK_R1_.fastq.gz data/LK682_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656442_1.fastq.gz sra/SRR5656442_2.fastq.gz to data/BI0097_matK_R1_.fastq.gz data/BI0097_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656451_1.fastq.gz sra/SRR5656451_2.fastq.gz to 16s/data/KR57_16s_R1_.fastq.gz 16s/data/KR57_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656462_1.fastq.gz sra/SRR5656462_2.fastq.gz to data/BI0055_matK_R1_.fastq.gz data/BI0055_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656463_1.fastq.gz sra/SRR5656463_2.fastq.gz to data/TL3842_rbcL_R1_.fastq.gz data/TL3842_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656464_1.fastq.gz sra/SRR5656464_2.fastq.gz to data/TL3842_matK_R1_.fastq.gz data/TL3842_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656465_1.fastq.gz sra/SRR5656465_2.fastq.gz to data/BI0055_rbcL_R1_.fastq.gz data/BI0055_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656483_1.fastq.gz sra/SRR5656483_2.fastq.gz to 16s/data/IK3513_16s_R1_.fastq.gz 16s/data/IK3513_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656488_1.fastq.gz sra/SRR5656488_2.fastq.gz to data/TL3889_matK_R1_.fastq.gz data/TL3889_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656489_1.fastq.gz sra/SRR5656489_2.fastq.gz to data/TL3889_rbcL_R1_.fastq.gz data/TL3889_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656490_1.fastq.gz sra/SRR5656490_2.fastq.gz to data/TL3856_matK_R1_.fastq.gz data/TL3856_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656491_1.fastq.gz sra/SRR5656491_2.fastq.gz to data/TL3856_rbcL_R1_.fastq.gz data/TL3856_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656492_1.fastq.gz sra/SRR5656492_2.fastq.gz to data/TL3862_matK_R1_.fastq.gz data/TL3862_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656493_1.fastq.gz sra/SRR5656493_2.fastq.gz to data/TL3862_rbcL_R1_.fastq.gz data/TL3862_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656497_1.fastq.gz sra/SRR5656497_2.fastq.gz to 16s/data/PA1065_16s_R1_.fastq.gz 16s/data/PA1065_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656498_1.fastq.gz sra/SRR5656498_2.fastq.gz to data/TL3905_matK_R1_.fastq.gz data/TL3905_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656499_1.fastq.gz sra/SRR5656499_2.fastq.gz to data/TL3905_rbcL_R1_.fastq.gz data/TL3905_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656500_1.fastq.gz sra/SRR5656500_2.fastq.gz to 16s/data/PA1044_16s_R1_.fastq.gz 16s/data/PA1044_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656501_1.fastq.gz sra/SRR5656501_2.fastq.gz to 16s/data/PA1039_16s_R1_.fastq.gz 16s/data/PA1039_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656512_1.fastq.gz sra/SRR5656512_2.fastq.gz to data/PA1038_rbcL_R1_.fastq.gz data/PA1038_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656513_1.fastq.gz sra/SRR5656513_2.fastq.gz to data/PA0367_matK_R1_.fastq.gz data/PA0367_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656524_1.fastq.gz sra/SRR5656524_2.fastq.gz to data/TL3910_rbcL_R1_.fastq.gz data/TL3910_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656525_1.fastq.gz sra/SRR5656525_2.fastq.gz to data/TL3910_matK_R1_.fastq.gz data/TL3910_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656528_1.fastq.gz sra/SRR5656528_2.fastq.gz to data/TL3915_rbcL_R1_.fastq.gz data/TL3915_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656529_1.fastq.gz sra/SRR5656529_2.fastq.gz to data/TL3915_matK_R1_.fastq.gz data/TL3915_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656535_1.fastq.gz sra/SRR5656535_2.fastq.gz to data/PA0367_rbcL_R1_.fastq.gz data/PA0367_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656544_1.fastq.gz sra/SRR5656544_2.fastq.gz to 16s/data/UB1435_16s_R1_.fastq.gz 16s/data/UB1435_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656545_1.fastq.gz sra/SRR5656545_2.fastq.gz to data/KR57_rbcL_R1_.fastq.gz data/KR57_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656546_1.fastq.gz sra/SRR5656546_2.fastq.gz to data/LK685_matK_R1_.fastq.gz data/LK685_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656548_1.fastq.gz sra/SRR5656548_2.fastq.gz to data/PA0368_matK_R1_.fastq.gz data/PA0368_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656566_1.fastq.gz sra/SRR5656566_2.fastq.gz to 16s/data/PA1049_16s_R1_.fastq.gz 16s/data/PA1049_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656568_1.fastq.gz sra/SRR5656568_2.fastq.gz to data/PA0368_rbcL_R1_.fastq.gz data/PA0368_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656579_1.fastq.gz sra/SRR5656579_2.fastq.gz to data/LK685_rbcL_R1_.fastq.gz data/LK685_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656601_1.fastq.gz sra/SRR5656601_2.fastq.gz to data/PA0370_matK_R1_.fastq.gz data/PA0370_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656602_1.fastq.gz sra/SRR5656602_2.fastq.gz to data/PA0456_rbcL_R1_.fastq.gz data/PA0456_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656603_1.fastq.gz sra/SRR5656603_2.fastq.gz to data/PA0370_rbcL_R1_.fastq.gz data/PA0370_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656612_1.fastq.gz sra/SRR5656612_2.fastq.gz to 16s/data/PA0456_16s_R1_.fastq.gz 16s/data/PA0456_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656619_1.fastq.gz sra/SRR5656619_2.fastq.gz to 16s/data/BI0097_16s_R1_.fastq.gz 16s/data/BI0097_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656629_1.fastq.gz sra/SRR5656629_2.fastq.gz to 16s/data/TL3862_16s_R1_.fastq.gz 16s/data/TL3862_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656630_1.fastq.gz sra/SRR5656630_2.fastq.gz to 16s/data/TL3856_16s_R1_.fastq.gz 16s/data/TL3856_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656631_1.fastq.gz sra/SRR5656631_2.fastq.gz to 16s/data/TL3842_16s_R1_.fastq.gz 16s/data/TL3842_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656632_1.fastq.gz sra/SRR5656632_2.fastq.gz to 16s/data/LK685_16s_R1_.fastq.gz 16s/data/LK685_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656633_1.fastq.gz sra/SRR5656633_2.fastq.gz to 16s/data/TL3910_16s_R1_.fastq.gz 16s/data/TL3910_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656634_1.fastq.gz sra/SRR5656634_2.fastq.gz to 16s/data/TL3905_16s_R1_.fastq.gz 16s/data/TL3905_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656635_1.fastq.gz sra/SRR5656635_2.fastq.gz to 16s/data/TL3889_16s_R1_.fastq.gz 16s/data/TL3889_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656636_1.fastq.gz sra/SRR5656636_2.fastq.gz to 16s/data/LK682_16s_R1_.fastq.gz 16s/data/LK682_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656637_1.fastq.gz sra/SRR5656637_2.fastq.gz to 16s/data/TL3915_16s_R1_.fastq.gz 16s/data/TL3915_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656638_1.fastq.gz sra/SRR5656638_2.fastq.gz to 16s/data/BI0093_16s_R1_.fastq.gz 16s/data/BI0093_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656639_1.fastq.gz sra/SRR5656639_2.fastq.gz to 16s/data/BI0055_16s_R1_.fastq.gz 16s/data/BI0055_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656640_1.fastq.gz sra/SRR5656640_2.fastq.gz to data/PA0456_matK_R1_.fastq.gz data/PA0456_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656641_1.fastq.gz sra/SRR5656641_2.fastq.gz to data/PA1038_matK_R1_.fastq.gz data/PA1038_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656650_1.fastq.gz sra/SRR5656650_2.fastq.gz to data/IK3513_rbcL_R1_.fastq.gz data/IK3513_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656651_1.fastq.gz sra/SRR5656651_2.fastq.gz to data/IK3513_matK_R1_.fastq.gz data/IK3513_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656653_1.fastq.gz sra/SRR5656653_2.fastq.gz to data/UB1452_matK_R1_.fastq.gz data/UB1452_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656656_1.fastq.gz sra/SRR5656656_2.fastq.gz to data/UB1435_rbcL_R1_.fastq.gz data/UB1435_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656657_1.fastq.gz sra/SRR5656657_2.fastq.gz to data/UB1435_matK_R1_.fastq.gz data/UB1435_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656658_1.fastq.gz sra/SRR5656658_2.fastq.gz to data/UB1446_rbcL_R1_.fastq.gz data/UB1446_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656659_1.fastq.gz sra/SRR5656659_2.fastq.gz to data/UB1446_matK_R1_.fastq.gz data/UB1446_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656660_1.fastq.gz sra/SRR5656660_2.fastq.gz to 16s/data/PA0368_16s_R1_.fastq.gz 16s/data/PA0368_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656661_1.fastq.gz sra/SRR5656661_2.fastq.gz to 16s/data/PA0370_16s_R1_.fastq.gz 16s/data/PA0370_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656663_1.fastq.gz sra/SRR5656663_2.fastq.gz to 16s/data/PA0367_16s_R1_.fastq.gz 16s/data/PA0367_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656668_1.fastq.gz sra/SRR5656668_2.fastq.gz to data/KR57_matK_R1_.fastq.gz data/KR57_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656669_1.fastq.gz sra/SRR5656669_2.fastq.gz to data/BI0097_rbcL_R1_.fastq.gz data/BI0097_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656670_1.fastq.gz sra/SRR5656670_2.fastq.gz to data/BI0093_matK_R1_.fastq.gz data/BI0093_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656671_1.fastq.gz sra/SRR5656671_2.fastq.gz to data/BI0093_rbcL_R1_.fastq.gz data/BI0093_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656672_1.fastq.gz sra/SRR5656672_2.fastq.gz to data/PA1065_rbcL_R1_.fastq.gz data/PA1065_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656673_1.fastq.gz sra/SRR5656673_2.fastq.gz to data/PA1065_matK_R1_.fastq.gz data/PA1065_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656676_1.fastq.gz sra/SRR5656676_2.fastq.gz to data/PA1049_rbcL_R1_.fastq.gz data/PA1049_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656677_1.fastq.gz sra/SRR5656677_2.fastq.gz to data/PA1049_matK_R1_.fastq.gz data/PA1049_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656678_1.fastq.gz sra/SRR5656678_2.fastq.gz to data/PA1059_rbcL_R1_.fastq.gz data/PA1059_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656679_1.fastq.gz sra/SRR5656679_2.fastq.gz to data/PA1059_matK_R1_.fastq.gz data/PA1059_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656680_1.fastq.gz sra/SRR5656680_2.fastq.gz to data/PA1039_rbcL_R1_.fastq.gz data/PA1039_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656681_1.fastq.gz sra/SRR5656681_2.fastq.gz to data/PA1039_matK_R1_.fastq.gz data/PA1039_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656682_1.fastq.gz sra/SRR5656682_2.fastq.gz to data/PA1044_rbcL_R1_.fastq.gz data/PA1044_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656683_1.fastq.gz sra/SRR5656683_2.fastq.gz to data/PA1044_matK_R1_.fastq.gz data/PA1044_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656684_1.fastq.gz sra/SRR5656684_2.fastq.gz to 16s/data/UB1446_16s_R1_.fastq.gz 16s/data/UB1446_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656685_1.fastq.gz sra/SRR5656685_2.fastq.gz to 16s/data/UB1452_16s_R1_.fastq.gz 16s/data/UB1452_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656435_1.fastq.gz sra/SRR5656435_2.fastq.gz to data/BI0246_rbcL_R1_.fastq.gz data/BI0246_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656436_1.fastq.gz sra/SRR5656436_2.fastq.gz to data/BI0246_matK_R1_.fastq.gz data/BI0246_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656437_1.fastq.gz sra/SRR5656437_2.fastq.gz to 16s/data/TL3916_16s_R1_.fastq.gz 16s/data/TL3916_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656438_1.fastq.gz sra/SRR5656438_2.fastq.gz to 16s/data/LG4314_16s_R1_.fastq.gz 16s/data/LG4314_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656439_1.fastq.gz sra/SRR5656439_2.fastq.gz to 16s/data/LG4322_16s_R1_.fastq.gz 16s/data/LG4322_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656443_1.fastq.gz sra/SRR5656443_2.fastq.gz to data/LK670_matK_R1_.fastq.gz data/LK670_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656444_1.fastq.gz sra/SRR5656444_2.fastq.gz to data/LK668_rbcL_R1_.fastq.gz data/LK668_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656445_1.fastq.gz sra/SRR5656445_2.fastq.gz to data/LK668_matK_R1_.fastq.gz data/LK668_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656446_1.fastq.gz sra/SRR5656446_2.fastq.gz to 16s/data/KR12_16s_R1_.fastq.gz 16s/data/KR12_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656447_1.fastq.gz sra/SRR5656447_2.fastq.gz to 16s/data/KR21_16s_R1_.fastq.gz 16s/data/KR21_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656448_1.fastq.gz sra/SRR5656448_2.fastq.gz to 16s/data/KR33_16s_R1_.fastq.gz 16s/data/KR33_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656449_1.fastq.gz sra/SRR5656449_2.fastq.gz to 16s/data/KR35_16s_R1_.fastq.gz 16s/data/KR35_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656450_1.fastq.gz sra/SRR5656450_2.fastq.gz to 16s/data/KR52_16s_R1_.fastq.gz 16s/data/KR52_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656452_1.fastq.gz sra/SRR5656452_2.fastq.gz to 16s/data/KR67_16s_R1_.fastq.gz 16s/data/KR67_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656453_1.fastq.gz sra/SRR5656453_2.fastq.gz to 16s/data/LG4300_16s_R1_.fastq.gz 16s/data/LG4300_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656454_1.fastq.gz sra/SRR5656454_2.fastq.gz to data/TL3838_rbcL_R1_.fastq.gz data/TL3838_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656455_1.fastq.gz sra/SRR5656455_2.fastq.gz to data/TL3838_matK_R1_.fastq.gz data/TL3838_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656456_1.fastq.gz sra/SRR5656456_2.fastq.gz to data/TL3826_rbcL_R1_.fastq.gz data/TL3826_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656457_1.fastq.gz sra/SRR5656457_2.fastq.gz to data/TL3826_matK_R1_.fastq.gz data/TL3826_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656458_1.fastq.gz sra/SRR5656458_2.fastq.gz to data/TL3824_rbcL_R1_.fastq.gz data/TL3824_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656459_1.fastq.gz sra/SRR5656459_2.fastq.gz to data/TL3824_matK_R1_.fastq.gz data/TL3824_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656460_1.fastq.gz sra/SRR5656460_2.fastq.gz to data/TL3821_rbcL_R1_.fastq.gz data/TL3821_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656461_1.fastq.gz sra/SRR5656461_2.fastq.gz to data/TL3821_matK_R1_.fastq.gz data/TL3821_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656466_1.fastq.gz sra/SRR5656466_2.fastq.gz to data/LK665_rbcL_R1_.fastq.gz data/LK665_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656467_1.fastq.gz sra/SRR5656467_2.fastq.gz to data/LK647_matK_R1_.fastq.gz data/LK647_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656468_1.fastq.gz sra/SRR5656468_2.fastq.gz to data/LK647_rbcL_R1_.fastq.gz data/LK647_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656469_1.fastq.gz sra/SRR5656469_2.fastq.gz to data/LK645_matK_R1_.fastq.gz data/LK645_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656470_1.fastq.gz sra/SRR5656470_2.fastq.gz to data/LK645_rbcL_R1_.fastq.gz data/LK645_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656471_1.fastq.gz sra/SRR5656471_2.fastq.gz to data/LK661_matK_R1_.fastq.gz data/LK661_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656472_1.fastq.gz sra/SRR5656472_2.fastq.gz to data/LK661_rbcL_R1_.fastq.gz data/LK661_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656473_1.fastq.gz sra/SRR5656473_2.fastq.gz to data/LK653_matK_R1_.fastq.gz data/LK653_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656474_1.fastq.gz sra/SRR5656474_2.fastq.gz to data/LK653_rbcL_R1_.fastq.gz data/LK653_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656475_1.fastq.gz sra/SRR5656475_2.fastq.gz to data/LK665_matK_R1_.fastq.gz data/LK665_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656476_1.fastq.gz sra/SRR5656476_2.fastq.gz to 16s/data/KR05_16s_R1_.fastq.gz 16s/data/KR05_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656477_1.fastq.gz sra/SRR5656477_2.fastq.gz to 16s/data/KR02_16s_R1_.fastq.gz 16s/data/KR02_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656478_1.fastq.gz sra/SRR5656478_2.fastq.gz to 16s/data/IK4214_16s_R1_.fastq.gz 16s/data/IK4214_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656480_1.fastq.gz sra/SRR5656480_2.fastq.gz to 16s/data/IK3777_16s_R1_.fastq.gz 16s/data/IK3777_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656481_1.fastq.gz sra/SRR5656481_2.fastq.gz to 16s/data/IK3701_16s_R1_.fastq.gz 16s/data/IK3701_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656482_1.fastq.gz sra/SRR5656482_2.fastq.gz to 16s/data/IK3650_16s_R1_.fastq.gz 16s/data/IK3650_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656484_1.fastq.gz sra/SRR5656484_2.fastq.gz to 16s/data/KR10_16s_R1_.fastq.gz 16s/data/KR10_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656485_1.fastq.gz sra/SRR5656485_2.fastq.gz to 16s/data/KR07_16s_R1_.fastq.gz 16s/data/KR07_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656486_1.fastq.gz sra/SRR5656486_2.fastq.gz to data/TL3882_matK_R1_.fastq.gz data/TL3882_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656487_1.fastq.gz sra/SRR5656487_2.fastq.gz to data/TL3882_rbcL_R1_.fastq.gz data/TL3882_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656494_1.fastq.gz sra/SRR5656494_2.fastq.gz to 16s/data/UB0599_16s_R1_.fastq.gz 16s/data/UB0599_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656495_1.fastq.gz sra/SRR5656495_2.fastq.gz to 16s/data/UB0445_16s_R1_.fastq.gz 16s/data/UB0445_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656496_1.fastq.gz sra/SRR5656496_2.fastq.gz to 16s/data/UB0439_16s_R1_.fastq.gz 16s/data/UB0439_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656502_1.fastq.gz sra/SRR5656502_2.fastq.gz to data/LG4314_rbcL_R1_.fastq.gz data/LG4314_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656503_1.fastq.gz sra/SRR5656503_2.fastq.gz to data/LG4314_matK_R1_.fastq.gz data/LG4314_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656504_1.fastq.gz sra/SRR5656504_2.fastq.gz to data/LG4322_rbcL_R1_.fastq.gz data/LG4322_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656505_1.fastq.gz sra/SRR5656505_2.fastq.gz to data/LG4322_matK_R1_.fastq.gz data/LG4322_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656506_1.fastq.gz sra/SRR5656506_2.fastq.gz to data/KR67_rbcL_R1_.fastq.gz data/KR67_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656507_1.fastq.gz sra/SRR5656507_2.fastq.gz to data/KR67_matK_R1_.fastq.gz data/KR67_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656508_1.fastq.gz sra/SRR5656508_2.fastq.gz to data/LG4300_rbcL_R1_.fastq.gz data/LG4300_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656509_1.fastq.gz sra/SRR5656509_2.fastq.gz to data/LG4300_matK_R1_.fastq.gz data/LG4300_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656510_1.fastq.gz sra/SRR5656510_2.fastq.gz to data/LG4327_rbcL_R1_.fastq.gz data/LG4327_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656511_1.fastq.gz sra/SRR5656511_2.fastq.gz to data/LG4327_matK_R1_.fastq.gz data/LG4327_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656516_1.fastq.gz sra/SRR5656516_2.fastq.gz to 16s/data/LK686_16s_R1_.fastq.gz 16s/data/LK686_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656518_1.fastq.gz sra/SRR5656518_2.fastq.gz to 16s/data/TL3821_16s_R1_.fastq.gz 16s/data/TL3821_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656519_1.fastq.gz sra/SRR5656519_2.fastq.gz to 16s/data/TL3824_16s_R1_.fastq.gz 16s/data/TL3824_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656522_1.fastq.gz sra/SRR5656522_2.fastq.gz to 16s/data/TL3826_16s_R1_.fastq.gz 16s/data/TL3826_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656523_1.fastq.gz sra/SRR5656523_2.fastq.gz to 16s/data/TL3838_16s_R1_.fastq.gz 16s/data/TL3838_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656526_1.fastq.gz sra/SRR5656526_2.fastq.gz to data/TL3911_rbcL_R1_.fastq.gz data/TL3911_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656527_1.fastq.gz sra/SRR5656527_2.fastq.gz to data/TL3911_matK_R1_.fastq.gz data/TL3911_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656530_1.fastq.gz sra/SRR5656530_2.fastq.gz to data/TL3916_rbcL_R1_.fastq.gz data/TL3916_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656531_1.fastq.gz sra/SRR5656531_2.fastq.gz to data/TL3916_matK_R1_.fastq.gz data/TL3916_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656532_1.fastq.gz sra/SRR5656532_2.fastq.gz to data/TL3918_rbcL_R1_.fastq.gz data/TL3918_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656533_1.fastq.gz sra/SRR5656533_2.fastq.gz to data/TL3918_matK_R1_.fastq.gz data/TL3918_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656534_1.fastq.gz sra/SRR5656534_2.fastq.gz to 16s/data/TL3911_16s_R1_.fastq.gz 16s/data/TL3911_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656536_1.fastq.gz sra/SRR5656536_2.fastq.gz to data/KR21_matK_R1_.fastq.gz data/KR21_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656537_1.fastq.gz sra/SRR5656537_2.fastq.gz to data/KR21_rbcL_R1_.fastq.gz data/KR21_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656538_1.fastq.gz sra/SRR5656538_2.fastq.gz to data/KR33_matK_R1_.fastq.gz data/KR33_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656539_1.fastq.gz sra/SRR5656539_2.fastq.gz to data/KR33_rbcL_R1_.fastq.gz data/KR33_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656540_1.fastq.gz sra/SRR5656540_2.fastq.gz to data/KR35_matK_R1_.fastq.gz data/KR35_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656541_1.fastq.gz sra/SRR5656541_2.fastq.gz to data/KR35_rbcL_R1_.fastq.gz data/KR35_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656542_1.fastq.gz sra/SRR5656542_2.fastq.gz to data/KR52_matK_R1_.fastq.gz data/KR52_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656543_1.fastq.gz sra/SRR5656543_2.fastq.gz to data/KR52_rbcL_R1_.fastq.gz data/KR52_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656547_1.fastq.gz sra/SRR5656547_2.fastq.gz to 16s/data/UB1430_16s_R1_.fastq.gz 16s/data/UB1430_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656549_1.fastq.gz sra/SRR5656549_2.fastq.gz to data/LK670_rbcL_R1_.fastq.gz data/LK670_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656550_1.fastq.gz sra/SRR5656550_2.fastq.gz to 16s/data/LK645_16s_R1_.fastq.gz 16s/data/LK645_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656551_1.fastq.gz sra/SRR5656551_2.fastq.gz to 16s/data/LG4327_16s_R1_.fastq.gz 16s/data/LG4327_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656552_1.fastq.gz sra/SRR5656552_2.fastq.gz to 16s/data/LK653_16s_R1_.fastq.gz 16s/data/LK653_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656553_1.fastq.gz sra/SRR5656553_2.fastq.gz to 16s/data/LK647_16s_R1_.fastq.gz 16s/data/LK647_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656554_1.fastq.gz sra/SRR5656554_2.fastq.gz to 16s/data/LK665_16s_R1_.fastq.gz 16s/data/LK665_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656555_1.fastq.gz sra/SRR5656555_2.fastq.gz to 16s/data/LK661_16s_R1_.fastq.gz 16s/data/LK661_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656556_1.fastq.gz sra/SRR5656556_2.fastq.gz to data/TL3932_matK_R1_.fastq.gz data/TL3932_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656557_1.fastq.gz sra/SRR5656557_2.fastq.gz to data/TL3932_rbcL_R1_.fastq.gz data/TL3932_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656558_1.fastq.gz sra/SRR5656558_2.fastq.gz to data/TL3929_matK_R1_.fastq.gz data/TL3929_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656559_1.fastq.gz sra/SRR5656559_2.fastq.gz to data/TL3929_rbcL_R1_.fastq.gz data/TL3929_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656560_1.fastq.gz sra/SRR5656560_2.fastq.gz to data/TL3927_matK_R1_.fastq.gz data/TL3927_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656561_1.fastq.gz sra/SRR5656561_2.fastq.gz to data/TL3927_rbcL_R1_.fastq.gz data/TL3927_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656562_1.fastq.gz sra/SRR5656562_2.fastq.gz to data/TL3926_matK_R1_.fastq.gz data/TL3926_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656563_1.fastq.gz sra/SRR5656563_2.fastq.gz to data/TL3926_rbcL_R1_.fastq.gz data/TL3926_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656564_1.fastq.gz sra/SRR5656564_2.fastq.gz to data/TL3925_matK_R1_.fastq.gz data/TL3925_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656565_1.fastq.gz sra/SRR5656565_2.fastq.gz to data/TL3925_rbcL_R1_.fastq.gz data/TL3925_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656567_1.fastq.gz sra/SRR5656567_2.fastq.gz to 16s/data/TL3939_16s_R1_.fastq.gz 16s/data/TL3939_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656569_1.fastq.gz sra/SRR5656569_2.fastq.gz to data/KR05_rbcL_R1_.fastq.gz data/KR05_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656570_1.fastq.gz sra/SRR5656570_2.fastq.gz to data/KR05_matK_R1_.fastq.gz data/KR05_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656571_1.fastq.gz sra/SRR5656571_2.fastq.gz to data/KR02_rbcL_R1_.fastq.gz data/KR02_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656572_1.fastq.gz sra/SRR5656572_2.fastq.gz to data/KR02_matK_R1_.fastq.gz data/KR02_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656573_1.fastq.gz sra/SRR5656573_2.fastq.gz to data/KR10_rbcL_R1_.fastq.gz data/KR10_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656574_1.fastq.gz sra/SRR5656574_2.fastq.gz to data/KR10_matK_R1_.fastq.gz data/KR10_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656575_1.fastq.gz sra/SRR5656575_2.fastq.gz to data/KR07_rbcL_R1_.fastq.gz data/KR07_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656576_1.fastq.gz sra/SRR5656576_2.fastq.gz to data/KR07_matK_R1_.fastq.gz data/KR07_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656577_1.fastq.gz sra/SRR5656577_2.fastq.gz to data/KR12_rbcL_R1_.fastq.gz data/KR12_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656578_1.fastq.gz sra/SRR5656578_2.fastq.gz to data/KR12_matK_R1_.fastq.gz data/KR12_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656580_1.fastq.gz sra/SRR5656580_2.fastq.gz to data/LK686_matK_R1_.fastq.gz data/LK686_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656583_1.fastq.gz sra/SRR5656583_2.fastq.gz to data/TL3939_rbcL_R1_.fastq.gz data/TL3939_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656584_1.fastq.gz sra/SRR5656584_2.fastq.gz to data/TL3939_matK_R1_.fastq.gz data/TL3939_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656585_1.fastq.gz sra/SRR5656585_2.fastq.gz to data/TL3936_rbcL_R1_.fastq.gz data/TL3936_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656586_1.fastq.gz sra/SRR5656586_2.fastq.gz to data/TL3936_matK_R1_.fastq.gz data/TL3936_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656587_1.fastq.gz sra/SRR5656587_2.fastq.gz to data/TL3942_rbcL_R1_.fastq.gz data/TL3942_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656588_1.fastq.gz sra/SRR5656588_2.fastq.gz to data/TL3942_matK_R1_.fastq.gz data/TL3942_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656589_1.fastq.gz sra/SRR5656589_2.fastq.gz to data/TL3940_rbcL_R1_.fastq.gz data/TL3940_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656590_1.fastq.gz sra/SRR5656590_2.fastq.gz to data/TL3940_matK_R1_.fastq.gz data/TL3940_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656591_1.fastq.gz sra/SRR5656591_2.fastq.gz to 16s/data/TL3929_16s_R1_.fastq.gz 16s/data/TL3929_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656592_1.fastq.gz sra/SRR5656592_2.fastq.gz to 16s/data/TL3932_16s_R1_.fastq.gz 16s/data/TL3932_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656593_1.fastq.gz sra/SRR5656593_2.fastq.gz to 16s/data/TL3936_16s_R1_.fastq.gz 16s/data/TL3936_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656594_1.fastq.gz sra/SRR5656594_2.fastq.gz to 16s/data/LK670_16s_R1_.fastq.gz 16s/data/LK670_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656595_1.fastq.gz sra/SRR5656595_2.fastq.gz to 16s/data/TL3918_16s_R1_.fastq.gz 16s/data/TL3918_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656596_1.fastq.gz sra/SRR5656596_2.fastq.gz to 16s/data/TL3925_16s_R1_.fastq.gz 16s/data/TL3925_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656597_1.fastq.gz sra/SRR5656597_2.fastq.gz to 16s/data/TL3926_16s_R1_.fastq.gz 16s/data/TL3926_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656598_1.fastq.gz sra/SRR5656598_2.fastq.gz to 16s/data/TL3927_16s_R1_.fastq.gz 16s/data/TL3927_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656599_1.fastq.gz sra/SRR5656599_2.fastq.gz to 16s/data/TL3940_16s_R1_.fastq.gz 16s/data/TL3940_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656600_1.fastq.gz sra/SRR5656600_2.fastq.gz to 16s/data/TL3942_16s_R1_.fastq.gz 16s/data/TL3942_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656605_1.fastq.gz sra/SRR5656605_2.fastq.gz to 16s/data/LK668_16s_R1_.fastq.gz 16s/data/LK668_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656606_1.fastq.gz sra/SRR5656606_2.fastq.gz to data/IK3777_matK_R1_.fastq.gz data/IK3777_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656607_1.fastq.gz sra/SRR5656607_2.fastq.gz to data/IK3777_rbcL_R1_.fastq.gz data/IK3777_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656608_1.fastq.gz sra/SRR5656608_2.fastq.gz to data/IK3701_matK_R1_.fastq.gz data/IK3701_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656609_1.fastq.gz sra/SRR5656609_2.fastq.gz to data/IK3701_rbcL_R1_.fastq.gz data/IK3701_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656610_1.fastq.gz sra/SRR5656610_2.fastq.gz to data/IK3650_matK_R1_.fastq.gz data/IK3650_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656611_1.fastq.gz sra/SRR5656611_2.fastq.gz to data/IK3650_rbcL_R1_.fastq.gz data/IK3650_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656613_1.fastq.gz sra/SRR5656613_2.fastq.gz to 16s/data/BI0246_16s_R1_.fastq.gz 16s/data/BI0246_16s_R2_.fastq.gz
```

```
## Moving sra/SRR5656614_1.fastq.gz sra/SRR5656614_2.fastq.gz to data/LK686_rbcL_R1_.fastq.gz data/LK686_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656615_1.fastq.gz sra/SRR5656615_2.fastq.gz to data/IK4214_matK_R1_.fastq.gz data/IK4214_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656616_1.fastq.gz sra/SRR5656616_2.fastq.gz to data/IK4214_rbcL_R1_.fastq.gz data/IK4214_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656617_1.fastq.gz sra/SRR5656617_2.fastq.gz to data/UB1430_matK_R1_.fastq.gz data/UB1430_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656618_1.fastq.gz sra/SRR5656618_2.fastq.gz to data/UB1430_rbcL_R1_.fastq.gz data/UB1430_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656622_1.fastq.gz sra/SRR5656622_2.fastq.gz to data/UB0439_matK_R1_.fastq.gz data/UB0439_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656623_1.fastq.gz sra/SRR5656623_2.fastq.gz to data/UB0439_rbcL_R1_.fastq.gz data/UB0439_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656624_1.fastq.gz sra/SRR5656624_2.fastq.gz to data/UB0445_matK_R1_.fastq.gz data/UB0445_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656625_1.fastq.gz sra/SRR5656625_2.fastq.gz to data/UB0445_rbcL_R1_.fastq.gz data/UB0445_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656626_1.fastq.gz sra/SRR5656626_2.fastq.gz to data/UB0599_matK_R1_.fastq.gz data/UB0599_matK_R2_.fastq.gz
```

```
## Moving sra/SRR5656627_1.fastq.gz sra/SRR5656627_2.fastq.gz to data/UB0599_rbcL_R1_.fastq.gz data/UB0599_rbcL_R2_.fastq.gz
```

```
## Moving sra/SRR5656628_1.fastq.gz sra/SRR5656628_2.fastq.gz to 16s/data/TL3882_16s_R1_.fastq.gz 16s/data/TL3882_16s_R2_.fastq.gz
```



