## PCoA and t-SNE of matK and rbcL abundances


```r
#set seed so reproducible
set.seed(12357)
#stop on errors
knitr::opts_chunk$set(error=FALSE,tidy=TRUE)
```

### Load libraries

```r
library(phyloseq)
packageVersion("phyloseq")
```

```
## [1] '1.20.0'
```

```r
library(ape)
packageVersion("ape")
```

```
## [1] '4.1'
```

```r
library(Rtsne)
packageVersion("Rtsne")
```

```
## [1] '0.13'
```

```r
library(vegan)
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.4-3
```

```r
packageVersion("vegan")
```

```
## [1] '2.4.3'
```

```r
library(ade4)
```

```
## 
## Attaching package: 'ade4'
```

```
## The following object is masked from 'package:vegan':
## 
##     cca
```

```r
packageVersion("ade4")
```

```
## [1] '1.7.6'
```

```r
source("functions.R")
```

### Load data

```r
source("loadData.R")
```

```
## matK
```

```
## Cache work/matK_rarefyOtus.Rdat does exist. Loading data
```

```
## rbcL
```

```
## Cache work/rbcL_rarefyOtus.Rdat does exist. Loading data
```

### Plot plant PCoA and t-SNE

```r
tlAdonis <- interactAdonis <- plantAdonis <- chimpAdonis <- bonoboAdonis <- list()
mantels <- list()
for (ii in names(swarmData)) {
    plotProp <- swarmData[[ii]][["props"]][swarmData[[ii]][["isEnough"]] & rownames(swarmData[[ii]][["props"]]) %in% 
        rownames(samples), ]
    plotProp2 <- swarmData[[ii]][["rare"]][swarmData[[ii]][["isEnough"]] & rownames(swarmData[[ii]][["rare"]]) %in% 
        rownames(samples), ]
    phyOtuW <- otu_table(plotProp, taxa_are_rows = FALSE)
    phyOtuU <- otu_table(plotProp2, taxa_are_rows = FALSE)
    qiimeDataW <- phyloseq(otu_table = phyOtuW, phy_tree = swarmData[[ii]][["tree"]])
    qiimeDataU <- phyloseq(otu_table = phyOtuU, phy_tree = swarmData[[ii]][["tree"]])
    brayDist <- distance(qiimeDataU, "bray", binary = TRUE)
    brayDistW <- distance(qiimeDataW, "bray", binary = FALSE)
    uniDist <- UniFrac(qiimeDataU, weighted = FALSE)
    uniDistW <- UniFrac(qiimeDataW, weighted = TRUE)
    mantels[[ii]] <- list(uniW = ade4::mantel.rtest(uniDist, uniDistW, nrepet = 10000), 
        brayW = ade4::mantel.rtest(uniDist, brayDistW, nrepet = 10000), brayUW = ade4::mantel.rtest(uniDist, 
            brayDist, nrepet = 10000))
    uniPca <- pcoa(uniDist)
    tsneUni <- Rtsne(uniDist, is_distance = TRUE, verbose = TRUE, perplexity = 15, 
        max_iter = 3000)
    selectDist <- uniDist
    
    importance <- uniPca$values$Relative_eig
    colnames(uniPca$vectors) <- sprintf("Principal coordinate %d (%d%% of variance)", 
        1:length(importance), round(importance * 100))[1:ncol(uniPca$vectors)]
    selectSamples <- samples[rownames(plotProp), ]
    selectSamples <- selectSamples[order(selectSamples$bonobo, selectSamples$area2, 
        selectSamples$malaria), ]
    colorBrew <- c("#e41a1cBB", "#377eb8BB", "#4daf4aBB", "#984ea3BB", "#ff7f00BB", 
        "#ffff33BB", "#a65628BB", "#f781bfBB", "#999999BB", "#88ddffBB")
    nArea <- length(unique(selectSamples$area2))
    if (nArea > length(colorBrew)) 
        stop("Need to adjust colors for more areas")
    areaCols <- colorBrew[1:nArea]
    names(areaCols) <- unique(selectSamples$area2[order(selectSamples$chimpBonobo)])
    areaPch <- sapply(names(areaCols), function(x) mostAbundant(selectSamples$chimpBonobo[selectSamples$area2 == 
        x]))
    malariaCols3 <- c("#00000022", "#000000CC")
    mediumMalariaCol <- "#00000077"
    malariaCols <- c("#00000022", "#000000CC")
    malariaCols2 <- rainbow.lab(2, alpha = 0.9, lightMultiple = 0.7)
    names(malariaCols3) <- names(malariaCols2) <- names(malariaCols) <- c("Laverania negative", 
        "Laverania positive")
    speciesPch <- 20 + 1:length(unique(selectSamples$Species))
    speciesCols <- rainbow.lab(length(unique(selectSamples$Species)), start = -2, 
        end = 1, alpha = 0.8, lightMultiple = 0.8)
    names(speciesCols) <- names(speciesPch) <- sort(unique(selectSamples$chimpBonobo))
    
    predictors <- model.matrix(~0 + Species + malaria + SIV + area, selectSamples)
    
    # pcoa
    pos <- my.biplot.pcoa(uniPca, predictors, plot.axes = 1:2, pch = 21, bg = speciesCols[selectSamples$chimpBonobo], 
        col = malariaCols3[selectSamples$malaria + 1], cex = 2.25, lwd = 4, 
        arrowsFilter = Inf, las = 1, mgp = c(2.75, 0.75, 0), sameAxis = FALSE, 
        bty = "l", type = "n")
    points(pos[!selectSamples$malaria, ], col = malariaCols3[1], cex = 2.25, 
        lwd = 4, bg = speciesCols[selectSamples[!selectSamples$malaria, "chimpBonobo"]], 
        pch = 21)
    points(pos[selectSamples$malaria, ], col = malariaCols3[2], cex = 2.25, 
        lwd = 4, bg = speciesCols[selectSamples[selectSamples$malaria, "chimpBonobo"]], 
        pch = 21)
    title(main = sprintf("%s", ii, 1, 2))
    
    # tsne
    par(mar = c(4, 4, 1.5, 10))
    plot(tsneUni$Y, pch = speciesPch[selectSamples$chimpBonobo], bg = areaCols[selectSamples$area2], 
        col = malariaCols[selectSamples$malaria + 1], cex = 2.5, lwd = 3, ylab = "t-SNE 2", 
        xlab = "t-SNE 1", main = sprintf("%s", ii), bty = "l", las = 1)
    legend(par("usr")[2] + 0.01 * diff(par("usr")[1:2]), mean(par("usr")[3:4]), 
        c(names(malariaCols), names(areaCols), names(speciesPch)), col = c(malariaCols, 
            rep(c(malariaCols[1], mediumMalariaCol), c(length(areaCols), length(speciesPch)))), 
        pch = c(rep(21, length(malariaCols)), speciesPch[areaPch], speciesPch), 
        pt.bg = c(rep(NA, length(malariaCols)), areaCols, rep(NA, length(speciesPch))), 
        inset = 0.01, pt.lwd = 3, pt.cex = 2.5, xjust = 0, xpd = NA, bty = "n")
    
    ss <- samples[labels(selectDist), ]
    plantAdonis[[ii]] <- adonis(selectDist ~ bonobo + area2 + malaria, data = ss, 
        permutations = 1e+07, parallel = 10)
    interactAdonis[[ii]] <- adonis(selectDist ~ bonobo + area2 + malaria * bonobo, 
        data = ss, permutations = 1e+07, parallel = 10)
    chimpDist <- as.matrix(selectDist)
    chimpDist <- as.dist(chimpDist[!ss$bonobo, !ss$bonobo])
    chimpAdonis[[ii]] <- adonis(chimpDist ~ area2 + malaria, data = ss[!ss$bonobo, 
        ], permutations = 1e+07, parallel = 10)
    bonoboDist <- as.matrix(selectDist)
    bonoboDist <- as.dist(bonoboDist[ss$bonobo, ss$bonobo])
    bonoboAdonis[[ii]] <- adonis(bonoboDist ~ area2 + malaria, data = ss[ss$bonobo, 
        ], permutations = 1e+07, parallel = 10)
    tlDist <- as.dist(as.matrix(selectDist)[ss$isTL, ss$isTL])
    tlAdonis[[ii]] <- adonis(tlDist ~ malaria, data = ss[ss$isTL, ], permutations = 1e+07, 
        parallel = 5)
}
```

```
## Read the 92 x 92 data matrix successfully!
## Using no_dims = 2, perplexity = 15.000000, and theta = 0.500000
## Computing input similarities...
## Building tree...
##  - point 0 of 92
## Done in 0.01 seconds (sparsity = 0.612949)!
## Learning embedding...
## Iteration 50: error is 59.308155 (50 iterations in 0.02 seconds)
## Iteration 100: error is 58.593163 (50 iterations in 0.02 seconds)
## Iteration 150: error is 60.708217 (50 iterations in 0.02 seconds)
## Iteration 200: error is 59.794104 (50 iterations in 0.02 seconds)
## Iteration 250: error is 58.457877 (50 iterations in 0.02 seconds)
## Iteration 300: error is 1.975426 (50 iterations in 0.02 seconds)
## Iteration 350: error is 1.578345 (50 iterations in 0.02 seconds)
## Iteration 400: error is 0.921924 (50 iterations in 0.02 seconds)
## Iteration 450: error is 0.490788 (50 iterations in 0.02 seconds)
## Iteration 500: error is 0.407768 (50 iterations in 0.02 seconds)
## Iteration 550: error is 0.388120 (50 iterations in 0.02 seconds)
## Iteration 600: error is 0.386493 (50 iterations in 0.02 seconds)
## Iteration 650: error is 0.387590 (50 iterations in 0.02 seconds)
## Iteration 700: error is 0.387036 (50 iterations in 0.02 seconds)
## Iteration 750: error is 0.389788 (50 iterations in 0.02 seconds)
## Iteration 800: error is 0.383201 (50 iterations in 0.02 seconds)
## Iteration 850: error is 0.384217 (50 iterations in 0.02 seconds)
## Iteration 900: error is 0.385083 (50 iterations in 0.02 seconds)
## Iteration 950: error is 0.384430 (50 iterations in 0.02 seconds)
## Iteration 1000: error is 0.383389 (50 iterations in 0.02 seconds)
## Iteration 1050: error is 0.385522 (50 iterations in 0.02 seconds)
## Iteration 1100: error is 0.385691 (50 iterations in 0.02 seconds)
## Iteration 1150: error is 0.386405 (50 iterations in 0.02 seconds)
## Iteration 1200: error is 0.383763 (50 iterations in 0.02 seconds)
## Iteration 1250: error is 0.384346 (50 iterations in 0.02 seconds)
## Iteration 1300: error is 0.384842 (50 iterations in 0.02 seconds)
## Iteration 1350: error is 0.383884 (50 iterations in 0.02 seconds)
## Iteration 1400: error is 0.386937 (50 iterations in 0.02 seconds)
## Iteration 1450: error is 0.386463 (50 iterations in 0.02 seconds)
## Iteration 1500: error is 0.385156 (50 iterations in 0.02 seconds)
## Iteration 1550: error is 0.384901 (50 iterations in 0.02 seconds)
## Iteration 1600: error is 0.386133 (50 iterations in 0.02 seconds)
## Iteration 1650: error is 0.384762 (50 iterations in 0.02 seconds)
## Iteration 1700: error is 0.384369 (50 iterations in 0.02 seconds)
## Iteration 1750: error is 0.386143 (50 iterations in 0.02 seconds)
## Iteration 1800: error is 0.384961 (50 iterations in 0.02 seconds)
## Iteration 1850: error is 0.388030 (50 iterations in 0.02 seconds)
## Iteration 1900: error is 0.386162 (50 iterations in 0.02 seconds)
## Iteration 1950: error is 0.388101 (50 iterations in 0.02 seconds)
## Iteration 2000: error is 0.386207 (50 iterations in 0.02 seconds)
## Iteration 2050: error is 0.383446 (50 iterations in 0.02 seconds)
## Iteration 2100: error is 0.387669 (50 iterations in 0.02 seconds)
## Iteration 2150: error is 0.385478 (50 iterations in 0.02 seconds)
## Iteration 2200: error is 0.384587 (50 iterations in 0.02 seconds)
## Iteration 2250: error is 0.384071 (50 iterations in 0.02 seconds)
## Iteration 2300: error is 0.383635 (50 iterations in 0.02 seconds)
## Iteration 2350: error is 0.385712 (50 iterations in 0.02 seconds)
## Iteration 2400: error is 0.382410 (50 iterations in 0.02 seconds)
## Iteration 2450: error is 0.382535 (50 iterations in 0.02 seconds)
## Iteration 2500: error is 0.381783 (50 iterations in 0.02 seconds)
## Iteration 2550: error is 0.385160 (50 iterations in 0.02 seconds)
## Iteration 2600: error is 0.384765 (50 iterations in 0.02 seconds)
## Iteration 2650: error is 0.383549 (50 iterations in 0.02 seconds)
## Iteration 2700: error is 0.383773 (50 iterations in 0.02 seconds)
## Iteration 2750: error is 0.386198 (50 iterations in 0.02 seconds)
## Iteration 2800: error is 0.388300 (50 iterations in 0.02 seconds)
## Iteration 2850: error is 0.384728 (50 iterations in 0.02 seconds)
## Iteration 2900: error is 0.386068 (50 iterations in 0.02 seconds)
## Iteration 2950: error is 0.386615 (50 iterations in 0.02 seconds)
## Iteration 3000: error is 0.383927 (50 iterations in 0.02 seconds)
## Fitting performed in 1.23 seconds.
```

![plot of chunk plantPcoaTsne](figure/plantPcoaTsne-1.png)![plot of chunk plantPcoaTsne](figure/plantPcoaTsne-2.png)

```
## Read the 93 x 93 data matrix successfully!
## Using no_dims = 2, perplexity = 15.000000, and theta = 0.500000
## Computing input similarities...
## Building tree...
##  - point 0 of 93
## Done in 0.01 seconds (sparsity = 0.622731)!
## Learning embedding...
## Iteration 50: error is 62.584210 (50 iterations in 0.03 seconds)
## Iteration 100: error is 58.409560 (50 iterations in 0.03 seconds)
## Iteration 150: error is 61.465872 (50 iterations in 0.03 seconds)
## Iteration 200: error is 60.004902 (50 iterations in 0.03 seconds)
## Iteration 250: error is 61.116887 (50 iterations in 0.03 seconds)
## Iteration 300: error is 2.227520 (50 iterations in 0.02 seconds)
## Iteration 350: error is 1.562859 (50 iterations in 0.02 seconds)
## Iteration 400: error is 1.179063 (50 iterations in 0.02 seconds)
## Iteration 450: error is 0.758984 (50 iterations in 0.02 seconds)
## Iteration 500: error is 0.559950 (50 iterations in 0.02 seconds)
## Iteration 550: error is 0.472228 (50 iterations in 0.02 seconds)
## Iteration 600: error is 0.410945 (50 iterations in 0.02 seconds)
## Iteration 650: error is 0.354830 (50 iterations in 0.02 seconds)
## Iteration 700: error is 0.330662 (50 iterations in 0.02 seconds)
## Iteration 750: error is 0.333520 (50 iterations in 0.02 seconds)
## Iteration 800: error is 0.334242 (50 iterations in 0.02 seconds)
## Iteration 850: error is 0.335824 (50 iterations in 0.02 seconds)
## Iteration 900: error is 0.333828 (50 iterations in 0.02 seconds)
## Iteration 950: error is 0.338101 (50 iterations in 0.02 seconds)
## Iteration 1000: error is 0.336269 (50 iterations in 0.02 seconds)
## Iteration 1050: error is 0.333786 (50 iterations in 0.03 seconds)
## Iteration 1100: error is 0.336763 (50 iterations in 0.02 seconds)
## Iteration 1150: error is 0.331937 (50 iterations in 0.02 seconds)
## Iteration 1200: error is 0.333717 (50 iterations in 0.02 seconds)
## Iteration 1250: error is 0.332186 (50 iterations in 0.02 seconds)
## Iteration 1300: error is 0.332396 (50 iterations in 0.02 seconds)
## Iteration 1350: error is 0.334494 (50 iterations in 0.02 seconds)
## Iteration 1400: error is 0.333192 (50 iterations in 0.02 seconds)
## Iteration 1450: error is 0.334094 (50 iterations in 0.02 seconds)
## Iteration 1500: error is 0.330746 (50 iterations in 0.02 seconds)
## Iteration 1550: error is 0.332237 (50 iterations in 0.02 seconds)
## Iteration 1600: error is 0.331692 (50 iterations in 0.02 seconds)
## Iteration 1650: error is 0.329536 (50 iterations in 0.02 seconds)
## Iteration 1700: error is 0.330279 (50 iterations in 0.02 seconds)
## Iteration 1750: error is 0.331773 (50 iterations in 0.02 seconds)
## Iteration 1800: error is 0.330683 (50 iterations in 0.02 seconds)
## Iteration 1850: error is 0.330146 (50 iterations in 0.02 seconds)
## Iteration 1900: error is 0.330949 (50 iterations in 0.02 seconds)
## Iteration 1950: error is 0.330708 (50 iterations in 0.02 seconds)
## Iteration 2000: error is 0.331051 (50 iterations in 0.02 seconds)
## Iteration 2050: error is 0.329961 (50 iterations in 0.02 seconds)
## Iteration 2100: error is 0.332392 (50 iterations in 0.02 seconds)
## Iteration 2150: error is 0.330381 (50 iterations in 0.02 seconds)
## Iteration 2200: error is 0.333374 (50 iterations in 0.02 seconds)
## Iteration 2250: error is 0.332611 (50 iterations in 0.02 seconds)
## Iteration 2300: error is 0.332196 (50 iterations in 0.02 seconds)
## Iteration 2350: error is 0.331385 (50 iterations in 0.02 seconds)
## Iteration 2400: error is 0.329833 (50 iterations in 0.03 seconds)
## Iteration 2450: error is 0.331202 (50 iterations in 0.02 seconds)
## Iteration 2500: error is 0.330622 (50 iterations in 0.02 seconds)
## Iteration 2550: error is 0.330168 (50 iterations in 0.02 seconds)
## Iteration 2600: error is 0.333865 (50 iterations in 0.03 seconds)
## Iteration 2650: error is 0.335342 (50 iterations in 0.03 seconds)
## Iteration 2700: error is 0.335566 (50 iterations in 0.03 seconds)
## Iteration 2750: error is 0.337193 (50 iterations in 0.02 seconds)
## Iteration 2800: error is 0.335639 (50 iterations in 0.02 seconds)
## Iteration 2850: error is 0.337042 (50 iterations in 0.02 seconds)
## Iteration 2900: error is 0.334436 (50 iterations in 0.02 seconds)
## Iteration 2950: error is 0.333870 (50 iterations in 0.02 seconds)
## Iteration 3000: error is 0.337300 (50 iterations in 0.02 seconds)
## Fitting performed in 1.40 seconds.
```

![plot of chunk plantPcoaTsne](figure/plantPcoaTsne-3.png)![plot of chunk plantPcoaTsne](figure/plantPcoaTsne-4.png)

### Mantel test of different distances

```r
mantels
```

```
## $matK
## $matK$uniW
## Monte-Carlo test
## Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
## 
## Observation: 0.3447672 
## 
## Based on 10000 replicates
## Simulated p-value: 9.999e-05 
## Alternative hypothesis: greater 
## 
##       Std.Obs   Expectation      Variance 
## 11.3227990708 -0.0003840247  0.0009292044 
## 
## $matK$brayW
## Monte-Carlo test
## Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
## 
## Observation: 0.3717721 
## 
## Based on 10000 replicates
## Simulated p-value: 9.999e-05 
## Alternative hypothesis: greater 
## 
##      Std.Obs  Expectation     Variance 
## 1.163445e+01 5.545013e-04 1.018040e-03 
## 
## $matK$brayUW
## Monte-Carlo test
## Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
## 
## Observation: 0.5529449 
## 
## Based on 10000 replicates
## Simulated p-value: 9.999e-05 
## Alternative hypothesis: greater 
## 
##       Std.Obs   Expectation      Variance 
## 18.3486620614 -0.0003570226  0.0009093174 
## 
## 
## $rbcL
## $rbcL$uniW
## Monte-Carlo test
## Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
## 
## Observation: 0.689035 
## 
## Based on 10000 replicates
## Simulated p-value: 9.999e-05 
## Alternative hypothesis: greater 
## 
##      Std.Obs  Expectation     Variance 
## 2.272761e+01 1.117186e-04 9.188278e-04 
## 
## $rbcL$brayW
## Monte-Carlo test
## Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
## 
## Observation: 0.6839004 
## 
## Based on 10000 replicates
## Simulated p-value: 9.999e-05 
## Alternative hypothesis: greater 
## 
##      Std.Obs  Expectation     Variance 
## 2.711865e+01 2.996594e-05 6.359333e-04 
## 
## $rbcL$brayUW
## Monte-Carlo test
## Call: ade4::mantel.rtest(m1 = uniDist, m2 = brayDist, nrepet = 10000)
## 
## Observation: 0.8692566 
## 
## Based on 10000 replicates
## Simulated p-value: 9.999e-05 
## Alternative hypothesis: greater 
## 
##       Std.Obs   Expectation      Variance 
##  3.332147e+01 -4.095156e-05  6.805949e-04
```

### PERMANOVA testing of covariates

```r
plantAdonis
```

```
## $matK
## 
## Call:
## adonis(formula = selectDist ~ bonobo + area2 + malaria, data = ss,      permutations = 1e+07, parallel = 10) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
## bonobo     1    0.5987 0.59874  2.1572 0.02049 0.002361 ** 
## area2      8    5.8053 0.72567  2.6146 0.19864    1e-07 ***
## malaria    1    0.3402 0.34024  1.2259 0.01164 0.179559    
## Residuals 81   22.4814 0.27755         0.76923             
## Total     91   29.2257                 1.00000             
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $rbcL
## 
## Call:
## adonis(formula = selectDist ~ bonobo + area2 + malaria, data = ss,      permutations = 1e+07, parallel = 10) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## bonobo     1     1.084 1.08356 2.95779 0.02801  1e-07 ***
## area2      8     7.229 0.90357 2.46646 0.18684  1e-07 ***
## malaria    1     0.336 0.33640 0.91826 0.00869 0.7085    
## Residuals 82    30.040 0.36634         0.77646           
## Total     92    38.689                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
interactAdonis
```

```
## $matK
## 
## Call:
## adonis(formula = selectDist ~ bonobo + area2 + malaria * bonobo,      data = ss, permutations = 1e+07, parallel = 10) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##                Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
## bonobo          1    0.5987 0.59874  2.1643 0.02049 0.002277 ** 
## area2           8    5.8053 0.72567  2.6231 0.19864    1e-07 ***
## malaria         1    0.3402 0.34024  1.2299 0.01164 0.176727    
## bonobo:malaria  1    0.3498 0.34978  1.2644 0.01197 0.152542    
## Residuals      80   22.1316 0.27664         0.75726             
## Total          91   29.2257                 1.00000             
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $rbcL
## 
## Call:
## adonis(formula = selectDist ~ bonobo + area2 + malaria * bonobo,      data = ss, permutations = 1e+07, parallel = 10) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## bonobo          1     1.084 1.08356 2.95440 0.02801  1e-07 ***
## area2           8     7.229 0.90357 2.46364 0.18684  1e-07 ***
## malaria         1     0.336 0.33640 0.91720 0.00869 0.7117    
## bonobo:malaria  1     0.332 0.33230 0.90604 0.00859 0.7487    
## Residuals      81    29.708 0.36676         0.76787           
## Total          92    38.689                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
chimpAdonis
```

```
## $matK
## 
## Call:
## adonis(formula = chimpDist ~ area2 + malaria, data = ss[!ss$bonobo,      ], permutations = 1e+07, parallel = 10) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
## area2      2    1.3470 0.67350 2.28831 0.16106 0.0009979 ***
## malaria    1    0.2469 0.24691 0.83891 0.02952 0.6179779    
## Residuals 23    6.7694 0.29432         0.80942              
## Total     26    8.3633                 1.00000              
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $rbcL
## 
## Call:
## adonis(formula = chimpDist ~ area2 + malaria, data = ss[!ss$bonobo,      ], permutations = 1e+07, parallel = 10) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## area2      2    1.6955 0.84775  2.1430 0.14237  9e-07 ***
## malaria    1    0.3238 0.32379  0.8185 0.02719  0.858    
## Residuals 25    9.8898 0.39559         0.83044           
## Total     28   11.9091                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
bonoboAdonis
```

```
## $matK
## 
## Call:
## adonis(formula = bonoboDist ~ area2 + malaria, data = ss[ss$bonobo,      ], permutations = 1e+07, parallel = 10) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)    
## area2      6    4.4583 0.74306  2.7571 0.22002   1e-07 ***
## malaria    1    0.4431 0.44311  1.6441 0.02187 0.02601 *  
## Residuals 57   15.3622 0.26951         0.75812            
## Total     64   20.2636                 1.00000            
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $rbcL
## 
## Call:
## adonis(formula = bonoboDist ~ area2 + malaria, data = ss[ss$bonobo,      ], permutations = 1e+07, parallel = 10) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## area2      6    5.5330 0.92217 2.60581 0.21533  1e-07 ***
## malaria    1    0.3449 0.34491 0.97461 0.01342 0.5051    
## Residuals 56   19.8180 0.35389         0.77125           
## Total     63   25.6959                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
tlAdonis
```

```
## $matK
## 
## Call:
## adonis(formula = tlDist ~ malaria, data = ss[ss$isTL, ], permutations = 1e+07,      parallel = 5) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)   
## malaria    1    0.5703 0.57031  1.9651 0.06148 0.00648 **
## Residuals 30    8.7065 0.29022         0.93852           
## Total     31    9.2769                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $rbcL
## 
## Call:
## adonis(formula = tlDist ~ malaria, data = ss[ss$isTL, ], permutations = 1e+07,      parallel = 5) 
## 
## Permutation: free
## Number of permutations: 1e+07
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
## malaria    1    0.4097 0.40966  1.0546 0.03396 0.3114
## Residuals 30   11.6533 0.38844         0.96604       
## Total     31   12.0630                 1.00000
```


