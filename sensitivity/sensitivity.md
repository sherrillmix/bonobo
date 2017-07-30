# Analyze sensitivity of plasmodium PCR detection

## Libraries

```r
library(actuar)
```

```
## 
## Attaching package: 'actuar'
```

```
## The following object is masked from 'package:grDevices':
## 
##     cm
```

```r
packageVersion('actuar')
```

```
## [1] '2.1.1'
```

## Load data
```{r]
pcr<-read.csv('pcrPos.csv',stringsAsFactors=FALSE)
```

## Find maximum likelihood estimate for sensitivity using a zero truncated binomial

```r
like<-function(theta,positives,replicates)-sum(dztbinom(positives,replicates,theta,log=TRUE))
mle<-optimize(like,interval=0:1,pcr$positive,pcr$replicates)
```

```
## Error in dztbinom(positives, replicates, theta, log = TRUE): object 'pcr' not found
```

```r
print(mle)
```

```
## Error in print(mle): object 'mle' not found
```

## Compare the observed data with that expected from the estimated sensitivity

```r
hist(pcr$positive[pcr$replicates==8],breaks=0:9-.5,xlab='Number of positive PCRs',main='8 replicate samples',las=1)
```

```
## Error in hist(pcr$positive[pcr$replicates == 8], breaks = 0:9 - 0.5, xlab = "Number of positive PCRs", : object 'pcr' not found
```

```r
ps<-dztbinom(0:8,8,mle$minimum)
```

```
## Error in dztbinom(0:8, 8, mle$minimum): object 'mle' not found
```

```r
preds<-ps*sum(pcr$replicates==8)
```

```
## Error in eval(expr, envir, enclos): object 'ps' not found
```

```r
segments(0:8-.5,preds,0:8+.5,preds,col='red')
```

```
## Error in segments(0:8 - 0.5, preds, 0:8 + 0.5, preds, col = "red"): object 'preds' not found
```

```r
hist(pcr$positive[pcr$replicates==10],breaks=0:10-.5,xlab='Number of positive PCRs',main='10 replicate samples',ylim=c(0,1.2),las=1)
```

```
## Error in hist(pcr$positive[pcr$replicates == 10], breaks = 0:10 - 0.5, : object 'pcr' not found
```

```r
ps<-dztbinom(0:10,10,mle$minimum)
```

```
## Error in dztbinom(0:10, 10, mle$minimum): object 'mle' not found
```

```r
preds<-ps*sum(pcr$replicates==10)
```

```
## Error in eval(expr, envir, enclos): object 'ps' not found
```

```r
segments(0:10-.5,preds,0:10+.5,preds,col='red')
```

```
## Error in segments(0:10 - 0.5, preds, 0:10 + 0.5, preds, col = "red"): object 'preds' not found
```

```r
dev.off()
```

```
## null device 
##           1
```
