# factorcpt

implements a two-stage methodology for consistent multiple change-point detection under factor modelling. It performs multiple change-point analysis on the common and idiosyncratic components separately, and thus automatically identifies their origins. The package also provides options to implement `Binary Segmentation`, `Wild Binary Segmentation`, `Sparsified Binary Segmentation`, `Wild Sparsified Binary Segmentation`, or `Double CUSUM Binary Segmentation` algorithms, which are proposed for multiple change-point detection in high-dimensional panel data with breaks. 

See

> M. Barigozzi, H. Cho and P. Fryzlewicz (2018) Simultaneous multiple change-point and factor analysis for high-dimensional time series.  

> Y.-N. Li, D. Li and P. Fryzlewicz (2022) Detection of multiple structural breaks in large covariance matrices.


## Installation


The developmental version can be installed from within R using the devtools-package:
```
library(devtools) install_github("markov10000/factorcpt")
```

## Usage

We can generate an example dataset
```
set.seed(1)
x1 <- matrix(rnorm(40*50),nrow=40,ncol=50)
x2 <- matrix(rnorm(40*50),nrow=40,ncol=50)*1.3
x <- cbind(x1,x2)
```

Fit 
```
model1.res <- factor.seg.alg(x, r=1, do.parallel=0, idio.diag=F)
```

Or use get.args() to get arguments
```
args2 <- get.args(x,ref='baseline') #ref='baseline', 'BCF18', or 'LLF22'
model2.res <- factor.seg.alg(x, args=args2)
#or
#model2.res <- factor.seg.alg(x, r=1, args=args2)
```

Print results of model1 and model2
```
model1.res$common.est.cps
model2.res$common.est.cps

model1.res$idio.est.cps
model2.res$idio.est.cps
```
