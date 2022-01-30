get.args <- function(x, ref = 'baseline'){
p <- dim(x)[1]
n <- dim(x)[2]
args <- list()
if(is.null(args[['dw']])) args[['dw']] <- round(min(log(n)^2, n^(6/7)/4))
if(is.null(args[['rule']])) args[['rule']] <- round(log(n, 2)/2)
if(is.null(args[['scales']])) args[['scales']] <- -(1:floor(log(log(n, 2), 2)))

if(ref == 'LLF22'){
  if(is.null(args[['common.bs.alg']])) args[['common.bs.alg']] <- 'WBS'
  if(is.null(args[['idio.bs.alg']])) args[['idio.bs.alg']] <- 'WSBS'
  if(is.null(args[['transform']])) args[['transform']] <- 'cov'
  if(is.null(args[['bn.op']])) args[['bn.op']] <- 2
  if(is.null(args[['r']])) args[['r']] <- NULL
  if(is.null(args[['gfm.norm']])) args[['gfm.norm']] <- FALSE
  if(is.null(args[['q.len']])) args[['q.len']] <- 1
  if(is.null(args[['factor.as.input']])) args[['factor.as.input']] <- TRUE
  if(is.null(args[['idio.diag']])) args[['idio.diag']] <- FALSE
  if(is.null(args[['do.parallel']])) args[['do.parallel']] <- FALSE
  if(is.null(args[['n.cores']])) args[['n.cores']] <- min(parallel::detectCores() - 1, 3)
  #WSBS
  if(is.null(args[['num.rndint']])) args[['num.rndint']] <- 200
  if(is.null(args[['agg.std.op']])) args[['agg.std.op']] <- 3
  if(is.null(args[['SSIC.pen.fun']]))args[['SSIC.pen.fun']] <- function(x,y)(x^0.5/2)
  if(is.null(args[['common.agg.norm']])) args[['common.agg.norm']] <- 2
  if(is.null(args[['common.agg.thr']])) args[['common.agg.thr']] <- Inf
  if(is.null(args[['common.thr.op']])) args[['common.thr.op']] <- FALSE
  if(is.null(args[['idio.agg.norm']])) args[['common.agg.norm']] <- 2
  if(is.null(args[['idio.agg.thr']])) args[['idio.agg.thr']] <- Inf
  if(is.null(args[['idio.thr.op']])) args[['idio.thr.op']] <- FALSE
  if(is.null(args[['sig.lev']])) args[['sig.lev']] <- .001
  if(is.null(args[['num.boots']])) args[['num.boots']] <- 10
}

if(ref == 'BCF18'){
  if(is.null(args[['common.bs.alg']])) args[['common.bs.alg']] <- 'DBS'
  if(is.null(args[['idio.bs.alg']])) args[['idio.bs.alg']] <- 'DBS'
  if(is.null(args[['transform']])) args[['transform']] <- 'wavelet'
  if(is.null(args[['bn.op']])) args[['bn.op']] <- 2
  if(is.null(args[['r']])) args[['r']] <- NULL
  if(is.null(args[['gfm.norm']])) args[['gfm.norm']] <- TRUE
  if(is.null(args[['max.q']])) args[['max.q']] <- max(round(20, sqrt(min(n, p))))
  if(is.null(args[['q.seq']])) args[['q.seq']] <- NULL
  if(is.null(args[['q.len']])) args[['q.len']] <- 5
  if(is.null(args[['factor.as.input']])) args[['factor.as.input']] <- FALSE
  if(is.null(args[['scales']])) args[['scales']] <- -(1:floor(log(log(n, 2), 2)))
  if(is.null(args[['idio.diag']])) args[['idio.diag']] <- FALSE
  if(is.null(args[['do.parallel']])) args[['do.parallel']] <- FALSE
  if(is.null(args[['n.cores']])) args[['n.cores']] <- min(parallel::detectCores() - 1, 3)
  #DBS
  if(is.null(args[['sig.lev']])) args[['sig.lev']] <- .01
  if(is.null(args[['num.boots']])) args[['num.boots']] <- 200
}

if(ref == 'baseline'){
  if(is.null(args[['common.bs.alg']])) args[['common.bs.alg']] <- 'BS'
  if(is.null(args[['idio.bs.alg']])) args[['idio.bs.alg']] <- 'BS'
  if(is.null(args[['transform']])) args[['transform']] <- 'wavelet'
  if(is.null(args[['bn.op']])) args[['bn.op']] <- 2
  if(is.null(args[['r']])) args[['r']] <- NULL
  if(is.null(args[['gfm.norm']])) args[['gfm.norm']] <- FALSE
  if(is.null(args[['q.len']])) args[['q.len']] <- 1
  if(is.null(args[['factor.as.input']])) args[['factor.as.input']] <- TRUE
  if(is.null(args[['idio.diag']])) args[['idio.diag']] <- FALSE
  if(is.null(args[['do.parallel']])) args[['do.parallel']] <- FALSE
  if(is.null(args[['n.cores']])) args[['n.cores']] <- min(parallel::detectCores() - 1, 3)
  #WSBS
  if(is.null(args[['num.rndint']])) args[['num.rndint']] <- 200
  if(is.null(args[['agg.std.op']])) args[['agg.std.op']] <- 3
  if(is.null(args[['SSIC.pen.fun']]))args[['SSIC.pen.fun']] <- function(x,y)(x^0.5/2)
  if(is.null(args[['common.agg.norm']])) args[['common.agg.norm']] <- 2
  if(is.null(args[['common.agg.thr']])) args[['common.agg.thr']] <- Inf
  if(is.null(args[['common.thr.op']])) args[['common.thr.op']] <- FALSE
  if(is.null(args[['idio.agg.norm']])) args[['common.agg.norm']] <- 2
  if(is.null(args[['idio.agg.thr']])) args[['idio.agg.thr']] <- Inf
  if(is.null(args[['idio.thr.op']])) args[['idio.thr.op']] <- FALSE
  if(is.null(args[['sig.lev']])) args[['sig.lev']] <- .001
  if(is.null(args[['num.boots']])) args[['num.boots']] <- 10
}
args
}