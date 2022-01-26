# r = NULL; bn.op = 2; sig.lev = .01; max.q = NULL; q.seq = NULL; scales = NULL; rule = NULL; B = 200; do.diag = idio.diag = FALSE; do.parallel = FALSE
library(doParallel)

factor.seg.alg <- function(x, r = NULL, bn.op = 2, sig.lev = .01, max.q = NULL, q.seq = NULL,
                           scales = NULL, rule = NULL, B = 200, idio.diag = FALSE,
                           do.parallel = FALSE, n.cores = min(parallel::detectCores() - 1, 3)){

	p <- dim(x)[1]
	n <- dim(x)[2]
	q.len <- 5
	dw <- round(min(log(n)^2, n^(6/7)/4))
	
	if(is.null(max.q)) max.q <- max(round(20, sqrt(min(n, p))))
	if(is.null(rule)) rule <- round(log(n, 2)/2)
	if(do.parallel){
	  cl <- parallel::makeCluster(n.cores); doParallel::registerDoParallel(cl)
	}
  if(is.null(scales)) scales <- -(1:floor(log(log(n, 2), 2)))
	
	gfm <- get.factor.model(x, max.q = max.q, q = r, bn.op = bn.op)
	if(is.null(r)) q.hat <- gfm$q.hat else q.hat <- r
	if(is.null(q.seq)) q.seq <- sort(unique(c(max(1, q.hat), max.q, 
	                                          round(seq(max(1, q.hat), max(min(p - 1, 2 * q.hat), max.q), length.out = q.len)))), 
	                                 decreasing = FALSE)

	est.cps <- cs.list <- list()
	for(qq in q.seq){
		cs <- common.seg.alg(gfm, q = qq, scales = scales, sig.lev = sig.lev, 
		                     rule = rule, dw = dw, B = B, do.parallel = do.parallel)
		est.cps <- c(est.cps, list(cs$est.cps))
		cs.list <- c(cs.list, list(cs))
	}
	qq <- max(which(unlist(lapply(est.cps, length)) == max(unlist(lapply(est.cps, length)))))
	q <- q.seq[qq]
	cs <- cs.list[[qq]]
  common.est.cps <- est.cps[[qq]]

	is <- idio.seg.alg(gfm, q = q, scales = scales, sig.lev = sig.lev, rule = rule, dw = dw, 
	                   do.diag = idio.diag, B = B, do.parallel = do.parallel)
	idio.est.cps <- is$est.cps
	if(do.parallel) parallel::stopCluster(cl)

	return(list(common.res = cs.list, r = q, common.est.cps = common.est.cps, 
	            idio.res = is, idio.est.cps = idio.est.cps, gfm = gfm, q.seq = q.seq))

}

#' @keywords internal
get.factor.model <- function(x, max.q = NULL, q = NULL, bn.op = 2, normalisation = TRUE){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
	cnt <- min(n, p)
	if(is.null(max.q)) max.q <- round(sqrt(cnt))

	if(normalisation){
	  mean.x <- apply(x, 1, mean)
	  sd.x <- apply(x, 1, sd)
	  x <- t(scale(t(x)))
	} else mean.x <- rep(0, p); sd.x <- rep(1, p)
	
	xx <- x %*% t(x)/n
	eig <- svd(xx, nu = max.q, nv = 0)
	lam <- eig$u * sqrt(p)
	f <- t(eig$u)%*%x / sqrt(p)

	if(is.null(q)){
		ic <- rep(0, 1 + max.q)
		ic[1] <- (bn.op <= 4) * log(mean(x^2)) + (bn.op == 5) * mean(x^2)
		l <- 1
		while(l <= max.q){
			hchi <- lam[, 1:l, drop = FALSE]%*%f[1:l,, drop = FALSE]
			ic[l + 1] <- (bn.op <= 4) * log(mean((x - hchi)^2)) +
			  (bn.op == 1) * l * (n + p)/(n * p) * log(n * p/(n + p)) +
				(bn.op == 2) * l * (n + p)/(n * p) * log(cnt) +
				(bn.op == 3) * l * log(cnt)/cnt +
				(bn.op == 4) * l * ((n + p - l) * log(n * p)/(n * p) + (n + p)/(n * p) * log(cnt))/2 +
				(bn.op == 5) * (mean((x - hchi)^2) + l * mean((x - hchi)^2)*(n + p - l) * log(n * p)/(n * p))
			l <- l + 1
		}
		q.hat <- which(ic == min(ic)) - 1
	} else{
		ic <- rep(0, max.q + 1)
		q.hat <- q
	}

	return(list(lam = lam, f = f, norm.x = x, mean.x = mean.x, sd.x = sd.x,
	            q.hat = q.hat, max.q = max.q, ic = ic))
}


post.cpts.analysis <- function(x, est.cps, cutoff.seq = seq(.5, .95, by = .05), do.plot = TRUE){

  p <- dim(x)[1]
  n <- dim(x)[2]
  B <- length(est.cps)
  brks <- c(0, est.cps, n)
  C <- length(cutoff.seq)

  heat.mat <- matrix(0, ncol = C, nrow = n)
  for(b in 1:(B + 1)){
    s <- brks[b] + 1; e <- brks[b + 1]
    len <- e - s + 1
    z <- x[, s:e, drop = FALSE]
    z <- z[!apply(z, 1, sd) == 0, ]
    z <- t(scale(t(z)))
    eig <- eigen(z %*% t(z) / len)
    z.eval <- eig$values[1:min(dim(z) - 1)]
    for(c in 1:C) heat.mat[s:e, c] <- min(which(cumsum(z.eval)/sum(z.eval) > cutoff.seq[c]))
  }
  
  if(do.plot){
    par(mar = c(4, 4.5, 2, 7))
    image(heat.mat, col = tim.colors(diff(range(heat.mat)) + 1), 
          breaks = min(heat.mat):(max(heat.mat) + 1),
          xlab = 'time', ylab = 'c', axes = FALSE)
    abline(v = est.cps / n, col = 8, lwd = 2)
    axis(side = 1, at = seq(0, 1, length.out = 10), labels = round(seq(1, n, length.out = 10)))
    axis(side = 2, at = seq(0, 1, length.out = 10), labels = seq(min(cutoff.seq), max(cutoff.seq), length.out = 10))
    fields::imagePlot(heat.mat, col = tim.colors(diff(range(heat.mat)) + 1), 
                      breaks = min(heat.mat):(max(heat.mat) + 1), legend.only = TRUE)
  }
  
  return(heat.mat)

}

