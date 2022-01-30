common.seg.alg <- function(gfm, q, bs.alg = 'DBS', ...){
  
  lam <- gfm$lam; f <- gfm$f; nx <- gfm$norm.x
  p <- dim(nx)[1]
  n <- dim(nx)[2]
  
  # -------------------- Input Checks -------------------
  
  args <- list(...)
  # ----- Fill in Defaults -----
  if(is.null(bs.alg)) bs.alg <- 'DBS'
  if(is.null(args[['transform']]))args[['transform']] <- 'wavelet'
  if(is.null(args[['factor.as.input']])) args[['factor.as.input']] <- FALSE
  if(is.null(args[['rule']])) args[['rule']] <- round(log(n, 2)/2)
  if(is.null(args[['dw']])) args[['dw']] <- round(min(log(n)^2, n^(6/7)/4))
  if(is.null(args[['do.parallel']])) args[['do.parallel']] <- FALSE
  if(is.null(args[['n.cores']])) args[['n.cores']] <- min(parallel::detectCores() - 1, 3)
  if(is.null(args[['num.boots']])) args[['num.boots']] <- 200
  if(is.null(args[['sig.lev']])) args[['sig.lev']] <- .01
  if(bs.alg=='DBS'){
  if(is.null(args[['scales']])) args[['scales']] <- -(1:floor(log(log(n, 2), 2)))
  }
  if(bs.alg%in%c('BS','WBS','SBS','WSBS')){
    if(is.null(args[['num.rndint']])) args[['num.rndint']] <- 200
    if(is.null(args[['SSIC.pen.fun']]))args[['SSIC.pen.fun']] <- function(x,y)(sqrt(x))
    if(is.null(args[['agg.std.op']])) args[['agg.std.op']] <- 0
    if(is.null(args[['agg.norm']])) args[['agg.norm']] <- 2
    if(is.null(args[['agg.thr']])) args[['agg.thr']] <- Inf
    if(is.null(args[['thr.op']])) args[['thr.op']] <- FALSE
  }
  # -------------------- Transform hat.chi (hat.f) to get y -------------------  
  if(args[['factor.as.input']]){
    hat.chi <- f[1:q,, drop = FALSE]
  }else{
    hat.chi <- lam[, 1:q, drop = FALSE] %*% f[1:q,, drop = FALSE]
  }
  dim.hat.chi <- dim(hat.chi)[1]
  y <- c()
  if (args[['transform']] == 'wavelet'){
    for(sc in args[['scales']]){
      cc <- func_coef(hat.chi, sc)
      if(sc > min(args[['scales']])) cc <- cc[, -(1:(2^(-min(args[['scales']]))-2^(-sc))), drop = FALSE]
      chi.input <- t(func_input_on(cc))
      y <- rbind(y, log(chi.input))
    }
  }
  if(args[['transform']] == 'cov'){ #covariance
    all.pairs <- combn(1:dim.hat.chi,2)
    y <- rbind(hat.chi^2, hat.chi[all.pairs[1,], ,drop=F] * hat.chi[all.pairs[2,], ,drop=F])
  }
  if(args[['transform']] == 'addmin'){ #add and minus
    all.pairs <- combn(1:dim.hat.chi,2)
    y <- rbind(4*hat.chi^2, (hat.chi[all.pairs[1,], ,drop=F] + hat.chi[all.pairs[2,], ,drop=F])^2, (hat.chi[all.pairs[1,], ,drop=F] - hat.chi[all.pairs[2,], ,drop=F])^2)
  }
  d <- nrow(y); len <- ncol(y)
  
  # -------------------- Make tree -------------------
  matrnd <- c()
  if (bs.alg == 'DBS'){ 
    mat <- DBS.make.tree(y, args[['dw']], args[['rule']])$mat
  }
  if (bs.alg == 'BS'){
    mt <- WBS.make.tree(y, args[['dw']], args[['rule']], agg.std.op=args[['agg.std.op']], WBS=0,SBS=0, M=args[['num.rndint']], norm0=args[['agg.norm']], norm.thr=args[['agg.thr']], thr.op=args[['thr.op']],
    do.parallel=args[['do.parallel']])
    mat <- mt$mat
  }
  if (bs.alg == 'WBS'){
    mt <- WBS.make.tree(y, args[['dw']], args[['rule']], agg.std.op=args[['agg.std.op']], WBS=1,SBS=0, M=args[['num.rndint']], norm0=args[['agg.norm']], norm.thr=args[['agg.thr']], thr.op=args[['thr.op']],
    do.parallel=args[['do.parallel']])
    mat <- mt$mat
    matrnd <- mt$matrnd
  }
  if (bs.alg == 'WSBS'){
    mt <- WBS.make.tree(y, args[['dw']], args[['rule']], agg.std.op=args[['agg.std.op']], WBS=1,SBS=1, M=args[['num.rndint']], norm0=args[['agg.norm']], norm.thr=args[['agg.thr']], thr.op=args[['thr.op']],
    do.parallel=args[['do.parallel']])
    mat <- mt$mat
    matrnd <- mt$matrnd
  }
  if (bs.alg == 'SBS'){
    mt <- WBS.make.tree(y, args[['dw']], args[['rule']], agg.std.op=args[['agg.std.op']], WBS=0,SBS=1, M=args[['num.rndint']], norm0=args[['agg.norm']], norm.thr=args[['agg.thr']], thr.op=args[['thr.op']],
    do.parallel=args[['do.parallel']])
    mat <- mt$mat
  }

  # -------------------- Prune tree -------------------
  if (bs.alg == 'DBS'){  
    mby <- ceiling(log(d))
    tby <- ceiling(log(len))
    prob.seq <- apply(f[1:q,, drop = FALSE], 1, 
                    function(z){g <- get.gg(z); min(.5, ((g[2]/g[1])^2)^(-1/3)*n^(-1/5))})
    ns <- common.null.stat(mat, lam, f, q, prob.seq, args[['scales']], mby, tby, args[['num.boots']], args[['do.parallel']])
  
    k <- 1; pos.seq <- c()
    while(k <= ncol(mat)){
      pos <- mean(mat[5, k] > ns[k, ])
      pos.seq <- c(pos.seq, pos)
      if(pos <= 1 - args[['sig.lev']]){
        l <- mat[6, k] + 1; rm.ind <- mat[1, k]
        while(l <= max(mat[6, ])){
          ind <- which(mat[6, ] == l & is.element(mat[1, ], c(rm.ind * 2 - 1, rm.ind * 2)))
          if(length(ind) > 0){
            rm.ind <- mat[1, ind]
            mat <- mat[, -ind, drop = FALSE]
            l <- l + 1
          } else break
        }
      }
      k <- k+1
    }
    mat <- rbind(mat, pos.seq); mat <- mat[, pos.seq > 1 - args[['sig.lev']], drop = FALSE]
  


    tree <- list(matrix(0, 6, 1))
    if(dim(mat)[2] > 0){
      for(l in 1:length(unique(mat[6, ]))){
        j <- unique(mat[6, ])[l]
        for(ncc in 1:sum(mat[6, ]==j)){
          k <- sum(mat[6, ] < j) + ncc
          if(length(tree) < j) tree <- c(tree, list(matrix(0, 6, 0)))
          tree[[l]] <- matrix(c(tree[[l]], matrix(0, 6, 1)), 6, ncc)
          tree[[l]][1, ncc] <- mat[1, k]
          tree[[l]][2, ncc] <- mat[2, k]
          tree[[l]][3, ncc] <- mat[3, k]
          tree[[l]][4, ncc] <- mat[4, k]
          tree[[l]][5, ncc] <- mat[8, k]
          tree[[l]][6, ncc] <- mat[5, k]
        }
      }
    }
  est.cps <- sort(mat[3, ]) + (args[['transform']]=='wavelet') * (2^(-min(args[['scales']])) - 1)
  ls <- list(tree = tree, est.cps = est.cps)
  }else{ #BS, WBS, SBS, WSBS
    if(is.null(mat)) {cpt.trylist<-c()}else{cpt.trylist <- mat[3,order(mat[5,], decreasing = TRUE)]}
    ct <- SSIC(y, cpt.trylist, SSIC.pen.fun=args[['SSIC.pen.fun']])
    cus.qtl <- select.thr(y - means.between.cpt(y, ct$cpt), args[['dw']], args[['agg.std.op']], args[['num.boots']], select= c(1-args[['sig.lev']],seq(0.99,1,0.001)))
    if (mt$flag[1]==1){
      cat('Try common.norm.thr=',cus.qtl[1])
    }
    est.cps <- ct$cpt + (args[['transform']]=='wavelet') * (2^(-min(args[['scales']])) - 1)
    ls <- list( BSmat=mat, matrnd=matrnd, SSIC=ct, est.cps=est.cps, norm=args[['agg.norm']],idio.norm.thr=args[['agg.thr']], cus.qtl=cus.qtl[-1], cus.thr=cus.qtl[1])    
  }
  return(ls)
  
}

common.null.stat <- function(mat, lam, f, q, prob.seq, scales, mby, tby, B, do.parallel){
  
  len <- dim(f)[2]
  if(do.parallel){
    null.stat <- foreach::foreach(ll = iterators::iter(1:B), .combine=cbind, .packages = c('Rcpp', 'RcppArmadillo', 'factorcpt')) %dopar% {
      set.seed(ll)
      
      boot.chi <- 0
      for(qq in 1:q){
        ind <- c()
        while(length(ind) < len){
          L <- max(rgeom(1, prob.seq[qq]), 1 + 2^(-min(scales))); I <- sample(len, 1)
          ind <- c(ind, rep(1:len, 1 + ceiling(L/len))[I:(I + L - 1)])
        }
        ind <- ind[1:len]
        boot.chi <- boot.chi + lam[, qq, drop = FALSE]%*%f[qq, ind, drop = FALSE]
      }
      by <- c()
      for(sc in scales){
        cc <- func_coef(boot.chi, sc)
        if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop = FALSE]
        z <- t(func_input_on(cc))
        z[z == 0] <- min(z[z > 0])
        by <- rbind(by, log(z))
      }
      tmp <- c()
      for(i in 1:ncol(mat)){
        s <- mat[2, i]; e <- mat[4, i]
        bfd <- func_dc_by(by[, s:e], .5, mby, tby)
        tmp <- c(tmp, max(bfd$res))
        rm(bfd)
      }
      rm(by)
      tmp
    }
  } else{
    null.stat <- foreach::foreach(ll = iterators::iter(1:B), .combine = cbind, .packages = c('Rcpp', 'RcppArmadillo')) %do% {
      set.seed(ll)
      
      boot.chi <- 0
      for(qq in 1:q){
        ind <- c()
        while(length(ind) < len){
          L <- max(rgeom(1, prob.seq[qq]), 1 + 2^(-min(scales))); I <- sample(len, 1)
          ind <- c(ind, rep(1:len, 1 + ceiling(L/len))[I:(I + L - 1)])
        }
        ind <- ind[1:len]
        boot.chi <- boot.chi + lam[, qq, drop = FALSE]%*%f[qq, ind, drop = FALSE]
      }
      by <- c()
      for(sc in scales){
        cc <- func_coef(boot.chi, sc)
        if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop = FALSE]
        cc[cc == 0] <- mean(cc)
        by <- rbind(by, log(t(func_input_on(cc))))
      }
      tmp <- c()
      for(i in 1:ncol(mat)){
        s <- mat[2, i]; e <- mat[4, i]
        bfd <- func_dc_by(by[, s:e], .5, mby, tby)
        tmp <- c(tmp, max(bfd$res))
        rm(bfd)
      }
      rm(by)
      tmp
    }
  }
  null.stat
  
}
