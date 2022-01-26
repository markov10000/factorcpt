common.seg.alg <- function(gfm, q, scales = -1, sig.lev = 0.01, rule, dw, B = 200, do.parallel = TRUE){
  
  lam <- gfm$lam; f <- gfm$f; nx <- gfm$norm.x
  n <- dim(nx)[2]
  prob.seq <- apply(f[1:q,, drop = FALSE], 1, 
                    function(z){g <- get.gg(z); min(.5, ((g[2]/g[1])^2)^(-1/3)*n^(-1/5))})
  
  y <- c()
  hat.chi <- lam[, 1:q, drop = FALSE] %*% f[1:q,, drop = FALSE]
  for(sc in scales){
    cc <- func_coef(hat.chi, sc)
    if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop = FALSE]
    chi.input <- t(func_input_on(cc))
    y <- rbind(y, log(chi.input))
  }
  
  d <- nrow(y); len <- ncol(y)
  mby <- ceiling(log(d))
  tby <- ceiling(log(len))
  
  mat <- make.tree(y, dw, rule)$mat
  ns <- common.null.stat(mat, lam, f, q, prob.seq, scales, mby, tby, B, do.parallel)
  
  k <- 1; pos.seq <- c()
  while(k <= ncol(mat)){
    pos <- mean(mat[5, k] > ns[k, ])
    pos.seq <- c(pos.seq, pos)
    if(pos <= 1 - sig.lev){
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
  mat <- rbind(mat, pos.seq); mat <- mat[, pos.seq > 1 - sig.lev, drop = FALSE]
  
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
  est.cps <- sort(mat[3, ])
  
  ls <- list(tree = tree, est.cps = est.cps + 2^(-min(scales)) - 1)
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
