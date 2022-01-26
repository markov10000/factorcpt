#' @keywords internal
make.tree <- function(y, dw, rule){
  
  len <- ncol(y)
  tree <- list(matrix(0, 6, 1))
  mat <- c()
  
  fd <- func_dc(y, .5)
  stat <- fd$res
  test.stat <- max(stat[-c((1:dw), (len - dw):len)])
  hat.chp <- min(which(stat == test.stat))
  
  tree[[1]][1, 1] <- 1
  tree[[1]][2, 1] <- 1
  tree[[1]][3, 1] <- hat.chp
  tree[[1]][4, 1] <- len
  tree[[1]][5, 1] <- 0
  tree[[1]][6, 1] <- test.stat
  mat <- cbind(mat, c(tree[[1]][-5, ], 1, 1))
  
  j <- 1
  while(length(tree) == j & j < rule){
    npc <- dim(tree[[j]])[2]
    if(sum(tree[[j]][4, ] - tree[[j]][2, ] - rep(4 * dw, npc) > 0)){
      ncc <- 0; i <- 1
      while(i <= npc){
        if(tree[[j]][3, i] - tree[[j]][2, i] + 1 > 4 * dw){
          s <- tree[[j]][2, i]; e <- tree[[j]][3, i]
          fd <- func_dc(y[, s:e], .5)
          stat <- fd$res
          test.stat <- max(stat[-c((1:dw), (e - s + 1 - dw):(e - s + 1))])
          hat.chp <- s + min(which(stat == test.stat)) - 1
          
          if(length(tree) == j) tree <- c(tree, list(matrix(0, 6, 0)))
          ncc <- ncc + 1
          tree[[j + 1]] <- matrix(c(tree[[j + 1]], matrix(0, 6, 1)), 6, ncc)
          tree[[j + 1]][1, ncc] <- 2 * tree[[j]][1, i] - 1
          tree[[j + 1]][2, ncc] <- s
          tree[[j + 1]][3, ncc] <- hat.chp
          tree[[j + 1]][4, ncc] <- e
          tree[[j + 1]][5, ncc] <- 0
          tree[[j + 1]][6, ncc] <- test.stat
          mat <- cbind(mat, c(tree[[j + 1]][-5, ncc], j + 1, ncc))
        }
        if(tree[[j]][4, i] - tree[[j]][3, i] > 4 * dw){
          s <- tree[[j]][3, i] + 1; e <- tree[[j]][4, i]
          fd <- func_dc(y[, s:e], .5)
          stat <- fd$res
          test.stat <- max(stat[-c((1:dw), (e - s + 1 - dw):(e - s + 1))])
          hat.chp <- s + min(which(stat == test.stat)) - 1
          
          if(length(tree) == j) tree <- c(tree, list(matrix(0, 6, 0)))
          ncc <- ncc + 1
          tree[[j + 1]] <- matrix(c(tree[[j + 1]], matrix(0, 6, 1)), 6, ncc)
          tree[[j + 1]][1, ncc] <- 2 * tree[[j]][1, i]
          tree[[j + 1]][2, ncc] <- s
          tree[[j + 1]][3, ncc] <- hat.chp
          tree[[j + 1]][4, ncc] <- e
          tree[[j + 1]][5, ncc] <- 0
          tree[[j + 1]][6, ncc] <- test.stat
          mat <- cbind(mat, c(tree[[j + 1]][-5, ncc], j + 1, ncc))
        }
        i <- i + 1
      }
      j <- j + 1
    } else break
  }
  
  list(tree = tree, mat = mat)
  
}

#' @keywords internal
tri.kern <- function(h){
  
  filter <- rep(0, h + 1)
  i <- 0
  while (i <= h) {
    u <- i/h
    if (u < 1/2)
      filter[i + 1] <- 1
    if (u >= 1/2 & u < 1)
      filter[i + 1] <- 2 * (1 - u)
    if (u > 1)
      break
    i <- i + 1
  }
  filter
  
}

#' @keywords internal
get.gg <- function(z, M = NULL, C = 2, max.K = 5){
  
  len <- length(z)
  max.K <- max(max.K, sqrt(log(len)))
  acv <- acf(z, type = "covariance", lag.max = len - 1, plot = FALSE)$acf[,, 1]
  if(is.null(M)){
    l <- 1; ind <- 0
    while(l < sqrt(len)){
      if(abs(acv[l + 1]) / acv[1] < C * sqrt(log(len)/len)){
        ind <- ind + 1
      } else if(ind > 0) ind <- 0
      if(ind == max.K) break
      l <- l + 1
    }
    lam <- max(1/2, l - max.K); M <- 2*lam
  }
  k <- tri.kern(M)

  c(acv[1] + 2 * sum(k[-1] * acv[2:(M + 1)]), 2 * sum(k[-1] * (1:M) * acv[2:(M + 1)]))

}
