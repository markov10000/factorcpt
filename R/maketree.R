#' @keywords internal
DBS.make.tree <- function(y, dw, rule){
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
WBS.make.tree <- function(y, dw, rule=NULL, agg.std.op=0, WBS=1, SBS=1, M=100, norm0=2, norm.thr=0, thr.op=FALSE, do.parallel=TRUE, setseed=0 ,WBS.op=0){
  p <- dim(y)[1]
  len <- dim(y)[2]

  set.seed(setseed)
  flag=c(0,0)
  if (WBS==0) M <- 0
  if (M==0) WBS <- 0

  if(is.null(rule)) rule <- round(log(len, 2)/2)
  if(missing(norm0)) norm0<-1
  if(is.infinite(norm0))  {
    normf <- function(x,norm1=norm0)return(apply(abs(x),2,max))
  }else if(norm0>=1){
      normf <- function(x,norm1=norm0){
        #x <- as.matrix(x)
        #if (dim(x)[2] == 1) x <- t(x) 
        return(colSums((abs(x))^norm1)^(1/norm1))
      } ##alter / colSums(abs(x)>=norm.thr)
  }else{
    normf <- function(x,norm1=norm0)return(colSums(abs(x)))
  }
  
  tree <- list(matrix(0, 8, 1))
  screen.tree <- list(list())
  mat <- c()

  matrnd <- matrix(0,4,1)
  matrnd[1,1] <- 1
  matrnd[2,1] <- len
  stat <- calcul.cusum(y,dw=dw,agg.std.op=agg.std.op)
  screen.id <- 1:p
    if (SBS) screen.id <- which(apply(abs(stat),1,max)>norm.thr) #FIRST SCREENING
    if(!length(screen.id)){
      cat('norm.thr may be too large. norm.thr=')
      cat(norm.thr, '\n')
        matrnd[3,1] <- 0
        matrnd[4,1] <- NA
        flag[1] <- 1
        return(list(mat=mat, matrnd=matrnd, flagMT=flag))
    }else{
        norm.stat <- normf(stat[screen.id,-c((1:dw), (len-dw+1):len),drop=F])
        matrnd[3,1] <- max(norm.stat)
    matrnd[4,1] <- dw + min(which(norm.stat==matrnd[3,1])) #hat.chp
    }
if(WBS) {
  matrnd0<- matrix(0,4,M) #random intervals # at least 4*dw
  rnd1 <- sample(1:(len-4*dw), M, replace = TRUE)
  rnd2 <- sample(1:(len-4*dw), M, replace = TRUE)
  matrnd0[1,1:M] <- pmin(rnd1, rnd2)     # random l_m
  matrnd0[2,1:M] <- pmax(rnd1, rnd2) + 4*dw   # random u_m

  rnd.ind<-rep(TRUE,M)
  matrnd0[3:4,1:M] <- calcul.matrnd(y,dw,screen.id,matrnd0[,1:M],normf,agg.std.op,do.parallel) 
  matrnd <- cbind(matrnd0,matrnd)
}   

#AFTER first change point finded
test.stat.ind<-which.max(matrnd[3,])
#Make a tree root
  tree[[1]][1, 1] <- 1
  tree[[1]][2, 1] <- 1
  tree[[1]][3, 1] <- matrnd[4,test.stat.ind]
  tree[[1]][4, 1] <- len
  tree[[1]][5, 1] <- Inf
  tree[[1]][6, 1] <- matrnd[3,test.stat.ind]
  tree[[1]][7, 1] <- length(screen.id)
  tree[[1]][8, 1] <- length(matrnd[3,])
  mat <- cbind(mat, c(tree[[1]][-5, ], 1, 1))
  screen.tree[[1]][[1]] <- screen.id

j <- 1
  while(length(tree)==j & j < rule){
    npc <- dim(tree[[j]])[2]
    if(sum(tree[[j]][4,]-tree[[j]][2,]-rep(4*dw, npc)>0)){
      ncc <- 0; i <- 1
      while(i <= npc){
        if(tree[[j]][3, i]-tree[[j]][2, i]>4*dw & length(screen.tree[[j]][[i]])){#binary I
          matrnd <- matrix(0,4,1)
          s <- tree[[j]][2, i]; e <- tree[[j]][3, i]-1
          stat <- calcul.cusum(y[screen.tree[[j]][[i]], s:e,drop=F],dw=dw,agg.std.op=agg.std.op)
          if(thr.op) {
            thr.op.const <- ((e-s+1)/len)^0.5 ###
          }else{
            thr.op.const <- 1
          }
          if(SBS) screen.id <- screen.tree[[j]][[i]][apply(abs(stat),1,max) > thr.op.const * norm.thr] #SCREENING
          if(length(screen.id)){
            if(WBS){
              rnd.ind <- (matrnd0[1,] >= s) & (matrnd0[2,] <= e)
              rnd.len <- sum(rnd.ind)
              if(rnd.len){
                matrnd <- matrnd0[,rnd.ind,drop=F]
                if(SBS) matrnd[3:4,] <- calcul.matrnd(y,dw,screen.id,matrnd,normf,agg.std.op,do.parallel)
              }else{
                flag[1]=1
                cat('num.randint may be too small: num.randint=')
                cat(M, '\n')
              }
            }else{#WBS==0
              matrnd[1:2,] <- c(s,e)
              len <- e - s + 1
              norm.stat <- normf(stat[,-c((1:dw), (len-dw+1):len),drop=F])
              matrnd[3,1] <- max(norm.stat)
              matrnd[4,1] <- s + dw + min(which(norm.stat==matrnd[3,1])) -1 #hat.chp
            }
          }else{flag[2]=1}#thr should be smaller
          test.stat.ind<-which.max(matrnd[3,])
          if(matrnd[3,test.stat.ind]){# make a new tree node            
            if(length(tree)==j) {
              tree <- c(tree, list(matrix(0, 8, 0)))
              screen.tree[[j+1]] <- list()
            }
            ncc <- ncc+1
            tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 8, 1)), 8, ncc)
            tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]-1
            tree[[j+1]][2, ncc] <- s
            tree[[j+1]][3, ncc] <- matrnd[4,test.stat.ind]
            tree[[j+1]][4, ncc] <- e
            tree[[j+1]][5, ncc] <- Inf
            tree[[j+1]][6, ncc] <- matrnd[3,test.stat.ind]
            tree[[j+1]][7, ncc] <- length(screen.id)
            tree[[j+1]][8, ncc] <- length(matrnd[3,])
            screen.tree[[j+1]] <- cbind(screen.tree[[j+1]], list())
            screen.tree[[j+1]][[ncc]] <- screen.id
            mat <- cbind(mat, c(tree[[j+1]][-5, ncc], j+1, ncc))
          }
        }
        if(tree[[j]][4, i]-tree[[j]][3, i]+1>4*dw & length(screen.tree[[j]][[i]])){#binary II
          matrnd <- matrix(0,4,1)
          s <- tree[[j]][3, i]; e <- tree[[j]][4, i]
          stat <- calcul.cusum(y[screen.tree[[j]][[i]], s:e,drop=F],dw=dw,agg.std.op=agg.std.op)
          if(thr.op) {
            thr.op.const <- ((e-s+1)/len)^0.5 ###
          }else{
            thr.op.const <- 1
          }
          if(SBS) screen.id <- screen.tree[[j]][[i]][apply(abs(stat),1,max) > thr.op.const * norm.thr] #SCREENING
          if(length(screen.id)){
            if(WBS){
              rnd.ind <- (matrnd0[1,] >= s) & (matrnd0[2,] <= e)
              rnd.len <- sum(rnd.ind)
              if(rnd.len){
                matrnd <- matrnd0[,rnd.ind,drop=F]
                if(SBS) matrnd[3:4,] <- calcul.matrnd(y,dw,screen.id,matrnd,normf,agg.std.op,do.parallel)
              }else{
              flag[1]=1
              cat('num.randint may be too small: num.randint=')
              cat(M, '\n')
              }
            }else{#WBS==0
              matrnd[1:2,] <- c(s,e)
              len <- e - s + 1
              norm.stat <- normf(stat[,-c((1:dw), (len-dw+1):len),drop=F])              
              matrnd[3,1] <- max(norm.stat)
              matrnd[4,1] <- s + dw + min(which(norm.stat==matrnd[3,1])) -1#hat.chp
            }
          }else{flag[2]=1}#thr should be smaller
          test.stat.ind<-which.max(matrnd[3,])
          if(matrnd[3,test.stat.ind]){
            if(length(tree)==j) {
              tree <- c(tree, list(matrix(0, 8, 0)))
              screen.tree[[j+1]] <- list()
            }
            ncc <- ncc+1
            tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 8, 1)), 8, ncc)
            tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]
            tree[[j+1]][2, ncc] <- s
            tree[[j+1]][3, ncc] <- matrnd[4,test.stat.ind]
            tree[[j+1]][4, ncc] <- e
            tree[[j+1]][5, ncc] <- 0
            tree[[j+1]][6, ncc] <- matrnd[3,test.stat.ind]
            tree[[j+1]][7, ncc] <- length(screen.id)
            tree[[j+1]][8, ncc] <- length(matrnd[3,])
            screen.tree[[j+1]] <- cbind(screen.tree[[j+1]],list())
            screen.tree[[j+1]][[ncc]] <- screen.id
            mat <- cbind(mat, c(tree[[j+1]][-5, ncc], j+1, ncc))
          }
        }
        i <- i+1
      }
      j <- j+1
    } else{
      break
    }
  }
if(WBS) {
  list(mat=mat, matrnd=matrnd0, flagMT=flag)
  }
  else {
  list(mat=mat, matrnd=NULL, flagMT=flag)
  } #flagMT for further development
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
