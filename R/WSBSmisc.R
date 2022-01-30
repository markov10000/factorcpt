#' @keywords internal
calcul.cusum <-
function(x,dw=0, alpha=0.5,agg.std.op=0){
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation

    leftsums <- apply(x,1,cumsum)


    t <- 1:(n-1)
    tnt <- t*(n-t)
    cusum<-c()
    # constructing CUSUM matrix
    #rightsums <- t(leftsums[n,]-t(leftsums))
    #cusum$cusum0 <- t((rightsums[t,]/(n-t) - leftsums[t,]/t)*sqrt(n)*(tnt/n^2)^alpha)
    cusum$cusum <- t((leftsums[t,]*n - t%o%leftsums[n,])*(tnt)^(alpha-1) * n^(-2*alpha+0.5))
    cusum$SN <- rep(1,p)
    if(agg.std.op==0){#return(cusum)
    }else if (agg.std.op==1){#Self Normalize
        if(dw*2>n) {dw = 0; warning("dw*2>n; set dw=0")}
        rightsums <- t(leftsums[n,]-t(leftsums))
        #cpt at cusum$SNk+1
        cusum$SNk <- dw + apply(cusum$cusum[,-c(1:dw,(n-dw+1):n)],1,which.max)
        #        std<-function(x)mad(diff(x)/sqrt(2))
        for (pii in 1:p){
            SNn <- cusum$SNk[pii]
            SNt <- 1:(SNn-1)
            SNcumsum1 <- (leftsums[SNt,pii]*SNn - SNt%o%leftsums[SNn,pii]) * SNn^(-1)
            SNt2 <- 1:(n-SNn-1)
            SNcumsum2 <- (rightsums[(n-1):(SNn+1),pii]*(n-SNn) - SNt2%o%rightsums[SNn,pii]) * (n-SNn)^(-1)
            cusum$SN[pii] <- sum(rbind(SNcumsum1^2,SNcumsum2^2))^0.5/n

        }
    }else if (agg.std.op==2){#chi-square
        cusum$SN <- rowMeans(abs(x))
    }else if (agg.std.op==3){#normal diff-mad
        std<-function(x) mad(abs(diff(x)-mad(diff(x))))
        cusum$SN <- apply(x,1,std)
    }else if (agg.std.op==4){#normal mad
        std<-function(x) mad(abs(x-mad(x)))
        cusum$SN <- apply(x,1,std)
    }else if (agg.std.op==5){#sd
        cusum$SN <- apply(x,1,sd)
        }
    return(cusum$cusum/cusum$SN)
}

#' @keywords internal
calcul.matrnd <- function(y,dw,screen.id,matrnd,normf,agg.std.op,do.parallel){
    n <- dim(y)[1]
    T <- dim(y)[2]
    M <- dim(matrnd)[2]


    if(do.parallel&&length(screen.id)>500){#do.parallel only if dimension is big
        calcul.cusum=calcul.cusum
        mat34 <- foreach::foreach(MM=(1:M), .combine=cbind, .packages=c()) %dopar% {
            len <- matrnd[2,MM]-matrnd[1,MM] + 1
            stat<-calcul.cusum(y[screen.id, matrnd[1,MM]:matrnd[2,MM],drop=F],dw=dw,agg.std.op=agg.std.op)
            norm.stat <- normf(stat[,-c((1:dw), (len-dw+1):len,drop=F)])
            max.stat <- max(norm.stat)
            hat.chp <- matrnd[1,MM] + dw + min(which(norm.stat==max.stat)) -1
            c(max.stat,hat.chp)
        }
    }else{
        mat34 <- matrix(0,2,M)
        for (MM in 1:M){
            len <- matrnd[2,MM]-matrnd[1,MM] + 1
            stat<- calcul.cusum(y[screen.id, matrnd[1,MM]:matrnd[2,MM],drop=F],dw=dw,agg.std.op=agg.std.op)
            norm.stat <- normf(stat[,-c((1:dw), (len-dw+1):len),drop=F])
            mat34[1,MM] <- max(norm.stat)
            mat34[2,MM] <- matrnd[1,MM] + dw + min(which(norm.stat==mat34[1,MM])) -1   #hat.chp
        }
    }
    return(mat34)
}

#' @keywords internal
select.thr <-
function(x,dw=NULL, agg.std.op=0,num.boots=1,select=seq(0.99,1,0.001)){
    len <- dim(x)[2]
    if(is.null(dw)) dw <- round(min(log(len)^2, len^(6/7)/4))
    cdist <- c()
    for (Mi in 1:num.boots){
        dist<- abs(calcul.cusum(x[,sample(len)],dw=dw,agg.std.op=agg.std.op))
        cdist <- rbind(cdist, dist)
        }
    return(quantile(apply(cdist[,-c((1:dw), (len-dw+1):len),drop=F],1,max),select))
}

means.between.cpt <- function(x, cpt=NULL) {
    k<-dim(x)[1]
    n<-dim(x)[2]
    len.cpt <- length(cpt)

    if (len.cpt) {
        if(min(cpt)<1 || max(cpt)>=n) stop("change-points should be between 1 and n-1")
        cpt <- sort(cpt)
    }

    s <- e <- rep(0, len.cpt+1)
    s[1] <- 1
    e[len.cpt+1] <- n
    if (len.cpt) {
        s[2:(len.cpt+1)] <- cpt+1
        e[1:len.cpt] <- cpt
    }

    means <- matrix(0,k, len.cpt+1)
    for (i in 1:(len.cpt+1)) means[,i] <- rowMeans(x[,s[i]:e[i],drop=F])

    return(t(apply(means, 1, rep, e-s+1)))
}


SSIC <- function(x,cpt.trylist,SSIC.pen.fun=function(y)(sqrt(y))){
    p <- dim(x)[1]
    n <- dim(x)[2]
    len.cpt <- length(cpt.trylist)
    SIC.pen.value <- SSIC.pen.fun(n,p)
    
    w.cpt<-c()
    if(!len.cpt) {
        w.cpt$cpt<- as.numeric()
        w.cpt$no.cpt <- 0
        w.cpt$max.diff.pen <- rep(0,p)
        return(w.cpt)
    }

    if(SIC.pen.value<=0){#all candidates are returned
        w.cpt$cpt<- cpt.trylist
        w.cpt$no.cpt <- len.cpt
        w.cpt$max.diff.pen <- rep(0,p)
        return(w.cpt)
    }


    #ic part

        w.cpt$ic.curve <- matrix(0,p,len.cpt+1)

        if(len.cpt) for(i in len.cpt:1){
            min.log.lik <- n/2 * log(rowSums((x -  means.between.cpt(x,cpt.trylist[1:i]))^2)/(n-1)) 
            w.cpt$ic.curve[,i+1] <- min.log.lik +  i * SIC.pen.value
        }
        w.cpt$ic.curve[,1] <- n/2 * log(apply(x,1,var))
        diff.penalty <- t(diff(t(w.cpt$ic.curve))<0)  ##decreasing
        w.cpt$dim.cpt <- colSums(matrix(diff.penalty, nrow(diff.penalty)))
        w.cpt$no.cpt <- which(c(w.cpt$dim.cpt,0)==0)[1] - 1
        if(w.cpt$no.cpt>0){w.cpt$cpt <- cpt.trylist[1:w.cpt$no.cpt]}else{w.cpt$cpt <- numeric(0)}

return(w.cpt)
}
