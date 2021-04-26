print.UPCM <- function(x, ...){
  
  cat("Output UPCM estimation:","\n")
  
  cat("---","\n")
  
  cat("\n")
  
  

  q <- ncol(x$delta)
  I <- nrow(x$delta)
  p.X <- length(x$xi)
  
  if(is.na(x$xi[1])){
    p.X <- 0
  }
  
  delta.mat <- cbind(c(t(x$delta)), c(t(x$se.delta)))
  delta.mat <- cbind(delta.mat, delta.mat[,1]/delta.mat[,2])
  delta.mat <- cbind(delta.mat, 2*pnorm(-abs(delta.mat[, 3])))
  colnames(delta.mat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  rownames(delta.mat) <- paste(rep(rownames(x$delta),each=q),rep(1:q,I),sep=":")
  
  
  if(p.X > 0){
    xi.mat <- cbind(c(t(x$xi)), c(t(x$se.xi)))
    xi.mat <- cbind(xi.mat, xi.mat[,1]/xi.mat[,2])
    xi.mat <- cbind(xi.mat, 2*pnorm(-abs(xi.mat[, 3])))
    colnames(xi.mat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
    rownames(xi.mat) <- names(x$xi)
    
    
    
    alpha.mat <- cbind(c(t(x$alpha)), c(t(x$se.alpha)))
    alpha.mat <- cbind(alpha.mat, alpha.mat[,1]/alpha.mat[,2])
    alpha.mat <- cbind(alpha.mat, 2*pnorm(-abs(alpha.mat[, 3])))
    colnames(alpha.mat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
    rownames(alpha.mat) <- names(x$alpha)
  }
  
  cat("Estimates of item parameters delta","\n")
  printCoefmat(delta.mat, ...)
  
  cat("\n")
  
  if(p.X > 0){
    cat("Estimates of trait effects","\n")
    printCoefmat(xi.mat, ...)
    
    cat("\n")
    
    cat("Estimates of uncertainty effects","\n")
    printCoefmat(alpha.mat, ...)
    
    cat("\n")
  }
  
  cat("Estimates of covariance matrix Sigma","\n")
  print(x$Sigma, ...)
  
  
  invisible(x)
}