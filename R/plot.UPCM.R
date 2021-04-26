#' Plot function for UPCM
#' 
#' Plot function for a \code{UPCM} or a \code{UGPCM} object. Plots show coefficient estimates together with
#' confidence intervals displayed as star plots. 
#' 
#' @usage \method{plot}{UPCM}(x, sig = 0.05, KIfactor = 0.9, xlim, ylim, \dots)
#' @param x \code{UPCM} object
#' @param sig Significance level for confidence intervals, default is \code{sig = 0.05}. 
#' @param KIfactor Parameter to regulate the shape of the resulting star.
#' @param xlim See \code{xlim} in \code{\link{plot.default}}.
#' @param ylim See \code{ylim} in \code{\link{plot.default}}.
#' @param ... Further plot arguments.
#' @return No return value, called for side effects
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://www.sg.tum.de/epidemiologie/team/schauberger/}
#' @seealso \code{\link{UPCM}}
#' @references Tutz, Gerhard and Schauberger, Gunther (2020): Uncertainty in Latent Trait Models, 
#' \emph{Applied Psychological Measurement}, \url{https://journals.sagepub.com/doi/abs/10.1177/0146621620920932?journalCode=apma}
#' @examples
#' \donttest{
#' data(tenseness)
#' 
#' Y <- data.matrix(tenseness[,1:4])
#' X <- model.matrix(~ Gender + Age, data = tenseness)[,-1]
#' 
#' m_upcm <- UPCM(Y = Y, X = X, cores = 2, GPCM = FALSE)
#' m_upcm
#' plot(m_upcm)
#' }
plot.UPCM <- function(x, sig = 0.05, KIfactor = 0.9, xlim, ylim, ...){
  
  quant <- qnorm(1-sig/2)
  
  if(is.na(x$xi[1])|is.na(x$alpha[1])){
    stop("Plotting is only possible if covariates are used (both for 
         the location and heterogeneity effect (i.e. if argument X is specified)!")
  }
  # browser()
  xi <- x$xi
  xi.KI <- cbind(xi-quant*x$se.xi,xi+quant*x$se.xi)
  xi <- xi
  
  alpha <- x$alpha
  alpha.KI <- exp(cbind(alpha-quant*x$se.alpha,alpha+quant*x$se.alpha))
  alpha <- exp(alpha)
  
  
  
  if(missing(ylim)){
    ylim <- range(c(1,alpha.KI))
  }
  
  if(missing(xlim)){
    xlim <- range(c(0,xi.KI))
  }
  
  plot(xi,alpha,pch=16,xlim=xlim,ylim=ylim,
       xlab=expression(xi),ylab=expression(exp(alpha)),...)
  
  p.X <- length(xi)
  
  label.x <- label.y <- c()

  for(i in 1:p.X){
    
    x <- c(xi.KI[i,1],xi.KI[i,1]+(xi[i]-xi.KI[i,1])*(KIfactor),xi[i],xi[i]+(xi[i]-xi.KI[i,1])*(1-KIfactor),
           xi.KI[i,2],xi[i]+(xi[i]-xi.KI[i,1])*(1-KIfactor),xi[i],xi.KI[i,1]+(xi[i]-xi.KI[i,1])*(KIfactor),
           xi.KI[i,1])
    
    y <- c(alpha[i],alpha.KI[i,1]+(alpha[i]-alpha.KI[i,1])*(KIfactor),alpha.KI[i,1],alpha.KI[i,1]+(alpha[i]-alpha.KI[i,1])*(KIfactor),alpha[i],
           alpha[i]+(alpha[i]-alpha.KI[i,1])*(1-KIfactor),alpha.KI[i,2],alpha[i]+(alpha[i]-alpha.KI[i,1])*(1-KIfactor),alpha[i])
    
    polygon(x,y,col=grey(0.9))
    label.x <- c(label.x,x[6])
    label.y <- c(label.y,y[6])
  }
  points(xi,alpha,pch=16)
  abline(h=1,lty=2,lwd=2,col="gray")
  abline(v=0,lty=2,lwd=2,col="gray")
  
  text(label.x,label.y,labels=names(xi),adj=c(-0.1,-0.1))
}
