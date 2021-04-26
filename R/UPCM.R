


## response function for acat
responseFun <- function(eta){
  q <- length(eta)
  eta.help <- matrix(rep(c(0,eta),each=q+1),ncol=q+1)
  eta.help[upper.tri(eta.help)] <- 0
  pi <- cumprod(c(1,exp(eta[-q])))/sum(apply(exp(eta.help),1,prod))
  pi
}

## create responses for acat from ordinal values
createResponse <- function(Y){
  c(t(model.matrix(~0+Y)[,-length(levels(Y))]))
}

#' Uncertainty in (Generalized) Partial Credit Models
#' 
#' Performs UPCM, a method to model uncertainty in (Generalized) Partial Credit Models
#' 
#' 
#' @param Y Matrix containing the ordinal item response data (as ordered factors), 
#' one row per observation, one column per item.
#' @param X Matrix containing explanatory variables which are used both for 
#' trait parameters and uncertainty parameters, one row per observation, one column per variable.
#' @param GPCM Specifies the baseline model. \code{GPCM = TRUE} results in a \code{UGPCM} while  \code{GPCM = FALSE} results in a \code{UPCM}.
#' @param Q Number of nodes to be used (per dimension) in two-dimensional Gauss-Hermite-Quadrature.
#' @param cores Number of cores to be used in parallelized computation
#' @param lambda Tuning parameter for ridge penalty on all coefficients except sigma/slope parameters. Should be small, only used to stabilize results.
#' @param se Should standard errors be computed? Standard errors are necessary for \code{\link{plot.UPCM}}. Computation is 
#' time-consuming because numerical optimization methods are used.
#' @param method Specifies optimization algorithm used , either \code{\link{nlminb}} or \code{L-BFGS-B (\link{optim})}.
#' @param ctrl.nlminb List of control arguments for optimization procedure \code{\link{nlminb}}.
#' @return 
#' \item{delta}{Matrix containing all item parameters for the UPCM pr UGPCM model, one row
#' per item, one column per category.} 
#' \item{Sigma}{2*2 covariance matrix for both random effects, namely the trait parameters theta and the
#' uncertainty parameters alpha.}
#' \item{xi}{Estimates for covariate effects on trait parameters.}
#' \item{alpha}{Estimates for covariate effects on uncertainty parameters.}
#' \item{slopes}{Estimates item slope parameters (only for \code{GPCM = TRUE}).}
#' \item{se.delta}{}
#' \item{se.xi}{Estimates of standard errors for covariate effects on trait parameters.}
#' \item{se.alpha}{Estimates of standard errors for covariate effects on uncertainty parameters.}
#' \item{se.sigma}{Estimates of standard errors for covariance parameters. Attention: First and third parameter are estimates of se for
#' both variances, the variance of theta and the variance of alpha. Second parameter is the estimate for correlation coefficient between theta and alpha, 
#' NOT of the corresponding covariance.}
#' \item{se.slopes}{Estimates of standard errors of item slope parameters (only for \code{GPCM = TRUE}).}
#' \item{delta.GPCM}{Estimates of item parameters theta in the PCM or GPCM model.}
#' \item{sigma.GPCM}{Estimate of variance of  trait parameters theta in the PCM or GPCM model.}
#' \item{slopes.GPCM}{Estimates of slope parameters in the GPCM (only for \code{GPCM = TRUE}).}
#' \item{Y}{Matrix containing the ordinal item response data, one row per obeservation, one column per item.} 
#' \item{loglik}{Marginal log-likelihood}
#' \item{coefs}{Complete vector of all estimated parameters (for internal use).}
#' \item{se.vec}{Complete vector of all estimated standard errors (for internal use).}
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://www.sg.tum.de/epidemiologie/team/schauberger/}
#' @seealso \code{\link{plot.UPCM}} \code{\link{UPCM-package}}
#' @keywords UPCM
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
#' \dontshow{
#' set.seed(1860)
#' n <- 50
#' I <- 2
#' Y <- matrix(sample(1:3, I*n, replace = TRUE), ncol = I)
#' m_upcm <- UPCM(Y = Y, cores = 1, GPCM = FALSE, se = FALSE, ctrl.nlminb = list(rel.tol = 1e-06))
#' m_upcm
#' }
UPCM <- function(Y, X = NULL, GPCM = TRUE, Q = 10, cores = 2, lambda = 1e-2, se = TRUE,
                 method = c("nlminb", "L-BFGS-B"),
                 ctrl.nlminb = list(eval.max = 200, iter.max = 150, abs.tol = 1e-8,
                 rel.tol = 1e-8, trace = 0,step.min = 0.1, x.tol = 1e-8, xf.tol = 1e-8)){
# browser()
  method <- match.arg(method)
  if(!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  
  if(!is.matrix(Y)){
    stop("Y has to be a matrix!")
  }
  
  if(!is.null(X)){
  if(!is.matrix(X)){
    stop("X has to be a matrix!")
  }}
  
  ## initalize response vector and other parameters
  y.vec <- as.factor(c(t(Y)))
  k <- length(levels(y.vec))
  q <- k-1
  n <- nrow(Y)
  I <- ncol(Y)
  if(!is.null(X)){
  p.X <- ncol(X)
  }else{
    p.X <- 0
    X <- matrix(0, 1)
  }
  
  ## create response for acat family
  response <- createResponse(y.vec)
  
  response2 <- matrix(response,nrow=q)
  response2 <- rbind(response2,1-colSums(response2))

  ## initialize starting values, 
  ## just use a partial credit model for delta
  ## and the sigma for both sigmas,
  if(GPCM){
    m.ltm <- gpcm(Y, constraint = "gpcm")
    if(m.ltm$convergence == 1){
      stop("Estimation failed for basic GPCM model (without uncertainty)")
    }

    coef.ltm <- coef(m.ltm)
    sigma1.start <- coef.ltm[1,q+1]
    slope.start <- coef.ltm[-1,q+1]/coef.ltm[1,q+1]
    slope.start[slope.start<0] <- -slope.start[slope.start<0] 
    delta.start <- c(t(coef.ltm[,-(q+1)]*coef.ltm[,q+1]))
    sigma.start <- c(sigma1.start^2*0.9,0,0.1)
    betaX.start <- betaU.start <- rep(0,p.X)
    alpha.start <- c(delta.start,betaX.start,betaU.start,sigma.start, slope.start)
    p.delta <- I*q
    p.slope <- I-1
    p.fix <- p.delta + 2*p.X + p.slope 
  }else{
    m.ltm <- gpcm(Y, constraint = "1PL")
    if(m.ltm$convergence == 1){
      stop("Estimation failed for basic PCM model (without uncertainty component)")
    }
    coef.ltm <- coef(m.ltm)
    sigma1.start <- coef.ltm[1,q+1]
    delta.start <- c(t(coef.ltm[,-(q+1)])*sigma1.start)
    sigma.start <- c(sigma1.start^2*0.9,0,0.1)
    betaX.start <- betaU.start <- rep(0,p.X)
    alpha.start <- c(delta.start,betaX.start,betaU.start,sigma.start)
    p.delta <- I*q
    p.slope <- 0
    slope.start <- c()
    p.fix <- p.delta + 2*p.X
  }
  
  ## get nodes and weights
  her_poly <- gauss.quad(Q, "hermite")
  nodes <- her_poly$nodes
  weights <- her_poly$weights * exp(nodes^2)
  weights <- (weights %*% t(weights)) * (dnorm(nodes)%*%t(dnorm(nodes)))
  # node.probs = dnorm(nodes)%*%t(dnorm(nodes))

## lower bounds for variance parameters
## upper and lower bound for correlation
  l.bound <- rep(-Inf,length(alpha.start))
  l.bound[(p.fix-p.slope +1):(p.fix-p.slope+3)] <- c(1e-2, -1+1e-2, 1e-2)
  
  u.bound <- rep(Inf,length(alpha.start))
  u.bound[p.fix-p.slope + 2] <- 1-1e-2

   # browser()

conv <- FALSE
try.grad <- TRUE
convhelp <- 0.1
count.conv <- 0

if(GPCM){
  score_exec <- scoreUGPCM2
  loglik_exec <- loglikUGPCM
  loglik_exec2 <- loglikUGPCM2
  if(p.X == 0){
    loglik_exec <- loglikUGPCMnoX
    loglik_exec2 <- loglikUGPCM2noX
  }
}else{
  score_exec <- scoreUPCM2
  loglik_exec <- loglikUPCM
  loglik_exec2 <- loglikUPCM2
  if(p.X == 0){
    loglik_exec <- loglikUPCMnoX
    loglik_exec2 <- loglikUPCM2noX
  }
}


  while(!conv & count.conv<5){
    alpha.start <- (0.9+convhelp)*alpha.start
    s.start <- try(score_exec(alpha.start,Q = Q, q = q, I = I, n = n, Y = response, p = length(alpha.start),
                          GHweights = weights, GHnodes = nodes, 
                          pX = p.X, X = X,
                          cores = cores,lambda=lambda),silent = TRUE)
    l.start <- try(loglik_exec(alpha.start,Q = Q, q = q, I = I, n = n, Y = response, p = length(alpha.start),
                          GHweights = weights, GHnodes = nodes, 
                          pX = p.X, X = X,
                          cores = cores,lambda=lambda), silent = TRUE)
    if(inherits(s.start, "try-error")|inherits(l.start, "try-error")){
      s.start <- NA
    }
    
    conv <- (!is.na(sum(s.start))) & (!is.na(l.start))
    convhelp <- 0
    count.conv <- count.conv + 1
    if(count.conv==5){
      try.grad <- FALSE
      
      if(GPCM){
        alpha.start <- c(delta.start,betaX.start,betaU.start,sigma.start, slope.start)
    }else{
        alpha.start <- c(delta.start,betaX.start,betaU.start,sigma.start)
      }
    }
    
  }

## change scaling for variance and correlation parameters
 par.scale.null <- rep(1,length(alpha.start))

 
 ### If method is nlminb
if(method=="nlminb"){
  if(try.grad){
  m.opt <- try(nlminb(start = alpha.start, objective = loglik_exec, gradient = score_exec, 
                  Q = Q, q = q, I = I, n = n, Y = response, pall = length(alpha.start),
                  GHweights = weights, GHnodes = nodes, 
                  pX = p.X, X = X,
                  cores = cores, lambda = lambda, lower = l.bound, upper = u.bound, 
                  control=ctrl.nlminb), silent = TRUE)
  
  if( (class(m.opt) == "try-error")){
  try.grad <- FALSE
  }else{
    print("Fit with analytic gradient was successful!")
  }
  }
  if(!try.grad){
    par.scale <- par.scale.null
    
    
    m.opt <- try(nlminb(start = alpha.start, objective = loglik_exec2, gradient = NULL, 
                        Q = Q, q = q, I = I, n = n, Y = response2, pall = length(alpha.start),
                        GHweights = weights, GHnodes = nodes, 
                        pX = p.X, X = X,
                        cores = cores, lambda = lambda, lower = l.bound, upper = u.bound,
                        control=ctrl.nlminb), silent = TRUE)
    
  }
}
 
 ### If method is Rvmmin
 if(method=="L-BFGS-B"){
   if(try.grad){
     m.opt <- try(optim(par = alpha.start, fn = loglik_exec, gr = score_exec, 
                         Q = Q, q = q, I = I, n = n, Y = response, pall = length(alpha.start),
                         GHweights = weights, GHnodes = nodes, 
                         pX = p.X, X = X,
                         cores = cores, lambda = lambda, lower = l.bound, upper = u.bound, 
                         method = "L-BFGS-B"), silent = TRUE)
     
     
     if( (class(m.opt) == "try-error")){
       try.grad <- FALSE
     }else{
       print("Fit with analytic gradient was successful!")
     }
   }
   if(!try.grad){
     
     m.opt <- try(optim(par = alpha.start, fn = loglik_exec2, gr = NULL, 
                        Q = Q, q = q, I = I, n = n, Y = response2, pall = length(alpha.start),
                        GHweights = weights, GHnodes = nodes, 
                        pX = p.X, X = X,
                        cores = cores, lambda = lambda, lower = l.bound, upper = u.bound, 
                        method = "L-BFGS-B"), silent = TRUE)
     
   }
 }
 
########################
  
delta <- Sigma <- loglik <- xi <- beta <- slopes <- coefs <- NA
se.delta <- se.beta <- se.xi <- se.slopes <- se.sigma <- se.vec <-  NA

if( (class(m.opt) != "try-error")){
  ## extract results and prepare return
  coefs <- m.opt$par


  loglik <- -loglik_exec2(coefs,  Q = Q, q = q, I = I, n = n, Y = response2, pall = length(alpha.start),
                         GHweights = weights, GHnodes = nodes, pX = p.X, X = X,  cores = cores, lambda = 0)

  
if(se){

    hess <- optimHess(par = coefs, fn = loglik_exec2, gr = NULL, 
                      Q = Q, q = q, I = I, n = n, Y = response2, pall = length(alpha.start),
                      GHweights = weights, GHnodes = nodes, 
                      pX = p.X, X = X,
                      cores = cores, lambda = lambda)


  hess.inv <- solve(hess)
  se.vec <- sqrt(diag(hess.inv))

  se.delta <- matrix(se.vec[1:p.delta],byrow=TRUE,ncol=q)
  if(p.X>0){
    se.beta <- se.vec[(p.delta+p.X+1):(p.delta+2*p.X)]
    se.xi <- se.vec[(p.delta+1):(p.delta+p.X)]
  }
  se.sigma <- se.vec[(p.delta+ 2*p.X + 1):(p.delta+2*p.X+3)]

  if(p.X>0){
    names(se.xi) <- names(se.beta) <- colnames(X)
  }
  colnames(se.delta) <- paste("Catgr",1:q,sep=".")
  rownames(se.delta) <- colnames(Y)
  names(se.sigma) <- c("sigma.theta", "cor", "sigma.alpha")
  
  if(GPCM){
    se.slopes <- se.vec[(p.delta+2*p.X+4):(p.delta+2*p.X+I+2)]
    names(se.slopes) <- paste0("Item",2:I)
  }
}
  
  delta <- matrix(coefs[1:p.delta],byrow=TRUE,ncol=q)
  
  
  colnames(delta) <-  paste("Catgr",1:q,sep=".")
  rownames(delta) <-  colnames(Y)
  
  if(p.X>0){
    xi <- coefs[(p.delta+1):(p.delta+p.X)]
    beta <- coefs[(p.delta+p.X+1):(p.delta+2*p.X)]
    names(xi) <- names(beta) <-  colnames(X)
  }
  
  if(GPCM){
    slopes <- tail(coefs, I-1)
    names(slopes) <- paste0("Item",2:I)
  }
  

  
  if(GPCM){
    sigma <- coefs[(length(coefs)-1-I):(length(coefs)-I+1)]
  }else{
    sigma <- coefs[(length(coefs)-2):(length(coefs))]
  }
  if(sigma[3]==0){
    Sigma <- matrix(c(sigma[1],0,0,0), ncol = 2)
  }else{
    Sigma <- matrix(c(sigma[1], sigma[2]*sqrt(sigma[1])*sqrt(sigma[3]),
                      sigma[2]*sqrt(sigma[1])*sqrt(sigma[3]), sigma[3]),ncol=2)
  }


colnames(Sigma) <- rownames(Sigma) <- c("theta","alpha")


  if(!GPCM){
    slopes <- se.slopes <- NA
  }
}else{
  stop("Estimation did not converge!")
}

delta.GPCM <- matrix(delta.start, ncol=q, byrow=TRUE)
colnames(delta.GPCM) <- paste("Catgr",1:q,sep=".")
rownames(delta.GPCM) <- colnames(Y)

sigma.GPCM <- sigma1.start^2

if(GPCM){
  slopes.GPCM <- slope.start
}else{
  slopes.GPCM <- NA
}


  ret.list <- list(delta = delta, Sigma = Sigma, xi = xi, alpha = beta, slopes = slopes,
                   se.delta = se.delta, se.xi = se.xi, se.alpha = se.beta, se.sigma = se.sigma, se.slopes = se.slopes,
                   delta.GPCM = delta.GPCM, sigma.GPCM = sigma.GPCM, slopes.GPCM = slopes.GPCM,
                   Y = Y, loglik = loglik, coefs = coefs, se.vec = se.vec)

  class(ret.list) <- "UPCM"
  if(GPCM){class(ret.list) <- c("UGPCM","UPCM")}
  
  return(ret.list)
}

scoreUPCM2 <- function(alpha, Y, Q, q, n, I, pall,
                       GHweights, GHnodes, pX, X, cores, lambda){
  # browser()
  sigma <- alpha[(length(alpha)-2):length(alpha)]
  alpha <- head(alpha,length(alpha)-3)
  
  if(!pX == 0){
  grad.fixed <- scoreUPCM(alpha, Y, Q, q, n, I, pall-3,
                          GHweights, GHnodes,  pX, X, cores, sigma, lambda)
  
  betaX <- alpha[(length(alpha)-2*pX+1):(length(alpha)-pX)]
  betaU <- alpha[(length(alpha)-pX+1):(length(alpha))]
  delta <- alpha[1:(I*q)]
  
  grad.sigma <- grad(loglikUPCM4, sigma, Q = Q, q = q, I = I, n = n, Y = Y, 
                     GHweights = GHweights, GHnodes = GHnodes, X = X,
                     cores = cores, betaX = betaX, betaU = betaU, delta = delta,
                     lambda = lambda)
  
  }else{

    grad.fixed <- scoreUPCMnoX(alpha, Y, Q, q, n, I, pall-3,
                            GHweights, GHnodes,  cores, sigma, lambda)
    
    delta <- alpha[1:(I*q)]
    # browser()
    grad.sigma <- grad(loglikUPCM4noX, sigma, Q = Q, q = q, I = I, n = n, Y = Y, 
                       GHweights = GHweights, GHnodes = GHnodes, 
                       cores = cores, delta = delta, lambda = lambda)
  }
  
  return(c(grad.fixed, grad.sigma))
}



scoreUGPCM2 <- function(alpha, Y, Q, q, n, I, pall,
                       GHweights, GHnodes, pX, X, cores, lambda){

  sigma <- alpha[(length(alpha)-1-I):(length(alpha)-I+1)]
  alpha <-  alpha[-((length(alpha)-1-I):(length(alpha)-I+1))]

  if(!pX == 0){
  grad.fixed <- scoreUGPCM(alpha, Y, Q, q, n, I, pall-3,
            GHweights, GHnodes,  pX, X, cores, sigma, lambda)
  
  betaX <- alpha[(length(alpha)-2*pX+2-I):(length(alpha)-pX-I+1)]
  betaU <- alpha[(length(alpha)-pX+2-I):(length(alpha)-I+1)]
  delta <- alpha[1:(I*q)]
  slopes <- c(1, tail(alpha, I-1))
  
  grad.sigma <- grad(loglikUGPCM4, sigma, Q = Q, q = q, I = I, n = n, Y = Y, 
                     GHweights = GHweights, GHnodes = GHnodes, 
                     X = X,
                     cores = cores, betaX = betaX, betaU = betaU, delta = delta, slopes = slopes,
                     lambda = lambda)
  
  }else{
    grad.fixed <- scoreUGPCMnoX(alpha, Y, Q, q, n, I, pall-3,
                             GHweights, GHnodes, cores, sigma, lambda)
    delta <- alpha[1:(I*q)]
    slopes <- c(1, tail(alpha, I-1))
    
    grad.sigma <- grad(loglikUGPCM4, sigma, Q = Q, q = q, I = I, n = n, Y = Y, 
                       GHweights = GHweights, GHnodes = GHnodes, 
                       cores = cores,  delta = delta, slopes = slopes,
                       lambda = lambda)
}


  return(c(grad.fixed[1:(q*I+2*pX)], grad.sigma,grad.fixed[(q*I+2*pX+1):length(grad.fixed)]))
}