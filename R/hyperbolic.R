
########################################
###  Hyperbolic Distributions for R  ###
########################################


# Function to fit a hyperbolic distribution to data.

# Different Starting Values (start.values) Include:
    # Approach based on Barndorff-Nielsen 1977, p.407 ("bn")
    # Fitted Normal distribution ("fn")
    # Skew Laplace ("sl")
    # All of the Above ("all")
    # User Specified ("us")


fit.hyperb <- function(x,freq=NULL,breaks=NULL,theta.start=NULL,
            start.values=c("all"),method=c("BFGS"),
            max.lh=FALSE,plots=TRUE,leave.on=FALSE,
            controlbfgs=list(maxit=100),
            controlnm=list(maxit=1000),
            maxitnlm=1500,...){
  
    xname <- deparse(substitute(x))
    if(!is.null(freq)){
        if(length(freq)!=length(x)) 
        stop("vectors x and freq are not of the same length")
        x <- rep(x,freq)
    }

    
    if(any(start.values=="all")) start.values=c("bn","fn","sl")
    if(any(method=="all")) method=c("BFGS","NLM","Nelder-Mead")

    if(!is.null(theta.start)&&all(start.values!="us")) 
        start.values <- c(start.values,"us")
    if(any(start.values=="us")&&is.null(theta.start)) 
        stop("theta.start not specified")
    
    if(!is.null(theta.start)){
        if(length(theta.start)!=4) 
           stop("theta.start must contain 4 values")
        if(theta.start[2]<=0) 
           stop("Zeta in theta.start must be greater than zero")
        if(theta.start[3]<=0) 
           stop("Delta in theta.start must be greater than zero")
    }

    start.values <- sort(start.values)
    method <- sort(method)

    sv.names <- vector("character",length(start.values))
    if(any(start.values=="bn")) 
         sv.names[start.values=="bn"] <- "Barndorff-Nielsen 1977"
    if(any(start.values=="fn")) 
         sv.names[start.values=="fn"] <- "Fitted Normal"
    if(any(start.values=="sl")) 
         sv.names[start.values=="sl"] <- "Skew Laplace"
    if(any(start.values=="us")) 
         sv.names[start.values=="us"] <- "User Specified"

    output.bfgs <- NULL
    output.nm <- NULL
    output.nlm <- NULL

    x <- as.numeric(na.omit(x))

    if(is.null(breaks)){
        breaks <- hist(x,plot=FALSE,right=FALSE)$breaks
        }

    midpoints <- hist(x,breaks,plot=FALSE,right=FALSE)$mids
    x.freq <- hist(x,breaks,plot=FALSE,right=FALSE)$counts
    emp.dens <- hist(x,breaks,plot=FALSE,right=FALSE)$density
    
    emp.dens <- ifelse(!is.finite(log(emp.dens)),NA,emp.dens)
    max.index <- order(emp.dens,na.last=FALSE)[length(emp.dens)]

    if(length(na.omit(emp.dens[1:max.index]))<2|
           length(na.omit(emp.dens[max.index:length(emp.dens)]))<2){
         if(any(start.values=="bn")|any(start.values=="sl")) 
        stop("not enough breaks to estimate asymptotes to log-density")
    }

    left.asymptote <- lm(log(emp.dens)[1:max.index] 
        ~ midpoints[1:max.index])$coef
    right.asymptote <- lm(log(emp.dens)[max.index:length(emp.dens)]
        ~ midpoints[max.index:length(emp.dens)])$coef 

    start.matrix <- matrix(nrow=length(start.values),ncol=4)
    i <- 1

    if(any(start.values=="bn")){
        hyp.phi <- as.numeric(left.asymptote[2])
        hyp.gamma <- as.numeric(-right.asymptote[2])
        hyp.mu <- as.numeric(-(left.asymptote[1]-
                    right.asymptote[1])/ 
                                (left.asymptote[2]-
                    right.asymptote[2])) 
        intersection.value <- as.numeric(left.asymptote[1]+
                    hyp.mu*left.asymptote[2]) 
        log.modal.dens <- log(max(emp.dens,na.rm=TRUE))
        hyp.zeta <- intersection.value - log.modal.dens
        if(hyp.zeta<=0){
            hyp.zeta <- 0.1 # This set arbitrarily
            }
        hyp.delta <- hyp.zeta/sqrt(hyp.phi*hyp.gamma) 
        hyp.pi <- as.numeric(hyperb.change.pars(3,1,
            c(hyp.phi,hyp.gamma,hyp.delta,hyp.mu))[1])
        theta.start.bn <- c(hyp.pi,log(hyp.zeta),
            log(hyp.delta),hyp.mu)
        start.matrix[i,1:4] <- theta.start.bn
        i <- i + 1
        }

    if(any(start.values=="fn")){
        hyp.nu <- as.numeric(midpoints[max.index])
        hyp.mu <- mean(x)
        hyp.delta <- sd(x)
        hyp.pi <- (hyp.nu - hyp.mu)/hyp.delta
        hyp.zeta <- 1 + hyp.pi^2
        theta.start.fn <- c(hyp.pi,log(hyp.zeta),
                    log(hyp.delta),hyp.mu)
        start.matrix[i,1:4] <- theta.start.fn
        i <- i + 1
        }

    if(any(start.values=="sl")){
        llsklp <- function(theta){
            -sum(log(dskewlap(x,theta)))
            }

        # SkewLP optimisation:
        skew.alpha <- as.numeric(1/left.asymptote[2]) 
        skew.beta <- as.numeric(-1/right.asymptote[2]) 
        skew.mu <- as.numeric(midpoints[max.index]) 
        theta.start.sl <- c(skew.alpha, skew.beta, skew.mu) 
        skewlp.optim <- optim(theta.start.sl,llsklp,NULL,
            method="Nelder-Mead",hessian=FALSE)
        
        hyp.phi <- skewlp.optim$par[1]
        hyp.gamma <- skewlp.optim$par[2]
        hyp.delta <- 0.1 # Taking delta to be small
        hyp.mu <- skewlp.optim$par[3]
        hyp.pi <- as.numeric(hyperb.change.pars(3,1,
            c(hyp.phi,hyp.gamma,hyp.delta,hyp.mu))[1])
        hyp.zeta <- as.numeric(hyperb.change.pars(3,1,
            c(hyp.phi,hyp.gamma,hyp.delta,hyp.mu))[2])
        theta.start.sl <- c(hyp.pi,log(hyp.zeta),
            log(hyp.delta),hyp.mu)
        start.matrix[i,1:4] <- theta.start.sl
        i <- i + 1
        }

    if(any(start.values=="us")){
        hyp.pi <- theta.start[1]
        hyp.zeta <- theta.start[2]
        hyp.delta <- theta.start[3]
        hyp.mu <- theta.start[4]
        theta.start.us <- c(hyp.pi,log(hyp.zeta),
            log(hyp.delta),hyp.mu)
        start.matrix[i,1:4] <- theta.start.us
        }

    llfunc <- function(theta){
        K.nu <- besselK(exp(theta[2]),nu=1)
        -sum(log(dhyperb(x,theta,K.nu=K.nu,log.pars=TRUE)))
        }
    
    if(any(method=="BFGS")){
        bfgs <- matrix(rep(NA,7*length(start.values)),
            nrow=length(start.values),ncol=7)
        i <- 1
        if(any(start.values=="bn")){
            optim.bn.bfgs <- optim(theta.start.bn,llfunc,
                NULL,method="BFGS",hessian=FALSE,
                control=controlbfgs)
            bfgs[i,1:4] <- optim.bn.bfgs$par
            bfgs[i,5] <- -optim.bn.bfgs$value
            bfgs[i,6] <- optim.bn.bfgs$convergence
            bfgs[i,7] <- optim.bn.bfgs$counts[2]
            i <- i + 1
            }
        if(any(start.values=="fn")){
            optim.fn.bfgs <- optim(theta.start.fn,llfunc,
                NULL,method="BFGS",hessian=FALSE,
                control=controlbfgs)
            bfgs[i,1:4] <- optim.fn.bfgs$par
            bfgs[i,5] <- -optim.fn.bfgs$value
            bfgs[i,6] <- optim.fn.bfgs$convergence
            bfgs[i,7] <- optim.fn.bfgs$counts[2]
            i <- i + 1
            }
        if(any(start.values=="sl")){
            optim.sl.bfgs <- optim(theta.start.sl,llfunc,
                NULL,method="BFGS",hessian=FALSE,
                control=controlbfgs)
            bfgs[i,1:4] <- optim.sl.bfgs$par
            bfgs[i,5] <- -optim.sl.bfgs$value
            bfgs[i,6] <- optim.sl.bfgs$convergence
            bfgs[i,7] <- optim.sl.bfgs$counts[2]
            i <- i + 1
            }
        if(any(start.values=="us")){
            optim.us.bfgs <- optim(theta.start.us,llfunc,
                NULL,method="BFGS",hessian=FALSE,
                control=controlbfgs)
            bfgs[i,1:4] <- optim.us.bfgs$par
            bfgs[i,5] <- -optim.us.bfgs$value
            bfgs[i,6] <- optim.us.bfgs$convergence
            bfgs[i,7] <- optim.us.bfgs$counts[2]
            }

        output.bfgs <- data.frame(Pi=bfgs[,1],
            Zeta=exp(bfgs[,2]),Delta=exp(bfgs[,3]),
            Mu=bfgs[,4],LogLikelihood=bfgs[,5],
            Convergence=paste("Optim",bfgs[,6]),
            Iterations=bfgs[,7],
            row.names=paste("BFGS",sv.names))
        }

    if(any(method=="Nelder-Mead")){
        nm <- matrix(rep(NA,7*length(start.values)),
            nrow=length(start.values),ncol=7)
        i <- 1
        if(any(start.values=="bn")){
            optim.bn.nm <- optim(theta.start.bn,llfunc,
                NULL,method="Nelder-Mead",
                hessian=FALSE,
                control=controlnm)
            nm[i,1:4] <- optim.bn.nm$par
            nm[i,5] <- -optim.bn.nm$value
            nm[i,6] <- optim.bn.nm$convergence
            nm[i,7] <- optim.bn.nm$counts[1]
            i <- i + 1
            }
        if(any(start.values=="fn")){
            optim.fn.nm <- optim(theta.start.fn,llfunc,
                NULL,method="Nelder-Mead",
                hessian=FALSE,
                control=controlnm)
            nm[i,1:4] <- optim.fn.nm$par
            nm[i,5] <- -optim.fn.nm$value
            nm[i,6] <- optim.fn.nm$convergence
            nm[i,7] <- optim.fn.nm$counts[1]
            i <- i + 1
            }
        if(any(start.values=="sl")){
            optim.sl.nm <- optim(theta.start.sl,llfunc,
                NULL,method="Nelder-Mead",
                hessian=FALSE,
                control=controlnm)
            nm[i,1:4] <- optim.sl.nm$par
            nm[i,5] <- -optim.sl.nm$value
            nm[i,6] <- optim.sl.nm$convergence
            nm[i,7] <- optim.sl.nm$counts[1]
            i <- i + 1
            }
        if(any(start.values=="us")){
            optim.us.nm <- optim(theta.start.us,llfunc,
                NULL,method="Nelder-Mead",
                hessian=FALSE,
                control=controlnm)
            nm[i,1:4] <- optim.us.nm$par
            nm[i,5] <- -optim.us.nm$value
            nm[i,6] <- optim.us.nm$convergence
            nm[i,7] <- optim.us.nm$counts[1]
            }
        output.nm <- data.frame(Pi=nm[,1],Zeta=exp(nm[,2]),
            Delta=exp(nm[,3]),Mu=nm[,4],
            LogLikelihood=nm[,5],
            Convergence=paste("Optim",nm[,6]),
            Iterations=nm[,7],
            row.names=paste("Nelder-Mead",sv.names))
        }

    if(any(method=="nlm")){
        nlm <- matrix(rep(NA,7*length(start.values)),
            nrow=length(start.values),ncol=7)
        i <- 1
        if(any(start.values=="bn")){
            optim.bn.nlm <- nlm(llfunc,theta.start.bn,
                hessian=TRUE,iterlim=maxitnlm)
            nlm[i,1:4] <- optim.bn.nlm$estimate
            nlm[i,5] <- -optim.bn.nlm$minimum
            nlm[i,6] <- optim.bn.nlm$code
            nlm[i,7] <- optim.bn.nlm$iterations
            i <- i + 1
            }
        if(any(start.values=="fn")){
            optim.fn.nlm <- nlm(llfunc,theta.start.fn,
                hessian=TRUE,iterlim=maxitnlm)
            nlm[i,1:4] <- optim.fn.nlm$estimate
            nlm[i,5] <- -optim.fn.nlm$minimum
            nlm[i,6] <- optim.fn.nlm$code
            nlm[i,7] <- optim.fn.nlm$iterations
            i <- i + 1
            }
        if(any(start.values=="sl")){
            optim.sl.nlm <- nlm(llfunc,theta.start.sl,
                hessian=TRUE,iterlim=maxitnlm)
            nlm[i,1:4] <- optim.sl.nlm$estimate
            nlm[i,5] <- -optim.sl.nlm$minimum
            nlm[i,6] <- optim.sl.nlm$code
            nlm[i,7] <- optim.sl.nlm$iterations
            i <- i + 1
            }
        if(any(start.values=="us")){
            optim.us.nlm <- nlm(llfunc,theta.start.us,
                hessian=TRUE,iterlim=maxitnlm)
            nlm[i,1:4] <- optim.us.nlm$estimate
            nlm[i,5] <- -optim.us.nlm$minimum
            nlm[i,6] <- optim.us.nlm$code
            nlm[i,7] <- optim.us.nlm$iterations
            }
        output.nlm <- data.frame(Pi=nlm[,1],
            Zeta=exp(nlm[,2]),Delta=exp(nlm[,3]),
            Mu=nlm[,4],LogLikelihood=nlm[,5],
            Convergence=paste("Nlm",nlm[,6]),
            Iterations=nlm[,7],
            row.names=paste("NLM",sv.names))
        }

    names.vector <- c(row.names(output.bfgs),
                row.names(output.nm),
                row.names(output.nlm))

    output.df <- data.frame(Pi=c(output.bfgs$Pi,output.nm$Pi,
            output.nlm$Pi),Zeta=c(output.bfgs$Zeta,
            output.nm$Zeta,output.nlm$Zeta),
            Delta=c(output.bfgs$Delta,output.nm$Delta,
            output.nlm$Delta),Mu=c(output.bfgs$Mu,
            output.nm$Mu,output.nlm$Mu),
            LogLikelihood=c(output.bfgs$LogLikelihood,
            output.nm$LogLikelihood,
            output.nlm$LogLikelihood),Convergence=
            c(as.character(output.bfgs$Convergence),
            as.character(output.nm$Convergence),
            as.character(output.nlm$Convergence)),
            Iterations=c(output.bfgs$Iterations,
            output.nm$Iterations,output.nlm$Iterations),
            row.names=names.vector)

    starting.values <- data.frame(Pi.start=start.matrix[,1],
            Zeta.start=exp(start.matrix[,2]),
            Delta.start=exp(start.matrix[,3]),
            Mu.start=start.matrix[,4],row.names=sv.names)

    output <- list(pars=output.df[
            rev(order(output.df$LogLikelihood)),],
            breaks=breaks,
            starting.values=starting.values)

    if(plots==TRUE){
        old.par <- par(no.readonly = TRUE)
        par(mfrow=c(1,2),...)

        theta <- as.numeric(output$pars[
            output$pars$LogLikelihood==
            max(output$pars$LogLikelihood),1:4])
        K.nu <- besselK(theta[2],nu=1)
        hyp.dens <- function(x) 
              dhyperb(x,theta,K.nu=K.nu,log.pars=FALSE)    
        log.hyp.dens <- function(x) 
              log(dhyperb(x,theta,K.nu=K.nu,log.pars=FALSE)) 
        ymax <- 1.06*max(hyp.dens(seq(min(breaks),
            max(breaks),0.1)),emp.dens,na.rm=TRUE)

        hist(x,breaks,right=FALSE,freq=FALSE,
            ylim=c(0,ymax),main = paste("Histogram of", xname),...)
        curve(hyp.dens,min(breaks)-1,max(breaks)+1,
            add=TRUE,ylab=NULL) 
        log.hist(x,breaks,include.lowest=TRUE,right=FALSE,
                 main = paste("Log-Histogram of", xname),...)
        curve(log.hyp.dens,min(breaks)-1,max(breaks)+1,
            add=TRUE,ylab=NULL,xlab=NULL) 

        if(getOption("device")==names(dev.cur())) 
            par(ask=TRUE)
        par(mfrow=c(1,1))
        alphas <- (1:length(x)-0.5)/length(x) 
        quantiles <- qhyperb(alphas,theta)
        qqplot(quantiles,x,
            main = "Q-Q Plot of Hyperbolic Distribution",
            xlab = " Hyperbolic quantiles", 
            ylab = paste("Sample quantiles of",xname),...) 
        abline(0,1) 

        par(old.par)
        if((getOption("device")!=names(dev.cur()))
           &(leave.on==FALSE)) 
            dev.off()
        }

    if(max.lh==TRUE) output$pars[output$pars$LogLikelihood==
        max(output$pars$LogLikelihood),] else output
} # End of fit.hyperb()



#  Function to draw a log histogram.
#  Based on hist.default - only the plotting part differs.

log.hist <- 
function (x, breaks = "Sturges",
    include.lowest = TRUE, right = TRUE, 
    main = paste("Log-Histogram of", xname), 
    xlim = range(breaks), ylim = NULL, xlab = xname, 
    ylab = "Log Density", nclass = NULL, ...)
{ 
    if (!is.numeric(x)) 
        stop("`x' must be numeric")
    xname <- deparse(substitute(x))
    n <- length(x <- x[!is.na(x)])
    use.br <- !missing(breaks)
    if (use.br) {
        if (!missing(nclass)) 
            warning("`nclass' not used when `breaks' specified")
    }
    else if (!is.null(nclass) && length(nclass) == 1) 
        breaks <- nclass
    use.br <- use.br && (nB <- length(breaks)) > 1
    if (use.br) 
        breaks <- sort(breaks)
    else {
        if (is.character(breaks)) {
            breaks <- match.arg(tolower(breaks), c("sturges", 
                "fd", "freedman-diaconis", "scott"))
            breaks <- switch(breaks,
                             sturges = nclass.Sturges(x), 
                             "freedman-diaconis" = ,
                             fd = nclass.FD(x),
                             scott = nclass.scott(x), 
                stop("Unknown breaks algorithm"))
        }
        else if (is.function(breaks)) {
            breaks <- breaks(x)
        }
        if (!is.numeric(breaks) || is.na(breaks) || breaks < 2) 
            stop("invalid number of breaks")
        breaks <- pretty(range(x), n = breaks, min.n = 1)
        nB <- length(breaks)
        if (nB <= 1) 
            stop(paste("hist.default: pretty() error, breaks=", 
                format(breaks)))
    }
    h <- diff(breaks)
    equidist <- !use.br || diff(range(h)) < 1e-07 * mean(h)
    if (!use.br && any(h <= 0)) 
        stop("not strictly increasing `breaks'.")
    #if (is.null(freq)) {
    #    freq <- if (!missing(probability)) 
    #       !as.logical(probability)
    #    else equidist
    #}
    #else if (!missing(probability) && any(probability == freq)) 
    #    stop("`probability' is an alias for `!freq', however they differ.")
    diddle <- 1e-07 * max(abs(range(breaks)))
    fuzz <- if (right){ 
        c(if (include.lowest) -diddle else diddle,
          rep(diddle,length(breaks) - 1))
    }
    else{
      c(rep(-diddle, length(breaks) - 1),
        if (include.lowest) diddle else -diddle)
    }
    breaks <- breaks + fuzz
    h <- diff(breaks)
    storage.mode(x) <- "double"
    storage.mode(breaks) <- "double"
    counts <- .C("bincount", x, n, breaks, nB, counts = integer(nB - 
        1), right = as.logical(right), include = as.logical(include.lowest), 
        naok = FALSE, NAOK = FALSE, DUP = FALSE, PACKAGE = "base")$counts
    if (any(counts < 0)) 
        stop("negative `counts'. Internal Error in C-code for \"bincount\"")
    if (sum(counts) < n) 
        stop("some `x' not counted; maybe `breaks' do not span range of `x'")
    dens <- counts/(n * h)
    mids <- 0.5 * (breaks[-1] + breaks[-nB])
    
    ### up to here same as hist.default 
    log.density <- log(dens) 
    height <- range(log.density,finite=TRUE)[2] - 
              range(log.density,finite=TRUE)[1] 
    base <- min(log.density[is.finite(log.density)]) - 0.25 *height 

    ### yMax is the max value of log.density plus another 25%. 
    yMax <- 0.25*abs(max(log.density))+ max(log.density) 
    ### plot the log-histogram 
    if(is.null(ylim)) ylim <- range(base,yMax)
    plot(mids, log.density, xlim = xlim, ylim = ylim,
    type="n", xlab=xlab, ylab=ylab, main=main, ...) 
    points(mids, log.density) 
    heights <- rep(0,nB) 

    for (j in  2:(nB-1)) { 
        if(is.finite(max(log.density[j-1],log.density[j]))){ 
        heights[j] <- max(log.density[j-1],log.density[j]) 
    } 
    else { 
        heights[j] <- NA } 
    } 
    heights[1] <- 
        ifelse(is.finite(log.density[1]),log.density[1],NA) 
    heights[nB] <- 
    ifelse(is.finite(log.density[nB-1]),log.density[nB-1],NA) 

    i <- 1:(nB) 
    segments(breaks[i],log.density[i],breaks[i+1],log.density[i]) 
    segments(breaks[i],heights[i],breaks[i],base,lty=2) 
    segments(breaks[nB],heights[nB],breaks[nB],base,lty=2) 

    r<-list(breaks = breaks, counts = counts, 
        log.density = log.density, mids = mids, 
        xname = xname, heights = heights, ylim = ylim) 
    invisible(r) 
} # End of log.hist()


#  Function giving the density of the skew laplace distribution.

dskewlap <- function(x,theta){

    if (length(theta)!=3)
        stop("parameter vector must contain 3 values")
    skew.alpha <- theta[1] 
    skew.beta <- theta[2] 
    skew.mu <- theta[3] 
    below.mu <- 1/(skew.alpha + skew.beta)*
        exp((x - skew.mu)/skew.alpha)
    above.mu <- 1/(skew.alpha + skew.beta)*
        exp((skew.mu - x)/skew.beta)
    skewlap.dens <- ifelse(x <= skew.mu, below.mu, above.mu) 
    skewlap.dens
} # End of dskewlap



#  Function giving the density of the hyperbolic distribution.

dhyperb <- function(x,theta,K.nu=NULL,log.pars=FALSE){ 

  if (length(theta)!=4)
    stop("parameter vector must contain 4 values")
  if(log.pars==TRUE){ 
    hyperb.pi <- theta[1] 
    L.hyperb.zeta <- theta[2] 
    L.hyperb.delta <- theta[3] 
    hyperb.mu <- theta[4] 

    if(is.null(K.nu)){ 
      K.nu <- besselK(exp(L.hyperb.zeta), nu=1)     
    } 
     
    hyperb.dens <- (2 * exp(L.hyperb.delta) * 
                    sqrt(1 + hyperb.pi^2) * K.nu)^(-1) * 
                    exp(-exp(L.hyperb.zeta) * 
            (sqrt(1 + hyperb.pi^2) * 
                    sqrt(1 + ((x - hyperb.mu)/
            exp(L.hyperb.delta))^2) - 
                    hyperb.pi*(x - hyperb.mu)/
            exp(L.hyperb.delta))) 
  } 
  else{ 
    hyperb.pi <- theta[1] 
    hyperb.zeta <- theta[2] 
    hyperb.delta <- theta[3] 
    hyperb.mu <- theta[4] 
     
    if(is.null(K.nu)){ 
      K.nu <- besselK(hyperb.zeta, nu=1)     
    }   

    hyperb.dens <- (2 * hyperb.delta * 
            sqrt(1 + hyperb.pi^2) * K.nu)^(-1) * 
                    exp(-hyperb.zeta * (sqrt(1 + hyperb.pi^2) * 
                    sqrt(1 + ((x - hyperb.mu)/hyperb.delta)^2) - 
                    hyperb.pi*(x - hyperb.mu)/hyperb.delta)) 
  } 
  hyperb.dens 
} # End of dhyperb()



#  The hyperbolic distribution function.

phyperb <- function(q,theta,tol=10^(-5),subdivisions=100){ 

   if (length(theta)!=4){
    stop("parameter vector must contain 4 values")
   }

   hyperb.pi <- theta[1]
   zeta <- theta[2] 
   delta <- theta[3] 
   mu <- theta[4]

   K.nu <- besselK(zeta,nu=1)
   phi <- as.numeric(hyperb.change.pars(1,3,theta)[1])
   gamma <- as.numeric(hyperb.change.pars(1,3,theta)[2])

   c <- 1/(2*delta*(1+hyperb.pi^2)^(1/2)*K.nu)
   
   y.lower <- 1/phi*log(tol*phi/c)+mu 
   y.upper <- -1/gamma*log(tol*gamma/c)+mu
   range <- c(y.lower,y.upper)

   new.upper <- sort(q)
  
   low.limits <- c(y.lower,new.upper) 
   high.limits <- c(new.upper,y.upper)

   dhyp.int <- function(q){ 
     dhyperb(q,theta,K.nu)
   }

   int.fun <- rep(NA,length(new.upper)) 
   
   for(i in 1:length(low.limits)){ 
     int.fun[i] <- integrate(dhyp.int,low.limits[i],high.limits[i],
            subdivisions)$value
   } 

   int.fun.left <- cumsum(int.fun) 
   int.fun.right <- (1 - rev(cumsum(rev(int.fun))))
   
   int.fun.ave <- (int.fun.left[-length(int.fun.left)] +
            int.fun.right[-1])/2

   return(int.fun.ave[rank(q)]) 
} # End of phyperb()



# Function giving the quantiles of the hyperbolic distribution.

qhyperb <- function(p,theta,tol=10^(-5),
            n.interpol=100,subdivisions=100,...){ 
   if(length(theta) < 4){
    stop("parameter vector must contain 4 values") 
   } 

   hyperb.pi <- theta[1]
   zeta <- theta[2] 
   delta <- theta[3] 
   mu <- theta[4]

   K.nu <- besselK(zeta,nu=1)
   phi <- as.numeric(hyperb.change.pars(1,3,theta)[1])
   gamma <- as.numeric(hyperb.change.pars(1,3,theta)[2])

   c <- 1/(2*delta*(1+hyperb.pi^2)^(1/2)*K.nu)
   
   lower <- 1/phi*log(tol*phi/c)+mu 
   upper <- -1/gamma*log(tol*gamma/c)+mu
   range <- c(lower,upper)

   x.values <- seq(lower,upper,length=n.interpol)
   phyperb.values <- phyperb(x.values,theta,
                subdivisions=subdivisions) 
   phyperb.spline <- splinefun(x.values,phyperb.values) 
   q <- rep(NA,length(p)) 
   for(i in 1:length(p)){ 
     zero.fun<-function(x){ 
       phyperb.spline(x)-p[i] 
     } 
     if((0<p[i]) & (p[i]<1)) 
    q[i] <- uniroot(zero.fun,interval=c(lower,upper),...)$root
     if(p[i]==0) q[i] <- -Inf
     if(p[i]==1) q[i] <- Inf
     if((p[i]<0)|(p[i]>1)) q[i] <- NA
   } 
   return(q) 
} # End of qhyperb()



# Function to generate pseudo-random observations from a
# hyperbolic distribution.
# The algorithm is as given in Atkinson, A.C.
# "The Simulation of Generalized Inverse Gaussian and
# Hyperbolic Random Variables" (1982).

# Note: "theta" as given in A.C.Atkinson (1982) is different 
# to theta as used here.  Atkinson's "theta" is called 
# "theta.start" in this program


rhyperb <- function(n,theta){
    hyp.pi <- theta[1]
    zeta <- theta[2]
    delta <- theta[3]
    mu <- theta[4]
    alpha <- as.numeric(hyperb.change.pars(1,2,theta))[1]*delta
    beta <- as.numeric(hyperb.change.pars(1,2,theta))[2]*delta
    phi <- as.numeric(hyperb.change.pars(1,3,theta))[1]*delta
    gamma <- as.numeric(hyperb.change.pars(1,3,theta))[2]*delta

    theta.start <- -sqrt(phi*gamma)

    t <- -sqrt(gamma/phi)
    w <- sqrt(phi/gamma)

    delta1 <- exp(theta.start)/phi
    delta2 <- (w-t)*exp(theta.start)
    delta3 <- exp(-gamma*w)/gamma

    k <- 1/(delta1+delta2+delta3)

    r <- k*delta1
    v <- 1-k*delta3

    output <- numeric(n)
    need.value <- TRUE

    for(i in 1:n){
        while(need.value==TRUE){
            # Generate U & E
            U <- runif(1)
            E <- rexp(1)
            if(U<=r){
                x <- 1/phi*log(phi*U/k)
                if(E>=alpha*(sqrt(1+x^2)+x)){
                    need.value <- FALSE
                    }
                }
            if((U>r)&(U<=v)){
                x <- t-1/phi+U*exp(-theta.start)/k
                if(E>=alpha*sqrt(1+x^2)-beta*x+
                        theta.start){
                    need.value <- FALSE
                    }
                }
            if(U>v){
                x <- 1/gamma*log(k/gamma)-
                    1/gamma*log(1-U)
                if(E>=alpha*(sqrt(1+x^2)-x)){
                    need.value <- FALSE
                    }
                }
            } # End of while loop
        output[i] <- delta*x+mu
        need.value <- TRUE
    } # End of for loop
    output
} # End of rhyperb()



# Function to interchange parameterisations of the hyperbolic
# distribution.

hyperb.change.pars <- function(from,to,theta,no.names=FALSE){ 

  if(length (theta)!=4){
    stop("parameter vector must contain 4 values") 
  } 
  if((from!=1) & (from!=2) & (from!=3)){
    stop("the argument 'from' must be either 1,2 or 3")
  }  
  if((to!=1) & (to!=2) & (to!=3)){
    stop("the argument 'to' must be either 1,2 or 3")
  }

  delta <- theta[3]
  if(delta<=0) stop("delta must be greater than zero")
  mu <- theta[4]

  if(from==1){
    hyperb.pi <- theta[1]
    zeta <- theta[2]
    if(zeta<=0) stop("zeta must be greater than zero")
  }
  if(from==2){
    alpha <- theta[1]
    beta <- theta[2]
    if(alpha<=0) stop("alpha must be greater than zero")
    if(abs(beta)>=alpha) 
    stop("absolute value of beta must be less than alpha")
  }
  if(from==3){
    phi <- theta[1]
    gamma <- theta[2]
    if(phi<=0) stop("phi must be greater than zero")
    if(gamma<=0) stop("gamma must be greater than zero")
  }

  if(from==1 && to==2){ 
    alpha <- zeta*sqrt(1+hyperb.pi^2)/delta 
    beta <- zeta*hyperb.pi/delta 
    output=c(alpha=alpha,beta=beta,delta=delta,mu=mu) 
  } 
  if(from==1 && to==3){ 
    phi <- zeta/delta*(sqrt(1+hyperb.pi^2)+hyperb.pi)
    gamma <- zeta/delta*(sqrt(1+hyperb.pi^2)-hyperb.pi)
    output=c(phi=phi,gamma=gamma,delta=delta,mu=mu) 
  } 
  if(from==2 && to==3){ 
    phi <- alpha+beta
    gamma <- alpha-beta 
    output=c(phi=phi,gamma=gamma,delta=delta,mu=mu) 
  } 
  if(from==2 && to==1){ 
    hyperb.pi <- beta/sqrt(alpha^2-beta^2)
    zeta <- delta*sqrt(alpha^2-beta^2)
    output=c(hyperb.pi=hyperb.pi,zeta=zeta,delta=delta,mu=mu) 
  } 
  if(from==3 && to==1){ 
    hyperb.pi <- (phi-gamma)/(2*sqrt(phi*gamma))
    zeta <- delta*sqrt(phi*gamma)
    output=c(hyperb.pi=hyperb.pi,zeta=zeta,delta=delta,mu=mu) 
  } 
  if(from==3 && to==2){ 
    alpha <- (phi+gamma)/2
    beta <- (phi-gamma)/2
    output=c(alpha=alpha,beta=beta,delta=delta,mu=mu) 
  } 
  if(from==to){
    if(from==1) 
    output=c(hyperb.pi=hyperb.pi,zeta=zeta,delta=delta,mu=mu)
    if(from==2) 
    output=c(alpha=alpha,beta=beta,delta=delta,mu=mu)
    if(from==3) 
    output=c(phi=phi,gamma=gamma,delta=delta,mu=mu)
  }
  if(no.names==TRUE)
    names(output) <- NULL
  output  
} # End of hyperb.change.pars()



# Function to calculate the theoretical mode point of a 
# hyperbolic distribution given its parameters.

hyperb.mode <- function(theta) {
    hyperb.pi <- theta[1]
    zeta <- theta[2]
    delta <- theta[3]
    mu <- theta[4]
    hyperb.nu <- mu + delta*hyperb.pi
    hyperb.nu
} # End of hyperb.mode()



# Function to calculate the range over which the density of a
# hyperbolic distribution is concentrated.

calculate.range <- function(theta,tol=10^(-5)){
    hyperb.pi <- theta[1]
    zeta <- theta[2]
    delta <- theta[3]
    mu <- theta[4]
    K.nu <- besselK(zeta,nu=1)
    phi <- as.numeric(hyperb.change.pars(1,3,theta)[1])
    gamma <- as.numeric(hyperb.change.pars(1,3,theta)[2])
    c <- 1/(2*delta*(1+hyperb.pi^2)^(1/2)*K.nu)
    y.lower <- 1/phi*log(tol*phi/c) + mu
    y.upper <- -1/gamma*log(tol*gamma/c) + mu
    range <- c(y.lower,y.upper)
    return(range)
}



# R.lambda & S.lambda - functions used in calculation of the
# theoretical mean and variance of a hyperbolic distribution 
# (as given in Barndorff-Nielsen & Blaesild).

R.lambda <- function(zeta,lambda=1){
    besselK(zeta,nu=lambda+1)/besselK(zeta,nu=lambda)
    }

S.lambda <- function(zeta,lambda=1){
    (besselK(zeta,nu=lambda+2)*besselK(zeta,nu=lambda)-
     besselK(zeta,nu=lambda+1)^2)/besselK(zeta,nu=lambda)^2
    }



# Functions to calculate the theoretical mean and variance of a
# hyperbolic distribution (as given in Barndorff-Nielsen & Blaesild).

hyperb.mean <- function(theta){
    hyp.pi <- theta[1]
    zeta <- theta[2]
    delta <- theta[3]
    mu <- theta[4]
    mu + delta*hyp.pi*R.lambda(zeta,lambda=1)
    }

hyperb.var <- function(theta){
    hyp.pi <- theta[1]
    zeta <- theta[2]
    delta <- theta[3]
    mu <- theta[4]
    delta^2*(1/zeta*R.lambda(zeta) + hyp.pi^2*S.lambda(zeta))
    } 
