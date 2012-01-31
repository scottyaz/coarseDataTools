## dic.fit() fits a log-normal model to doubly interval censored incubation period data
## Parameters:
##    dat = a matrix with columns named "EL", "ER", "SL", "SR"
##    type = either "dic" or "sic"
##    [because optim() can be sensitive to starting parameters, some of them may
##     be changed in the options of dic.fit().]
##    start.log.sigma = the log-log-scale starting value for the dispersion
##    mu.int = the log-scale interval of possible median values (in days)
##    log.sigma.int = the log-log-scale interval of possible dispersion values
##    ... = additional options passed to optim


dic.fit <- function(dat,
		    start.log.sigma=log(log(2)),
		    opt.method="L-BFGS-B",
		    mu.int=c(log(.5), log(13)),
		    log.sigma.int=c(log(log(1.01)), log(log(5))),
		    ptiles=c(.05, .95, .99),
                    dist="L",...) {

	## check format of dat
	cnames <- colnames(dat)
	if(!("EL" %in% cnames)) stop("dat must have column named EL")
	if(!("ER" %in% cnames)) stop("dat must have column named ER")
	if(!("SL" %in% cnames)) stop("dat must have column named SL")
	if(!("SR" %in% cnames)) stop("dat must have column named SR")

	if(!("type" %in% cnames)) stop("dat must have column named type")
	if(!all(dat[,"type"] %in% c(0,1,2)))
		stop("values in type column must be either 0, 1 or 2.")

        if(!dist %in% c("G","W","L")) stop("Please use either L,W, or G for the distributions")
           ## fix sample size
           n <- nrow(dat)

           ## make sure dat is a matrix
           dat <- as.matrix(dat[,c("EL", "ER", "SL", "SR", "type")])
           if(class(dat)=="data.frame") stop("dat should be a matrix.")

           ## find starting values for DIC analysis using profile likelihoods
           start.mu <- optimize(f=pl.mu, interval=mu.int,
                                log.sigma=start.log.sigma, dat=dat,dist=dist)$min
           start.log.sigma <- optimize(f=pl.sigma, interval=log.sigma.int, mu=start.mu,
                                       dat=dat,dist=dist)$min

           ## find MLEs for doubly censored data using optim
           tmp <- list(convergence=1)
           msg <- NULL
           fail <- FALSE
           tryCatch(tmp <- optim(par=c(start.mu, start.log.sigma),
                                 method=opt.method, hessian=TRUE,
                                 lower=c(log(0.5), log(log(1.04))),
                                 fn=loglik, dat=dat,dist=dist, ...),
		 error = function(e) {
			 msg <<- e$message
			 fail <<- TRUE
		 },
		 warning = function(w){
			 msg <<- w$message
			 fail <<- TRUE
		 })
	## also, to catch a few more errors
	if(tmp$convergence!=0 | all(tmp$hessian==0) ){
		msg <- tmp$message
		if(all(tmp$hessian==0)) msg <- paste(msg, "& hessian is singular")
		fail <- TRUE
	}

	if(!fail ){
                cat(sprintf("Loglikelihood = %f \n",tmp$value))
                cat(sprintf("Parameter Values = %f \n ",tmp$par))
                if (dist == "L"){
                    med <- exp(tmp$par[1])
                    disp <- exp(exp(tmp$par[2]))
                    norm.quants <- qnorm(ptiles)
                    ests <- c(med,
                              disp,
                              med*disp^norm.quants)
                    Sig <- solve(tmp$hessian)
                    ses <- dic.getSE(log(med), log(log(disp)), Sig, ptiles,dist=dist)
                    cil <- ests - qt(.975, n-1)*ses
                    cih <- ests + qt(.975, n-1)*ses
                } else {
                    ses <- dic.getSE(dat=dat,mu=tmp$par[1],log.s=tmp$par[2], Sig=NULL, ptiles=ptiles,dist=dist,opt.method=opt.method)
                }
		## save the quantile estimates
		quant.matrix <- matrix(c(ests, cil, cih, ses),
				       nrow=2+length(ptiles), byrow=FALSE)
		ptiles.names <- paste("p", 100*ptiles, sep="")
		rownames(quant.matrix) <- c("p50", "disp", ptiles.names)
		colnames(quant.matrix) <- c("est", "CIlow", "CIhigh", "StdErr")

		return(list(ests=round(quant.matrix, 3),
			    conv=1,
			    MSG=NULL,
			    Sig.log.scale=Sig,
                            loglik=tmp$value))
	}
	else {
		return(list(ests=matrix(NA, nrow=5, ncol=4),
			    conv=0,
			    MSG=msg,
			    Sig.log.scale=NULL,
                            loglik=NULL))
	}


}


## profile likelihood for mu -- used by dic.fit() to get starting values
pl.mu <- function(mu, log.sigma, dat, dist){
	loglik(pars=c(mu, log.sigma),dist=dist,dat=dat)
}


## profile likelihood for sigma -- used by dic.fit() to get starting values
pl.sigma <- function(log.sigma, mu, dat, dist){
	loglik(pars=c(mu, log.sigma), dist=dist, dat=dat)
}

## functions that manipulate/calculate the likelihood for the censored data

## the functions coded here are taken directly from the
## doubly interval censored likelihood notes.

fw1 <- function(t, EL, ER, SL, SR, mu, sigma, dist){
	## function that calculates the first function for the DIC integral
    if (dist=="W"){
        (ER-SL+t) * (ER-SL+t) * dweibull(x=t,shape=mu,scale=sigma)
    } else if (dist=="G") {
        (ER-SL+t) * dgamma(x=t, shape=mu, scale=sigma)
    } else {
        (ER-SL+t) * dlnorm(x=t, meanlog=mu, sdlog=sigma)
    }
}

fw3 <- function(t, EL, ER, SL, SR, mu, sigma, dist){
	## function that calculates the third function for the DIC integral
    if (dist == "W"){
	(SR-EL-t) * dweibull(x=t, shape=mu, scale=sigma)
    } else if (dist == "G"){
    	(SR-EL-t) * dgamma(x=t, shape=mu, scale=sigma)
    } else {
        (SR-EL-t) * dlnorm(x=t, meanlog=mu, sdlog=sigma)
    }
}


lik <- function(mu, sigma, EL, ER, SL, SR, type, dist){
	## returns the right likelihood for the type of data
	## 0 = DIC, 1=SIC, 2=exact
	if(type==0) return(diclik2(mu, sigma, EL, ER, SL, SR, dist))
	if(type==1) return(siclik(mu, sigma, EL, ER, SL, SR, dist))
	if(type==2) return(exactlik(mu, sigma, EL, ER, SL, SR, dist))
}


diclik <- function(mu, sigma, EL, ER, SL, SR, dist){
	## calculates the DIC likelihood by integration

	## if symptom window is bigger than exposure window
	if(SR-SL>ER-EL){
		dic1 <- integrate(fw1, lower=SL-ER, upper=SL-EL,
				  subdivisions=10,
				  mu=mu, sigma=sigma,
				  EL=EL, ER=ER, SL=SL, SR=SR,
                                  dist=dist)$value
                if (dist == "W"){
                    dic2 <- (ER-EL)*
                        (pweibull(SR-ER, shape=mu, scale=sigma) - pweibull(SL-EL, shape=mu, scale=sigma))
                } else if (dist == "G"){
                    dic2 <- (ER-EL)*
                        (pgamma(SR-ER, shape=mu, scale=sigma) - pgamma(SL-EL, shape=mu, scale=sigma))
                } else {
                    dic2 <- (ER-EL)*
                        (plnorm(SR-ER, mu, sigma) - plnorm(SL-EL, mu, sigma))
                }
		dic3 <- integrate(fw3, lower=SR-ER, upper=SR-EL,
				  subdivisions=10,
				  mu=mu, sigma=sigma,
				  EL=EL, ER=ER, SL=SL, SR=SR,
                                  dist=dist)$value
		return(dic1 + dic2 + dic3)
	}

	## if exposure window is bigger than symptom window
	else{
            dic1 <- integrate(fw1, lower=SL-ER, upper=SR-ER,
                              subdivisions=10,
                              mu=mu, sigma=sigma,
                              EL=EL, ER=ER, SL=SL, SR=SR,
                              dist=dist)$value
            if (dist == "W"){
		dic2 <- (SR-SL)*
                    (pweibull(SL-EL, shape=mu, scale=sigma) - pweibull(SR-ER, shape=mu, scale=sigma))
            } else if (dist == "G"){
                dic2 <- (SR-SL)*
                    (pgamma(SL-EL, shape=mu, scale=sigma) - pgamma(SR-ER, shape=mu, scale=sigma))
            } else {
                dic2 <- (SR-SL)*
                    (plnorm(SL-EL, mu, sigma) - plnorm(SR-ER, mu, sigma))
            }
		dic3 <- integrate(fw3, lower=SL-EL, upper=SR-EL,
				  subdivisions=10,
				  mu=mu, sigma=sigma,
				  EL=EL, ER=ER, SL=SL, SR=SR,
                                  dist=dist)$value
		return(dic1 + dic2 + dic3)
	}
    }

## this dic likelihood is designed for data that has overlapping intervals
diclik2 <- function(mu, sigma, EL, ER, SL, SR, dist){
	if(SL>ER) {
		return(diclik(mu, sigma, EL, ER, SL, SR, dist))
	} else {
		lik1 <- integrate(diclik2.helper1, lower=EL, upper=SL,
				  SL=SL, SR=SR, mu=mu, sigma=sigma, dist=dist)$value
		lik2 <- integrate(diclik2.helper2, lower=SL, upper=ER,
				  SR=SR, mu=mu, sigma=sigma, dist=dist)$value
		return(lik1+lik2)
	}
}

## likelihood functions for diclik2
diclik2.helper1 <- function(x, SL, SR, mu, sigma, dist){
    if (dist =="W"){
        pweibull(SR-x, shape=mu, scale=sigma) - pweibull(SL-x, shape=mu, scale=sigma)
    } else if (dist =="G") {
        pgamma(SR-x, shape=mu, scale=sigma) - pgamma(SL-x, shape=mu, scale=sigma)
    } else {
        plnorm(SR-x, mu, sigma) - plnorm(SL-x, mu, sigma)
    }
}

diclik2.helper2 <- function(x, SR, mu, sigma, dist){
    if (dist =="W"){
        pweibull(SR-x, shape=mu, scale=sigma)
    } else if (dist =="G") {
        pgamma(SR-x, shape=mu, scale=sigma)
    } else {
	plnorm(SR-x, mu, sigma)
    }
}

siclik <- function(mu, sigma, EL, ER, SL, SR, dist){
	## calculates the SIC likelihood as the difference in CDFs
    if (dist =="W"){
        pweibull(SR-EL, shape=mu, scale=sigma) - pweibull(SL-EL, shape=mu, scale=sigma)
    } else if (dist =="G") {
        pgamma(SR-EL, shape=mu, scale=sigma) - pgamma(SL-EL, shape=mu, scale=sigma)
    } else {
        plnorm(SR-EL, mu, sigma) - plnorm(SL-ER, mu, sigma)
    }
}

exactlik <- function(mu, sigma, EL, ER, SL, SR, dist){
	## calculates the likelihood for an exact observation

	## NB: the two Ss should be equal and the two Es should be equal
	##     so it doesn't matter which pair we use in the formula below.
    if (dist =="W"){
        dweibull(SR-EL, shape=mu, scale=sigma)
    } else if (dist =="G") {
        dgamma(SR-EL, shape=mu, scale=sigma)
    } else {
        dlnorm(SR-EL, mu, sigma)
    }
}

loglik <- function(pars, dat, dist) {
	## calculates the log-likelihood of DIC data
	## dat must have EL, ER, SL, SR and type columns
	mu <- pars[1]
	sigma <- exp(pars[2])
	sprintf("mu = %.2f, sigma = %.2f",mu, sigma)  ## for debugging
	n <- nrow(dat)
	totlik <- 0
	for(i in 1:n){
		totlik <- totlik +
			log(lik(mu, sigma, type=dat[i,"type"],
				EL=dat[i,"EL"], ER=dat[i,"ER"],
				SL=dat[i,"SL"], SR=dat[i,"SR"],
                                dist=dist))

	}
	return(-totlik) ## NB: NEEDS TO BE -totlik IF WE ARE MAXIMIZING USING OPTIM!
}


## calculates the standard errors for estimates from dic.fit() using delta method
        dic.getSE <- function(mu, log.s, Sig, ptiles, dist,dat=dat,opt.method){
            n.boots <- 1000
            boots <- vector("list",n.boots)
            if (dist == "L") {
                cat(sprintf("Asymtotic Confidence Intervals for Log Normal"))
            s <- exp(log.s)
            qnorms <- qnorm(ptiles)
            df <- matrix(c(exp(mu), 0, exp(mu+qnorms*s),
                           0, exp(s+log.s), qnorms * exp(mu + qnorms*s + log.s)),
                         nrow=2, ncol=2+length(ptiles), byrow=TRUE)
            ses <- sqrt(diag(t(df)%*%Sig%*%df))
           } else {
               cat(sprintf("Bootstrapping (n=1000) Standard Errors for %s \n",dist))
               #draw line numbers for 1000 new data sets
               line.nums <- matrix(sample(1:nrow(dat),nrow(dat)*n.boots,replace=T),nrow=nrow(dat),ncol=n.boots)
               for (i in 1:n.boots){
                   if (i %% 5 == 0) cat(".")
                   boots[[i]] <-
                       single.boot(mu.s=mu,log.s.s=log.s,opt.method=opt.method,dat.tmp=dat[line.nums[,i],],dist=dist)
               }
           }
            return(ses)
        }

#estimates one set of parameters for bootstraps
single.boot <- function(mu.s,log.s.s,opt.method,dat.tmp,dist,...){
    ## tmp <- list(convergence=1)
    ## msg <- NULL
    ## fail <- FALSE
#    tryCatch(
             tmp <- optim(par=c(mu.s,log.s.s),
                          method=opt.method, hessian=FALSE,
                          lower=c(-10,-10),
                          fn=loglik, dat=dat.tmp,dist=dist, ...)
    ## ,
    ##          error = function(e) {
    ##              msg <- e$message
    ##              fail <- TRUE
    ##          },
    ##          warning = function(w){
    ##              msg <- w$message
    ##              fail <- TRUE
    ##          })
    ## also, to catch a few more errors
                    ##     if(tmp$convergence!=0 | all(tmp$hessian==0) ){
                    ##         msg <- tmp$message
                    ##         if(all(tmp$hessian==0)) msg <- paste(msg, "& hessian is singular")
                    ##         fail <- TRUE
                    ## }
    return(tmp)
}

        get.obs.type <- function(dat) {
            type <- rep(0, nrow(dat))

                                        #get the single interval censored
            type[dat[,"EL"]==dat[,"ER"]]<-1
            type[dat[,"SL"]==dat[,"SR"]]<-1
            type[dat[,"ER"]>=dat[,"SL"]]<-1

                                        #some of those are actually exact!
            type[(dat[,"EL"]==dat[,"ER"]) & (dat[,"SL"]==dat[,"SR"])]<- 2

            return(type)
        }
