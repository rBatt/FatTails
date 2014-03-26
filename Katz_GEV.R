gev.fit <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, 

    mulink = identity, siglink = identity, shlink = identity, 

    show = TRUE, method = "Nelder-Mead", maxit = 10000, ...) 

{

    z <- list() #these are going to be aspects of the output, more later

    npmu <- length(mul) + 1 #number of parameters to be fitted for mu.  if mul is NULL, then length(mul) == 0, and npmu is 1

    npsc <- length(sigl) + 1 #number of params for sigma

    npsh <- length(shl) + 1 # n param for shape

    z$trans <- FALSE #transformation?

    in2 <- sqrt(6 * var(xdat))/pi #the initial guess for sigma to be used in the optimization

    in1 <- mean(xdat) - 0.57722 * in2 #initial guess for the location

    if (is.null(mul)) {

        mumat <- as.matrix(rep(1, length(xdat))) #if mul is NULL, then only 1 parameter, the intercept, is needed for mu, and the predictor variable is 1 

        muinit <- in1

    }

    else {

        z$trans <- TRUE #if we do have non-NULL mul, then the transformation is TRUE

        mumat <- cbind(rep(1, length(xdat)), ydat[, mul]) #the predictor variable matrix will include a column of 1's (for the intercept), and then all of the columns in ydat that correspond to other desired predictor variables.  E.g., mul could be equal to "Year", in which case our predictor matrix, mumat, would have a column of 1's of length xdat, and a column of the Years, also of length xdat.  The model would be similar to Y[i] = 1*B0 + Year[i]*B1, or Y = mumat%*%params, where Y is actually xdat, xdat has N rows and 1 column, mumat has N rows and 2 columns, and params is a matrix with 2 rows and 1 column

        muinit <- c(in1, rep(0, length(mul))) #now we have to make more guesses for the initial values --- we calculate the intercept as the "default" estimate of the parameter, and 0 (or no effect) as the default estimate of the other parameters

    }

    if (is.null(sigl)) { #same idea as for mul, and subsequent estiamtes of initials and matrices

        sigmat <- as.matrix(rep(1, length(xdat)))

        siginit <- in2

    }

    else {

        z$trans <- TRUE

        sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])

        siginit <- c(in2, rep(0, length(sigl)))

    }

    if (is.null(shl)) { #same idea as for mu and sigma

        shmat <- as.matrix(rep(1, length(xdat)))

        shinit <- 0.1

    }

    else {

        z$trans <- TRUE

        shmat <- cbind(rep(1, length(xdat)), ydat[, shl])

        shinit <- c(0.1, rep(0, length(shl)))

    }

    z$model <- list(mul, sigl, shl) #this is a list of 3 elements, where each element is a character vector whose elements in turn specify a column (name) of ydat to be used as a predictor variables for the respective GEV parameter

    z$link <- deparse(substitute(c(mulink, siglink, shlink))) #turns the link functions used into character strings

    init <- c(muinit, siginit, shinit)

    gev.lik <- function(a){
        mu <- mulink(mumat %*% (a[1:npmu]))
        sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
        xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
        y <- (xdat - mu)/sc #this line is part of t(x)
        y <- 1 + xi * y #this line is more of t(x), and is now only missing the ^(-1/Xi) part
        if (any(y <= 0) || any(sc <= 0)){
            return(10^6) #if you get an impossible value for y or sc, just return a huge nll, which will discourage the optimizer from picking the current parameters
			#also, this is a shortcut for the support that X is an element of (mu-sigma/Xi to INF) when Xi >0, and X is an element of (-INF to mu-sigma/Xi) when Xi < 0
				#there is only support for those values of X because t_x must be defined... 
				#in the case where Xi!=0, (1+((x-mu)/sigma)*xi), which is defined as "y" above, ends up being in the denominator when Y^(-1/xi)
				# Y^(-1/xi) is t(x)
				#while Y can't be zero, b/c denom, it can't be negative either b/c we take the log of y in the likelihood equation
		}
        sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y)*(1/xi + 1)) #negative log likelihood (from the pdf of the gev)
			
			#sum(log(sc)) is the -log() of 1/sigma
			#sum(y^(-1/xi)) is the -log() of exp(-t(x))
			#sum(log(y) * (1/xi + 1)) is the -log() of t(x)^(xi+1) 
				#because had (y^-1/xi)^*(xi+1)  
				#==> y^(-1 + -1/xi) 
				#==> log(y) * -(1 + 1/xi) 
				#==> -log(y)*(1/xi + 1) (and then take neg for nll)
    }

    x <- optim(init, gev.lik, hessian = TRUE, method = method, 

        control = list(maxit = maxit, ...))

    z$conv <- x$convergence

    mu <- mulink(mumat %*% (x$par[1:npmu]))

    sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))

    xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))

    z$nllh <- x$value

    z$data <- xdat

    if (z$trans) {

        z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))

    }

    z$mle <- x$par

    z$cov <- solve(x$hessian)

    z$se <- sqrt(diag(z$cov))

    z$vals <- cbind(mu, sc, xi)

    if (show) {

        if (z$trans) 

            print(z[c(2, 3, 4)])

        else print(z[4])

        if (!z$conv) 

            print(z[c(5, 7, 9)])

    }

    invisible(z)

}



gpd.fit

function (xdat, threshold, npy = 365, ydat = NULL, sigl = NULL, 

    shl = NULL, siglink = identity, shlink = identity, show = TRUE, 

    method = "Nelder-Mead", maxit = 10000, ...) 

{

    z <- list()

    npsc <- length(sigl) + 1

    npsh <- length(shl) + 1

    n <- length(xdat)

    z$trans <- FALSE

    if (is.function(threshold)) 

        stop("`threshold' cannot be a function")

    u <- rep(threshold, length.out = n)

    if (length(unique(u)) > 1) 

        z$trans <- TRUE

    xdatu <- xdat[xdat > u]

    xind <- (1:n)[xdat > u]

    u <- u[xind]

    in2 <- sqrt(6 * var(xdat))/pi

    in1 <- mean(xdat, na.rm = TRUE) - 0.57722 * in2

    if (is.null(sigl)) {

        sigmat <- as.matrix(rep(1, length(xdatu)))

        siginit <- in2

    }

    else {

        z$trans <- TRUE

        sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])

        siginit <- c(in2, rep(0, length(sigl)))

    }

    if (is.null(shl)) {

        shmat <- as.matrix(rep(1, length(xdatu)))

        shinit <- 0.1

    }

    else {

        z$trans <- TRUE

        shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])

        shinit <- c(0.1, rep(0, length(shl)))

    }

    init <- c(siginit, shinit)

    z$model <- list(sigl, shl)

    z$link <- deparse(substitute(c(siglink, shlink)))

    z$threshold <- threshold

    z$nexc <- length(xdatu)

    z$data <- xdatu

    gpd.lik <- function(a) {

        sc <- siglink(sigmat %*% (a[seq(1, length = npsc)]))

        xi <- shlink(shmat %*% (a[seq(npsc + 1, length = npsh)]))

        y <- (xdatu - u)/sc

        y <- 1 + xi * y

        if (min(sc) <= 0) 

            l <- 10^6

        else {

            if (min(y) <= 0) 

                l <- 10^6

            else {

                l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))

            }

        }

        l

    }

    x <- optim(init, gpd.lik, hessian = TRUE, method = method, 

        control = list(maxit = maxit, ...))

    sc <- siglink(sigmat %*% (x$par[seq(1, length = npsc)]))

    xi <- shlink(shmat %*% (x$par[seq(npsc + 1, length = npsh)]))

    z$conv <- x$convergence

    z$nllh <- x$value

    z$vals <- cbind(sc, xi, u)

    if (z$trans) {

        z$data <- -log(as.vector((1 + (xi * (xdatu - u))/sc)^(-1/xi)))

    }

    z$mle <- x$par

    z$rate <- length(xdatu)/n

    z$cov <- solve(x$hessian)

    z$se <- sqrt(diag(z$cov))

    z$n <- n

    z$npy <- npy

    z$xdata <- xdat

    if (show) {

        if (z$trans) 

            print(z[c(2, 3)])

        if (length(z[[4]]) == 1) 

            print(z[4])

        print(z[c(5, 7)])

        if (!z$conv) 

            print(z[c(8, 10, 11, 13)])

    }

    invisible(z)

}



pp.fit

function (xdat, threshold, npy = 365, ydat = NULL, mul = NULL, 

    sigl = NULL, shl = NULL, mulink = identity, siglink = identity, 

    shlink = identity, show = TRUE, method = "Nelder-Mead", maxit = 10000, 

    ...) 

{

    z <- list()

    npmu <- length(mul) + 1

    npsc <- length(sigl) + 1

    npsh <- length(shl) + 1

    n <- length(xdat)

    z$trans <- FALSE

    if (is.function(threshold)) 

        stop("`threshold' cannot be a function")

    u <- rep(threshold, length.out = n)

    if (length(unique(u)) > 1) 

        z$trans <- TRUE

    xdatu <- xdat[xdat > u]

    xind <- (1:n)[xdat > u]

    u <- u[xind]

    in2 <- sqrt(6 * var(xdat))/pi

    in1 <- mean(xdat) - 0.57722 * in2

    if (is.null(mul)) {

        mumat <- as.matrix(rep(1, length(xdatu)))

        muinit <- in1

    }

    else {

        z$trans <- TRUE

        mumat <- cbind(rep(1, length(xdatu)), ydat[xind, mul])

        muinit <- c(in1, rep(0, length(mul)))

    }

    if (is.null(sigl)) {

        sigmat <- as.matrix(rep(1, length(xdatu)))

        siginit <- in2

    }

    else {

        z$trans <- TRUE

        sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])

        siginit <- c(in2, rep(0, length(sigl)))

    }

    if (is.null(shl)) {

        shmat <- as.matrix(rep(1, length(xdatu)))

        shinit <- 0.1

    }

    else {

        z$trans <- TRUE

        shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])

        shinit <- c(0.1, rep(0, length(shl)))

    }

    init <- c(muinit, siginit, shinit)

    z$model <- list(mul, sigl, shl)

    z$link <- deparse(substitute(c(mulink, siglink, shlink)))

    z$threshold <- threshold

    z$npy <- npy

    z$nexc <- length(xdatu)

    z$data <- xdatu

    pp.lik <- function(a) {

        mu <- mulink(mumat %*% (a[1:npmu]))

        sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))

        xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))

        if (any(sc <= 0)) 

            return(10^6)

        if (min(1 + ((xi * (u - mu))/sc)) < 0) {

            l <- 10^6

        }

        else {

            y <- (xdatu - mu)/sc

            y <- 1 + xi * y

            if (min(y) <= 0) 

                l <- 10^6

            else l <- sum(log(sc)) + sum(log(y) * (1/xi + 1)) + 

                n/npy * mean((1 + (xi * (u - mu))/sc)^(-1/xi))

        }

        l

    }

    x <- optim(init, pp.lik, hessian = TRUE, method = method, 

        control = list(maxit = maxit, ...))

    mu <- mulink(mumat %*% (x$par[1:npmu]))

    sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))

    xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))

    z$conv <- x$convergence

    z$nllh <- x$value

    z$vals <- cbind(mu, sc, xi, u)

    z$gpd <- apply(z$vals, 1, ppp, npy)

    if (z$trans) {

        z$data <- as.vector((1 + (xi * (xdatu - u))/z$gpd[2, 

            ])^(-1/xi))

    }

    z$mle <- x$par

    z$cov <- solve(x$hessian)

    z$se <- sqrt(diag(z$cov))

    if (show) {

        if (z$trans) 

            print(z[c(2, 3)])

        if (length(z[[4]]) == 1) 

            print(z[4])

        print(z[c(5, 6, 8)])

        if (!z$conv) 

            print(z[c(9, 12, 14)])

    }

    invisible(z)

}
