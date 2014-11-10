# RDB
# 2014-06-28
# random number generator based on accumulating coin flip errors
# used in Fat Tails

rcoin <- function(N, nCoins=4, acc=FALSE){
	# N is the number of random numbers to generate
	# nCoins is the number of coin flips ... it behaves somewhat like a variance parameter. The higher it is, the bigger the shocks (in either direction)
	# If acc is set to TRUE, the print-out will be a matrix where the first row of each column contains output that would be equivalent to acc=FALSE. The last row is the shock (startShock) that initializes the process, and the values above the bottom row show how the shocks accmulate after each coin flip. If N=100 and nCoins=4, you'll get a 4x100 matrix. Note that if you look from the last to first row in a column, successive rows that have the same value indicate that the coin flip was a 0.
	Funcall <- function(f, ...){
		f(...)	
	} 
	Iterate <- function(f, n = 1){
	    function(x){
			Reduce(Funcall, rep.int(list(f), n), x, right=TRUE, accumulate=acc)
		} 
	}
	startShock <- rnorm(N)
	if(acc){
		matrix(c(unlist(
			Iterate(f=function(x){coin <- rbinom(n=N, size=1, prob=0.5); shock <- (rnorm(N, sd=1)); x + coin*shock*x}, n=nCoins)(startShock)
			)), ncol=N, byrow=T)
	}else{
		Iterate(f=function(x){coin <- rbinom(n=N, size=1, prob=0.5); shock <- (rnorm(N, sd=1)); x + coin*shock*x}, n=nCoins)(startShock)	
	}
}
