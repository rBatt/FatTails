test <- function(){
	tryCatch(
		{
			o <- if(rbinom(n=1,1,0.1)==0){asdf}else{1}
			print(o)
			t <- if(rbinom(n=1,1,0.75)==0){asdf}else{2} 
			print(t)
			th <- if(rbinom(n=1,1,0.75)==0){asdf}else{3}
			print(th)
		
			return(o+t+th)
		},
		error=function(cond){
			return(Inf)
		}
	)
}