#FatColor_v1.0
#RDB, 9-Dec-2010
#Inspired by: Vasseur & Yodzis. 2004. The Color of Environmental Noise. Ecology
rm(list=ls())
graphics.off()

#My function for generating serially correlated observations
nObs=1000 #Number of simulated observations
#This is a function to simulate serially correlated data.  A description of the arguments:
	# "phi" is the AR(1) Coefficient.
	# "theta" is the MA(1) coefficient, where a positive value for theta allows the series to remain stationary.
	# "AbsShock" is my cheap way of mimicking 'resilience'--- what portion of random shocks actually make it to the response variable?
	# "rMu" is the mean of the random normal errors to be generated
	# "rSd" is the standard deviation of the random normal errors to be generated
	# "C" is some random constant from the Ecology paper.  
	#I honestly don't understand why you do the whole "C*(1-phi^2)^0.5" thing in front of Epsilon.
	#See the caption for Figure 5c.
MyArma <- function(nObs=1000, phi=0.7, theta=0, AbsShock=1, rMu=0, rSd=1, C=20){
	Eps <- rnorm(nObs, rMu, rSd) #Create all the random noise
	Y0 <- rnorm(1) + Eps[1] #The first data point.
	Y <- c(Y0) #Define the vector of observations
	for(i in 2:nObs){#Begin loop to create data
		Y[i] <- ifelse(phi!=0, phi*Y[i-1], 0) - ifelse(theta!=0, theta*Eps[i-1], 0)  + C*(1-phi^2)^0.5 *AbsShock*Eps[i] #+ AbsShock*Eps[i]
	}
	return(Y)#Return the data
}

#My function for conducting spectral analysis-- yields a spectral density plot and a Beta, and fft values... or w/e combo
#I basically got these calculations from <http://knol.google.com/k/suresh-emre/fast-fourier-transform-and-power/35vsnxisjn2mw/142#>
#The calculations seemed reasonable, although something way funky is going on with my x-axis.
MySpectral <- function(Data, LOG=TRUE, Return=Beta, PlotB=TRUE){
	Raw <- fft(Data) #get the fast Fourier Transform from the simluated data
	PhaseRaw <- Arg(Raw); PhaseRaw <- PhaseRaw[1:length(PhaseRaw)/2] #get the phase from the first half by dividing the imaginary elements by the real elements
	PowerRaw <- Mod(Raw); PowerRaw <- PowerRaw[1:length(PowerRaw)/2] #calculate the power values by: adding the reals squared to the imaginaries squared, then taking sqrt of that sum.
	xRaw <- 1:length(PowerRaw)/(2*length(PowerRaw)) #this is how I get the x-axis, with the min being 0 and max being .5 ... which is kinda what it should be.....

	# Power=1/(Frequency^Beta)
	#log(P) = -Beta*log(f)
	Model <- lm(log(PowerRaw) ~ log(xRaw)) #Best fit line to the PowerRaw vs. xRaw
	Inter <- Model$coef[[1]] #the intercept
	Beta <- -Model$coef[[2]] #Beta, the spectral exponent.
	
	#What plots do I want drawn? Log? Best fit line?
	#dev.new()
	if(LOG==FALSE){
		plot(xRaw, PowerRaw, type="l") #no log scale, no best fit line
	}else{
		if(PlotB==TRUE){plot(log(xRaw), log(PowerRaw), type="l", main="Power Spectrum (?)", sub=paste("Beta=", round(Beta, digits=3),sep=""), cex.lab=.75, col.sub="red") #If we want to plot on a log scale, but don't want the BFL
						lines(log(xRaw), fitted(Model), col="red") #If we want to plot on a log scale, and we want the Beta slope drawn
		}else{
			plot(log(xRaw), log(PowerRaw), type="l")
		}
	}
	qqnorm(Model$resid)#qq plot of the residuals of the model used to find Beta
	qqline(Model$resid)#draw the line through the qq plot
	#What values do I want returned? Beta? fft?  Both?
	if(Return==Beta){
		return(Beta)
	}else{
		if(Return==Power){
			return(PowerRaw)
		}else{
			return(list(PowerRaw,Beta))
		}
	}	
}#End MySpectral function


Data <- MyArma(nObs=nObs)#Simulate data

dev.new()
par(mfrow=c(2,2))
plot(1:length(Data),Data, type="l")#plot the raw simulated data

ylimMax <- max(c(max(curve(nObs*dnorm(x,mean=mean(Data),sd=sd(Data)), type="n", add=TRUE)$y), max(hist(Data,plot=F)$counts))) #figure out what I want the upper bound of ylim to be for next plot
#area under the plot == nObs
hist(Data, ylim=c(0,ylimMax), main="Hist of Data; Curve of Nrml Apprx.") #histogram of the simulated data
curve(nObs*dnorm(x,mean=mean(Data),sd=sd(Data)), add=TRUE, col="blue"); #draw a Normal PDF*nObs that is based on the mean and sd of the raw simulated data
#I did the above 2 steps in hopes of getting some basis of comparison for the fatness of the simulated data vs. the fatness of the normal distribution.  Kinda works?

#Density under the plot== 1
# DataDist <- hist(Data,plot=F)$counts/1000
# plot(seq(min(Data),max(Data), length.out=length(DataDist)), DataDist, type="h", ylim=c(0,ylimMax/1000))
# curve(dnorm(x,mean=mean(Data),sd=sd(Data)), add=TRUE, col="green")


MySpectral(Data)
#plot(seq(1,4*pi, length.out=1000), 4*sin(.5*seq(1,4*pi, length.out=1000)))