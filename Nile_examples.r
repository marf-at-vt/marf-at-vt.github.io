# Nile: Measurements of the annual flow (in 10^8 cubic meters) 
#       of the Nile River at Ashwan 1871-1970
# Data available (and loaded by default) in R, taken from Durbin 
#    & Koopman (2001), Time Series Analysis of State Space models. 
# The data are analyzed using the dlm package in Petris, Petrone 
#    & Campagnoli (2009)

# Updated by Marco A. R. Ferreira (2025)

library(dlm)

############################
# Model 1: {1,1,15100,755} #
############################

plot(Nile,col='darkblue',xlab="time",ylab="",lwd=2,main="Model 1")
points(Nile,pch='*',col='darkblue')

model_1=dlmModPoly(order=1,dV=15100,dW=755) # Build model
NileShort=ts(Nile[1:95],start=1871,freq=1)
NileFilt_1=dlmFilter(NileShort,model_1) # Filtering

postmean_1=ts(c(dropFirst(NileFilt_1$m),rep(NA,5)),start=1871,freq=1)
postsd_1=ts(c(NileFilt_1$D.C,rep(NA,5)),start=1871,freq=1)
postsd_1_forecast=postsd_1+sqrt(model_1$V)
lines(postmean_1,col='red',lwd=2)
lines(postmean_1+1.96*postsd_1,lty=2,lwd=2,col='red')
lines(postmean_1-1.96*postsd_1,lty=2,lwd=2,col='red')

############################################
# Model 3: {1,1,v,w} with MLEs for v and w #
############################################

# Build model
build=function(param){
dlmModPoly(order=1,dV=exp(param[1]),dW=exp(param[2]))}

# Estimate variances using MLE
fit=dlmMLE(NileShort,rep(0,2),build)

# Here are the MLEs of the variances
exp(fit$par)

# Build model with variances equal to MLEs (Empirical Bayes)
model_3 = build(fit$par)

NileFilt_3=dlmFilter(NileShort,model_3) # Filtering

# FORECASTING

# Forecasting under Model 1
NileFor_1=dlmForecast(NileFilt_1,nAhead=5)
forecastmean_1=ts(c(rep(NA,94),postmean_1[95],NileFor_1$a),start=1871,freq=1)
forecastsd_1=ts(c(rep(NA,95),sqrt(as.double(NileFor_1$R))),start=1871,freq=1)

plot(Nile,col='white',xlab="time",ylab="",lwd=2,main="",ylim=c(450,1500))
points(Nile,pch='*',col='darkblue',cex=2)
lines(postmean_1,lwd=2,col='red')
lines(postmean_1-1.96*postsd_1,lwd=2,col='red',lty=2)
lines(postmean_1+1.96*postsd_1,lwd=2,col='red',lty=2)
lines(forecastmean_1,lwd=2,col='darkgreen')
lines(forecastmean_1+1.96*forecastsd_1,lwd=2,col='darkgreen',lty=2)
lines(forecastmean_1-1.96*forecastsd_1,lwd=2,col='darkgreen',lty=2)

# SMOOTHING...

# Model 1 smoothing
NileSmooth_1=dlmSmooth(NileFilt_1)
postmeansmooth_1=NileSmooth_1$s

plot(NileShort,ylab="",main="",lwd=1,col='gray')
points(NileShort,pch="*",cex=2,col='black')
lines(postmeansmooth_1,col='blue',lty=1,lwd=3)
lowersmooth_1=postmeansmooth_1-1.96*sqrt(unlist(dlmSvd2var(NileSmooth_1$U.S,NileSmooth_1$D.S)))
uppersmooth_1=postmeansmooth_1+1.96*sqrt(unlist(dlmSvd2var(NileSmooth_1$U.S,NileSmooth_1$D.S)))
lines(lowersmooth_1,lty=2,lwd=2,col='blue')
lines(uppersmooth_1,lty=2,lwd=2,col='blue')
