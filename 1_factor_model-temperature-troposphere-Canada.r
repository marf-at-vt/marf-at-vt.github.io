# Analysis of troposphere air daily average temperature at pressure levels 1000, 925, 850, 700, 600, 500, 400, 300, and 250 over Canada from June 1 to August 31, 2015.

le = read.csv("temperature-troposphere.csv")
tseq = as.Date(le[,1])
nc.a = t(le[,2:10])

# Plot data through time
ll = min(nc.a)
ul = max(nc.a)
plot(as.Date(tseq,"%D"),nc.a[1,],type="l",ylim=c(ll,ul),cex=2,ylab="Temperature (Kelvin)",xlab="Date")
for(i in 2:9) lines(as.Date(tseq,"%D"),nc.a[i,])

m = apply(nc.a,1,mean)
nc.a.centered = matrix(NA,nrow=nrow(nc.a),ncol=ncol(nc.a))
for(i in 1:9) nc.a.centered[i,] = nc.a[i,] - m[i]
ll = min(nc.a.centered)
ul = max(nc.a.centered)
plot(nc.a.centered[1,],type="l",ylim=c(ll,ul),cex=2,ylab="")
for(i in 1:9) lines(nc.a.centered[i,])

nc.a.standardized = matrix(NA,nrow=nrow(nc.a),ncol=ncol(nc.a))
for(i in 1:9)
{
  s=sqrt(var(nc.a[i,]))
  nc.a.standardized[i,] = nc.a.centered[i,] / s
}
ll = min(nc.a.standardized)
ul = max(nc.a.standardized)
plot(nc.a.standardized[1,],type="l",ylim=c(ll,ul),cex=2,ylab="")
for(i in 1:9) lines(nc.a.standardized[i,])


#######################################################################
#                                                                     #
# Gibbs Sampler for factor model with one factor.                     #
#                                                                     #
# Marco A. R. Ferreira, January 2017.                                 #
# Model (with k=1):                                                   #
# Y_i = B x_i + v_i                                                   #
# where                                                               #
# Y_i: i-th observation of r-dimensional vector of interest           #
# B: r by k factor loadings matrix                                    #
# x_i: i-th realization of the k-dimensional latent factor process    #
# x_i independent N(0, H), i = 1, ..., n                              #
# H = diag(exp(lambda_1), ..., exp(lambda_k))                         #
# v_i: r-dimensional error vector                                     #
# v_i i.i.d. N(0,V)                                                   #
# V = diag(sigma2_1 , ..., sigma2_r)                                  #
#                                                                     #
# Alternatively, we can write in matrix notation                      #
# Y = X B' + v                                                        #
# where                                                               #
# Y = (Y_1', ..., Y_n')'is n by r matrix of observations              #
# X = (x_1', ..., x_n')' is n by k matrix of latent factor process    #
# v = (v_1', ..., v_n')' is n by r matrix of errors                   #
#                                                                     #
#######################################################################


##########################################
##########################################
##                                      ##
##    Assign temperature data to Y      ##
##                                      ##
##########################################
##########################################

set.seed(12345)  # set pseudo-random generator seed for reproducibility of 
                 # numerical results and plots

Y = t(nc.a.standardized)

n = nrow(Y)      # Sample size
k =    1         # Number of factors
r = ncol(Y)      # Dimensional of each observational vector

# Checking the eigenvalues of the sample correlation matrix
sample.correlation.matrix = cor(Y)
plot(eigen(sample.correlation.matrix)$values)
# The first component explains about 75.8% of the variability.

##########################################
##########################################
##                                      ##
##             Gibbs Sampler            ##
##                                      ##
##########################################
##########################################

##############################
#                            #
# Set priors hyperparameters #
#                            #
##############################

n.b = 1
n.s2.b = 1
n.H = 1
n.s2.H = 1

# Priors for sigma^2_j recommended by Lopes and West (2004) for modeling
# standardized variables (and therefore the correlation matrix). This choice
# implies prior means equal to 0.5 for each sigma^2_j.
n.sigma =  2.2
n.s2.sigma = 0.1


######################
#                    #
# Set initial values #
#                    #
######################

Gbig = 5000

b.j.current = c(1,rep(0,r-1))
tau2b.current = 1
x.i.current = apply(Y,1,mean)
H.current = 100
sigma2.j.current = rep(1,r)


#######################################
#                                     #
# Define arrays to store Gibbs sample #
#                                     #
#######################################

b.j.Gibbs = matrix(NA,nrow=Gbig,ncol=r)
tau2b.Gibbs = rep(NA,Gbig)  
x.i.Gibbs = matrix(NA,nrow=Gbig,ncol=n)
H.Gibbs = rep(NA,Gbig)
sigma2.j.Gibbs = matrix(NA,nrow=Gbig,ncol=r)

b.j.Gibbs[1,] = b.j.current
tau2b.Gibbs[1] = tau2b.current
x.i.Gibbs[1,] = x.i.current
H.Gibbs[1] = H.current
sigma2.j.Gibbs[1,] = sigma2.j.current


######################
#                    #
# Gibbs Sampler loop #
#                    #
######################

for(g in 2:Gbig)
{

  cat("g=",g,"\n")

  # Simulate x_i, i = 1, ..., n
  for(i in 1:n)
  {
    sum.bj2.precision.j = sum(b.j.current^2 / sigma2.j.current)
    sum.bj.yij.precision.j = sum(b.j.current * Y[i,] / sigma2.j.current)
    mean.aux = sum.bj.yij.precision.j / (sum.bj2.precision.j + 1/H.current)
    var.aux = 1 / (sum.bj2.precision.j + 1/H.current)
    x.i.current[i] = rnorm(1,mean.aux,sqrt(var.aux))
  }  
  
  # Simulate sigma_j^2, j = 1, ..., r
  for(j in 1:r)
  {
    n.aux = n.sigma + n
    n.s2.aux = n.s2.sigma + sum((Y[,j] - b.j.current[j] * x.i.current)^2) 
    sigma2.j.current[j] = 1 / rgamma(1,shape=n.aux/2,rate=n.s2.aux/2)
  }
  
  # Simulate b_j, j = 2, ..., r
  for(j in 2:r)
  {
    sum.x2 = sum(x.i.current^2)
    sum.xy = sum(x.i.current * Y[,j])
    mean.aux = sum.xy/(sum.x2 + sigma2.j.current[j]/tau2b.current)
    var.aux = sigma2.j.current[j] /(sum.x2 + sigma2.j.current[j]/tau2b.current)
    b.j.current[j] = rnorm(1,mean.aux,sqrt(var.aux))
  }
    
  # Simulate tau_b^2
  n.aux = n.b + r - 1
  n.s2.aux = n.s2.b + sum(b.j.current^2)
  tau2b.current = 1 / rgamma(1,shape=n.aux/2,rate=n.s2.aux/2)
  
  # Simulate H
  n.aux = n.H + n
  n.s2.aux = n.s2.H + sum(x.i.current^2)
  H.current = 1 / rgamma(1,shape=n.aux/2,rate=n.s2.aux/2)
  
  # Copy current values into storage arrays
  b.j.Gibbs[g,] = b.j.current
  tau2b.Gibbs[g] = tau2b.current
  x.i.Gibbs[g,] = x.i.current
  H.Gibbs[g] = H.current
  sigma2.j.Gibbs[g,] = sigma2.j.current

}

######################
#                    #
#      Results       #
#                    #
######################

# Trace plots
plot(H.Gibbs,type="l",ylab="H",xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(tau2b.Gibbs,type="l",ylab=expression(tau^2),xlab="iteration",cex=4.0,main="",lwd=2.0)

burnin = 1000

plot(H.Gibbs[burnin:Gbig],type="l",ylab="H",xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(tau2b.Gibbs[burnin:Gbig],type="l",ylab=expression(tau^2),xlab="iteration",cex=4.0,main="",lwd=2.0)

hist(H.Gibbs[burnin:Gbig],xlab="H",cex=4.0,main="",lwd=2.0)

hist(tau2b.Gibbs[burnin:Gbig],xlab=expression(tau^2),cex=4.0,main="",lwd=2.0)

ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlab="Date",ylab="",main="Common factor")
lines(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,],2,quantile,p=0.025),lty=2,lwd=2)  
lines(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,],2,quantile,p=0.975),lty=2,lwd=2)



# Boxplots of posterior draws for different b_j's
pressure = c(1000, 925, 850, 700, 600, 500, 400, 300, 250)
aux = cbind(b.j.Gibbs[burnin:Gbig, 1],rep(pressure[1], Gbig - burnin + 1))  
for(j in 2:9) aux = rbind(aux,cbind(b.j.Gibbs[burnin:Gbig, j],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(b[j]),cex=2,main="Factor loadings")
abline(h=1,lty=2)

# Boxplots of posterior draws for different sigma2_j's
pressure = c(1000, 925, 850, 700, 600, 500, 400, 300, 250)
aux = cbind(sigma2.j.Gibbs[burnin:Gbig, 1],rep(pressure[1], Gbig - burnin + 1))  
for(j in 2:9) aux = rbind(aux,cbind(sigma2.j.Gibbs[burnin:Gbig, j],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(sigma[j]^2),cex=2,ylim=c(0,1.7),main="Idiosyncratic variances")
abline(h=1,lty=2)


percentage.unexplained.variability = apply(sigma2.j.Gibbs[burnin:Gbig,],1,sum)/9
hist(percentage.unexplained.variability,xlab="",ylab="",xlim=c(0,0.4),main="Percentage unexplained variability")


#####################################################################
#####################################################################
##                                                                                                                                  ##
##                                                                                                                                  ##
## For model selection, compute a Laplace-Metropolis approximation                     ##
## to the marginal data density.                                                                                   ##
##                                                                                                                                 ##
##                                                                                                                                 ##
#####################################################################
#####################################################################

######################################################
# Compute the logarithm of the integrated likelihood #
# at each iteration of the MCMC                      #
######################################################

log.integrated.likelihood = function(Y,B,H,V)
{
  n = nrow(Y)
  cov.mat = H * B %*% t(B) + V
  result = -(n/2) * determinant(cov.mat,logarithm=TRUE)$modulus
  cov.mat.inv = solve(cov.mat)
  for(i in 1:n) result = result - 0.5 * t(Y[i,]) %*% cov.mat.inv %*% Y[i,]
  result
}

log.int.lik = rep(NA,Gbig)
for(g in 1:Gbig)
{
  B = b.j.Gibbs[g,]
  H = H.Gibbs[g]
  V = diag(sigma2.j.Gibbs[g,])
  log.int.lik[g] = log.integrated.likelihood(Y,B,H,V)
}

burnin = 2000

max.log.int.lik = max(log.int.lik)
index.max = which(log.int.lik == max.log.int.lik)
B.hat = b.j.Gibbs[index.max,]
H.hat = H.Gibbs[index.max]
V.hat = diag(sigma2.j.Gibbs[index.max,])
tau2b.hat = tau2b.Gibbs[index.max]

d = (2*r*k - k^2 + k + 2*r + 2)/2 # number of unknown parameters in (B,H,V,tau2b)

# Compute the log of the determinant of the posterior covariance matrix of the unknown parameters. Note that we are applying the logarithm transformation to the variance parameters to improve the Laplace-Metropolis approximation.
theta = matrix(NA,nrow=Gbig,ncol=d)
for(g in 1:Gbig)
{
  aux = c()
  for(j in (k+1):r) aux = c(aux,b.j.Gibbs[g,j])
  aux = c(aux,log(H.Gibbs[g]))
  aux = c(aux,log(sigma2.j.Gibbs[g,]))
  aux = c(aux,log(tau2b.Gibbs[g]))    
  theta[g,] = aux
}
cov.theta = cov(theta[burnin:Gbig,])
log.det.cov.theta = determinant(cov.theta,logarithm=TRUE)$modulus

# Compute the contribution of the prior densities. Note that because we are applying the logarithm transformation to the variance parameters, we need to consider the Jacobian of the transformation.
log.prior = 0
# Contribution of b.j
aux = theta[index.max,1:(k*(2*r-k-1)/2)]
log.prior = log.prior - 0.5*(k*(2*r-k-1)/2)*(log(2*pi) + log(tau2b.hat)) - 0.5 * sum(aux^2) / tau2b.hat

# Contribution of log(h). Note: because of the Jacobian, we add log(h) to log(prior).
log.prior = log.prior + k*(0.5*n.H*log(n.s2.H/2)-lgamma(0.5*n.H)) + sum(-0.5*n.H*log(H.Gibbs[index.max]) - n.s2.H/(2*H.Gibbs[index.max]))

# Contribution of log(sigma2.j)
log.prior = log.prior + r*(0.5*n.sigma*log(n.s2.sigma/2)-lgamma(0.5*n.sigma)) + sum(-0.5*n.sigma*log(sigma2.j.Gibbs[index.max,]) - n.s2.sigma/(2*sigma2.j.Gibbs[index.max,]))

# Contribution of log(tau2b)
log.prior = log.prior + 0.5*n.b*log(n.s2.b/2)-lgamma(0.5*n.b) + sum(-0.5*n.b*log(tau2b.Gibbs[index.max]) - n.s2.b/(2*tau2b.Gibbs[index.max]))

# Compute the logarithm of the Metropolis-Laplace approximation to the marginal data density
log.marginal.density.ML = 0.5 * d * log(2*pi) + 0.5 * log.det.cov.theta + max.log.int.lik + log.prior

log.marginal.density.ML
# [1] 89.90676
