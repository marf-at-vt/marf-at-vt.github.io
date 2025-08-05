# Analysis of troposphere air daily average temperature at pressure levels 1000, 925, 850, 700, 600, 500, 400, 300, and 250 over Oklahoma from June 1 to August 31, 2015.

# Using a factor model with k factors. Here k = 3.

le = read.csv("temperature-troposphere.csv")
tseq = as.Date(le[,1])
nc.a = t(le[,2:10])

m = apply(nc.a,1,mean)
nc.a.centered = matrix(NA,nrow=nrow(nc.a),ncol=ncol(nc.a))
for(i in 1:9) nc.a.centered[i,] = nc.a[i,] - m[i]

nc.a.standardized = matrix(NA,nrow=nrow(nc.a),ncol=ncol(nc.a))
for(i in 1:9)
{
  s=sqrt(var(nc.a[i,]))
  nc.a.standardized[i,] = nc.a.centered[i,] / s
}

#######################################################################
#                                                                     #
# Gibbs Sampler for factor model with k factors.                      #
#                                                                     #
# Marco A. R. Ferreira, March 29, 2017.                               #
# Model:                                                              #
# Y_i = B x_i + v_i                                                   #
# where                                                               #
# Y_i: i-th observation of r-dimensional vector of interest           #
# B: r by k factor loadings matrix                                    #
# x_i: i-th realization of the k-dimensional latent factor process    #
# x_i independent N(0, H), i = 1, ..., n                              #
# H = diag(h_1, ..., h_k)                                             #
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

# Here we are reordering the variables so that the first variable is the
# temperature at pressure level 600, the second variable is the 
# temperature at pressure level 1000, and the third variable is the temperature at pressure level 250.

permutation = c(5,1,9,2:4,6:8)

Y = t(nc.a.standardized)[,permutation]

n = nrow(Y)      # Sample size
k =    3         # Number of factors
r = ncol(Y)      # Dimensional of each observational vector

# Checking the eigenvalues of the sample correlation matrix
sample.correlation.matrix = cor(Y)
plot(eigen(sample.correlation.matrix)$values)
# The first component explains about 78% of the variability.

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

b.j.current = matrix(0,nrow=r,ncol=k)
for(l in 1:k) b.j.current[l,l] = 1

x.i.current = matrix(0,nrow=n,ncol=k)
x.i.current[,1] = apply(Y,1,mean)

tau2b.current = 1

h.current = rep(100,k)

sigma2.j.current = rep(1,r)


#######################################
#                                     #
# Define arrays to store Gibbs sample #
#                                     #
#######################################

b.j.Gibbs = array(NA,dim=c(Gbig,r,k))
tau2b.Gibbs = rep(NA,Gbig)  
x.i.Gibbs = array(NA,dim=c(Gbig,n,k))
h.Gibbs = matrix(NA,nrow=Gbig,ncol=k)
sigma2.j.Gibbs = matrix(NA,nrow=Gbig,ncol=r)

b.j.Gibbs[1,,] = b.j.current
tau2b.Gibbs[1] = tau2b.current
x.i.Gibbs[1,,] = x.i.current
h.Gibbs[1,] = h.current
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
  var.aux = solve(diag(1/h.current) + t(b.j.current) %*% diag(1/sigma2.j.current) %*% b.j.current)
  aux.aux =  var.aux %*% t(b.j.current) %*% diag(1/sigma2.j.current)
  chol.var.aux = chol(var.aux)
  for(i in 1:n)
  {
    mean.aux = aux.aux %*% Y[i,]
    x.i.current[i,] = mean.aux + t(chol.var.aux) %*% rnorm(k)
  }  

  # Simulate sigma_j^2, j = 1, ..., r
  n.aux = n.sigma + n
  for(j in 1:r)
  {
    n.s2.aux = n.s2.sigma 
    for (i in 1:n) n.s2.aux = n.s2.aux + (Y[i,j] - b.j.current[j,] %*% x.i.current[i,])^2
    sigma2.j.current[j] = 1 / rgamma(1,shape=n.aux/2,rate=n.s2.aux/2)
  }
  
  # Simulate b_j, j = k+1, ..., r
  XtX = t(x.i.current) %*% x.i.current
  XtY = t(x.i.current) %*% Y
  identity.k = diag(k)  
  for(j in (k+1):r)
  {
    var.aux = solve( XtX / sigma2.j.current[j] + identity.k / tau2b.current)
    mean.aux = var.aux %*% XtY[,j] / sigma2.j.current[j]
    chol.var.aux = chol(var.aux)
    b.j.current[j,] = mean.aux + t(chol.var.aux) %*% rnorm(k)
  }
  # Simulate b_j, j = 2, ..., k
  for(j in 2:k)
  {
    XtX = t(x.i.current[,1:(j-1)]) %*% x.i.current[,1:(j-1)]
    XtY = t(x.i.current[,1:(j-1)]) %*% (Y[,j]-x.i.current[,j])
    var.aux = solve( XtX / sigma2.j.current[j] + diag(j-1) / tau2b.current)
    mean.aux = var.aux %*% XtY / sigma2.j.current[j]
    chol.var.aux = chol(var.aux)  
    b.j.current[j,1:(j-1)] = mean.aux + t(chol.var.aux) %*% rnorm(j-1)
  }

  # Simulate tau_b^2
  n.aux = n.b + k*(2*r-k-1)/2
  n.s2.aux = n.s2.b 
  for(l in 1:k) for(j in (l+1):r) n.s2.aux = n.s2.aux + b.j.current[j,l]^2
  tau2b.current = 1 / rgamma(1,shape=n.aux/2,rate=n.s2.aux/2)
  
  # Simulate h[l], l=1,...,k
  n.aux = n.H + n
  for(l in 1:k)
  {
    n.s2.aux = n.s2.H + sum(x.i.current[,l]^2)
    h.current[l] = 1 / rgamma(1,shape=n.aux/2,rate=n.s2.aux/2)
  }
    
  # Copy current values into storage arrays
  b.j.Gibbs[g,,] = b.j.current
  tau2b.Gibbs[g] = tau2b.current
  x.i.Gibbs[g,,] = x.i.current
  h.Gibbs[g,] = h.current
  sigma2.j.Gibbs[g,] = sigma2.j.current

}

######################
#                    #
#      Results       #
#                    #
######################

# Trace plots
plot(h.Gibbs[,1],type="l",ylab="H",xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(h.Gibbs[,2],type="l",ylab="H",xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(h.Gibbs[,3],type="l",ylab="H",xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(tau2b.Gibbs,type="l",ylab=expression(tau^2),xlab="iteration",cex=4.0,main="",lwd=2.0)

burnin = 1000

hist(h.Gibbs[burnin:Gbig,1],xlab="H",cex=4.0,main="",lwd=2.0)

hist(h.Gibbs[burnin:Gbig,2],xlab="H",cex=4.0,main="",lwd=2.0)

hist(h.Gibbs[burnin:Gbig,3],xlab="H",cex=4.0,main="",lwd=2.0)

hist(tau2b.Gibbs[burnin:Gbig],xlab=expression(tau^2),cex=4.0,main="",lwd=2.0)


pressure = c(1000, 925, 850, 700, 600, 500, 400, 300, 250)[permutation]

aux = cbind(b.j.Gibbs[burnin:Gbig, 1, 1],rep(pressure[1], Gbig - burnin + 1)) 
for(j in 2:9) aux = rbind(aux,cbind(b.j.Gibbs[burnin:Gbig, j,1],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(b[j]),cex=2,main="Loadings, factor 1")
abline(h=1,lty=2)

aux = cbind(b.j.Gibbs[burnin:Gbig, 1, 2],rep(pressure[1], Gbig - burnin + 1)) 
for(j in 2:9) aux = rbind(aux,cbind(b.j.Gibbs[burnin:Gbig, j,2],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(b[j]),cex=2,main="Loadings, factor 2")
abline(h=1,lty=2)

aux = cbind(b.j.Gibbs[burnin:Gbig, 1, 3],rep(pressure[1], Gbig - burnin + 1)) 
for(j in 2:9) aux = rbind(aux,cbind(b.j.Gibbs[burnin:Gbig, j,3],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(b[j]),cex=2,main="Loadings, factor 3")
abline(h=1,lty=2)

aux = cbind(sigma2.j.Gibbs[burnin:Gbig, 1],rep(pressure[1], Gbig - burnin + 1)) 
for(j in 2:9) aux = rbind(aux,cbind(sigma2.j.Gibbs[burnin:Gbig, j],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(sigma[j]^2),cex=2,ylim=c(0,1.7),main="Idiosyncratic variances")
abline(h=1,lty=2)

ll = min(Y)
ul = max(Y)  
plot(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,1],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlab="Date",ylab="",main="Common factor 1")
lines(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,1],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,1],2,quantile,p=0.975),lty=2,lwd=2)


ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,2],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlab="Date",ylab="",main="Common factor 2")
lines(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,2],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,2],2,quantile,p=0.975),lty=2,lwd=2)

ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,3],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlab="Date",ylab="",main="Common factor 3")
lines(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,3],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.i.Gibbs[burnin:Gbig,,3],2,quantile,p=0.975),lty=2,lwd=2)

percentage.unexplained.variability = apply(sigma2.j.Gibbs[burnin:Gbig,],1,sum)/9
hist(percentage.unexplained.variability,xlab="",ylab="",xlim=c(0,0.4),main="Percentage unexplained variability")


#####################################################################
#####################################################################
##                                                                 ##
##                                                                 ##
## For model selection, compute a Laplace-Metropolis approximation ##
## to the marginal data density.                                   ##
##                                                                 ##
##                                                                 ##
#####################################################################
#####################################################################

######################################################
# Compute the logarithm of the integrated likelihood #
# at each iteration of the MCMC                      #
######################################################

log.integrated.likelihood = function(Y,B,H,V)
{
  n = nrow(Y)
  cov.mat = B %*% H %*% t(B) + V
  result = -(n/2) * determinant(cov.mat,logarithm=TRUE)$modulus
  cov.mat.inv = solve(cov.mat)
  for(i in 1:n) result = result - 0.5 * t(Y[i,]) %*% cov.mat.inv %*% Y[i,]
  result
}

log.int.lik = rep(NA,Gbig)
for(g in 1:Gbig)
{
  B = b.j.Gibbs[g,,]
  H = diag(h.Gibbs[g,])
  V = diag(sigma2.j.Gibbs[g,])
  log.int.lik[g] = log.integrated.likelihood(Y,B,H,V)
}

burnin = 2000

max.log.int.lik = max(log.int.lik)
index.max = which(log.int.lik == max.log.int.lik)
B.hat = b.j.Gibbs[index.max,,]
H.hat = diag(h.Gibbs[index.max,])
V.hat = diag(sigma2.j.Gibbs[index.max,])
tau2b.hat = tau2b.Gibbs[index.max]

d = (2*r*k - k^2 + k + 2*r + 2)/2 # number of unknown parameters in (B,H,V,tau2b)

# Compute the log of the determinant of the posterior covariance matrix of the unknown parameters. Note that we are applying the logarithm transformation to the variance parameters to improve the Laplace-Metropolis approximation.
theta = matrix(NA,nrow=Gbig,ncol=d)
for(g in 1:Gbig)
{
  aux = c()
  for(j in 2:k) aux = c(aux,b.j.Gibbs[g,j,1:(j-1)])
  for(j in (k+1):r) aux = c(aux,b.j.Gibbs[g,j,])
  aux = c(aux,log(h.Gibbs[g,]))
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
log.prior = log.prior + k*(0.5*n.H*log(n.s2.H/2)-lgamma(0.5*n.H)) + sum(-0.5*n.H*log(h.Gibbs[index.max,]) - n.s2.H/(2*h.Gibbs[index.max,]))

# Contribution of log(sigma2.j)
log.prior = log.prior + r*(0.5*n.sigma*log(n.s2.sigma/2)-lgamma(0.5*n.sigma)) + sum(-0.5*n.sigma*log(sigma2.j.Gibbs[index.max,]) - n.s2.sigma/(2*sigma2.j.Gibbs[index.max,]))

# Contribution of log(tau2b)
log.prior = log.prior + 0.5*n.b*log(n.s2.b/2)-lgamma(0.5*n.b) + sum(-0.5*n.b*log(tau2b.Gibbs[index.max]) - n.s2.b/(2*tau2b.Gibbs[index.max]))

# Compute the logarithm of the Metropolis-Laplace approximation to the marginal data density
log.marginal.density.ML = 0.5 * d * log(2*pi) + 0.5 * log.det.cov.theta + max.log.int.lik + log.prior

log.marginal.density.ML
#[1] 326.498
