# Dynamic factor model (with 3 factors)

# Analysis of troposphere air daily average temperature at pressure levels 1000, 925, 850, 700, 600, 500, 400, 300, and 250 over Oklahoma from June 1 to August 31, 2015.

library(Matrix)

le = read.csv("temperature-troposphere.csv")
tseq = as.Date(le[,1])
nc.a = t(le[,2:10])

m = apply(nc.a,1,mean)
nc.a.centered = matrix(NA,nrow=nrow(nc.a),ncol=ncol(nc.a))
for(i in 1:9) nc.a.centered[i,] = nc.a[i,] - m[i]

s = rep(NA,9)
nc.a.standardized = matrix(NA,nrow=nrow(nc.a),ncol=ncol(nc.a))
for(i in 1:9)
{
  s[i]=sqrt(var(nc.a[i,]))
  nc.a.standardized[i,] = nc.a.centered[i,] / s[i]
}



#######################################################################
#                                                                     #
# Gibbs Sampler for factor model with k factors.                      #
#                                                                     #
# Marco A. R. Ferreira, March 29, 2017.                               #
# Model:                                                              #
# Y_i = B x_i + v_i                                                   #
# x_i = G x_{i-t} + w_i
# where                                                               #
# Y_i: i-th observation of r-dimensional vector of interest           #
# B: r by k factor loadings matrix                                    #
# x_i: i-th realization of the k-dimensional latent factor process    #
# G = diag(phi_1,...,phi_k)
# W_i = diag(w_1, ..., w_k)                                             #
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
phi.order = c(2,2,1) # AR order for the latent factors
phi.max.dim = max(phi.order)  # Maximum AR order for the latent factors
theta.length = sum(phi.order)


FF = c(1,rep(0,phi.order[1]-1)) # observation matrix
if(k>1) for(l in 2:k) FF = bdiag(FF,c(1,rep(0,phi.order[l]-1)))
FF = t(FF)
FF = matrix(FF,nrow=nrow(FF))


##########################################
##########################################
##                                      ##
##             Gibbs Sampler            ##
##                                      ##
##########################################
##########################################

##############################
#                            #
# Define auxiliary functions #
#                            #
##############################

matrix.sqrt <- function(H)
{
  # Computes square root of nonnegative definite symmetric matrix using spectral decomposition
  
  if(nrow(H)==1) {H.sqrt = matrix(sqrt(H),nrow=1,ncol=1)} else
  {
    H.eigen = eigen(H)
    H.eigen.values = H.eigen$values    
    H.eigen.values[abs(H.eigen$values) < 10^(-10)] = 0
    H.sqrt = H.eigen$vectors %*% diag(sqrt(H.eigen.values)) %*% t(H.eigen$vectors)
  }  

  H.sqrt
}

matrix.Moore.Penrose <- function(H)
{
  # Computes Moore Penrose generalized inverse of symmetric matrix using spectral decomposition
  
  H.eigen = eigen(H)
  inverse.values = rep(0,nrow(H))
  inverse.values[abs(H.eigen$values) > 10^(-10)] = 1/H.eigen$values[abs(H.eigen$values) > 10^(-10)]
  H.MP = H.eigen$vectors %*% diag(inverse.values) %*% t(H.eigen$vectors)
  
  H.MP
}

theta.simulation <- function(TT,y,sigma2,m0,C0,FF,GG,W)
{
  # Simulate the latent process

  p = ncol(GG)
  r = length(sigma2)

  # Define auxiliary matrices (T by p matrices):
  a = matrix(0,nrow=TT,ncol=p)
  m = matrix(0,nrow=TT,ncol=p)
  A = array(0,dim=c(TT,p,r))
  f = matrix(0,nrow=TT,ncol=r)
  e = matrix(0,nrow=TT,ncol=r)

  # Define auxiliary arrays (T by p by arrays):
  R = array(0,dim=c(TT,p,p))
  C = array(0,dim=c(TT,p,p))
  Q = array(0,dim=c(TT,r,r))

  # This matrix will store the simulated theta_t, t=1,...,T:
  theta <- matrix(0,nrow=TT,ncol=p)

  # Kalman filter:
  a[1,] <- as.vector(GG %*% m0)
  R[1,,] <- GG %*% C0 %*% t(GG) + W
  f[1,] <- FF %*% a[1,]
  Q[1,,] <- FF %*% R[1,,] %*% t(FF) + diag(sigma2)
  A[1,,] <- R[1,,] %*% t(FF) %*% solve(Q[1,,])
  e[1,] <- y[1,] - f[1,]
  m[1,] <- a[1,] + A[1,,] %*% e[1,]
  C[1,,] <- R[1,,] - A[1,,] %*% Q[1,,] %*% t(A[1,,])

  for(t in 2:TT)
  {
    a[t,] <- as.vector(GG %*% m[t-1,])
    R[t,,] <- GG %*% C[t-1,,] %*% t(GG) + W
    f[t,] <- FF %*% a[t,]
    Q[t,,] <- FF %*% R[t,,] %*% t(FF) + diag(sigma2)
    A[t,,] <- R[t,,] %*% t(FF) %*% solve(Q[t,,])
    e[t,] <- y[t,] - f[t,]
    m[t,] <- a[t,] + A[t,,] %*% e[t,]
    C[t,,] <- R[t,,] - A[t,,] %*% Q[t,,] %*% t(A[t,,])

  }

  # Now, the backward sampler:
  C.sqrt = matrix.sqrt(C[TT,,])
  theta[TT,] <- m[TT,] + C.sqrt %*% rnorm(p,mean=0,sd=1)

  for(t in (TT-1):1)
  {
    B <- C[t,,] %*% t(GG) %*% matrix.Moore.Penrose(R[t+1,,])
    h <- m[t,] + B %*% (theta[t+1,]-a[t+1,])
    H <- C[t,,] - B %*% R[t+1,,] %*% t(B)   

    H.sqrt = matrix.sqrt(H)

    theta[t,] <- h + H.sqrt %*% rnorm(p,mean=0,sd=1)
  }

  theta
}

phi.simulation <- function(TT,p,x,w)
{
  # Simulate the vector of autoregressive coefficients
  # This assumes a noninformative prior p(phi) propto 1

  n = TT - p # number of terms in the conditional likelihood
  X = matrix(NA,nrow=n,ncol=p)
  for(i in 1:n) X[i,] = x[(TT-i):(TT-i-p+1)]

  cov.matrix = solve(t(X) %*% X) 
  mean.vector = solve(t(X) %*% X) %*% t(X) %*% x[TT:(p+1)]

  cov.matrix.sqrt = matrix.sqrt(cov.matrix)
  phi = mean.vector + sqrt(w) * cov.matrix.sqrt %*% rnorm(p)
  phi
}

w.simulation <- function(Tbig,theta,GG,a.prior,b.prior)
{ # This is a Gibbs step.

  # Dimension of state-space vector:
  p <- ncol(theta)

  a <- a.prior + 0.5 * (Tbig - p) 
  b <- b.prior
  
  for(t in (p+1):Tbig)
  {
    aux = theta[t,] - GG %*% theta[t-1,]
    b <- b + 0.5 * aux[1]^2
  }
  w <- 1.0 / rgamma(1,a,b)
  w
}


##############################
#                            #
# Set priors hyperparameters #
#                            #
##############################

n.b = 1
n.s2.b = 1
n.w = 1
n.s2.w = 1
m0 <- rep(0,theta.length)
C0 <- 10^5 * diag(theta.length)


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

b.current = matrix(0,nrow=r,ncol=k)
for(l in 1:k) b.current[l,l] = 1

x.current = matrix(0,nrow=n,ncol=k)
x.current[,1] = apply(Y,1,mean)

theta.current = matrix(0,nrow=n,ncol=theta.length)
theta.current[,1] = x.current[,1]

tau2.current = 1

sigma2.current = rep(1,r)

phi.max.dim = 2
phi.current = matrix(NA,nrow=k,ncol=phi.max.dim)
phi.current[1,] = c(1.05,-0.29)
phi.current[2,] = c(0.73,-0.24)
phi.current[3,] = c(0.28,0)

w.current = c(0.30,0.22,0.2)

#######################################
#                                     #
# Define arrays to store Gibbs sample #
#                                     #
#######################################

b.Gibbs = array(NA,dim=c(Gbig,r,k))
sigma2.Gibbs = matrix(NA,nrow=Gbig,ncol=r)
tau2.Gibbs = rep(NA,Gbig)  
x.Gibbs = array(NA,dim=c(Gbig,n,k))
theta.Gibbs = array(NA,dim=c(Gbig,n,theta.length))
phi.Gibbs = array(NA,dim=c(Gbig,k,phi.max.dim))
w.Gibbs = matrix(NA,nrow=Gbig,ncol=k)

b.Gibbs[1,,] = b.current
sigma2.Gibbs[1,] = sigma2.current
tau2.Gibbs[1] = tau2.current
x.Gibbs[1,,] = x.current
theta.Gibbs[1,,] = theta.current
phi.Gibbs[1,,] = phi.current
w.Gibbs[1,] = w.current

# Compute evolution matrix GG
GG = rbind(phi.current[1,1:phi.order[1]],cbind(diag(phi.order[1]-1),rep(0,phi.order[1]-1)))  
if(k>1) for(l in 2:k) GG = bdiag(GG,rbind(phi.current[l,1:phi.order[l]],cbind(diag(phi.order[l]-1),rep(0,phi.order[l]-1))))
GG = matrix(GG,nrow=nrow(GG))

    
# Compute matrix W
aux = matrix(0,nrow=phi.order[1],ncol=phi.order[1])
aux[1,1] = w.current[1]
W = aux
if(k>1) for(l in 2:k)
{
  aux = matrix(0,nrow=phi.order[l],ncol=phi.order[l])
  aux[1,1] = w.current[l]
  W = bdiag(W,aux)    
}
W = matrix(W,nrow=nrow(W))

######################
#                    #
# Gibbs Sampler loop #
#                    #
######################

for(g in 2:Gbig)
{

  cat("g=",g,"\n")

  # Simulate theta_t, t = 1, ..., n
  theta.current = theta.simulation(n,Y,sigma2.current,m0,C0,b.current%*%FF,GG,W)


  # Simulate x_t, t = 1, ..., n
  for(t in 1:n) x.current[t,] = FF %*% theta.current[t,]

  # Simulate sigma_j^2, j = 1, ..., r
  n.aux = n.sigma + n
  for(j in 1:r)
  {
    n.s2.aux = n.s2.sigma 
    for (i in 1:n) n.s2.aux = n.s2.aux + (Y[i,j] - b.current[j,] %*% x.current[i,])^2
    sigma2.current[j] = 1 / rgamma(1,shape=n.aux/2,rate=n.s2.aux/2)
  }

  # Simulate b_j, j = k+1, ..., r
  XtX = t(x.current) %*% x.current
  XtY = t(x.current) %*% Y
  identity.k = diag(k)  
  for(j in (k+1):r)
  {
    var.aux = solve( XtX / sigma2.current[j] + identity.k / tau2.current)
    mean.aux = var.aux %*% XtY[,j] / sigma2.current[j]
    chol.var.aux = chol(var.aux)
    b.current[j,] = mean.aux + t(chol.var.aux) %*% rnorm(k)
  }
  # Simulate b_j, j = 2, ..., k
  for(j in 2:k)
  {
    XtX = t(x.current[,1:(j-1)]) %*% x.current[,1:(j-1)]
    XtY = t(x.current[,1:(j-1)]) %*% (Y[,j]-x.current[,j])
    var.aux = solve( XtX / sigma2.current[j] + diag(j-1) / tau2.current)
    mean.aux = var.aux %*% XtY / sigma2.current[j]
    chol.var.aux = chol(var.aux)  
    b.current[j,1:(j-1)] = mean.aux + t(chol.var.aux) %*% rnorm(j-1)
  }

  # Simulate tau_b^2
  n.aux = n.b + k*(2*r-k-1)/2
  n.s2.aux = n.s2.b 
  for(l in 1:k) for(j in (l+1):r) n.s2.aux = n.s2.aux + b.current[j,l]^2
  tau2.current = 1 / rgamma(1,shape=n.aux/2,rate=n.s2.aux/2)

  # Simulate phi_l, l = 1, ..., k
  for(l in 1:k) 
  {
    if(l==1){start.index = 1} else {start.index = 1 + sum(phi.order[1:(l-1)])}
    aux = theta.current[,start.index]
    phi.current[l,1:phi.order[l]] = phi.simulation(n,phi.order[l],aux,w.current[l])
  }
  
  # Simulate w_l, l = 1, ..., k
  for(l in 1:k) 
  {
    if(l==1){start.index = 1} else {start.index = 1 + sum(phi.order[1:(l-1)])}
    end.index = sum(phi.order[1:l])
    thetap = matrix(theta.current[,start.index:end.index],ncol=phi.order[l])   # theta.current[,start.index:end.index]
    GGp = rbind(phi.current[l,1:phi.order[l]],cbind(diag(phi.order[l]-1),rep(0,phi.order[l]-1)))
    w.current[l] = w.simulation(n,thetap,GGp,n.w/2,n.s2.w/2)
  }
  
  # Compute evolution matrix GG
  GG = rbind(phi.current[1,1:phi.order[1]],cbind(diag(phi.order[1]-1),rep(0,phi.order[1]-1)))  
  if(k>1) for(l in 2:k) GG = bdiag(GG,rbind(phi.current[l,1:phi.order[l]],cbind(diag(phi.order[l]-1),rep(0,phi.order[l]-1))))
  GG = matrix(GG,nrow=nrow(GG))

  # Compute matrix W
  aux = matrix(0,nrow=phi.order[1],ncol=phi.order[1])
  aux[1,1] = w.current[1]
  W = aux
  if(k>1) for(l in 2:k)
  {
    aux = matrix(0,nrow=phi.order[l],ncol=phi.order[l])
    aux[1,1] = w.current[l]
    W = bdiag(W,aux)    
  }
  W = matrix(W,nrow=nrow(W))

  # Copy current values into storage arrays
  b.Gibbs[g,,] = b.current
  sigma2.Gibbs[g,] = sigma2.current
  tau2.Gibbs[g] = tau2.current
  x.Gibbs[g,,] = x.current
  theta.Gibbs[g,,] = theta.current
  phi.Gibbs[g,,] = phi.current
  w.Gibbs[g,] = w.current

}
    

######################
#                    #
#      Results       #
#                    #
######################

# Trace plots
plot(w.Gibbs[,1],type="l",ylab="w1",xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(w.Gibbs[,2],type="l",ylab="w2",xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(w.Gibbs[,3],type="l",ylab="w3",xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(tau2.Gibbs,type="l",ylab=expression(tau^2),xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(phi.Gibbs[,1,1],type="l",ylab=expression(phi[11]),xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(phi.Gibbs[,1,2],type="l",ylab=expression(phi[12]),xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(phi.Gibbs[,2,1],type="l",ylab=expression(phi[21]),xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(phi.Gibbs[,2,2],type="l",ylab=expression(phi[22]),xlab="iteration",cex=4.0,main="",lwd=2.0)

plot(phi.Gibbs[,3,1],type="l",ylab=expression(phi[31]),xlab="iteration",cex=4.0,main="",lwd=2.0)


burnin = 1000

hist(w.Gibbs[burnin:Gbig,1],xlab="w1",cex=4.0,main="",lwd=2.0)

hist(w.Gibbs[burnin:Gbig,2],xlab="w2",cex=4.0,main="",lwd=2.0)

hist(w.Gibbs[burnin:Gbig,3],xlab="w3",cex=4.0,main="",lwd=2.0)

hist(tau2.Gibbs[burnin:Gbig],xlab=expression(tau^2),cex=4.0,main="",lwd=2.0)

hist(phi.Gibbs[burnin:Gbig,1,1],ylab="",xlab=expression(phi[11]),cex=4.0,main="",lwd=2.0)

hist(phi.Gibbs[burnin:Gbig,1,2],xlab=expression(phi[12]),ylab="",cex=4.0,main="",lwd=2.0)

hist(phi.Gibbs[burnin:Gbig,2,1],xlab=expression(phi[21]),ylab="",cex=4.0,main="",lwd=2.0)

hist(phi.Gibbs[burnin:Gbig,2,2],xlab=expression(phi[22]),ylab="",cex=4.0,main="",lwd=2.0)

hist(phi.Gibbs[burnin:Gbig,3,1],xlab=expression(phi[31]),ylab="",cex=4.0,main="",lwd=2.0)


pressure = c(1000, 925, 850, 700, 600, 500, 400, 300, 250)[permutation]

# Boxplots of posterior draws for different b_j's
aux = cbind(b.Gibbs[burnin:Gbig, 1, 1],rep(pressure[1], Gbig - burnin + 1)) 
for(j in 2:9) aux = rbind(aux,cbind(b.Gibbs[burnin:Gbig, j,1],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(b[j]),cex=2,main="Loadings, factor 1")
abline(h=1,lty=2)

aux = cbind(b.Gibbs[burnin:Gbig, 1, 2],rep(pressure[1], Gbig - burnin + 1)) 
for(j in 2:9) aux = rbind(aux,cbind(b.Gibbs[burnin:Gbig, j,2],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(b[j]),cex=2,main="Loadings, factor 2")
abline(h=1,lty=2)

aux = cbind(b.Gibbs[burnin:Gbig, 1, 3],rep(pressure[1], Gbig - burnin + 1)) 
for(j in 2:9) aux = rbind(aux,cbind(b.Gibbs[burnin:Gbig, j,3],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(b[j]),cex=2,main="Loadings, factor 3")
abline(h=1,lty=2)

aux = cbind(sigma2.Gibbs[burnin:Gbig, 1],rep(pressure[1], Gbig - burnin + 1)) 
for(j in 2:9) aux = rbind(aux,cbind(sigma2.Gibbs[burnin:Gbig, j],rep(pressure[j], Gbig - burnin + 1)))
boxplot(aux[,1] ~ aux[,2],xlab="Pressure level",ylab=expression(sigma[j]^2),cex=2,ylim=c(0,1.7),main="Idiosyncratic variances")
abline(h=1,lty=2)


ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,1],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlab="Date",ylab="",main="Common factor 1")
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,1],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,1],2,quantile,p=0.975),lty=2,lwd=2)

ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,2],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlab="Date",ylab="",main="Common factor 2")
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,2],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,2],2,quantile,p=0.975),lty=2,lwd=2)

ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,3],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlab="Date",ylab="",main="Common factor 3")
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,3],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,3],2,quantile,p=0.975),lty=2,lwd=2)

percentage.unexplained.variability = apply(sigma2.Gibbs[burnin:Gbig,],1,sum)/9
hist(percentage.unexplained.variability,xlab="",ylab="",xlim=c(0,0.4),main="Percentage unexplained variability")

# Analyzing the autoregressive coefficients of the dynamic common factors:
mean(phi.Gibbs[burnin:Gbig,1,1])
mean(phi.Gibbs[burnin:Gbig,1,2])
mean(phi.Gibbs[burnin:Gbig,2,1])
mean(phi.Gibbs[burnin:Gbig,2,2])
mean(phi.Gibbs[burnin:Gbig,3,1])

# Computing the wavelength or period for first and second common factors
library(polynom)

# For common factor 1:
phi <- c(mean(phi.Gibbs[burnin:Gbig,1,1]),mean(phi.Gibbs[burnin:Gbig,1,2]))
charac.polyn.coeff <- c(1,-phi)
charac.polyn <- polynomial(charac.polyn.coeff)  # Characteristic polynomial
charac.polyn.roots <- solve(charac.polyn) # Compute roots of characteristic polynomial
modulus.roots <- Mod(charac.polyn.roots) # Compute moduli of the roots 
modulus.recipr.roots <- 1/Mod(charac.polyn.roots) # Compute moduli of reciprocal roots
omega <- Arg(charac.polyn.roots) # Compute argument of the roots and reciprocal roots
wavelength <- 2 * pi / abs(omega) # Compute wavelength or period
wavelength[1]

phi <- c(mean(phi.Gibbs[burnin:Gbig,2,1]),mean(phi.Gibbs[burnin:Gbig,2,2]))
charac.polyn.coeff <- c(1,-phi)
charac.polyn <- polynomial(charac.polyn.coeff)  # Characteristic polynomial
charac.polyn.roots <- solve(charac.polyn) # Compute roots of characteristic polynomial
modulus.roots <- Mod(charac.polyn.roots) # Compute moduli of the roots 
modulus.recipr.roots <- 1/Mod(charac.polyn.roots) # Compute moduli of reciprocal roots
omega <- Arg(charac.polyn.roots) # Compute argument of the roots and reciprocal roots
wavelength <- 2 * pi / abs(omega) # Compute wavelength or period
wavelength[1]

###############################################
# Forecasting common factors and temperatures #
###############################################

forecast.window = 10
theta.forecast.Gibbs = array(NA,dim=c(Gbig,forecast.window,theta.length))
y.forecast.Gibbs = array(NA,dim=c(Gbig,forecast.window,r))

for(g in 1:Gbig)
{
  b.current = b.Gibbs[g,,]
  sigma2.current = sigma2.Gibbs[g,]
  tau2.current = tau2.Gibbs[g]
  x.current = x.Gibbs[g,,]
  theta.current = theta.Gibbs[g,,]
  phi.current = phi.Gibbs[g,,]
  w.current = w.Gibbs[g,]

  # Compute evolution matrix GG
  GG = rbind(phi.current[1,1:phi.order[1]],cbind(diag(phi.order[1]-1),rep(0,phi.order[1]-1)))  
  if(k>1) for(l in 2:k) GG = bdiag(GG,rbind(phi.current[l,1:phi.order[l]],cbind(diag(phi.order[l]-1),rep(0,phi.order[l]-1))))
  GG = matrix(GG,nrow=nrow(GG))

  # Compute matrix W
  aux = matrix(0,nrow=phi.order[1],ncol=phi.order[1])
  aux[1,1] = w.current[1]
  W = aux
  if(k>1) for(l in 2:k)
  {
    aux = matrix(0,nrow=phi.order[l],ncol=phi.order[l])
    aux[1,1] = w.current[l]
    W = bdiag(W,aux)    
  }
  W = matrix(W,nrow=nrow(W))

  t=1
  theta.forecast.Gibbs[g,t,] = GG %*% theta.current[n,] +  matrix.sqrt(W) %*% rnorm(theta.length,mean=0,sd=1)
  y.forecast.Gibbs[g,t,] = b.current %*% FF %*% theta.forecast.Gibbs[g,t,] + rnorm(r,mean=0,sd=sqrt(sigma2.current))

  for(t in 2: forecast.window)
  {
    theta.forecast.Gibbs[g,t,] = GG %*% theta.forecast.Gibbs[g,t-1,] +  matrix.sqrt(W) %*% rnorm(theta.length,mean=0,sd=1)
    y.forecast.Gibbs[g,t,] = b.current %*% FF %*% theta.forecast.Gibbs[g,t,] + rnorm(r,mean=0,sd=sqrt(sigma2.current))
  }
}

time.s=as.POSIXct('2015-06-01',tz='UTC')
time.e=as.POSIXct('2015-08-31',tz='UTC')
tseq=seq(from=time.s, to=time.e, by="24 hours")
aux = tseq[length(tseq)] + 60*60*24*(1:forecast.window)

ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,1],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlim=c(as.Date(tseq[1],"%D"),as.Date(aux[forecast.window],"%D")),xlab="Date",ylab="",main="Forecast - Common factor 1")
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,1],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,1],2,quantile,p=0.975),lty=2,lwd=2)
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,1],2,mean),lty=1,lwd=2,col="red")
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,1],2,quantile,p=0.025),lty=2,lwd=2,col="red")
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,1],2,quantile,p=0.975),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,1],2,mean)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,1],2,mean)[1]),lty=1,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,1],2,quantile,p=0.025)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,1],2,quantile,p=0.025)[1]),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,1],2,quantile,p=0.975)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,1],2,quantile,p=0.975)[1]),lty=2,lwd=2,col="red")
abline(h=0)

ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,2],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlim=c(as.Date(tseq[1],"%D"),as.Date(aux[forecast.window],"%D")),xlab="Date",ylab="",main="Forecast - Common factor 2")
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,2],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,2],2,quantile,p=0.975),lty=2,lwd=2)
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,3],2,mean),lty=1,lwd=2,col="red")
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,3],2,quantile,p=0.025),lty=2,lwd=2,col="red")
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,3],2,quantile,p=0.975),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,2],2,mean)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,3],2,mean)[1]),lty=1,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,2],2,quantile,p=0.025)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,3],2,quantile,p=0.025)[1]),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,2],2,quantile,p=0.975)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,3],2,quantile,p=0.975)[1]),lty=2,lwd=2,col="red")
abline(h=0)

ll = min(Y)
ul = max(Y)
plot(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,3],2,mean),type="l",lwd=2,ylim=c(ll,ul),xlim=c(as.Date(tseq[1],"%D"),as.Date(aux[forecast.window],"%D")),xlab="Date",ylab="",main="Forecast - Common factor 3")
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,3],2,quantile,p=0.025),lty=2,lwd=2)
lines(as.Date(tseq,"%D"),apply(x.Gibbs[burnin:Gbig,,3],2,quantile,p=0.975),lty=2,lwd=2)
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,5],2,mean),lty=1,lwd=2,col="red")
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,5],2,quantile,p=0.025),lty=2,lwd=2,col="red")
lines(as.Date(aux,"%D"),apply(theta.forecast.Gibbs[burnin:Gbig,,5],2,quantile,p=0.975),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,3],2,mean)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,5],2,mean)[1]),lty=1,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,3],2,quantile,p=0.025)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,5],2,quantile,p=0.025)[1]),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux[1]),"%D"),c(apply(x.Gibbs[burnin:Gbig,,3],2,quantile,p=0.975)[length(tseq)],apply(theta.forecast.Gibbs[burnin:Gbig,,5],2,quantile,p=0.975)[1]),lty=2,lwd=2,col="red")
abline(h=0)

ll = min(nc.a)
ul = max(nc.a)
  
# Altitude pressure level 600 milibars
plot(as.Date(tseq,"%D"),m[5]+s[5]*Y[,1],type="l",lwd=2,ylim=c(ll,ul),xlim=c(as.Date(tseq[1],"%D"),as.Date(aux[forecast.window],"%D")),xlab="Date",ylab="Kelvin",main="Forecast - temperature at 1000 600 400 and 250 milibars")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[5]+s[5]*c(Y[length(tseq),1],apply(y.forecast.Gibbs[burnin:Gbig,,1],2,mean)),lty=1,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[5]+s[5]*c(Y[length(tseq),1],apply(y.forecast.Gibbs[burnin:Gbig,,1],2,quantile,p=0.025)),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[5]+s[5]*c(Y[length(tseq),1],apply(y.forecast.Gibbs[burnin:Gbig,,1],2,quantile,p=0.975)),lty=2,lwd=2,col="red")
  
# Altitude pressure level 1000 milibars  
  lines(as.Date(tseq,"%D"),m[1]+s[1]*Y[,2],lty=1,lwd=2)
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[1]+s[1]*c(Y[length(tseq),2],apply(y.forecast.Gibbs[burnin:Gbig,,2],2,mean)),lty=1,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[1]+s[1]*c(Y[length(tseq),2],apply(y.forecast.Gibbs[burnin:Gbig,,2],2,quantile,p=0.025)),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[1]+s[1]*c(Y[length(tseq),2],apply(y.forecast.Gibbs[burnin:Gbig,,2],2,quantile,p=0.975)),lty=2,lwd=2,col="red")

# Altitude pressure level 250 milibars  
lines(as.Date(tseq,"%D"),m[9]+s[9]*Y[,3],lty=1,lwd=2)
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[9]+s[9]*c(Y[length(tseq),3],apply(y.forecast.Gibbs[burnin:Gbig,,3],2,mean)),lty=1,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[9]+s[9]*c(Y[length(tseq),3],apply(y.forecast.Gibbs[burnin:Gbig,,3],2,quantile,p=0.025)),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[9]+s[9]*c(Y[length(tseq),3],apply(y.forecast.Gibbs[burnin:Gbig,,3],2,quantile,p=0.975)),lty=2,lwd=2,col="red")

# Altitude pressure level 400 milibars  
lines(as.Date(tseq,"%D"),m[7]+s[7]*Y[,8],lty=1,lwd=2)
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[7]+s[7]*c(Y[length(tseq),8],apply(y.forecast.Gibbs[burnin:Gbig,,8],2,mean)),lty=1,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[7]+s[7]*c(Y[length(tseq),8],apply(y.forecast.Gibbs[burnin:Gbig,,8],2,quantile,p=0.025)),lty=2,lwd=2,col="red")
lines(as.Date(c(tseq[length(tseq)],aux),"%D"),m[7]+s[7]*c(Y[length(tseq),8],apply(y.forecast.Gibbs[burnin:Gbig,,8],2,quantile,p=0.975)),lty=2,lwd=2,col="red")
  

