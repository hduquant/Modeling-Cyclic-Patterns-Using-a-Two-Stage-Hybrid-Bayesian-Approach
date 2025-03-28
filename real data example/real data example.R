library(MASS)
library(rjags)
# Posterior mode
mode.v <- function(v){
  dens.y <- density(v)
  mode <- dens.y$x[order(dens.y$y,decreasing=T)][1]
  return(mode)}

#  HPD: Monte Carlo method
emp.hpd <- function (theta, alpha = 0.95) {
  alpha <- min(alpha, 1 - alpha) 
  n <- length(theta) 
  L.start <- round(n * alpha) 
  theta <- sort(theta) 
  e <- theta[(n - L.start + 1):n] - theta[1:L.start] 
  ind <- which(e == min(e))[1] 
  return(c(theta[ind], theta[n - L.start + ind])) 
}


data <- read.csv("data.csv",header = TRUE)
data=data[!is.na(data[2]),]
data=data[!is.na(data[3]),]

dim(data)
N=nrow(data)
#w <- 2*pi/rep(7,N)
w <- 2*pi/(data[,3]-1)
y <- data[,5]

ind <- as.integer(factor(data[,1]))
length(ind)
n=max(ind)

x <- as.numeric(data[,2]=="in_relationship")
#x <- data[,14]
x <- tapply(x, ind, function(y) y[1])

t <-  data[,4]


# Step 1: Model
model = "
model {

for (i in 1:n){ # n individuals

betai[i,1:3] ~ dmnorm(mu.beta[i,1:3],omega[1:3,1:3])
mu.beta[i,1]<- mu.D 
mu.beta[i,2]<-r1+beta1*x[i]
mu.beta[i,3]<-r2+beta2*x[i]

Ri[i] <- sqrt((betai[i,2])^2+(betai[i,3])^2)

phii[i]<-ifelse( (betai[i,3]) >= 0 && (betai[i,2]) > 0, arctan((betai[i,3])/(betai[i,2])),
     ifelse( (betai[i,3]) <= 0 && (betai[i,2]) < 0, arctan((betai[i,3])/(betai[i,2])) + pi,
     ifelse( (betai[i,3]) >= 0 && (betai[i,2]) < 0, arctan((betai[i,3])/(betai[i,2])) + pi,
     ifelse( (betai[i,3]) <= 0 && (betai[i,2]) > 0, arctan((betai[i,3])/(betai[i,2])),
     ifelse( (betai[i,3]) >= 0 && (betai[i,2]) == 0, pi / 2,
     ifelse( (betai[i,3]) <= 0 && (betai[i,2]) == 0, 3 * pi / 2,  NAvar))))))

}

 for (j in 1:N){ # t time points
 mu[j]<-betai[ind[j],1] + betai[ind[j],2]*cos(w[j]*t[j]) - betai[ind[j],3]*sin(w[j]*t[j]) 
 y[j] ~ dnorm(mu[j],sigma.pre) 
  }

#### define parameters of interest 
beta <- sum((Ri-mean(Ri))*(x-mean(x)))/sum((x-mean(x))*(x-mean(x))) 
R <- mean(Ri)-beta*mean(x)
Rih <- R+beta*x
rse<- sqrt(sum((Ri-Rih)^2)/(n-2))
se.beta<-(rse/sqrt(sum((x-mean(x))*(x-mean(x)))))^2
se.R<-(rse*sqrt(1/n+mean(x)^2/(sum((x-mean(x))*(x-mean(x))))))^2

phi<- mean(phii)
se.phi<- sd(phii)^2/n

mu.M<-mean(betai[,1])
se.mu<-sd(betai[,1])^2/n

###pirors
sigma.pre ~ dgamma(0.001,0.001)
sigma2 <- 1/sigma.pre 
r1 ~ dnorm (0,0.001) 
r2 ~ dnorm (0,0.001) 
mu.D ~ dnorm (0,0.001) 
beta1 ~ dnorm (0,0.001)
beta2 ~ dnorm (0,0.001) 

omega[1:3,1:3]~dwish(S[1:3,1:3],3) 
Cov[1:3,1:3]<-inverse(omega[1:3,1:3]) 

}
"
# Save model
writeLines( model , con="model.txt" ) #creates a text file object
#------------------------------------------------------------------------------
# Load data 
dataList = list(    # Put the information into a list.
  t = t,
  y = y ,
  x = x ,
  N=N,
  n=n,
  w=w,
  S =diag(3),
  pi=pi,
  NAvar=-999999,
  ind=ind
)

#------------
# Step 2: Specifying a starting value
r1 <- 0.1
# You need to set a seed; otherwise, you will get a different result every time you run the code
initsList1 = list(  r1=r1,
                   .RNG.name="base::Super-Duper", .RNG.seed=1)

initsList2 = list(  r1=r1+0.1,
                    .RNG.name="base::Super-Duper", .RNG.seed=2)
#------------------------------------------------------------------------------ 
# Step 3: Adaptation and burn-in 
parameters = c( "mu.D","beta1","beta2",
                "r1","r2","sigma2","Cov","R","se.R","phi","se.phi","mu.M","se.mu",
                "beta","se.beta")   # Specify the estimated parameters 
adaptSteps =100           # Adaptive period
burnInSteps = 3000      # Burn-in period
nChains = 2
nIter =5000  # The number of kept iterations
jagsModel = jags.model( "model.txt" , data=dataList ,
                        inits=list(initsList1,initsList2), 
                        n.chains=nChains , n.adapt=adaptSteps )
update( jagsModel , n.iter=burnInSteps)
#------------------------------------------------------------------------------ 
# Step 4: Sample from posterior distributions
codaSamples = coda.samples( jagsModel , variable.names=parameters, 
                            n.iter=nIter , thin=1)
convergence=gelman.diag(codaSamples,multivariate=FALSE)[[1]][,1]
convergence

#------------------------------------------------------------------------------ 
# Step 5: Summarize posterior distributions
mcmcChain = as.matrix( codaSamples)
mcmcChain<-mcmc(  mcmcChain)
mean <- apply(mcmcChain, 2, mean)
median <- apply(mcmcChain, 2, median)
mode <- apply(mcmcChain, 2,mode.v)
sd <- apply(mcmcChain, 2,sd)
qbp <- apply(mcmcChain, 2, function (x) quantile(x, c(0.025,0.975)))
hpd <- apply(mcmcChain, 2,function (x) emp.hpd(x, c(0.025,0.975)))

vw.R=mean(mcmcChain[,19])
vb.R=var(mcmcChain[,10])
se.R=sqrt(vw.R+vb.R+vb.R/length(mcmcChain[,10]))
CI.R<-c(mean[10]-1.96*se.R,mean[10]+1.96*se.R)

vw.beta=mean(mcmcChain[,20])
vb.beta=var(mcmcChain[,11])
se.beta=sqrt(vw.beta+vb.beta+vb.beta/length(mcmcChain[,11]))
CI.beta<-c(mean[11]-1.96*se.beta,mean[11]+1.96*se.beta)

vw.mu=mean(mcmcChain[,21])
vb.mu=var(mcmcChain[,15])
se.mu=sqrt(vw.mu+vb.mu+vb.mu/length(mcmcChain[,15]))
CI.mu<-c(mean[15]-1.96*se.mu,mean[15]+1.96*se.mu)

vw.phi=mean(mcmcChain[,22])
vb.phi=var(mcmcChain[,16])
se.phi=sqrt(vw.phi+vb.phi+vb.phi/length(mcmcChain[,16]))
CI.phi<-c(mean[16]-1.96*se.phi,mean[16]+1.96*se.phi)

# Now we can obtain a better summary
Bayes_sum <- rbind(mean,median, mode, sd, qbp, hpd,convergence1)
Bayes_sum <- Bayes_sum[,-c(19:22)]
CI<-cbind(matrix(NA,nrow=2,ncol=9),CI.R,CI.beta,matrix(NA,nrow=2,ncol=3),CI.mu,CI.phi,matrix(NA,nrow=2,ncol=3))
Bayes_sum<-rbind(Bayes_sum,CI)
rownames(Bayes_sum) <- c("mean", "median", "mode", "sd", "qbp.lower", "qbp.upper",
                         "hpd.lower", "hpd.upper","psrf","ci.lower", "ci.upper")

Bayes_sum
