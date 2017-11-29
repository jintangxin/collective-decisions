########################simulate to get final states for all neurons###########################################
require(plotly) # for nice plot at the end

# TRACK COMPUTING TIME OF THIS SIMULATION
ptm_total <- proc.time()

# PARAMETERS
total_time <- 810 #1620 #162 # in ms # 1000 # (BCD final value should be 1620 ms)
timesteps <- 500 #500 #100 # 1000
dt <- total_time/timesteps
time_silence <- 810 #81 # in ms # 500 # (BCD final value should be 810 ms)
timestep_silence <- time_silence/total_time*timesteps

num_sims <- 100 # 2 #30 # number of simulations for each 'pixel'

num_s <- 10 # number of different input signals we are looking at

s <- seq(-0.4,0.4,length.out=num_s) # input signals 
#s <- 0
M <- 10 # number of individual neurons considered in this analysis
tau <- rep(10,M) # in ms

# BCD create list of times for plots
times <- array(1:timesteps)*dt

# connectivity matrix (!fixed across simulations)
gamma <- 0 # .03, fixed gamma, used in building cij (connectivity heterogeneity / noise)
cijnoise <- matrix(rnorm(M^2,0,gamma),nrow=M,ncol=M) #matrix(c_avg+rnorm(M^2,0,gamma),nrow=M,ncol=M) # are these row/columns correct?
# Make the matrix symmetrical, so that cij[a,b] == cij[b,a], and cij[a,a]==0:
cijnoise[lower.tri(cijnoise)] = t(cijnoise)[lower.tri(cijnoise)]
diag(cijnoise)<-0


# FUNCTIONS
g <- function(x){ # firing rate function - I changed this from (1-tanh(x)).
  (tanh(x))
}
L <- function(r){ # LDA conversion function - right now, we're just using the average instead, and it seems fine.
  #  -.5*log(v%*%Cpos%*%v) - (((r-mupos)%*%v)^2)/(2*v%*%Cpos%*%v) + .5*log(v%*%Cneg%*%v) + (((r-muneg)%*%v)^2)/(2*v%*%Cneg%*%v)
  mean(r)
}

xf <- array(0,dim=c(num_s,num_sims))

for (k in 1:num_s){
  #
  # simulate over different input signals
  #
  
  c_range <- c(1.,1.2) # ! Put min and max c here 0,1.4
  cap_gamma_range <- c(1.,2.0) # ! Put min and max capital gamma values here 0,5
  resolution <- 1 # how many pixels wide/tall do you want the resulting graphic to be? 
  mat7a <- matrix(nrow=resolution, ncol=resolution)

  
  ptm <- proc.time() # track the time for each iteration
  c_avg <- 10 #! here just use c=1.2
  cbase <- matrix(rep(c_avg,M^2),nrow=M,ncol=M)
  diag(cbase)<- 0
  cij <- cijnoise + cbase
      
  cap_gamma <- cap_gamma_range[2] #! here just use cap_gamma=2
      
  # Initialize for this combination of c, capital gamma:
  lr <- matrix(nrow=num_sims, ncol=M)
  x <- array(0,dim=c(timesteps,M,num_sims))
  dx_array <- array(0,dim=c(timesteps,M,num_sims))
      
  # Loop through and compute data for each timestep, for (num_sims) trials:
  for(sim in 1:num_sims) {
    for(t in 1:timesteps) {
      lastx <- if(t>1) x[(t-1),,sim] else rep(0,M) # set the INITIAL STATE of the neurons here, rep(0,M) normally
      # NOTE that I'm using 1s here now, not 0s.
      for(i in 1:M) {
        dx <- ((ifelse(t<timestep_silence,s[k],0) - lastx[i] + cij[i,]%*%g(lastx)/M)*dt + rnorm(1,0,sqrt(dt*cap_gamma)))/tau[i]  #!!!! 2* for troubleshooting!
        #                 /sum(cij[i,]) as a divisor for the sigma cij term (last term inside parentheses)
        #       right now trying M*c_avg..
        x[t,i,sim]<- lastx[i] + dx
        dx_array[t,i,sim]<- dx 
      }
    }
  }
      
  # average x at final state over all neurons;
   xf[k,] <- apply(x[timesteps,,],c(2),mean)
      
}


############################plot##################
require(ggplot2)
plot(range(xf),c(0,5),type='n',main='distribution plot for different s', xlab='xf', ylab='density')
lines(density(xf[5,]), type='l', col='red') # s=0
lines(density(xf[1,]), type='l', col='black') # s= -0.4
lines(density(xf[2,]), type='l', col='yellow') # s= -0.3
lines(density(xf[3,]), type='l', col='blue') # s= -0.2
lines(density(xf[4,]), type='l', col='purple') # s= -0.1
lines(density(xf[6,]), type='l', col='green') # s= 0.1
lines(density(xf[7,]), type='l', col='cyan') # s= 0.2
lines(density(xf[8,]), type='l', col='coral') # s= 0.3
lines(density(xf[9,]), type='l', col='brown') # s= 0.4

###############################
plot(NULL,xlim=c(min(times),max(times)),ylim=c(-10,10), main = c("c_avg: ",c_avg,"cGAM: ",cap_gamma),xlab='time (ms)',ylab='mean neuron state')
#    title(main = c("c_avg: ",c_avg,"cGAM: ",cap_gamma), ps=2)
for(i in 1:num_sims)
  {points(times,rowMeans(x[,,i]),type="o", pch=".")} # prints out results of ea simulation, over time, for each pixel (the first simulation)


#############KL divergence using FNN
require(FNN)
KLMatrix <- array(0,dim=c(num_s,num_s))
for (i in 1:num_s){
  for (j in 1:num_s){
    KLMatrix[i,j] <- KL.divergence(X=xf[i,],Y=xf[j,], k=5)[5]
  }
}

###############KL divergence using kernel density estimator
density <- function(x, data){
  #density esimator from data using sd=h at new points x
  h <- 1
  n <- length(x) 
  y <- array(0,dim=n)
  for (i in 1:n){
    y[i] <- 1/h * mean(dnorm(x[i], mean=data, sd=h))
  }
  return(y)
}

kl_divergence <- function(data1, data2){
  # compute kl-divergence between data1 and data2
  # using their density estimation
  x = seq(-10,10,length.out = 1000) #initial points when compute kl-divergence
  p <- density(x, data1)
  q <- density(x, data2)
  return(sum(p*log(p/q)))
}

# initialize a matrix to store kl_divergence between different s
KLMatrix <- array(0,dim=c(num_s,num_s))
for (i in 1:num_s){
  for (j in 1:num_s){
    KLMatrix[i,j] <- kl_divergence(xf[i,],xf[j,])
  }
}



