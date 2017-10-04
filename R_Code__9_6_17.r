#require(MASS) # for LDA - not currently needed
require(plotly) # for nice plot at the end

# TRACK COMPUTING TIME OF THIS SIMULATION
ptm_total <- proc.time()

# PARAMETERS
total_time <- 162 # in ms # 1000
timesteps <- 100 # 1000
dt <- total_time/timesteps
time_silence <- 81 # in ms # 500
timestep_silence <- time_silence/total_time*timesteps

num_sims <- 30 # number of simulations for each 'pixel', 30

s <- sample(c(-1,1),num_sims,replace=TRUE) # input signals
M <- 100 # number of individual neurons considered in this analysis, 100
tau <- rep(10,M)

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

#
# ITERATE OVER TWO VARIABLES TO CREATE GRAPHICS:
#

c_range <- c(0,1.4) # ! Put min and max c here
cap_gamma_range <- c(0,5) # ! Put min and max capital gamma values here
resolution <- 20 # how many pixels wide/tall do you want the resulting graphic to be?, 20
mat7a <- matrix(nrow=resolution, ncol=resolution)
matFIM <- matrix(nrow=resolution, ncol=resolution)

# LOOP FOR EACH COMBINATION of c_avg, capital gamma:

for (w in 1:resolution){
  for (z in 1:resolution){

    ptm <- proc.time() # track the time for each iteration
    c_avg <- c_range[1]+w*((c_range[2]-c_range[1])/resolution)
    cbase <- matrix(rep(c_avg,M^2),nrow=M,ncol=M)
    diag(cbase)<- 0
    cij <- cijnoise + cbase
    
    cap_gamma <- cap_gamma_range[1]+z*((cap_gamma_range[2]-cap_gamma_range[1])/resolution)
    
    # Initialize for this combination of c, capital gamma:
    lr <- matrix(nrow=num_sims, ncol=M)
    x <- array(0,dim=c(timesteps,M,num_sims))
    dx_array <- array(0,dim=c(timesteps,M,num_sims))
    
    # Loop through and compute data for each timestep, for (num_sims) trials:
    for(sim in 1:num_sims) {
      for(t in 1:timesteps) {
        lastx <- if(t>1) x[(t-1),,sim] else rep(1,M) # set the INITIAL STATE of the neurons here, rep(0,M) normally
                                                     # NOTE that I'm using 1s here now, not 0s.
        for(i in 1:M) {
          dx <- (ifelse(t<timestep_silence,s[sim],0) - lastx[i] + rnorm(1,0,cap_gamma) + cij[i,]%*%g(lastx)/M)*dt/tau[i]  #!!!! 2* for troubleshooting!
#                 /sum(cij[i,]) as a divisor for the sigma cij term (last term inside parentheses)
          #       right now trying M*c_avg..
          x[t,i,sim]<- lastx[i] + dx
          dx_array[t,i,sim]<- dx
        }
      }
    }
    
    # My code part, compute FIM #
    x_avg = apply(x, c(1,2), mean) # average x over all simulations
    FIM <- rep(0,num_sims) #initial FIM
    for(sim in 1:num_sims){
      fim <- 0
      for(t in 1:(timesteps-1)){
        sum <- 0
        for(i in 1:M){
          sum <- sum + tau[i]^2/(dt^2*cap_gamma^2) * (x[t+1,i,sim]-(x[t,i,sim]+s[sim]/tau[i]-1/tau[i]*x[t,i,sim]+sum(cij[i,]*g(x[t,,sim]))/M)*dt)*dt/tau[i]
        }
        fim = fim + sum^2
      }
      FIM[sim] <- fim/(timesteps-1)
    }
    FIM <- mean(FIM)
    matFIM[w,z] <- FIM
    
    
    ## 7A ##
    
    lr <- t(g(x[timesteps,,])) # save final x values for each simulation to lr
    la <- data.frame(cbind(lr,s)) # combine row of x values with initial signal
    colnames(la)[M+1]<-'Signal'
    
    # LDA: !!! Not currently in use, because we're using the mean instead. So this section is commented out:
#    ld_obj <- lda(formula = Signal ~ ., data = la) # perform the LDA. Gets better with more training data = more simulations
#    
#    v <- as.vector(unname(ld_obj$scaling)) # pull data out of the LDA result object
#    muneg <- as.vector(ld_obj$means[1,])
#    mupos <- as.vector(ld_obj$means[2,])
#    Cpos <- unname(cov(la[la$Signal==1,1:M]))
#    Cneg <- unname(cov(la[la$Signal==-1,1:M]))
    
    
    # Apply formula for 7a and tabulate:
    llh <- rep(0,num_sims)
    for(i in 1:num_sims){
      llh[i] <- L(lr[i,])
    }
    results <-cbind(la[,M+1],llh)
    results <-cbind(results,(results>0)[,1]==(results>0)[,2])
    mat7a[w,z] <- sum(results[,3]/nrow(results)) # record accuracy to correct location in 7a matrix (mat7a object holds the values for fig 7a)
    print(results)
    
    # Display info tracking for sim - just to track the status of the simulation, and see how efficient it is
    cat("\nc_avg: ",c_avg,"  cap GAMMA: ",cap_gamma)
    plot(rowMeans(x[,,1])) # Prints the trajectory of one simulation, over time, for each pixel (the first simulation)
    cat(" Probability correct: ",mat7a[w,z],"\n\nTime for this iteration: ")# this is the probability we want to print for each pixel in Fig 7a!
    print(proc.time()[3] - ptm[3])
    #image(mat7a) # this updates a preview of the 7a figure as each 'pixel' is computed.
    
    #my code:
    image(log(matFIM))

  }  # End of loops iterating over c, cap gamma values:
}
cat("\nTotal Simulation Time: ",proc.time()[3] - ptm_total[3])

# Make a nice version of Plot 7A:
plot_ly(z = t(log(matFIM)), type = "heatmap", yaxis = "GAMMA", xaxis = "c")

