# Within-macrophage script
omega <- 0.19242
psi <- 177.953

n_max <- 300
t_max <- 300
dt <- 1e-02

# Q2[i, j] means from state j into state i
Q <- matrix(0, nrow = n_max, ncol = n_max)

# births
Q[2, 1] <- omega * (1 - 1 / psi)
for (i in 2:(n_max - 1)){
  Q[i + 1, i] <- omega * i
}

Q[n_max - 1, n_max] <- omega * n_max * (n_max / psi - 1)
for (i in 2:(n_max - 1)){
  Q[i - 1, i] <- omega * i ^ 2 / psi
}

for (i in 1:n_max){
  Q[i, i] <- -colSums(Q)[i]
}

# Initial probability distribution (start with n0 individuals)
n0 <- 1  # Starting population size
P0 <- rep(0, n_max)
P0[n0] <- 1  # P(n0, t=0) = 1

# Time evolution
time <- seq(0, t_max, by = dt)
P_over_time <- matrix(0, nrow = length(time), ncol = n_max)
P_over_time[1, ] <- P0  # Initial condition

P <- P0
threshold <- 1e-15
# runge kutta method for solving master equations
for (t in 2:length(time)) {
  # Compute k1, k2, k3, k4
  k1 <- Q %*% P
  k2 <- Q %*% (P + 0.5 * dt * k1)
  k3 <- Q %*% (P + 0.5 * dt * k2)
  k4 <- Q %*% (P + dt * k3)
  
  # Update P using the Runge-Kutta formula
  P <- P + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
  
  # Threshold small probabilities to zero
  P[P < threshold] <- 0
  
  # Normalize to sum to one
  P <- P / sum(P)
  
  # Store the result
  P_over_time[t, ] <- P
}

# Now estimate G code
# burr parameters are a=0.6327, b=11.0200, T=21.7964
dburrNEW <- function(x,a,b,T){(x/b+a)*(T/x)^(a)*exp((T-x)/b)/(x*(1+(T/x)^(a)*exp((T-x)/b))^2)}
pburrNEW <- function(x,a,b,T){(1+(T/x)^(a)*exp((T-x)/b))^(-1)}

dt <- 1e-02
t_max <- 300
time <- seq(0.01,t_max,by=dt)
rupture_pdf_discretized <- append(c(0),dt*dburrNEW(x=time,a=0.6327, b=11.0200, T=21.7964))
rupture_pdf_discretized <- rupture_pdf_discretized/sum(rupture_pdf_discretized)
time <- seq(0,t_max,by=dt)
rupture_dist <- c()
for (i in 1:300){
  rupture_dist <- append(rupture_dist,sum(rupture_pdf_discretized*P_over_time[,i]))
}

size <- seq(1,300,by=1)
G <- sum(size*rupture_dist)
cdf <- cumsum(rupture_dist)

sample_G <- function(sample_size) {
  
  # Generate uniform random values
  u <- runif(n = sample_size)
  
  # Initialize vector to store sampled G values
  G <- numeric(sample_size)
  
  for (i in 1:sample_size) {
    # Find the index where u[i] fits in the CDF
    lower_idx <- length(which(cdf < u[i]))
    upper_idx <- lower_idx + 1
    
    # Handle edge cases
    if (lower_idx == 0) {
      G[i] <- 1  # u[i] is smaller than the first CDF value
    } else if (upper_idx > length(cdf)) {
      G[i] <- length(cdf)  # u[i] is larger than the last CDF value
    } else {
      # Determine whether to pick lower or upper index
      lower_diff <- abs(cdf[lower_idx] - u[i])
      upper_diff <- abs(cdf[upper_idx] - u[i])
      
      if (lower_diff <= upper_diff) {
        G[i] <- lower_idx
      } else {
        G[i] <- upper_idx
      }
    }
  }
  
  return(G)
}

sample_G_time <- function(sample_size,time){
  
  time_vec <- seq(0,150,by=0.01)
  t_row <- which(abs(time_vec - time) < 1e-6)
  
  cdf_time <- cumsum(P_over_time[t_row,])
  
  # Generate uniform random values
  u <- runif(n = sample_size)
  
  # Initialize vector to store sampled G values
  G <- numeric(sample_size)
  
  for (i in 1:sample_size) {
    # Find the index where u[i] fits in the CDF
    lower_idx <- length(which(cdf_time < u[i]))
    upper_idx <- lower_idx + 1
    
    # Handle edge cases
    if (lower_idx == 0) {
      G[i] <- 1  # u[i] is smaller than the first CDF value
    } else if (upper_idx > length(cdf_time)) {
      G[i] <- length(cdf_time)  # u[i] is larger than the last CDF value
    } else {
      # Determine whether to pick lower or upper index
      lower_diff <- abs(cdf_time[lower_idx] - u[i])
      upper_diff <- abs(cdf_time[upper_idx] - u[i])
      
      if (lower_diff <= upper_diff) {
        G[i] <- lower_idx
      } else {
        G[i] <- upper_idx
      }
    }
  }
  
  return(G)
}

# Now to the actual code

# Burr distribution and its derivative
pdburr <- function(x, a, b, T) {
  (1 + (T / x)^(a) * exp((T - x) / b))^(-1)
}

load("ker_T.RData")

cdf_T <- cumsum(ker_T$y)/sum(ker_T$y)
sample_from_dens_T <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_T,ker_T$x,xout=u)$y
  return(sampled_values)
}

params <- c(alpha=0.089,beta=1.088,lambda=0.126,G=59.58144,theta=0.8554645)
results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(results) <- c("threshold","time")

xstart <- c(time=0,L=9,M1=0,M2=0,M3=0)

for (i in 1:1000){
  
  thresh <- sample_from_dens_T(1)
  
  # Burr distribution and its derivative
  pdburr <- function(x, a, b, T) {
    (1 + (T / x)^(a) * exp((T - x) / b))^(-1)
  }
  
  simulate_nMtau_R <- function(y){
    
    dose <- xstart["L"]
    
    # Parameters
    alpha <- 0.089
    beta <- 1.088
    G <- 62
    threshold <- thresh
    
    ill <- 0
    cleared <- 0
    
    dt <- 0.01
    t <- seq(0, 250, by = dt)
    x <- c("L" = dose, "M" = 0)
    ages <- c()
    
    for (i in t) {
      
      L <- x[1]
      M <- x[2]
      
      if (M > 0) {
        # Extract ages
        ages <- x[3:length(x)]
        u <- runif(n = length(ages))  # Uniform random numbers for rupture decision
        
        # Compute the CDF values and the survival function
        F_t <- pgamma(ages, shape = 3, rate = 0.126)  # CDF at time t
        F_t_plus_dt <- pgamma(dt + ages, shape = 3, rate = 0.126)  # CDF at t + dt
        survival_t <- 1 - F_t  # Survival function
        
        # Compute the probability of rupture occurring in the next time step
        probs <- (F_t_plus_dt - F_t) / survival_t
        
        # Rupture decision: If u < probs, rupture occurs
        rupture <- which(u < probs)  # This triggers when the uniform random number is less than the probability
        
        num_rupture <- length(rupture)
        M <- M - num_rupture  # Decrease the number of individuals
        
        if (num_rupture > 0) {
          ages <- ages[-rupture]  # Remove ruptured individuals
        }
        
        # Increment ages
        ages <- ages + dt
        
        # Add new individuals due to rupture events
        L <- L + sum(sample_G(sample_size = num_rupture))
      } else {
        ages <- c()  # If no individuals, reset ages
      }
      
      # Phagocytosis process (Poisson and Binomial processes)
      L_phago <- rpois(n = 1, lambda = (alpha + beta) * L * dt)
      
      if (L_phago > L){
        while (L_phago > L){
          L_phago <- rpois(n = 1, lambda = (alpha + beta) * L * dt)
        }
      }
      phago_survive <- rbinom(n = 1, size = L_phago, p = alpha / (alpha + beta))
      
      L <- L - L_phago  # Update L based on phagocytosis
      M <- M + phago_survive  # Update M based on surviving phagos
      if (exists("ages")) {
        ages <- append(ages, rep(0, phago_survive))  # Add new individuals with age 0
      } else {
        ages <- rep(0, phago_survive)  # If ages doesn't exist, create a new vector
      }
      
      # Update state vector
      x <- c("L" = L, "M" = M, ages)
      names(x)[1:2] <- c("L", "M")
      
      if (L > threshold){
        time <- i
        ill <- 1
        break
      }
      
      if (L==0 & M==0){
        cleared <- 1
        break
      }
    }
    
    if (ill == 1){
      output <- c(threshold,time)
      names(output) <- c("T"=threshold,"time")
      return(output)
    }
  }
  results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,simulate_nMtau_R,mc.cores=128)),na.action="omit")))
}

write.csv(results, "SA_Erlang.Intermediate.T.csv")