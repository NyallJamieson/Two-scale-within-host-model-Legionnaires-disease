# Now to the actual code

params <- c(alpha=0.089,beta=1.088,lambda=0.126,G=59.58144,theta=0.8554645)
results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(results) <- c("G","time")

thresh <- 50661

xstart <- c(time=0,L=9,M1=0,M2=0,M3=0)

for (i in 1:1000){

  print(i)
  
  # Within-macrophage script
  # Set up normal distribution for parameters
  mean <- c(177.9529526, 0.1924158)
  vcov <- matrix(c(82.49019097, -0.0254690939, -0.0254690939, 0.0001933546), nrow = 2)
  
  # Sample estimates
  sample <- MASS::mvrnorm(n = 1, mu = mean, Sigma = vcov)
  
  # Set parameter values
  psi <- sample[1]
  omega <- sample[2]
  
  param_means <- c(0.6327181, 11.0200301, 21.7964397)
  param_cov <- matrix(c(0.002182265, 0.012146420, -0.001563877, 0.012146420, 0.075391171, -0.009600021, -0.001563877, -0.009600021, 0.004187794), nrow = 3)
  vals <- MASS::mvrnorm(n = 1, mu = param_means, Sigma = param_cov)
  
  a_val <- vals[1]
  b_val <- vals[2]
  T_val <- vals[3]
  
  # Burr distribution and its derivative
  pdburr <- function(x, a, b, T) {
    (1 + (T / x)^(a) * exp((T - x) / b))^(-1)
  }
  
  # Initialise model run
  n_max <- 300
  t_max <- 300
  dt <- 1e-02
  
  # Q2[i, j] means from state j into state i
  Q <- matrix(0, nrow = n_max, ncol = n_max)
  
  # Births
  Q[2, 1] <- omega * (1 - 1 / psi)
  for (i in 2:(n_max - 1)){
    Q[i + 1, i] <- omega * i
  }
  
  # Deaths
  Q[n_max - 1, n_max] <- omega * n_max * (n_max / psi - 1)
  for (i in 2:(n_max - 1)){
    Q[i - 1, i] <- omega * i ^ 2 / psi
  }
  
  # Set diagonal values
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
  # Runge Kutta method for solving master equations
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
  rupture_pdf_discretized <- append(c(0),dt*dburrNEW(x=time,a=a_val, b=b_val, T=T_val))
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
  
  rate_val <- rnorm(n = 1, mean = 0.1260891, sd = 0.0007486)
  
  mean <- mean(sample_G(10000))
  
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
        F_t <- pgamma(ages, shape = 3, rate = rate_val) # CDF at time t
        F_t_plus_dt <- pgamma(dt + ages, shape = 3, rate = rate_val)  # CDF at t + dt
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
      
      if (L > 50661){
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
      output <- c(mean,time)
      names(output) <- c("G","time")
      return(output)
    }
  }
  
  results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,simulate_nMtau_R,mc.cores=150)),na.action="omit")))
}

write.csv(results, "SA_Erlang_Intermediate_G.csv")