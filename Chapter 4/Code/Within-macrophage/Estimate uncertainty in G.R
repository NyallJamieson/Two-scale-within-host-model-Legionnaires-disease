results <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(results) <- c("G")

G_est <- function(y){
  
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
  
  # Burr parameters are a=0.6327, b=11.0200, T=21.7964
  dburrNEW <- function(x,a,b,T){(x/b+a)*(T/x)^(a)*exp((T-x)/b)/(x*(1+(T/x)^(a)*exp((T-x)/b))^2)}
  pburrNEW <- function(x,a,b,T){(1+(T/x)^(a)*exp((T-x)/b))^(-1)}
  
  # Set up initialisation
  dt <- 1e-02
  t_max <- 300
  time <- seq(0.01,t_max,by=dt)
  
  # Calculate the rupture pdf
  rupture_pdf_discretized <- append(c(0),dt*dburrNEW(x=time,a=a_val, b=b_val, T=T_val))
  rupture_pdf_discretized <- rupture_pdf_discretized/sum(rupture_pdf_discretized)
  
  # Initialise again
  time <- seq(0,t_max,by=dt)
  
  # Finish calculation
  rupture_dist <- c()
  for (i in 1:300){
    rupture_dist <- append(rupture_dist,sum(rupture_pdf_discretized*P_over_time[,i]))
  }
  
  # Estimate rupture size
  size <- seq(1,300,by=1)
  G <- sum(size*rupture_dist)
  
  vec <- c("G" = G)
  return(vec)
}

# Run code in parallel for 1000 simulations
results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,G_est,mc.cores=10)),na.action="omit")))

# Save data for analysis
write.csv(results, "~/PhD-work/Chapter 4/Data/Within-macrophage/G_dist.csv")

# Fit Poisson distribution
fit <- fitdistrplus::fitdist(results[,1], "pois", method = "mle")
