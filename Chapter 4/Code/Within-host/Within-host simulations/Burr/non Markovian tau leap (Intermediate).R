source("~/PhD-work/Chapter 4/Code/Within-macrophage/Estimate G for Pois in Erlang within-macrophage model.R")

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
  threshold <- 50661
  
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
      F_t <- pdburr(ages, a = 0.6327, b = 11.0200, T = 21.7964)  # CDF at time t
      F_t_plus_dt <- pdburr(dt + ages, a = 0.6327, b = 11.0200, T = 21.7964)  # CDF at t + dt
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
    output <- c(dose,time)
    names(output) <- c("L","time")
    return(output)
  }
  
  if (cleared ==1){
    return(c(NA,NA))
  }
}

results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(results) <- c("L","time")

for (j in 255:500){
  print(j)
  xstart <- c("L"=j)
  results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,simulate_nMtau_R,mc.cores=15)),na.action="omit")))
}
