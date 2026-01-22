# Load in the kernel densities
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_alpha.RData")
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_beta.RData")
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_G.RData")
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_lambda.RData")
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_A.RData")
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_w.RData")
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_C.RData")
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_phi.RData")
load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_T.RData")

# Sampling functions from each parameter's kernel density
cdf_alpha <- cumsum(ker_alpha$y)/sum(ker_alpha$y)
sample_from_dens_alpha <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_alpha,ker_alpha$x,xout=u)$y
  return(sampled_values)
}

cdf_beta <- cumsum(ker_beta$y)/sum(ker_beta$y)
sample_from_dens_beta <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_beta,ker_beta$x,xout=u)$y
  return(sampled_values)
}

cdf_G <- cumsum(ker_G$y)/sum(ker_G$y)
sample_from_dens_G <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_G,ker_G$x,xout=u)$y
  return(sampled_values)
}

cdf_lambda <- cumsum(ker_lambda$y)/sum(ker_lambda$y)
sample_from_dens_lambda <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_lambda,ker_lambda$x,xout=u)$y
  return(sampled_values)
}

cdf_A <- cumsum(ker_A$y)/sum(ker_A$y)
sample_from_dens_A <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_A,ker_A$x,xout=u)$y
  return(sampled_values)
}

cdf_w <- cumsum(ker_w$y)/sum(ker_w$y)
sample_from_dens_w <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_w,ker_w$x,xout=u)$y
  return(sampled_values)
}

cdf_C <- cumsum(ker_C$y)/sum(ker_C$y)
sample_from_dens_C <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_C,ker_C$x,xout=u)$y
  return(sampled_values)
}

cdf_phi <- cumsum(ker_phi$y)/sum(ker_phi$y)
sample_from_dens_phi <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_phi,ker_phi$x,xout=u)$y
  return(sampled_values)
}

cdf_T <- cumsum(ker_T$y)/sum(ker_T$y)
sample_from_dens_T <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_T,ker_T$x,xout=u)$y
  return(sampled_values)
}
