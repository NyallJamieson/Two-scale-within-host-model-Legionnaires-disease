# Burr parameters are a=0.6327, b=11.0200, T=21.7964
dburrNEW <- function(x,a,b,T){(x/b+a)*(T/x)^(a)*exp((T-x)/b)/(x*(1+(T/x)^(a)*exp((T-x)/b))^2)}
pburrNEW <- function(x,a,b,T){(1+(T/x)^(a)*exp((T-x)/b))^(-1)}

P_over_time <- readRDS("~/PhD-work/Chapter 4/Data/Within-macrophage/Within_macro_stochastic_master_results.rds")
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

AAA <- read.csv("~/PhD-work/Chapter 4/Data/Within-macrophage/G_dist.csv")
P4_G_intermediate <- density(AAA[,2], from = 0, to = 300)

cols <- RColorBrewer::brewer.pal(9, "Set1")
cols

load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_G.RData")

# Save plot with 600 DPI resolution
png("RupturePmf.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(size,rupture_dist,xlab = "Rupture size (number of Legionella)",ylab="Probability mass",main="Rupture size distribution", col = cols[1], lwd=2,type="l")
lines(c(59.58144, 59.58144), c(0, 1), col = cols[2], lwd = 2)
lines(seq(0, 300, by = 1), dpois(seq(0, 300, by = 1), lambda = 59.58144), col = cols[3], lwd = 2)
lines(P4_G_intermediate, col = cols[4], lwd = 2)
lines(c(62.326, 62.326), c(0, 1), col = cols[5], lwd = 2, lty = 2)
lines(seq(0,300,by=1),dnbinom(seq(0,300,by=1),mu=62.792,size=32.072),lwd=2, col = cols[6], lty = 2)
lines(ker_G$x, ker_G$y, col = cols[7], lwd = 2, lty = 2)

# Legend
legend("topright", col = c(cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8], cols[9]), legend = c("Full G distribution", "G averaged", "G Poisson", "G intermediate", "[7] model A", "[7] model B", "[7] model C"), lty = c(1, 1, 1, 1, 2, 2, 2), cex = 1.5, lwd = 2)

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("RuptureCdf.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(size,cdf,xlab="Rupture size",ylab="Probability",main="Cumulative rupture size distribution",col="red",lwd=2,type="l")
lines(seq(0,300,by=1),pnbinom(seq(0,300,by=1),mu=62.792,size=32.072),lwd=2,col="blue")

# Legend
legend("bottomright",col=c("red","blue"),legend=c("Within-macrophage model","Markov model"),lty=1,cex=1.5)

# Close the device
dev.off()


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

# your existing table
tbl <- table(round(AAA[, 2]))

# define full support
support <- 1:300

# create full PMF counts with zeros for missing values
counts_full <- integer(length(support))
names(counts_full) <- support
counts_full[names(tbl)] <- as.integer(tbl)
pmf <- as.numeric(counts_full / length(AAA[,2]))

100 * (1 - sum(pmf * (beta / (alpha + beta))^(1:300)))
