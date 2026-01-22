# Set parameter estimates
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

saveRDS(P_over_time,"~/PhD-work/Chapter 4/Data/Within-macrophage/Within_macro_stochastic_master_results.rds")

# Save plot with 600 DPI resolution
tiff("MasterEq1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(1:n_max, P_over_time[1001, ], type = "h", col = "blue",
     xlab = "Population size", ylab = "Probability",
     main = paste("Probability distribution at t =", time[1001]))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("MasterEq2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(1:n_max, P_over_time[2501, ], type = "h", col = "blue",
     xlab = "Population size", ylab = "Probability",
     main = paste("Probability distribution at t =", time[2501]))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("MasterEq3.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(1:n_max, P_over_time[length(time), ], type = "h", col = "blue",
     xlab = "Population size", ylab = "Probability",
     main = paste("Probability distribution at t =", t_max))

# Close the device
dev.off()


png("MasterEq.png", width = 8, height = 6, units = "in", res = 600)

# Row-wise normalization
P_norm <- t(apply(P_over_time, 1, function(x) x / max(x)))

# Convert to long format
df <- as.data.frame.table(P_norm)
colnames(df) <- c("Row", "Col", "Value")

ggplot(df, aes(x = as.numeric(Row)*100, y = as.numeric(Col), fill = Value)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = colorRampPalette(c("white", "darkred"))(100),  # White → darkred
    name = "Normalised\nprobability mass"
  ) +
  scale_x_continuous(
    name = "Time post-infection in hours",
    labels = function(x) x / 10000,
    limits = c(0, 100 * 10000)
  ) +
  scale_y_continuous(
    name = "Intracellular Legionella population",
    limits = c(0, 250)
  ) +
  ggtitle("Master-equation solution of within-macrophage model")
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.margin = margin(5, 5, 4, 2)
  )

dev.off()