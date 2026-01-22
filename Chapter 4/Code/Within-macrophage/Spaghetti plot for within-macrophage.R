# One step of birth death CTMC
MODEL.onestep <- function(x,params){
  L <- x[2]
  beta <- params["beta"]
  C <- params["C"]
  
  if (L==1){
    rates <- c(
      birth <- beta*(1-1/C),
      death <- 0
    )
  } else {
  rates <- c(
    birth <- beta*L,
    death <- beta*L^2/C
  )
  }
  total.rates <- sum(rates)
  if (total.rates==0)
    tau <- -Inf
  else
    tau <- rexp(n=1,rate=total.rates)
  transitions <- list(
    Birth <- c(1),
    Death <- c(-1)
  )
  event <- sample.int(n=2,size=1,prob=rates/total.rates)
  x+c(tau,transitions[[event]])
}

# Run BD simulation
MODEL.simul <- function(y,maxstep=10000000,threshold=150){
  x <- c(time=1,L=1)
  params <- c(beta=0.19242,C=177.95295)
  output <- array(dim=c(maxstep+1,2))
  colnames(output) <- names(x)
  output[1,] <- x
  k <- 1
  while ((k <= maxstep) && (x["L"] > 0) && (x["time"] < threshold)){
    k <- k+1
    output[k,] <- x <- MODEL.onestep(x,params)}
  A <- as.data.frame(output[1:k,])
  
  DATA <- data.frame(matrix(, nrow=14901, ncol=0))
  t <- seq(1,150,by=0.01)
  
  col <- c(1)
  for (j in t){
    col <- append(col,A[length(which(A[,1]<j)),2])
  }
  DATA <- data.frame(DATA,col)
  return(DATA)
}


results <- as.data.frame(do.call(cbind,parallel::mclapply(1:1000,MODEL.simul,mc.cores=15)))
colnames(results) <- paste("col", seq_len(ncol(results)), sep = " ")

# Convert to long format for easier plotting
results_long <- results |>
  dplyr::mutate(time = seq(1, 150, by = 0.01)) |>
  tidyr::pivot_longer(-time, names_to = "trajectory", values_to = "L")

# Calculate standard deviation at each time point (or another measure of "commonness")
trajectory_stats <- results_long |>
  dplyr::group_by(trajectory) |>
  dplyr::summarize(weight = 1 / (1 + sd(L, na.rm = TRUE)))  # Example: Use 1 / (1 + SD)

# Merge weights back into the dataset
results_long <- results_long |>
  dplyr::left_join(trajectory_stats, by = "trajectory") |>
  dplyr::mutate(alpha = scales::rescale(weight, to = c(0.1, 1)))  # Normalize for alpha

# Define a color palette for trajectory fading
color_palette <- colorRampPalette(c("mistyrose","darkred"))

# Map weights to the color palette
results_long$color <- color_palette(100)[as.numeric(cut(results_long$weight, breaks = 100))]

# Unique trajectories
unique_trajectories <- unique(results_long$trajectory)

# Calculate summary statistics for filtering and weighting
trajectory_stats <- results_long |>
  dplyr::group_by(trajectory) |>
  dplyr::summarize(
    weight = 1 / (1 + sd(L, na.rm = TRUE)),  # Weight for fading
    mean_L = mean(L, na.rm = TRUE),         # Mean for filtering
    sd_L = sd(L, na.rm = TRUE)              # SD for filtering
  )

# Filter out near-constant trajectories
filtered_trajectories <- trajectory_stats |>
  dplyr::filter(mean_L > 10, sd_L > 0.01)  # Adjust thresholds as needed

# Merge back filtered trajectories into the dataset
results_long <- results_long |>
  dplyr::inner_join(filtered_trajectories, by = "trajectory") |>
  dplyr::mutate(alpha = scales::rescale(weight.y, to = c(0.1, 1)))  # Rescale for fading

# Define a threshold for "going towards zero"
zero_threshold <- 1e-2  # Adjust as needed

# Save plot with 600 DPI resolution
png("SpaghettiPlot.png", width = 8, height = 6, units = "in", res = 600)

library(Hmisc)  # for minor.tick

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins

# Produce plot
plot(
  NA, 
  xlim = c(0, 100), 
  ylim = range(results_long$L, na.rm = TRUE),
  xlab = "Time post-infection (hours)", 
  ylab = "Intracellular Legionella population", 
  main = "Spaghetti plot of within-macrophage model"
)

# Make spaces smaller
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

# Use a single base color and low alpha for all trajectories
base_col <- "firebrick"

for (trajectory in unique(results_long$trajectory)) {
  traj_data <- results_long[results_long$trajectory == trajectory, ]
  
  if (all(!is.na(traj_data$L))) {
    lines(traj_data$time, traj_data$L,
          col = adjustcolor(base_col, alpha.f = 0.05),
          lwd = 1)
  }
}

# Close the device
dev.off()

# Define a threshold for "going towards zero"
zero_threshold <- 1e-2  # Adjust as needed

# Save plot with 600 DPI resolution
png("SpaghettiPlot2.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(
  NA, 
  xlim = range(results_long$time), 
  ylim = range(results_long$L, na.rm = TRUE),
  xlab = "Time post-infection (hours)", 
  ylab = "Intracellular Legionella population\n(bacteria)", 
  main = "Spaghetti plot of within-macrophage model"
)

# Loop over unique trajectories
for (trajectory in unique_trajectories[1:100]) {
  traj_data <- results_long[results_long$trajectory == trajectory, ]
  
  # Check if the minimum L value is above the threshold
  if (min(traj_data$L, na.rm = TRUE) > zero_threshold) {
    
    # Set a default alpha value to reduce transparency (e.g., 0.5 for 50% opacity)
    alpha_value <- 0.5  # You can adjust this value between 0 and 1 as needed
    
    # Plot only if the trajectory does not go to zero
    lines(traj_data$time, traj_data$L, col = "red", lwd = 1)
  }
}

# Close the device
dev.off()

library(ggplot2)

# Sample 250 trajectories
set.seed(1)
sample_traj <- sample(unique(results_long$trajectory), 250)
results_sub <- results_long[results_long$trajectory %in% sample_traj, ]

# Save plot with same size & DPI as before
ggplot(results_sub, aes(x = time, y = L, group = trajectory)) +
  geom_line(alpha = 0.05, color = "firebrick", linewidth = 0.3) +
  labs(
    x = "Time post-infection (hours)",
    y = "Intracellular Legionella population",
    title = "Spaghetti plot of within-macrophage model"
  ) +
  coord_cartesian(xlim = c(0, 100)) +  # Limit x-axis to 0–100
  theme_minimal(base_size = 22) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )

# Save at 8x6 inches, 600 DPI — matches your base R PNG
ggsave(
  filename = "SpaghettiPlot.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
