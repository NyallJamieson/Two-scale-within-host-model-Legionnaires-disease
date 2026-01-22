get_sk_ku <- function(x) {
  data.frame(
    skewness = PerformanceAnalytics::skewness(x),
    kurtosis = PerformanceAnalytics::kurtosis(x, method = "moment")
  )
}

CF_E_ave_1 <- Erlang_averaged[which(Erlang_averaged[,1] > 100), 2]
CF_E_ave_2 <- Erlang_averaged[which(Erlang_averaged[,1] > 200), 2]
CF_E_ave_3 <- Erlang_averaged[which(Erlang_averaged[,1] > 300), 2]
CF_E_ave_4 <- Erlang_averaged[which(Erlang_averaged[,1] > 400), 2]
CF_E_ave_5 <- Erlang_averaged[which(Erlang_averaged[,1] == 500), 2]
CF_E_pois_1 <- Erlang_pois[which(Erlang_pois[,1] > 100), 2]
CF_E_pois_2 <- Erlang_pois[which(Erlang_pois[,1] > 200), 2]
CF_E_pois_3 <- Erlang_pois[which(Erlang_pois[,1] > 300), 2]
CF_E_pois_4 <- Erlang_pois[which(Erlang_pois[,1] > 400), 2]
CF_E_pois_5 <- Erlang_pois[which(Erlang_pois[,1] == 500), 2]
CF_E_int_1 <- Erlang_intermediate[which(Erlang_intermediate[,1] > 100), 2]
CF_E_int_2 <- Erlang_intermediate[which(Erlang_intermediate[,1] > 200), 2]
CF_E_int_3 <- Erlang_intermediate[which(Erlang_intermediate[,1] > 300), 2]
CF_E_int_4 <- Erlang_intermediate[which(Erlang_intermediate[,1] > 400), 2]
CF_E_int_5 <- Erlang_intermediate[which(Erlang_intermediate[,1] == 500), 2]

CF_B_ave_1 <- Burr_averaged[which(Burr_averaged[,1] > 100), 2]
CF_B_ave_2 <- Burr_averaged[which(Burr_averaged[,1] > 200), 2]
CF_B_ave_3 <- Burr_averaged[which(Burr_averaged[,1] > 300), 2]
CF_B_ave_4 <- Burr_averaged[which(Burr_averaged[,1] > 400), 2]
CF_B_ave_5 <- Burr_averaged[which(Burr_averaged[,1] == 500), 2]
CF_B_pois_1 <- Burr_pois[which(Burr_pois[,1] > 100), 2]
CF_B_pois_2 <- Burr_pois[which(Burr_pois[,1] > 200), 2]
CF_B_pois_3 <- Burr_pois[which(Burr_pois[,1] > 300), 2]
CF_B_pois_4 <- Burr_pois[which(Burr_pois[,1] > 400), 2]
CF_B_pois_5 <- Burr_pois[which(Burr_pois[,1] == 500), 2]
CF_B_int_1 <- Burr_intermediate[which(Burr_intermediate[,1] > 100), 2]
CF_B_int_2 <- Burr_intermediate[which(Burr_intermediate[,1] > 200), 2]
CF_B_int_3 <- Burr_intermediate[which(Burr_intermediate[,1] > 300), 2]
CF_B_int_4 <- Burr_intermediate[which(Burr_intermediate[,1] > 400), 2]
CF_B_int_5 <- Burr_intermediate[which(Burr_intermediate[,1] == 500), 2]


# Apply to each dataset separately and combine
results <- do.call(rbind, list(
  CF_E_ave_1 = get_sk_ku(CF_E_ave_1),
  CF_E_ave_2 = get_sk_ku(CF_E_ave_2),
  CF_E_ave_3 = get_sk_ku(CF_E_ave_3),
  CF_E_ave_4 = get_sk_ku(CF_E_ave_4),
  CF_E_ave_5 = get_sk_ku(CF_E_ave_5),
  CF_E_pois_1 = get_sk_ku(CF_E_pois_1),
  CF_E_pois_2 = get_sk_ku(CF_E_pois_2),
  CF_E_pois_3 = get_sk_ku(CF_E_pois_3),
  CF_E_pois_4 = get_sk_ku(CF_E_pois_4),
  CF_E_pois_5 = get_sk_ku(CF_E_pois_5),
  CF_E_int_1 = get_sk_ku(CF_E_int_1),
  CF_E_int_2 = get_sk_ku(CF_E_int_2),
  CF_E_int_3 = get_sk_ku(CF_E_int_3),
  CF_E_int_4 = get_sk_ku(CF_E_int_4),
  CF_E_int_5 = get_sk_ku(CF_E_int_5),
  CF_B_ave_1 = get_sk_ku(CF_B_ave_1),
  CF_B_ave_2 = get_sk_ku(CF_B_ave_2),
  CF_B_ave_3 = get_sk_ku(CF_B_ave_3),
  CF_B_ave_4 = get_sk_ku(CF_B_ave_4),
  CF_B_ave_5 = get_sk_ku(CF_B_ave_5),
  CF_B_pois_1 = get_sk_ku(CF_B_pois_1),
  CF_B_pois_2 = get_sk_ku(CF_B_pois_2),
  CF_B_pois_3 = get_sk_ku(CF_B_pois_3),
  CF_B_pois_4 = get_sk_ku(CF_B_pois_4),
  CF_B_pois_5 = get_sk_ku(CF_B_pois_5),
  CF_B_int_1 = get_sk_ku(CF_B_int_1),
  CF_B_int_2 = get_sk_ku(CF_B_int_2),
  CF_B_int_3 = get_sk_ku(CF_B_int_3),
  CF_B_int_4 = get_sk_ku(CF_B_int_4),
  CF_B_int_5 = get_sk_ku(CF_B_int_5)
))

print(results)




library(PerformanceAnalytics)
library(moments)
library(ggplot2)

# Your function
get_sk_ku <- function(x) {
  data.frame(
    skewness = PerformanceAnalytics::skewness(x),
    kurtosis = PerformanceAnalytics::kurtosis(x, method = "moment")
  )
}

# Put all datasets into a named list
datasets <- list(
  CF_E_ave_1 = CF_E_ave_1,
  CF_E_ave_2 = CF_E_ave_2,
  CF_E_ave_3 = CF_E_ave_3,
  CF_E_ave_4 = CF_E_ave_4,
  CF_E_ave_5 = CF_E_ave_5,
  CF_E_pois_1 = CF_E_pois_1,
  CF_E_pois_2 = CF_E_pois_2,
  CF_E_pois_3 = CF_E_pois_3,
  CF_E_pois_4 = CF_E_pois_4,
  CF_E_pois_5 = CF_E_pois_5,
  CF_E_int_1 = CF_E_int_1,
  CF_E_int_2 = CF_E_int_2,
  CF_E_int_3 = CF_E_int_3,
  CF_E_int_4 = CF_E_int_4,
  CF_E_int_5 = CF_E_int_5
)

# Apply function to all and bind into one df
results <- do.call(rbind, lapply(names(datasets), function(nm) {
  cbind(dataset = nm, get_sk_ku(datasets[[nm]]))
}))

print(results)

# Plot Cullen–Frey style scatter
ggplot(results, aes(x = skewness, y = kurtosis, color = dataset)) +
  geom_point(size = 3) +
  geom_text(aes(label = dataset), hjust = 1.2, vjust = 1.2, size = 3) +
  labs(
    title = "Cullen–Frey style plot",
    x = "Skewness",
    y = "Kurtosis"
  ) +
  theme_minimal()





# New plots
library(ggplot2)
library(PerformanceAnalytics)

# List of your datasets
datasets <- list(
  CF_E_ave_1, CF_E_ave_2, CF_E_ave_3, CF_E_ave_4, CF_E_ave_5,
  CF_E_pois_1, CF_E_pois_2, CF_E_pois_3, CF_E_pois_4, CF_E_pois_5,
  CF_E_int_1, CF_E_int_2, CF_E_int_3, CF_E_int_4, CF_E_int_5
)

# Names of the datasets
dataset_names <- c(
  "Erlang average 1", "Erlang average 2", "Erlang average 3", "Erlang average 4", "Erlang average 5",
  "Erlang Poisson 1", "Erlang Poisson 2", "Erlang Poisson 3", "Erlang Poisson 4", "Erlang Poisson 5",
  "Erlang intermediate 1", "Erlang intermediate 2", "Erlang intermediate 3", "Erlang intermediate 4", "Erlang intermediate 5"
)

# Compute skewness and kurtosis for each dataset
get_sk_ku <- function(x) {
  data.frame(
    skewness = PerformanceAnalytics::skewness(x),
    kurtosis = PerformanceAnalytics::kurtosis(x, method = "moment")
  )
}

results <- do.call(rbind, lapply(datasets, get_sk_ku))
results$dataset <- dataset_names

# Define hierarchical colors - COLOUR SAFE VERSION

library(colorspace)

# Base colors (colorblind-safe)
base_colors <- c("#E69F00", "#56B4E9", "#CC79A7")  # orange, blue, purple

# Generate 5 shades per base color
colors <- unlist(lapply(base_colors, function(col) {
  lighten(col, seq(0.4, 0, length.out = 5))  # from light to base color
}))

print(colors)

results$color <- colors

# Plot with hierarchical coloring
ggplot(results, aes(x = skewness, y = kurtosis)) +
  geom_point(aes(color = dataset), size = 4) +
  scale_color_manual(values = setNames(colors, results$dataset)) +
  labs(
    title = "Cullen-Frey Points with Hierarchical Colors",
    x = "Skewness",
    y = "Kurtosis",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 10)
  )

# --- Define theoretical curves ---
# Gamma curve: kurtosis = skewness^2 + 3
skew_seq <- seq(0, 1.1, length.out = 200)
gamma_curve <- data.frame(
  skewness = skew_seq,
  kurtosis = skew_seq^2 + 3
)

png("CF1.png", width = 8, height = 6, units = "in", res = 600)

# --- Plot ---
ggplot() +
  # Gamma and Lognormal curves
  geom_line(data = gamma_curve, aes(x = skewness, y = kurtosis),
            color = "grey", size = 1) +
  
  # Your hierarchical points
  geom_point(data = results, aes(x = skewness, y = kurtosis, color = dataset), size = 4) +
  
  # Manual color mapping
  scale_color_manual(values = setNames(colors, results$dataset)) +
  
  labs(
    title = "Cullen-Frey diagram (Erlang)",
    x = "Skewness",
    y = "Kurtosis",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    # Title size
    plot.title = element_text(size = 20, face = "bold"),
    # Axis labels size
    axis.title = element_text(size = 20),
    # Axis tick labels size
    axis.text = element_text(size = 14),
    # Legend title size
    legend.title = element_text(size = 16),
    # Legend text size
    legend.text = element_text(size = 14)
  ) +
  coord_cartesian(xlim = c(0, 1.1), ylim = c(2.5, 5))

dev.off()




# Now Burr plot
# List of your datasets
datasets <- list(
  CF_B_ave_1, CF_B_ave_2, CF_B_ave_3, CF_B_ave_4, CF_B_ave_5,
  CF_B_pois_1, CF_B_pois_2, CF_B_pois_3, CF_B_pois_4, CF_B_pois_5,
  CF_B_int_1, CF_B_int_2, CF_B_int_3, CF_B_int_4, CF_B_int_5
)

# Names of the datasets
dataset_names <- c(
  "Burr average 1", "Burr average 2", "Burr average 3", "Burr average 4", "Burr average 5",
  "Burr Poisson 1", "Burr Poisson 2", "Burr Poisson 3", "Burr Poisson 4", "Burr Poisson 5",
  "Burr intermediate 1", "Burr intermediate 2", "Burr intermediate 3", "Burr intermediate 4", "Burr intermediate 5"
)

# Compute skewness and kurtosis for each dataset
get_sk_ku <- function(x) {
  data.frame(
    skewness = PerformanceAnalytics::skewness(x),
    kurtosis = PerformanceAnalytics::kurtosis(x, method = "moment")
  )
}

results <- do.call(rbind, lapply(datasets, get_sk_ku))
results$dataset <- dataset_names

results$color <- colors

# Plot with hierarchical coloring
ggplot(results, aes(x = skewness, y = kurtosis)) +
  geom_point(aes(color = dataset), size = 4) +
  scale_color_manual(values = setNames(colors, results$dataset)) +
  labs(
    title = "Cullen-Frey diagram (Burr)",
    x = "Skewness",
    y = "Kurtosis",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 10)
  )

# --- Define theoretical curves ---
# Gamma curve: kurtosis = skewness^2 + 3
skew_seq <- seq(0, 1.1, length.out = 200)
gamma_curve <- data.frame(
  skewness = skew_seq,
  kurtosis = skew_seq^2 + 3
)

png("CF2.png", width = 8, height = 6, units = "in", res = 600)

# --- Plot ---
ggplot() +
  # Gamma and Lognormal curves
  geom_line(data = gamma_curve, aes(x = skewness, y = kurtosis),
            color = "grey", size = 1) +
  
  # Your hierarchical points
  geom_point(data = results, aes(x = skewness, y = kurtosis, color = dataset), size = 4) +
  
  # Manual color mapping
  scale_color_manual(values = setNames(colors, results$dataset)) +
  
  labs(
    title = "Cullen-Frey diagram (Burr)",
    x = "Skewness",
    y = "Kurtosis",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    # Title size
    plot.title = element_text(size = 20, face = "bold"),
    # Axis labels size
    axis.title = element_text(size = 20),
    # Axis tick labels size
    axis.text = element_text(size = 14),
    # Legend title size
    legend.title = element_text(size = 16),
    # Legend text size
    legend.text = element_text(size = 14)
  ) +
  coord_cartesian(xlim = c(0, 1.1), ylim = c(2.5, 5))

dev.off()
