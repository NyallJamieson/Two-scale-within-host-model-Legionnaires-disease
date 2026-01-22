# Read in results
Burr_constant_ab <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G constant/SA_Burr.Constant.ab.csv")
Burr_constant_G <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G constant/SA_Burr.Constant.G.csv")
Burr_constant_lambda <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G constant/SA_Burr.Constant.lambda.csv")
Burr_constant_T <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G constant/SA_Burr.Constant.T.csv")
Burr_pois_ab <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G pois/SA_Burr.Pois.ab.csv")
Burr_pois_G <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G pois/SA_Burr.Pois.G.csv")
Burr_pois_lambda <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G pois/SA_Burr.Pois.lambda.csv")
Burr_pois_T <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G pois/SA_Burr.Pois.T.csv")
Burr_intermediate_ab <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G intermediate/SA_Burr.Intermediate.ab.csv")
Burr_intermediate_G <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G intermediate/SA_Burr.Intermediate.G.csv")
Burr_intermediate_lambda <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G intermediate/SA_Burr.Intermediate.lambda.csv")
Burr_intermediate_T <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Burr/G intermediate/SA_Burr.Intermediate.T.csv")
Erlang_constant_ab <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G constant/SA_Erlang.Constant.AlphaBeta.csv")
Erlang_constant_G <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G constant/SA_Erlang.Constant.G.csv")
Erlang_constant_lambda <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G constant/SA_Erlang.Constant.Lambda.csv")
Erlang_constant_T <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G constant/SA_Erlang.Constant.T.csv")
Erlang_pois_ab <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G pois/SA_Erlang.Pois.ab.csv")
Erlang_pois_G <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G pois/SA_Erlang.Pois.G.csv")
Erlang_pois_lambda <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G pois/SA_Erlang.Pois.lambda.csv")
Erlang_pois_T <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G pois/SA_Erlang.Pois.T.csv")
Erlang_intermediate_ab <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G intermediate/SA_Erlang.Intermediate.ab.csv")
Erlang_intermediate_G <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G intermediate/SA_Erlang_Intermediate_G.csv")
Erlang_intermediate_lambda <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G intermediate/SA_Erlang.Intermediate.lambda.csv")
Erlang_intermediate_T <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Sensitivity analysis/Erlang/G intermediate/SA_Erlang.Intermediate.T.csv")

# Calculate correlations for incubation period
inc_Burr_constant_lambda <- cor(Burr_constant_lambda$lambda,Burr_constant_lambda$time,method="spearman",use = "complete.obs")
inc_Burr_constant_G <- cor(Burr_constant_G$G,Burr_constant_G$time,method="spearman",use = "complete.obs")
inc_Burr_constant_alpha <- cor(Burr_constant_ab$alpha,Burr_constant_ab$time,method="spearman",use = "complete.obs")
inc_Burr_constant_beta <- cor(Burr_constant_ab$beta,Burr_constant_ab$time,method="spearman",use = "complete.obs")
inc_Burr_constant_T <- cor(Burr_constant_T$threshold,Burr_constant_T$time,method="spearman",use = "complete.obs")
inc_Burr_pois_lambda <- cor(Burr_pois_lambda$lambda,Burr_pois_lambda$time,method="spearman",use = "complete.obs")
inc_Burr_pois_G <- cor(Burr_pois_G$G,Burr_pois_G$time,method="spearman",use = "complete.obs")
inc_Burr_pois_alpha <- cor(Burr_pois_ab$alpha,Burr_pois_ab$time,method="spearman",use = "complete.obs")
inc_Burr_pois_beta <- cor(Burr_pois_ab$beta,Burr_pois_ab$time,method="spearman",use = "complete.obs")
inc_Burr_pois_T <- cor(Burr_pois_T$threshold,Burr_pois_T$time,method="spearman",use = "complete.obs")
inc_Burr_intermediate_lambda <- cor(Burr_intermediate_lambda$lambda,Burr_intermediate_lambda$time,method="spearman",use = "complete.obs")
inc_Burr_intermediate_G <- cor(Burr_intermediate_G$G,Burr_intermediate_G$time,method="spearman",use = "complete.obs")
inc_Burr_intermediate_alpha <- cor(Burr_intermediate_ab$alpha,Burr_intermediate_ab$time,method="spearman",use = "complete.obs")
inc_Burr_intermediate_beta <- cor(Burr_intermediate_ab$beta,Burr_intermediate_ab$time,method="spearman",use = "complete.obs")
inc_Burr_intermediate_T <- cor(Burr_intermediate_T$threshold,Burr_intermediate_T$time,method="spearman",use = "complete.obs")
inc_Erlang_constant_lambda <- cor(Erlang_constant_lambda$lambda,Erlang_constant_lambda$time,method="spearman",use = "complete.obs")
inc_Erlang_constant_G <- cor(Erlang_constant_G$G,Erlang_constant_G$time,method="spearman",use = "complete.obs")
inc_Erlang_constant_alpha <- cor(Erlang_constant_ab$alpha,Erlang_constant_ab$time,method="spearman",use = "complete.obs")
inc_Erlang_constant_beta <- cor(Erlang_constant_ab$beta,Erlang_constant_ab$time,method="spearman",use = "complete.obs")
inc_Erlang_constant_T <- cor(Erlang_constant_T$threshold,Erlang_constant_T$time,method="spearman",use = "complete.obs")
inc_Erlang_pois_lambda <- cor(Erlang_pois_lambda$lambda,Erlang_pois_lambda$time,method="spearman",use = "complete.obs")
inc_Erlang_pois_G <- cor(Erlang_pois_G$G,Erlang_pois_G$time,method="spearman",use = "complete.obs")
inc_Erlang_pois_alpha <- cor(Erlang_pois_ab$alpha,Erlang_pois_ab$time,method="spearman",use = "complete.obs")
inc_Erlang_pois_beta <- cor(Erlang_pois_ab$beta,Erlang_pois_ab$time,method="spearman",use = "complete.obs")
inc_Erlang_pois_T <- cor(Erlang_pois_T$threshold,Erlang_pois_T$time,method="spearman",use = "complete.obs")
inc_Erlang_intermediate_lambda <- cor(Erlang_intermediate_lambda$lambda,Erlang_intermediate_lambda$time,method="spearman",use = "complete.obs")
inc_Erlang_intermediate_G <- cor(Erlang_intermediate_G$G,Erlang_intermediate_G$time,method="spearman",use = "complete.obs")
inc_Erlang_intermediate_alpha <- cor(Erlang_intermediate_ab$alpha,Erlang_intermediate_ab$time,method="spearman",use = "complete.obs")
inc_Erlang_intermediate_beta <- cor(Erlang_intermediate_ab$beta,Erlang_intermediate_ab$time,method="spearman",use = "complete.obs")
inc_Erlang_intermediate_T <- cor(Erlang_intermediate_T$threshold,Erlang_intermediate_T$time,method="spearman",use = "complete.obs")

# Calculate correlations for dose response
DR_Burr_constant_lambda <- cor(as.numeric(names(table(Burr_constant_lambda$lambda))), table(Burr_constant_lambda$lambda) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_constant_G <- cor(as.numeric(names(table(Burr_constant_G$G))), table(Burr_constant_G$G) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_constant_alpha <- cor(as.numeric(names(table(Burr_constant_ab$alpha))), table(Burr_constant_ab$alpha) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_constant_beta <- cor(as.numeric(names(table(Burr_constant_ab$beta))), table(Burr_constant_ab$beta) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_constant_T <- cor(as.numeric(names(table(Burr_constant_T$threshold))), table(Burr_constant_T$threshold) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_pois_lambda <- cor(as.numeric(names(table(Burr_pois_lambda$lambda))), table(Burr_pois_lambda$lambda) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_pois_G <- cor(as.numeric(names(table(Burr_pois_G$G))), table(Burr_pois_G$G) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_pois_alpha <- cor(as.numeric(names(table(Burr_pois_ab$alpha))), table(Burr_pois_ab$alpha) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_pois_beta <- cor(as.numeric(names(table(Burr_pois_ab$beta))), table(Burr_pois_ab$beta) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_pois_T <- cor(as.numeric(names(table(Burr_pois_T$threshold))), table(Burr_pois_T$threshold) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_intermediate_lambda <- cor(as.numeric(names(table(Burr_intermediate_lambda$lambda))), table(Burr_intermediate_lambda$lambda) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_intermediate_G <- cor(as.numeric(names(table(Burr_intermediate_G$G))), table(Burr_intermediate_G$G) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_intermediate_alpha <- cor(as.numeric(names(table(Burr_intermediate_ab$alpha))), table(Burr_intermediate_ab$alpha) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_intermediate_beta <- cor(as.numeric(names(table(Burr_intermediate_ab$beta))), table(Burr_intermediate_ab$beta) / 1000, method = "spearman",use = "complete.obs")
DR_Burr_intermediate_T <- cor(as.numeric(names(table(Burr_intermediate_T$threshold))), table(Burr_intermediate_T$threshold) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_constant_lambda <- cor(as.numeric(names(table(Erlang_constant_lambda$lambda))), table(Erlang_constant_lambda$lambda) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_constant_G <- cor(as.numeric(names(table(Erlang_constant_G$G))), table(Erlang_constant_G$G) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_constant_alpha <- cor(as.numeric(names(table(Erlang_constant_ab$alpha))), table(Erlang_constant_ab$alpha) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_constant_beta <- cor(as.numeric(names(table(Erlang_constant_ab$beta))), table(Erlang_constant_ab$beta) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_constant_T <- cor(as.numeric(names(table(Erlang_constant_T$threshold))), table(Erlang_constant_T$threshold) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_pois_lambda <- cor(as.numeric(names(table(Erlang_pois_lambda$lambda))), table(Erlang_pois_lambda$lambda) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_pois_G <- cor(as.numeric(names(table(Erlang_pois_G$G))), table(Erlang_pois_G$G) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_pois_alpha <- cor(as.numeric(names(table(Erlang_pois_ab$alpha))), table(Erlang_pois_ab$alpha) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_pois_beta <- cor(as.numeric(names(table(Erlang_pois_ab$beta))), table(Erlang_pois_ab$beta) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_pois_T <- cor(as.numeric(names(table(Erlang_pois_T$threshold))), table(Erlang_pois_T$threshold) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_intermediate_lambda <- cor(as.numeric(names(table(Erlang_intermediate_lambda$lambda))), table(Erlang_intermediate_lambda$lambda) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_intermediate_G <- cor(as.numeric(names(table(Erlang_intermediate_G$G))), table(Erlang_intermediate_G$G) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_intermediate_alpha <- cor(as.numeric(names(table(Erlang_intermediate_ab$alpha))), table(Erlang_intermediate_ab$alpha) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_intermediate_beta <- cor(as.numeric(names(table(Erlang_intermediate_ab$beta))), table(Erlang_intermediate_ab$beta) / 1000, method = "spearman",use = "complete.obs")
DR_Erlang_intermediate_T <- cor(as.numeric(names(table(Erlang_intermediate_T$threshold))), table(Erlang_intermediate_T$threshold) / 1000, method = "spearman",use = "complete.obs")

library(ragg)

# Incubation tornado plot - BURR CONSTANT
data <- data.frame(
  Variable = c("lambda", "G", "alpha", "beta", "T"),
  cor1 = c(inc_Burr_constant_lambda,inc_Burr_constant_G,inc_Burr_constant_alpha,inc_Burr_constant_beta,inc_Burr_constant_T),
  cor2 = c(DR_Burr_constant_lambda,DR_Burr_constant_G,DR_Burr_constant_alpha,DR_Burr_constant_beta,DR_Burr_constant_T)
)

# ragg device
ragg::agg_png("SA_Burr_constant_inc.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]

# Horizontal barplot, suppress names
bp <- barplot(
  data1$cor1,
  names.arg = rep("", nrow(data1)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data1$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "Incubation-period SA of Burr-average model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_Burr_constant_dr.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]

# Horizontal barplot, suppress names
bp <- barplot(
  data2$cor2,
  names.arg = rep("", nrow(data2)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data2$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "DR SA of Burr-average model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


# Incubation tornado plot - BURR Poisson
data <- data.frame(
  Variable = c("lambda", "G", "alpha", "beta", "T"),
  cor1 = c(inc_Burr_pois_lambda,inc_Burr_pois_G,inc_Burr_pois_alpha,inc_Burr_pois_beta,inc_Burr_pois_T),
  cor2 = c(DR_Burr_pois_lambda,DR_Burr_pois_G,DR_Burr_pois_alpha,DR_Burr_pois_beta,DR_Burr_pois_T)
)

# ragg device
ragg::agg_png("SA_Burr_pois_inc.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]

# Horizontal barplot, suppress names
bp <- barplot(
  data1$cor1,
  names.arg = rep("", nrow(data1)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data1$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "Incubation-period SA of Burr-Poisson model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_Burr_pois_dr.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]


# Horizontal barplot, suppress names
bp <- barplot(
  data2$cor2,
  names.arg = rep("", nrow(data2)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data2$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "DR SA of Burr-Poisson model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# Incubation tornado plot - BURR intermediate
data <- data.frame(
  Variable = c("lambda", "G", "alpha", "beta", "T"),
  cor1 = c(inc_Burr_intermediate_lambda,inc_Burr_intermediate_G,inc_Burr_intermediate_alpha,inc_Burr_intermediate_beta,inc_Burr_intermediate_T),
  cor2 = c(DR_Burr_intermediate_lambda,DR_Burr_intermediate_G,DR_Burr_intermediate_alpha,DR_Burr_intermediate_beta,DR_Burr_intermediate_T)
)

# ragg device
ragg::agg_png("SA_Burr_intermediate_inc.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]

# Horizontal barplot, suppress names
bp <- barplot(
  data1$cor1,
  names.arg = rep("", nrow(data1)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data1$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "Incubation-period SA of Burr-Intermediate model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_Burr_intermediate_dr.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]


# Horizontal barplot, suppress names
bp <- barplot(
  data2$cor2,
  names.arg = rep("", nrow(data2)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data2$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "DR SA of Burr-Intermediate model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


# Incubation tornado plot - Erlang CONSTANT
data <- data.frame(
  Variable = c("lambda", "G", "alpha", "beta", "T"),
  cor1 = c(inc_Erlang_constant_lambda,inc_Erlang_constant_G,inc_Erlang_constant_alpha,inc_Erlang_constant_beta,inc_Erlang_constant_T),
  cor2 = c(DR_Erlang_constant_lambda,DR_Erlang_constant_G,DR_Erlang_constant_alpha,DR_Erlang_constant_beta,DR_Erlang_constant_T)
)

# ragg device
ragg::agg_png("SA_Erlang_constant_inc.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]

# Horizontal barplot, suppress names
bp <- barplot(
  data1$cor1,
  names.arg = rep("", nrow(data1)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data1$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "Incubation-period SA of Erlang-average model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_Erlang_constant_dr.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]


# Horizontal barplot, suppress names
bp <- barplot(
  data2$cor2,
  names.arg = rep("", nrow(data2)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data2$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "DR SA of Erlang-average model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


# Incubation tornado plot - Erlang Poisson
data <- data.frame(
  Variable = c("lambda", "G", "alpha", "beta", "T"),
  cor1 = c(inc_Erlang_pois_lambda,inc_Erlang_pois_G,inc_Erlang_pois_alpha,inc_Erlang_pois_beta,inc_Erlang_pois_T),
  cor2 = c(DR_Erlang_pois_lambda,DR_Erlang_pois_G,DR_Erlang_pois_alpha,DR_Erlang_pois_beta,DR_Erlang_pois_T)
)

# ragg device
ragg::agg_png("SA_Erlang_pois_inc.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]

# Horizontal barplot, suppress names
bp <- barplot(
  data1$cor1,
  names.arg = rep("", nrow(data1)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data1$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "Incubation-period SA of Erlang-Poisson model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_Erlang_pois_dr.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]


# Horizontal barplot, suppress names
bp <- barplot(
  data2$cor2,
  names.arg = rep("", nrow(data2)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data2$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "DR SA of Erlang-Poisson model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# Incubation tornado plot - Erlang intermediate
data <- data.frame(
  Variable = c("lambda", "G", "alpha", "beta", "T"),
  cor1 = c(inc_Erlang_intermediate_lambda,inc_Erlang_intermediate_G,inc_Erlang_intermediate_alpha,inc_Erlang_intermediate_beta,inc_Erlang_intermediate_T),
  cor2 = c(DR_Erlang_intermediate_lambda,DR_Erlang_intermediate_G,DR_Erlang_intermediate_alpha,DR_Erlang_intermediate_beta,DR_Erlang_intermediate_T)
)

# ragg device
ragg::agg_png("SA_Erlang_intermediate_inc.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]

# Horizontal barplot, suppress names
bp <- barplot(
  data1$cor1,
  names.arg = rep("", nrow(data1)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data1$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "Incubation-period SA of Erlang-Intermediate model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_Erlang_intermediate_dr.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Sort your data
data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]


# Horizontal barplot, suppress names
bp <- barplot(
  data2$cor2,
  names.arg = rep("", nrow(data2)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data2$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "DR SA of Erlang-Intermediate model",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c(
  expression(lambda),
  expression(beta),
  expression(alpha),
  expression(T[L]),
  expression(G)
)

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()











# Now and try and group by parameter


# Incubation tornado plot - G
data <- data.frame(
  Variable = c("EA", "EP", "EI", "BA", "BP", "BI"),
  cor1 = c(inc_Erlang_constant_G, inc_Erlang_pois_G, inc_Erlang_intermediate_G, inc_Burr_constant_G, inc_Burr_pois_G, inc_Burr_intermediate_G),
  cor2 = c(DR_Erlang_constant_G, DR_Erlang_pois_G, DR_Erlang_intermediate_G, DR_Burr_constant_G, DR_Burr_pois_G, DR_Burr_intermediate_G)
)

# ragg device
ragg::agg_png("SA_G_INC.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor1,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "G's effect on the incubation period",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_G_DR.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor2,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = "G's effect on the dose-response",
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


# Incubation tornado plot - T
data <- data.frame(
  Variable = c("EA", "EP", "EI", "BA", "BP", "BI"),
  cor1 = c(inc_Erlang_constant_T, inc_Erlang_pois_T, inc_Erlang_intermediate_T, inc_Burr_constant_T, inc_Burr_pois_T, inc_Burr_intermediate_T),
  cor2 = c(DR_Erlang_constant_T, DR_Erlang_pois_T, DR_Erlang_intermediate_T, DR_Burr_constant_T, DR_Burr_pois_T, DR_Burr_intermediate_T)
)

# ragg device
ragg::agg_png("SA_TL_INC.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor1,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = expression(bold(T[L]*"'s effect on the incubation period")),
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_TL_DR.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor2,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = expression(bold(T[L]*"'s effect on the incubation period")),
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


# Incubation tornado plot - alpha
data <- data.frame(
  Variable = c("EA", "EP", "EI", "BA", "BP", "BI"),
  cor1 = c(inc_Erlang_constant_alpha, inc_Erlang_pois_alpha, inc_Erlang_intermediate_alpha, inc_Burr_constant_alpha, inc_Burr_pois_alpha, inc_Burr_intermediate_alpha),
  cor2 = c(DR_Erlang_constant_alpha, DR_Erlang_pois_alpha, DR_Erlang_intermediate_alpha, DR_Burr_constant_alpha, DR_Burr_pois_alpha, DR_Burr_intermediate_alpha)
)

# ragg device
ragg::agg_png("SA_alpha_INC.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor1,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = expression(bold(alpha*" 's effect on the incubation period")),
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_alpha_DR.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor2,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = expression(bold(alpha*" 's effect on the incubation period")),
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


# Incubation tornado plot - beta
data <- data.frame(
  Variable = c("EA", "EP", "EI", "BA", "BP", "BI"),
  cor1 = c(inc_Erlang_constant_beta, inc_Erlang_pois_beta, inc_Erlang_intermediate_beta, inc_Burr_constant_beta, inc_Burr_pois_beta, inc_Burr_intermediate_beta),
  cor2 = c(DR_Erlang_constant_beta, DR_Erlang_pois_beta, DR_Erlang_intermediate_beta, DR_Burr_constant_beta, DR_Burr_pois_beta, DR_Burr_intermediate_beta)
)

# ragg device
ragg::agg_png("SA_beta_INC.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor1,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = expression(bold(beta*" 's effect on the incubation period")),
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_beta_DR.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor2,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = expression(bold(beta*" 's effect on the incubation period")),
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


# Incubation tornado plot - lambdaa
data <- data.frame(
  Variable = c("EA", "EP", "EI", "BA", "BP", "BI"),
  cor1 = c(inc_Erlang_constant_lambda, inc_Erlang_pois_lambda, inc_Erlang_intermediate_lambda, inc_Burr_constant_lambda, inc_Burr_pois_lambda, inc_Burr_intermediate_lambda),
  cor2 = c(DR_Erlang_constant_lambda, DR_Erlang_pois_lambda, DR_Erlang_intermediate_lambda, DR_Burr_constant_lambda, DR_Burr_pois_lambda, DR_Burr_intermediate_lambda)
)

# ragg device
ragg::agg_png("SA_lambda_INC.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor1,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor1 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = expression(bold(lambda*" 's effect on the incubation period")),
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# ragg device
ragg::agg_png("SA_lambda_DR.png", width = 12, height = 8, units = "in", res = 150)

# Increase margin for Greek labels
par(mar = c(5, 8, 4, 2))  # bottom, left, top, right

# Horizontal barplot, suppress names
bp <- barplot(
  data$cor2,
  names.arg = rep("", nrow(data)),  # no default labels
  horiz = TRUE,
  las = 1,
  col = ifelse(data$cor2 > 0, "blue", "red"),
  xlim = c(-1, 1),
  xlab = "Spearman correlation coefficient",
  cex.lab = 2.5,   # axis labels bigger
  cex.axis = 2.5,  # tick labels bigger
  main = expression(bold(lambda*" 's effect on the incubation period")),
  cex.main = 2.5  # title bigger
)

# Greek/math labels vector
expr_labels <- c("EA", "EP", "EI", "BA", "BP", "BI")

# Place Greek labels, scaled up
for(i in seq_along(bp)) {
  text(x = -1.2, y = bp[i], labels = expr_labels[[i]], xpd = TRUE, adj = 1, cex = 2.5)
}

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()