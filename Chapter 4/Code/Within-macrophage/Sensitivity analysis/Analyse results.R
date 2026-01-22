# Load in sensitivity results
sens_omega_psi <- read.csv("~/PhD-work/Chapter 4/Data/Within-macrophage/Sensitivity analysis/WithinMacro_SA_omega_psi.csv")
sens_lambda <- read.csv("~/PhD-work/Chapter 4/Data/Within-macrophage/Sensitivity analysis/WithinMacro_SA_lambda.csv")

# Calculate correlations
cor_psi <- cor(sens_omega_psi$psi,sens_omega_psi$G,method="spearman")
cor_omega <- cor(sens_omega_psi$omega,sens_omega_psi$G,method="spearman")
cor_lambda <- cor(sens_lambda$lambda,sens_lambda$G,method="spearman")

# Incubation tornado plot
data <- data.frame(
  Variable = c("psi", "omega", "lambda"),
  cor1 = c(cor_psi, cor_omega, cor_lambda)
)

# Order data by the absolute values of the Result column
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]
data1$Variable <- expression(psi, omega, lambda)
barplot(
  data1$cor1, 
  names.arg = data1$Variable, 
  horiz = TRUE, 
  las = 1, # Rotate axis labels
  col = ifelse(data1$cor1 > 0, "blue", "red"), # Positive in blue, negative in red,
  xlim = c(-1, 1))
