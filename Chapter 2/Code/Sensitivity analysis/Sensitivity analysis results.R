# Load in sensitivity results
sens_lambda <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_lambda.csv")[,2:4]
sens_A <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_A.csv")[,2:4]
sens_w <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_w.csv")[,2:4]
sens_G <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_G.csv")[,2:4]
sens_alpha <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_alpha.csv")[,2:4]
sens_beta <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_beta.csv")[,2:4]
sens_eta <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_eta.csv")[,2:4]
sens_phi <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_phi.csv")[,2:4]
sens_T <- read.csv("~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_T.csv")[,2:4]

# Calculate correlations for incubation period
inc_lambda <- cor(sens_lambda$lambda,sens_lambda$mean_inc,method="spearman")
inc_A <- cor(sens_A$A,sens_A$mean_inc,method="spearman")
inc_w <- cor(sens_w$w,sens_w$mean_inc,method="spearman")
inc_G <- cor(sens_G$G,sens_G$mean_inc,method="spearman")
inc_alpha <- cor(sens_alpha$alpha,sens_alpha$mean_inc,method="spearman")
inc_beta <- cor(sens_beta$beta,sens_beta$mean_inc,method="spearman")
inc_eta <- cor(sens_eta$eta,sens_eta$mean_inc,method="spearman")
inc_phi <- cor(sens_phi$phi,sens_phi$mean_inc,method="spearman")
inc_T <- cor(sens_T$threshold,sens_T$mean_inc,method="spearman")

# Calculate correlations for dose response
prob_lambda <- cor(sens_lambda$lambda,sens_lambda$occurrences,method="spearman")
prob_A <- cor(sens_A$A,sens_A$occurrences,method="spearman")
prob_w <- cor(sens_w$w,sens_w$occurrences,method="spearman")
prob_G <- cor(sens_G$G,sens_G$occurrences,method="spearman")
prob_alpha <- cor(sens_alpha$alpha,sens_alpha$occurrences,method="spearman")
prob_beta <- cor(sens_beta$beta,sens_beta$occurrences,method="spearman")
prob_eta <- cor(sens_eta$eta,sens_eta$occurrences,method="spearman")
prob_phi <- cor(sens_phi$phi,sens_phi$occurrences,method="spearman")
prob_T <- cor(sens_T$threshold,sens_T$occurrences,method="spearman")

# Incubation tornado plot
data <- data.frame(
    Variable = c("lambda", "A", "w", "G", "alpha", "beta", "eta", "phi", "T"),
    cor1 = c(inc_lambda,inc_A,inc_w,inc_G,inc_alpha,inc_beta,inc_eta,inc_phi,inc_T),
    cor2 = c(prob_lambda,prob_A,prob_w,prob_G,prob_alpha,prob_beta,prob_eta,prob_phi,prob_T)
)

# Order data by the absolute values of the Result column
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]
data1$Variable <- expression(psi,lambda,omega,"G",phi,beta,alpha,"T"["L"],eta)
barplot(
    data1$cor1, 
    names.arg = data1$Variable, 
    horiz = TRUE, 
    las = 1, # Rotate axis labels
    col = ifelse(data1$cor1 > 0, "blue", "red"), # Positive in blue, negative in red,
    xlim = c(-1, 1))

data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]
data2$Variable <- expression(phi,"T"["L"],eta,alpha,beta,psi,"G",lambda,omega)
barplot(
    data2$cor2, 
    names.arg = data2$Variable, 
    horiz = TRUE, 
    las = 1, # Rotate axis labels
    col = ifelse(data2$cor2 > 0, "blue", "red"), # Positive in blue, negative in red,
    xlim = c(-1, 1))
