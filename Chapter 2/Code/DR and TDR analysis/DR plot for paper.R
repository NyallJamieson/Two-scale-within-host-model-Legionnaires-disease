dose_M <- c(1,5,50,100)
response_M <- c(0/4,1/4,8/8,8/8)
muller_fit <- minpack.lm::nlsLM(response_M~1-exp(-a*dose_M),start=list(a=0.1))

doses <- seq(0,100,by=1)
dr_A <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.Constant.csv")[,2:3]
dr_B <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.NegBin.csv")[,2:3]
dr_C <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.Hetero.Exact.csv")[,2:3]

R_A <- append(0,table(dr_A[,1])[1:100]/1000)
R_B <- append(0,table(dr_B[,1])[1:100]/1000)
R_C <- append(0,table(dr_C[,1])[1:100]/1000)

data <- data.frame(dose=c(1,5,50,100),ill=c(0,1,8,8),total=c(4,4,8,8))

dose <- c(1,5,50,100)
success <- c(0,1,8,8)
trial <- c(4,4,8,8)
neg_ll <- function(log_coeff) {
coeff <- exp(log_coeff)  
probs <- 1 - exp(-coeff * dose)
probs <- pmax(pmin(probs, 1 - 1e-10), 1e-10)
log_lik <- success * log(probs) + (trial - success) * log(1 - probs)
return(-sum(log_lik))
}
fit <- optim(par = log(0.1), fn = neg_ll, method = "Brent", lower = log(0.0001), upper = log(10), hessian = TRUE)

log_beta_hat <- fit$par
log_beta_var <- 1 / fit$hessian
log_beta_se <- sqrt(log_beta_var)
beta_hat <- exp(log_beta_hat)
ci_lower <- exp(log_beta_hat - 1.96 * log_beta_se)[,1]  # Lower bound of 95% CI
ci_upper <- exp(log_beta_hat + 1.96 * log_beta_se)[,1]  # Upper bound of 95% CI

point_lower <- c()
point_upper <- c()
for (i in 1:4){
ci <- binom::binom.confint(success[i], trial[i], conf.level = 0.95, methods = "wilson")
point_lower <- append(point_lower,ci$lower)
point_upper <- append(point_upper,ci$upper)
}

# Save plot with 600 DPI resolution
png("Fig3.png", width = 8, height = 6, units = "in", res = 300)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(dose_M,response_M,main="Dose-response for Legionnaires' disease",xlab="Deposited dose",ylab="Probability of symptom onset",lwd=2,ylim=c(0,1),xlim=c(0,100))

lines(doses,R_A,col="red",lwd=2)
lines(doses,R_B,col="green",lwd=2)
lines(doses,R_C,col="blue",lwd=2)
lines(doses,1-exp(-beta_hat*doses),lwd=2,col="black")
lines(doses,1-exp(-ci_upper*doses),col="black",lty=2)
lines(doses,1-exp(-ci_lower*doses),col="black",lty=2)
for (i in 1:4){
arrows(x0=dose_M[i], y0=point_lower[i], x1=dose_M[i], y1=point_upper[i], length=0.05, angle=90, code=3)
}
legend("bottomright",legend=c("Muller data [66]","Exponential DR model fitted to [66]","Exponential DR CI","Within-host model A simulation","Within-host model B simulation","Within-host model C simulation"),col=c("black","black","black","red","green","blue"),lty=c(NA,1,2,1,1,1),pch=c(1,NA,NA,NA,NA,NA),lwd=c(2,2,1,2,2,2),bty="n")

# Close the device
dev.off()