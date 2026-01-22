# Read in rupture data
l <- read.csv("~/PhD-work/Chapter 2/Data/Rupture data/lambda_data.csv",header=FALSE)
time <- l[,2]
prop <- l[,1]
pred_frame <- seq(1,72,by=1)
pred_frame <- data.frame(pred_frame)
death_prop <- prop/100
death_cdf <- death_prop/death_prop[72]
pburrNEW <- function(x,a,b,T){(1+(T/x)^(a)*exp((T-x)/b))^(-1)}

# Lambda
####################################################
# ci for lambda is obtained by method called non-parametric bootsrapping
part1 <- function(y){
u <- runif(n=1000)
samp <- c()
for (i in u){
val <- uniroot(function(x) pburrNEW(x,a=0.6327,b=11.0200,T=21.7964)-i,lower=0,upper=500)$root
samp <- append(samp,val)
}
return(unname(1/quantile(samp,probs=0.5/death_prop[72])))
}
lambda_dist <- as.numeric(parallel::mclapply(1:10000,part1,mc.cores=5))
####################################################

# Include T to finalise
####################################################
time <- c(2,24,48,72)
dose <- c(13584.15,588828.58,11207102.05,82907236.06)
L0 <- 10^6
yourdata2 <- data.frame(time,dose)
# Now stick in other parameters (also calculate theta, x1, C, T)

results <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(results) <- c("G")

part2 <- function(y){
h <- sample(size=1,lambda_dist,replace=TRUE)
A_w <- MASS::mvrnorm(n=1,mu=c(177.95295,0.19242),Sigma=matrix(c(82.49019097,-0.0254690939,-0.0254690939,0.0001933546),2,2))
A <- A_w[1]
w <- A_w[2]
G <- A/(1+(A-1)*exp(-w*(1/h-1)))

return(data.frame(G))
}

results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:10000,part2,mc.cores=5)),na.action="omit")))
fit_p <- fitdistrplus::fitdist(round(results[,1]),"pois",method="mle")
fit_nb <- fitdistrplus::fitdist(round(results[,1]),"nbinom",method="mle")

# Save plot with 600 DPI resolution
tiff("Fig2c.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

hist(round(results[,1]),breaks=50,freq=FALSE,xlab="Estimate of G",main="Distribution of G",ylim=c(0,0.05),xlim=c(20,100))
lines(seq(20,120,by=1),dpois(seq(20,120,by=1),lambda=62.862),lwd=2,col="blue")
lines(seq(20,120,by=1),dnbinom(seq(20,120,by=1),mu=62.86806,size=32.00728),lwd=2,col="red")

# Close the device
dev.off()

# Save plot with 600 DPI resolution
png("Fig2c.png", width = 8, height = 8, units = "in", res = 300)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

hist(round(results[,1]),breaks=50,freq=FALSE,xlab="Estimate of G",main="Distribution of G",ylim=c(0,0.05),xlim=c(20,100))
lines(seq(20,120,by=1),dpois(seq(20,120,by=1),lambda=62.862),lwd=2,col="blue")
lines(seq(20,120,by=1),dnbinom(seq(20,120,by=1),mu=62.86806,size=32.00728),lwd=2,col="red")

# Close the device
dev.off()
