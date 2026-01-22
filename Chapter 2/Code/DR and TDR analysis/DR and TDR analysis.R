source("~/PhD-work/Chapter 1/Code/1. Burr distributions.R")

library("ggplot2")

# Constant
R1 <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.Constant.csv")[,2:3]
dose <- seq(1,500,by=1)
prob <- unname(table(R1[,1]))/1000
plot(dose,prob,xlim=c(0,60),xlab="Dose",ylab="Probability of illness",main="DR curve for model A",yaxt="n")
axis(2, at = c(0,0.5,1), labels = c(0,0.5,1.0))

# dr1 is exponential dr model
dr1 <- minpack.lm::nlsLM(prob~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))
# dr2 is exact beta-poisson dr model
dr2 <- minpack.lm::nlsLM(prob~1-gsl::hyperg_1F1(a,a+b,-dose),start=list(a=1,b=1),control = nls.control(maxiter = 1000))
# dr3 is hill dr model
dr3 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)),start=list(a=1,d50=1),control = nls.control(maxiter = 1000))
# dr4 is burr 1 dr model
dr4 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)*exp(-b*(dose-d50))),start=list(a=1,b=1,d50=12),control = nls.control(maxiter = 1000))
# dr5 is burr 2 dr model
dr5 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)*exp(-b*(dose-d50))*(exp(y*d50)-1)/(exp(y*dose)-1)),start=list(a=0.5,b=0.5,y=0.5,d50=10))

code <- function(dose,a,b,y,d50){1/(1+(d50/dose)^(a)*exp(-b*(dose-d50))*(exp(y*d50)-1)/(exp(y*dose)-1))}
lines(seq(1,70,by=1),code(dose=seq(1,70,by=1),a=coef(dr5)[1],b=coef(dr5)[2],y=coef(dr5)[3],d50=coef(dr5)[4]))

# Calculate AIC of models
aic1 <- AIC(dr1)
aic2 <- AIC(dr2)
aic3 <- AIC(dr3)
aic4 <- AIC(dr4)
aic5 <- AIC(dr5)

# Get TDR data in correct form
d <- seq(1,500,by=1)
t <- seq(0,200,by=0.5)
val <- c()
for (i in t){
print(i)
for (j in d){
val <- append(val,length(R1[which(R1[,1]==j&R1[,2]<=i),2])/1000)}}
data1 <- expand.grid(X=d,Y=t)
data1$Z <- val
data1[length(d)*length(t)+1,] <- c(1000,1000,1)

# Save plot with 600 DPI resolution
png("Fig4.1.png", width = 8, height = 6, units = "in", res = 300)

# TDR heatmap
ggplot(data1, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
  )

# Close device
dev.off()

# Fit the TDR models
tdr1 <- minpack.lm::nlsLM(Z~1-exp(-exp(-k0/Y+k1)*X),start=list(k0=1,k1=1),data=data1,control = nls.control(maxiter = 1000))
tdr2 <- minpack.lm::nlsLM(Z~1-exp(-exp(-k0*Y^(-k2)+k1)*X),start=list(k0=1,k1=1,k2=1),data=data1,control = nls.control(maxiter = 1000))
tdr3 <- minpack.lm::nlsLM(Z~1-exp(-exp(-k0*(tau-Y))*X),start=list(k0=1,tau=1),data=data1,control = nls.control(maxiter = 1000))
tdr4 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(j0/Y+j1))^(-a),start=list(a=1,j0=1,j1=1),data=data1,control = nls.control(maxiter = 1000))
tdr5 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(j0*Y^(-j2)+j1))^(-a),start=list(a=1,j0=1,j1=1,j2=1),data=data1,control = nls.control(maxiter = 1000))
tdr6 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(-j0*Y+j1))^(-a),start=list(a=1,j0=1,j1=1),data=data1,control = nls.control(maxiter = 1000))
tdr7 <- minpack.lm::nlsLM(Z~code(dose=X,a=q1,b=q2,y=q3,d50=q4)*pburrNEW(x=Y,a=q5,b=q6,T=q7*d^(-q8)+q9),start=list(q1=1,q2=1,q3=1,q4=10,q5=1,q6=1,q7=1,q8=1,q9=1),data=data1,control = nls.control(maxiter = 1000))

# Moment calculations
m1 <- c()
m2 <- c()
m3 <- c()
m4 <- c()
for (i in 1:500){
times <- R1[which(R1[,1]==i),2]
m1 <- append(m1,mean(times))
m2 <- append(m2,var(times))
m3 <- append(m3,PerformanceAnalytics::skewness(times))
m4 <- append(m4,PerformanceAnalytics::kurtosis(times,method="moment"))
}
m1 <- m1/24
sd <- (m2^0.5)/24

# Mean at ID50 (9 bacteria)
m <- mean(R1[which(R1[,1]==9),2])
sd <- sd(R1[which(R1[,1]==9),2])
ci <- c(m-1.96*sd,m+1.96*sd)

plot(m3,col="blue",lwd=2,ylim=c(0,6),type="l",xlab="Deposited dose",ylab="Moments for model A",main="First four moments of the dose-dependent incubation period")
lines(log(m4),col="green",lwd=2)
lines(sd,type="l",lwd=2,col="red")
lines(m1,lwd=2)
legend("topright",col=c("black","red","blue","green"),legend=c("mean","standard deviation","skewness","log(kurtosis)"),lwd=2)

# Elbow method

summary_matrix <- as.matrix(data.frame(m1,m2,m3,m4))
# Install and load necessary packages
library(cluster)

# Compute hierarchical clustering
hc <- hclust(dist(summary_matrix), method = "ward.D2")
# Function to compute WCSS
calculate_wcss <- function(k) {
  clusters_k <- cutree(hc, k)
  clusters_list <- split(summary_matrix, clusters_k)
  wcss <- sum(sapply(clusters_list, function(cluster) {
    sum(dist(cluster)^2)
  }))
  return(wcss)
}

# Compute WCSS for a range of k
k_values <- 2:10
wcss_values <- sapply(k_values, calculate_wcss)

# Plot the Elbow Method
plot(k_values, wcss_values, type = "b", pch = 19, xlab = "Number of clusters", ylab = "Within-cluster sum of squares (WCSS)", main = "Elbow method for optimal k")

clusters <- cutree(hc, k = 3)
low <- c()
for (i in which(clusters==1)){
    low <- append(low,R1[which(R1[,1]==i),2])
}
medium <- c()
for (i in which(clusters==2)){
    medium <- append(medium,R1[which(R1[,1]==i),2])
}
high <- c()
for (i in which(clusters==3)){
    high <- append(high,R1[which(R1[,1]==i),2])
}
dens1 <- density(low)
dens2 <- density(medium)
dens3 <- density(high)
plot(dens1,lwd=2,type="l",xlim=c(40,200),ylim=c(0,0.1),xlab="Time in hours",main="Clustered dose-dependent incubation-period distributions")
lines(dens2,lwd=2,col="red")
lines(dens3,lwd=2,col="green")
legend("topright",lwd=2,col=c("black","red","green"),legend=c("low-dose","medium-dose","high-dose"))

plot(clusters,type="l",lwd=2,xlab="Deposited dose",main="Clustered deposited doses",yaxt="n")
# Suppress y-axis and add custom ticks
axis(2, at = c(1.0, 2.0, 3.0), labels = FALSE) # No labels yet
# Add rotated text at each tick position
text(x = par("usr")[1] - 30, y = c(1.0, 2.0, 3.0), 
     labels = c("low", "medium", "high"), 
     srt = 90, xpd = TRUE, adj = 0.5)

# Test of clustered distributions
wilcox.test(low, medium, alternative = "two.sided", paired = FALSE, exact = NULL)
wilcox.test(medium, high, alternative = "two.sided", paired = FALSE, exact = NULL)

# Nb
R2 <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.NegBin.csv")[,2:3]
dose <- seq(1,500,by=1)
prob <- unname(table(R2[,1]))/1000
plot(dose,prob,xlim=c(0,60),xlab="Dose",ylab="Probability of illness",main="DR curve for model B",yaxt="n")
axis(2, at = c(0,0.5,1), labels = c(0,0.5,1.0))

# Fit DR models
dr1 <- minpack.lm::nlsLM(prob~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))
dr2 <- minpack.lm::nlsLM(prob~1-gsl::hyperg_1F1(a,a+b,-dose),start=list(a=1,b=1),control = nls.control(maxiter = 1000))
dr3 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)),start=list(a=1,d50=1),control = nls.control(maxiter = 1000))
dr4 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)*exp(-b*(dose-d50))),start=list(a=1,b=1,d50=12),control = nls.control(maxiter = 1000))
dr5 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)*exp(-b*(dose-d50))*(exp(y*d50)-1)/(exp(y*dose)-1)),start=list(a=0.5,b=0.5,y=0.5,d50=10))

# Calculate AIC for models
aic1 <- AIC(dr1)
aic2 <- AIC(dr2)
aic3 <- AIC(dr3)
aic4 <- AIC(dr4)
aic5 <- AIC(dr5)

lines(seq(1,70,by=1),1-exp(-log(2)*seq(1,70,by=1)/coef(dr1)),lwd=2)
lines(seq(1,70,by=1),1-gsl::hyperg_1F1(coef(dr2)[1],coef(dr2)[1]+coef(dr2)[2],-seq(1,70,by=1)),col="red",lwd=2)
lines(seq(1,70,by=1),1/(1+(coef(dr3)[2]/seq(1,70,by=1))^(coef(dr3)[1])),col="blue")
lines(seq(1,70,by=1),1/(1+(coef(dr4)[3]/seq(1,70,by=1))^(coef(dr4)[1])*exp(-coef(dr4)[2]*(seq(1,70,by=1)-coef(dr4)[3]))),col="green",lwd=2)

# Get TDR data in form we want
d <- seq(1,500,by=1)
t <- seq(0,200,by=0.5)
val <- c()
for (i in t){
print(i)
for (j in d){
val <- append(val,length(R2[which(R2[,1]==j&R2[,2]<=i),2])/1000)}}
data2 <- expand.grid(X=d,Y=t)
data2$Z <- val
data2[length(d)*length(t)+1,] <- c(1000,1000,1)

# Save plot with 600 DPI resolution
png("Fig4.2.png", width = 8, height = 6, units = "in", res = 300)

# TDR heatmap
ggplot(data2, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
  )

# Close device
dev.off()

# Fit TDR models
tdr1 <- minpack.lm::nlsLM(Z~1-exp(-exp(-k0/Y+k1)*X),start=list(k0=1,k1=1),data=data2,control = nls.control(maxiter = 1000))
tdr2 <- minpack.lm::nlsLM(Z~1-exp(-exp(-k0*Y^(-k2)+k1)*X),start=list(k0=1,k1=0.1,k2=1),data=data2,control = nls.control(maxiter = 1000))
tdr3 <- minpack.lm::nlsLM(Z~1-exp(-exp(-k0*(tau-Y))*X),start=list(k0=1,tau=1),data=data2,control = nls.control(maxiter = 1000))
tdr4 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(j0/Y+j1))^(-a),start=list(a=1,j0=1,j1=1),data=data2,control = nls.control(maxiter = 1000))
tdr5 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(j0*Y^(-j2)+j1))^(-a),start=list(a=1,j0=1,j1=1,j2=1),data=data2,control = nls.control(maxiter = 1000))
tdr6 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(-j0*Y+j1))^(-a),start=list(a=1,j0=1,j1=1),data=data2,control = nls.control(maxiter = 1000))
tdr7 <- minpack.lm::nlsLM(Z~code(dose=X,a=q1,b=q2,y=q3,d50=q4)*pburrNEW(x=Y,a=q5,b=q6,T=q7*d^(-q8)+q9),start=list(q1=1,q2=1,q3=1,q4=10,q5=1,q6=1,q7=1,q8=1,q9=1),data=data2,control = nls.control(maxiter = 1000))

# Moment calculations
m1 <- c()
m2 <- c()
m3 <- c()
m4 <- c()
for (i in 1:500){
times <- R2[which(R2[,1]==i),2]
m1 <- append(m1,mean(times))
m2 <- append(m2,var(times))
m3 <- append(m3,PerformanceAnalytics::skewness(times))
m4 <- append(m4,PerformanceAnalytics::kurtosis(times,method="moment"))
}
m1 <- m1/24
sd <- (m2^0.5)/24

# Mean at ID50 (9 bacteria)
m <- mean(R2[which(R2[,1]==9),2])
sd <- sd(R2[which(R2[,1]==9),2])
ci <- c(m-1.96*sd,m+1.96*sd)/24

plot(m3,col="blue",lwd=2,ylim=c(0,6),type="l",xlab="Deposited dose",ylab="Moments for model B",main="First four moments of the dose-dependent incubation period")
lines(log(m4),col="green",lwd=2)
lines(sd,type="l",lwd=2,col="red")
lines(m1,lwd=2)
legend("topright",col=c("black","red","blue","green"),legend=c("mean","standard deviation","skewness","log(kurtosis)"),lwd=2)


# Elbow method

summary_matrix <- as.matrix(data.frame(m1,m2,m3,m4))
# Install and load necessary packages

# Compute hierarchical clustering
hc <- hclust(dist(summary_matrix), method = "ward.D2")
# Function to compute WCSS
calculate_wcss <- function(k) {
  clusters_k <- cutree(hc, k)
  clusters_list <- split(summary_matrix, clusters_k)
  wcss <- sum(sapply(clusters_list, function(cluster) {
    sum(dist(cluster)^2)
  }))
  return(wcss)
}

# Compute WCSS for a range of k
k_values <- 2:10
wcss_values <- sapply(k_values, calculate_wcss)

# Plot the Elbow Method
plot(k_values, wcss_values, type = "b", pch = 19, xlab = "Number of clusters", ylab = "Within-cluster sum of squares (WCSS)", main = "Elbow method for optimal k")

clusters <- cutree(hc, k = 4)
low <- c()
for (i in which(clusters==1)){
low <- append(low,R2[which(R2[,1]==i),2])
}
medium_low <- c()
for (i in which(clusters==2)){
medium_low <- append(medium_low,R2[which(R2[,1]==i),2])
}
medium_high <- c()
for (i in which(clusters==3)){
medium_high <- append(medium_high,R2[which(R2[,1]==i),2])
}
high <- c()
for (i in which(clusters==4)){
high <- append(high,R2[which(R2[,1]==i),2])
}
dens1 <- density(low)
dens2 <- density(medium_low)
dens3 <- density(medium_high)
dens4 <- density(high)
plot(dens1,lwd=2,type="l",xlim=c(40,200),ylim=c(0,0.1),xlab="Time in hours",main="Clustered dose-dependent incubation-period distributions")
lines(dens2,lwd=2,col="red")
lines(dens3,lwd=2,col="blue")
lines(dens4,lwd=2,col="green")
legend("topright",lwd=2,col=c("black","red","blue","green"),legend=c("low-dose","low-medium-dose","high-medium-dose","high-dose"))

plot(clusters,type="l",lwd=2,xlab="Deposited dose",main="Clustered deposited doses",yaxt="n")
# Suppress y-axis and add custom ticks
axis(2, at = c(1.0, 2.0, 3.0, 4.0), labels = FALSE) # No labels yet
# Add rotated text at each tick position
text(x = par("usr")[1] - 30, y = c(1.0, 2.0, 3.0, 4.0), 
     labels = c("low", "low\nmedium", "high\nmedium", "high"), 
     srt = 90, xpd = TRUE, adj = 0.5)


# Test of clustered distributions
wilcox.test(low, medium_low, alternative = "two.sided", paired = FALSE, exact = NULL)
wilcox.test(medium_low, medium_high, alternative = "two.sided", paired = FALSE, exact = NULL)
wilcox.test(medium, high, alternative = "two.sided", paired = FALSE, exact = NULL)


# Hetero exact
R3 <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.Hetero.Exact.1to500.csv")[,2:3]
dose <- seq(1,500,by=1)
prob <- unname(table(R3[,1]))/1000
plot(dose,prob,xlim=c(0,60),xlab="Dose",ylab="Probability of illness",main="DR curve for model C",yaxt="n")
axis(2, at = c(0,0.5,1), labels = c(0,0.5,1.0))

# Fit DR models
dr1 <- minpack.lm::nlsLM(prob~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))
dr2 <- minpack.lm::nlsLM(prob~1-gsl::hyperg_1F1(a,a+b,-dose),start=list(a=1,b=1),control = nls.control(maxiter = 1000))
dr3 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)),start=list(a=1,d50=1),control = nls.control(maxiter = 1000))
dr4 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)*exp(-b*(dose-d50))),start=list(a=1,b=1,d50=12),control = nls.control(maxiter = 1000))
dr5 <- minpack.lm::nlsLM(prob~1/(1+(d50/dose)^(a)*exp(-b*(dose-d50))*(exp(y*d50)-1)/(exp(y*dose)-1)),start=list(a=0.5,b=0.5,y=0.5,d50=10))

# Calculate AIC for fitted models
aic1 <- AIC(dr1)
aic2 <- AIC(dr2)
aic3 <- AIC(dr3)
aic4 <- AIC(dr4)
aic5 <- AIC(dr5)

# Get TDR data in form we want
d <- seq(1,500,by=1)
t <- seq(0,200,by=0.5)
val <- c()
for (i in t){
print(i)
for (j in d){
val <- append(val,length(R3[which(R3[,1]==j&R3[,2]<=i),2])/1000)}}
data3 <- expand.grid(X=d,Y=t)
data3$Z <- val
data3[length(d)*length(t)+1,] <- c(1000,1000,1)

# Save plot with 600 DPI resolution
png("Fig4.3.png", width = 8, height = 6, units = "in", res = 300)

# TDR heatmap
ggplot(data3, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
  )

# Close the dev
dev.off()

tdr1 <- minpack.lm::nlsLM(Z~1-exp(-exp(k0/Y+k1)*X),start=list(k0=1,k1=1),data=data3,control = nls.control(maxiter = 1000))
tdr2 <- minpack.lm::nlsLM(Z~1-exp(-exp(k0*Y^(-k2)+k1)*X),start=list(k0=1,k1=1,k2=1),data=data3,control = nls.control(maxiter = 1000))
tdr3 <- minpack.lm::nlsLM(Z~1-exp(-exp(k0*Y+k1)*X),start=list(k0=1,k1=1),data=data3,control = nls.control(maxiter = 1000))
tdr4 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(j0/Y+j1))^(-a),start=list(a=1,j0=1,j1=1),data=data3,control = nls.control(maxiter = 1000))
tdr5 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(j0*Y^(-j2)+j1))^(-a),start=list(a=1,j0=1,j1=1,j2=1),data=data3,control = nls.control(maxiter = 1000))
tdr6 <- minpack.lm::nlsLM(Z~1-(1+X*(2^(1/a)-1)/exp(j0*Y+j1))^(-a),start=list(a=1,j0=1,j1=1),data=data3,control = nls.control(maxiter = 1000))
tdr7 <- minpack.lm::nlsLM(Z~code(dose=X,a=q1,b=q2,y=q3,d50=q4)*pburrNEW(x=Y,a=q5,b=q6,T=q7*d^(-q8)+q9),start=list(q1=1,q2=1,q3=1,q4=10,q5=1,q6=1,q7=1,q8=1,q9=1),data=data3,control = nls.control(maxiter = 1000))

# moment calculations
m1 <- c()
m2 <- c()
m3 <- c()
m4 <- c()
for (i in 1:500){
times <- R3[which(R3[,1]==i),2]
m1 <- append(m1,mean(times))
m2 <- append(m2,var(times))
m3 <- append(m3,PerformanceAnalytics::skewness(times))
m4 <- append(m4,PerformanceAnalytics::kurtosis(times,method="moment"))
}
m1 <- m1/24
sd <- (m2^0.5)/24

plot(log(m3),col="blue",lwd=2,ylim=c(0,8),type="l",xlab="Deposited dose",ylab="Moments for model C",main="First four moments of the dose-dependent incubation period")
lines(log(m4),col="green",lwd=2)
lines(log(sd),type="l",lwd=2,col="red")
lines(m1,lwd=2)
legend("topright",col=c("black","red","blue","green"),legend=c("mean","log(standard deviation)","skewness","log(kurtosis)"),lwd=2)

# Elbow method
summary_matrix <- as.matrix(data.frame(m1,m2,m3,m4))
# Install and load necessary packages
library(cluster)

# Compute hierarchical clustering
hc <- hclust(dist(summary_matrix), method = "ward.D2")
# Function to compute WCSS
calculate_wcss <- function(k) {
  clusters_k <- cutree(hc, k)
  clusters_list <- split(summary_matrix, clusters_k)
  wcss <- sum(sapply(clusters_list, function(cluster) {
    sum(dist(cluster)^2)
  }))
  return(wcss)
}

# Compute WCSS for a range of k
k_values <- 2:10
wcss_values <- sapply(k_values, calculate_wcss)

# Plot the Elbow Method
plot(k_values, wcss_values, type = "b", pch = 19, xlab = "Number of clusters", ylab = "Within-cluster sum of squares (WCSS)", main = "Elbow method for optimal k")

clusters <- cutree(hc, k = 4)
low <- c()
for (i in which(clusters==1)){
low <- append(low,R3[which(R3[,1]==i),2])
}
medium_low <- c()
for (i in which(clusters==2)){
medium_low <- append(medium_low,R3[which(R3[,1]==i),2])
}
medium_high <- c()
for (i in which(clusters==3)){
medium_high <- append(medium_high,R3[which(R3[,1]==i),2])
}
high <- c()
for (i in which(clusters==4)){
high <- append(high,R3[which(R3[,1]==i),2])
}
dens1 <- density(low)
dens2 <- density(medium_low)
dens3 <- density(medium_high)
dens4 <- density(high)
plot(dens1,lwd=2,type="l",xlim=c(0,200),ylim=c(0,0.1),xlab="Time in hours",main="Clustered dose-dependent incubation-period distributions")
lines(dens2,lwd=2,col="red")
lines(dens3,lwd=2,col="blue")
lines(dens4,lwd=2,col="green")
legend("topright",lwd=2,col=c("black","red","blue","green"),legend=c("low-dose","low-medium-dose","high-medium-dose","high-dose"))

plot(clusters,type="l",lwd=2,xlab="Deposited dose",main="Clustered deposited doses",yaxt="n")
# Suppress y-axis and add custom ticks
axis(2, at = c(1.0, 2.0, 3.0, 4.0), labels = FALSE) # No labels yet
# Add rotated text at each tick position
text(x = par("usr")[1] - 30, y = c(1.0, 2.0, 3.0, 4.0), 
     labels = c("low", "low\nmedium", "high\nmedium", "high"), 
     srt = 90, xpd = TRUE, adj = 0.5)


# test of clustered distributions
wilcox.test(low, medium_low, alternative = "two.sided", paired = FALSE, exact = NULL)
wilcox.test(medium_low, medium_high, alternative = "two.sided", paired = FALSE, exact = NULL)
wilcox.test(medium, high, alternative = "two.sided", paired = FALSE, exact = NULL)

# Now calculate the mean (and standard deviation in the mean) for each of the three datasets
m1 <- mean(R1[which(R1[,1]==9),2])
m2 <- mean(R2[which(R2[,1]==9),2])
m3 <- mean(R3[which(R3[,1]==9),2])
sd1 <- sd(R1[which(R1[,1]==9),2])
sd2 <- sd(R2[which(R2[,1]==9),2])
sd3 <- sd(R3[which(R3[,1]==9),2])

# Confidence level 
conf_level <- 0.95
alpha <- 1 - conf_level

# t critical value
t_crit_1 <- qt(1 - alpha/2, df = length(which(R1[,1]==9)) - 1)
t_crit_2 <- qt(1 - alpha/2, df = length(which(R2[,1]==9)) - 1)
t_crit_3 <- qt(1 - alpha/2, df = length(which(R3[,1]==9)) - 1)

# Standard error
se_1 <- sd1 / sqrt(length(which(R1[,1]==9)))
se_2 <- sd2 / sqrt(length(which(R2[,1]==9)))
se_3 <- sd3 / sqrt(length(which(R3[,1]==9)))

# Confidence intervals
ci_1 <- c(m1 - t_crit_1 * se_1, m1 + t_crit_1 * se_1)
ci_2 <- c(m2 - t_crit_2 * se_2, m2 + t_crit_2 * se_2)
ci_3 <- c(m3 - t_crit_3 * se_3, m3 + t_crit_3 * se_3)
