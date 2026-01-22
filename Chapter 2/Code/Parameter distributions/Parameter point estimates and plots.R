# A and omega
t <- seq(1,72,by=1)
dep <- c(1,59.41733,161.69467,189.78458)
ind <- c(1,24,48,72)
yourdata <- data.frame(dep,ind)
a <- nls(dep~C/(1+(C-1)*exp(-w*(ind-1))),data=yourdata,start=list(C=100,w=0.2))
g <- function(t,pars){(t>=0 & t<=1)*1+(t>1)*pars[1]/(1+(pars[1]-1)*exp(-pars[2]*(t-1)))}

# Save plot with 600 DPI resolution
tiff("Fig2b.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(ind,dep,xlab="Time in hours",ylab="Intracellular Legionella per macrophage",main="Intracellular growth of Legionella",ylim=c(0,200))
lines(t,g(t,c(177.9530,0.1924)),col="red",lwd=2)
lines(t,g(t,c(177.9530+1.96*9.08241,0.19242+1.96*0.01391)),lwd=2)
lines(t,g(t,c(177.9530-1.96*9.08241,0.19242-1.96*0.01391)),lwd=2)

# Close the device
dev.off()

# Save plot with 600 DPI resolution
png("Fig2b.png", width = 8, height = 8, unit = "in", res = 300)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(ind,dep,xlab="Time in hours",ylab="Intracellular Legionella per macrophage",main="Intracellular growth of Legionella",ylim=c(0,200))
lines(t,g(t,c(177.9530,0.1924)),col="red",lwd=2)
lines(t,g(t,c(177.9530+1.96*9.08241,0.19242+1.96*0.01391)),lwd=2)
lines(t,g(t,c(177.9530-1.96*9.08241,0.19242-1.96*0.01391)),lwd=2)

# Close the device
dev.off()

A <- coef(a)[1]
w <- coef(a)[2]

# lambda and G
l <- read.csv("~/PhD-work/Chapter 2/Data/Rupture data/lambda_data.csv",header=FALSE)
t2 <- l[,2]
prop <- l[,1]
pred_frame <- seq(1,72,by=1)
pred_frame <- data.frame(pred_frame)
death_prop <- prop/100
death_cdf <- death_prop/death_prop[72]
pburrNEW <- function(x,a,b,T){(1+(T/x)^(a)*exp((T-x)/b))^(-1)}
model_fit <- nls(death_cdf~pburrNEW(t2,a=q1,b=q2,T=q3),start=list(q1=1,q2=10,q3=20))
lambda <- uniroot(function(t) pburrNEW(t,a=coef(model_fit)[1],b=coef(model_fit)[2],T=coef(model_fit)[3])-0.5/death_prop[72],interval=c(0,100))$root^(-1)
model_fit_exp <- nls(death_cdf~pexp(t2,rate=q1),start=list(q1=0.05))
model_fit_gamma <- nls(death_cdf~pgamma(t2,shape=q1,rate=q2),start=list(q1=0.1,q2=0.005))
G <- g(1/lambda,c(A,w))

# Save plot with 600 DPI resolution
tiff("Fig2a.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(t2,death_prop,xlab="Time in hours",ylab="Percentage of macrophages dead",main="Infected macrophages rupture")
lines(seq(0,72,by=1),death_prop[72]*pexp(seq(0,72,by=1),rate=coef(model_fit_exp)),lwd=2)
lines(seq(0,72,by=1),death_prop[72]*pgamma(seq(0,72,by=1),shape=coef(model_fit_gamma)[1],rate=coef(model_fit_gamma)[2]),lwd=2,col="green")
lines(seq(0,72,by=1),death_prop[72]*pburrNEW(seq(0,72,by=1),a=coef(model_fit)[1],b=coef(model_fit)[2],T=coef(model_fit)[3]),lwd=2,col="red")
legend("bottomright",lwd=2,legend=c("Exponential","Gamma","Burr"),col=c("black","green","red"), cex = 2)

# Close the device
dev.off()

# Save plot with 600 DPI resolution
png("FigHI.png", width = 8, height = 8, units = "in", res = 300)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(t2,death_prop,xlab="Time in hours",ylab="Percentage of macrophages dead",main="Infected macrophages rupture")
lines(seq(0,72,by=1),death_prop[72]*pexp(seq(0,72,by=1),rate=coef(model_fit_exp)),lwd=2,lty=1)
lines(seq(0,72,by=1),death_prop[72]*pgamma(seq(0,72,by=1),shape=coef(model_fit_gamma)[1],rate=coef(model_fit_gamma)[2]),lwd=2,col="green",lty=2)
lines(seq(0,72,by=1),death_prop[72]*pburrNEW(seq(0,72,by=1),a=coef(model_fit)[1],b=coef(model_fit)[2],T=coef(model_fit)[3]),lwd=2,col="red",lty=3)
legend("bottomright",lwd=2,legend=c("Exponential","Gamma","Burr"),col=c("black","green","red"), cex = 2,lty=c(1,2,3))

# Close the device
dev.off()

model_fit_erlang <- model_fit_gamma <- nls(death_cdf~pgamma(t2,shape=3,rate=q2),start=list(q2=0.05))
  
# Save plot with 600 DPI resolution
png("FigHI.png", width = 8, height = 8, units = "in", res = 300)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(t2,death_prop,xlab="Time (hours)",ylab="Proportion of macrophages dead",main="Infected macrophage rupture")
lines(seq(0,72,by=1),death_prop[72]*pexp(seq(0,72,by=1),rate=coef(model_fit_exp)),lwd=2,lty=1)
lines(seq(0,72,by=1),death_prop[72]*pgamma(seq(0,72,by=1),shape=3,rate=coef(model_fit_erlang)[1]),lwd=2,col="green",lty=2)
lines(seq(0,72,by=1),death_prop[72]*pburrNEW(seq(0,72,by=1),a=coef(model_fit)[1],b=coef(model_fit)[2],T=coef(model_fit)[3]),lwd=2,col="red",lty=2)
legend("bottomright",lwd=2,legend=c("Exponential","Erlang","Burr"),col=c("black","green","red"), cex = 2, lty = c(1, 2, 2))

# Close the device
dev.off()

# alpha and beta
time <- c(2,24,48,72)
logdose <- log(c(13584.15,588828.58,11207102.05,82907236.06))
L0 <- 10^5
yourdata2 <- data.frame(time,logdose)

ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
L <- function(time,gamma=gamma,Pi=Pi){log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)))}
alpha <- unname(coef(ab_model)[1]*coef(ab_model)[2])
beta <- unname(coef(ab_model)[1]*(1-coef(ab_model)[2]))

reps <- 10000
ab <- data.frame(MASS::mvrnorm(n=reps,mu=c(coef(ab_model)[1],coef(ab_model)[2]),Sigma=matrix(c(vcov(ab_model)[1],vcov(ab_model)[2],vcov(ab_model)[3],vcov(ab_model)[4]),2,2)))
data <- data.frame(matrix(NA,nrow=72,ncol=0))

ab2 <- ab |> dplyr::filter(ab[,1]>0 & ab[,2]>0)
print(nrow(ab2))
for (i in 1:nrow(ab2)){
print(i)
data <- data.frame(data,unname(L(t=seq(1,72,by=1),gamma=ab2[i,1],Pi=ab2[i,2])))
}
lwr <- c()
upr <- c()
for (i in 1:72){
print(i)
lwr <- append(lwr,quantile(na.omit(as.numeric(data[i,])),p=0.025))
upr <- append(upr,quantile(na.omit(as.numeric(data[i,])),p=0.975))
}

# Save plot with 600 DPI resolution
tiff("Fig2d.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(time,logdose,ylim=c(7,20),xlab="Time in hours",ylab="Log extracellular Legionella population",main="Log extracellular Legionella population")
lines(L(seq(1,72,by=1),gamma=coef(ab_model)[1],Pi=coef(ab_model)[2]),col="red",lwd=2)
lines(lwr,col="black",lwd=2)
lines(upr,col="black",lwd=2)

# Close the device
dev.off()

# Save plot with 600 DPI resolution
png("Fig2d.png", width = 8, height = 8, units = "in", res = 300)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(time,logdose,ylim=c(7,20),xlab="Time in hours",ylab="Log extracellular Legionella population",main="Log extracellular Legionella population")
lines(L(seq(1,72,by=1),gamma=coef(ab_model)[1],Pi=coef(ab_model)[2]),col="red",lwd=2)
lines(lwr,col="black",lwd=2)
lines(upr,col="black",lwd=2)

# Close the device
dev.off()

# Model based on only the first data
yourdata2_short <- yourdata2[1:3,]
ab_model2 <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0),data=yourdata2_short)
lines(L(seq(1,72,by=1),gamma=coef(ab_model2)[1],Pi=coef(ab_model2)[2]),col="blue",lwd=2)

# C
d3 <- c(4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,5*10^4,5*10^4,5*10^4,5*10^4,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2)
t3 <- 24*c(1,1,1,1,1,1,1,1,1,2,2,2,1,2,2,2,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4)
t_L <- t3-24
t_R <- t3
theta <- sqrt((lambda-alpha-beta)^2+4*lambda*alpha*G)
x1 <- 0.5*(-lambda-alpha-beta+theta)
rhs <- log(d3)/x1
Data <- data.frame(t_L,t_R,rhs)
vals <- log(d3)/x1

# Create an interval-censored survival object
surv_obj <- survival::Surv(time = t_L+vals, time2 = t_R+vals, type = "interval2")
model <- survival::survreg(surv_obj ~ 1, data = Data, dist = "gaussian")
C <- unname(coef(model))
se <- unname(predict(model,se.fit=TRUE)$se.fit[1])
sigma <- model$scale
C_lower <- C-1.96*sqrt(se^2+sigma^2)
C_upper <- C+1.96*sqrt(se^2+sigma^2)

x3 <- seq(5,13,by=1)

# Save plot with 600 DPI resolution
tiff("Fig2f.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(log(d3),t3,xlab="Log dose",ylab="Time to illness",main="Time for illness in guinea pigs",pch=19,col="black",cex=2.4)
lines(x3,-x3/x1+C,col="red",lwd=2)
lines(x3,-x3/x1+C_lower,col="black",lwd=2)
lines(x3,-x3/x1+C_upper,col="black",lwd=2)
text(8.29404964010203+0.01,24+0.1,"2",col="yellow")
text(10.8197782844103+0.01,24+0.1,"1",col="yellow")
text(12.8992198260901+0.03,24+0.1,"9",col="yellow")
text(8.29404964010203+0.01,48+0.1,"7",col="yellow")
text(10.8197782844103+0.01,48+0.1,"3",col="yellow")
text(12.8992198260901+0.03,48+0.1,"3",col="yellow")
text(5.29831736654804+0.02,72+0.1,"5",col="yellow")
text(5.29831736654804+0.02,96+0.1,"5",col="yellow")
text(8.29404964010203+0.01,72+0.1,"1",col="yellow")

# Close the device
dev.off()

# Save plot with 600 DPI resolution
png("Fig2f.png", width = 8, height = 8, units = "in", res = 300)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(log(d3),t3,xlab="Log dose",ylab="Time to illness (hours)",main="Time for illness in guinea pigs",pch=19,col="black",cex=2.4)
lines(x3,-x3/x1+C,col="red",lwd=2)
lines(x3,-x3/x1+C_lower,col="black",lwd=2)
lines(x3,-x3/x1+C_upper,col="black",lwd=2)
text(8.29404964010203+0.01,24+0.1,"2",col="yellow")
text(10.8197782844103+0.01,24+0.1,"1",col="yellow")
text(12.8992198260901+0.03,24+0.1,"9",col="yellow")
text(8.29404964010203+0.01,48+0.1,"7",col="yellow")
text(10.8197782844103+0.01,48+0.1,"3",col="yellow")
text(12.8992198260901+0.03,48+0.1,"3",col="yellow")
text(5.29831736654804+0.02,72+0.1,"5",col="yellow")
text(5.29831736654804+0.02,96+0.1,"5",col="yellow")
text(8.29404964010203+0.01,72+0.1,"1",col="yellow")

# Close the device
dev.off()

# phi
MPPDdepos <- c(0.2389,0.2284,0.2515,0.2226,0.2606,0.2391,0.2562,0.2657,0.2384,0.252,0.2523,0.2464,0.2531,0.2502,0.2433,0.2587,0.2504,0.2435,0.2767,0.2492,0.2585,0.2384,0.2688,0.2658,0.236,0.2677,0.2349,0.25,0.2658,0.2555,0.2452,0.2535,0.2543,0.2686,0.2218,0.2569,0.2594,0.2504,0.2781,0.2564,0.2638,0.2525,0.2221,0.2147,0.2607,0.2263,0.2396,0.2518,0.2558,0.2644,0.2118,0.2497,0.2601,0.2686,0.2662,0.249,0.2402,0.2538,0.2567,0.2565,0.25,0.2306,0.2533,0.2665,0.2164,0.2467,0.2501,0.235,0.2633,0.2404,0.2585,0.2655,0.2575,0.2674,0.2234,0.263,0.2401,0.2429,0.2384,0.2418,0.2266,0.2496,0.2735,0.2414,0.2473,0.2349,0.2444,0.2448,0.2564,0.2571,0.2641,0.2383,0.1927,0.2544,0.2653,0.2162,0.2665,0.2349,0.2476,0.2456)
phi_dist <- fitdistrplus::fitdist(MPPDdepos,"beta",method="mle")
phi <- unname(coef(phi_dist)[1]/(coef(phi_dist)[1]+coef(phi_dist)[2]))

# Save plot with 600 DPI resolution
tiff("Fig2e.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

hist(MPPDdepos,freq=FALSE,xlim=c(0.16,0.30),xlab="Probability of deposition",main="Probability of deposition")
lines(seq(0.16,0.30,by=0.001),dbeta(seq(0.16,0.30,by=0.001),shape1=coef(phi_dist)[1],shape2=coef(phi_dist)[2]),lwd=2,col="red")

# Close the device
dev.off()

# Save plot with 600 DPI resolution
png("Fig2e.png", width = 8, height = 8, units = "in", res = 300)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

hist(MPPDdepos,freq=FALSE,xlim=c(0.16,0.30),xlab="Probability of deposition",main="Probability of deposition")
lines(seq(0.16,0.30,by=0.001),dbeta(seq(0.16,0.30,by=0.001),shape1=coef(phi_dist)[1],shape2=coef(phi_dist)[2]),lwd=2,col="red")

# Close the device
dev.off()

# T
T <- phi*(x1+lambda)*exp(C*x1)/theta

# Now do plotting for t_T

# t_T
t <- function(d){
  return(C - log(d) / x1)
}

# Get dose range
dose_vals <- seq(1, 100, by = 1)

# Get point estimates
t_vals <- t(seq(1, 100, by = 1))

# Now get in form we want
line_data <- data.frame(A = dose_vals, B = t_vals)

# Now get 
# Read in stochastic results
t_list <- readRDS("~/PhD-work/Chapter 2/Data/Parameter distributions/t_T_dist.rds")

# Get TDR data in form we want
d <- seq(1,100,by=1)
t <- seq(0,150,by=0.5)
val <- c()
for (i in t){
  print(i)
  for (j in d){
    val <- append(val,length(which(t_list[[j]]<i))/10000)
  }
}
data <- expand.grid(X=d,Y=t)
data$Z <- val
data[length(d)*length(t)+1,] <- c(1000,1000,1)

# Save plot with 600 DPI resolution
tiff("Fig_tT.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Now plot the data (both point estimates and distributional estimates)
ggplot(data, aes(X, Y, fill = Z)) +
  geom_tile() +
  geom_line(
    data = line_data,
    aes(x = A, y = B),     # Replace with your actual column names for the line
    color = "black",
    size = 1.2,
    inherit.aes = FALSE    # Prevent it from trying to use `fill = Z`
  ) +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 150), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20, margin = margin(b = 10)),
    legend.text = element_text(size = 20)
  )

# Close the dev
dev.off()

# Now get some confidence intervals
c(quantile(t_list[[1]], p = 0.025), quantile(t_list[[1]], p = 0.975))
c(quantile(t_list[[9]], p = 0.025), quantile(t_list[[9]], p = 0.975))
c(quantile(t_list[[100]], p = 0.025), quantile(t_list[[100]], p = 0.975))
