# Code for obtaining distributions for each parameter to be used within the CTMC model.

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
logdose <- log(c(13584.15,588828.58,11207102.05,82907236.06))
L0 <- 10^5
yourdata2 <- data.frame(time,dose)
# Now stick in other parameters (also calculate theta, x1, C, T)

results <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(results) <- c("alpha_dist","beta_dist","G_dist","h_dist","A_dist","w_dist")


part2 <- function(y){
h <- sample(size=1,lambda_dist,replace=TRUE)
A_w <- MASS::mvrnorm(n=1,mu=c(177.95295,0.19242),Sigma=matrix(c(82.49019097,-0.0254690939,-0.0254690939,0.0001933546),2,2))
A <- A_w[1]
w <- A_w[2]
G <- A/(1+(A-1)*exp(-w*(1/h-1)))

reps <- 500
ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((h-gamma)^2+4*h*Pi*gamma*G)+log(0.5*(h-gamma+((h-gamma)^2+4*h*gamma*Pi*G)^0.5)*exp(0.5*time*(-h-gamma+((h-gamma)^2+4*h*gamma*Pi*G)^0.5))-0.5*(h-gamma-((h-gamma)^2+4*h*gamma*Pi*G)^0.5)*exp(0.5*time*(-h-gamma-((h-gamma)^2+4*h*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
ab <- data.frame(MASS::mvrnorm(n=reps,mu=c(coef(ab_model)[1],coef(ab_model)[2]),Sigma=matrix(c(vcov(ab_model)[1],vcov(ab_model)[2],vcov(ab_model)[3],vcov(ab_model)[4]),2,2)))


alpha_dist <- ab[,1]*ab[,2]
beta_dist <- ab[,1]*(1-ab[,2])
G_dist <- rep(G,reps)
h_dist <- rep(h,reps)
A_dist <- rep(A,reps)
w_dist <- rep(w,reps)
return(data.frame(alpha_dist,beta_dist,G_dist,h_dist,A_dist,w_dist))
}

results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:10000,part2,mc.cores=5)),na.action="omit")))
######################################################
# remove dodgy results and save
alpha_dist <- results[,1]
beta_dist <- results[,2]
int1 <- which(alpha_dist>0)
int2 <- which(beta_dist>0)
int3 <- intersect(int1,int2)

results <- results[int3,]
######################################################
# this is to repeat all rows so that we can sample C 100 times from normal distribution to incorporate uncertainty
Repeated <- results |> dplyr::slice(rep(1:nrow(results), each = 100))

alpha <- results[,1]
beta <- results[,2]
G <- results[,3]
h <- results[,4]
theta <- sqrt((h-alpha-beta)^2+4*h*alpha*G)
x1 <- 0.5*(-h-alpha-beta+theta)


d3 <- c(4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,5*10^4,5*10^4,5*10^4,5*10^4,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2)
t3 <- 24*c(1,1,1,1,1,1,1,1,1,2,2,2,1,2,2,2,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4)
t_L <- t3-24
t_R <- t3

C_vec <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(C_vec)="C_dist"

part3 <- function(y){
x2 <- x1[y]
rhs <- log(d3)/x2
left <- t_L+rhs
right <- t_R+rhs
Data <- data.frame(left,right)
# Create an interval-censored survival object
surv_obj <- survival::Surv(time = left, time2 = right, type = "interval2")
model <- survival::survreg(surv_obj ~ 1, data = Data, dist = "gaussian")
return(rnorm(n=100,mean=unname(coef(model)),sd=unname(predict(model,se.fit=TRUE)$se.fit[1])))
}
C <- unlist(parallel::mclapply(1:length(x1),part3,mc.cores=5))

Repeated <- Repeated |> tidyr::drop_na()
Repeated <- Repeated[which(!is.na(C)),]
results2 <- Repeated
results3 <- data.frame(results2,na.omit(C))

# Remove other variables now for ram
rm(Repeated)
rm(results)
rm(results2)
gc()

theta <- sqrt((results3[,4]-results3[,1]-results3[,2])^2+4*results3[,4]*results3[,1]*results3[,3])
x1 <- 0.5*(-results3[,4]-results3[,1]-results3[,2]+theta)
phi <- rbeta(n=nrow(results3),shape=191.18,shape2=577.26)
T <- phi*(x1+results3[,4])*exp(results3[,7]*x1)/theta

results3 <- results3 |> dplyr::mutate("phi"=phi)
results3 <- results3 |> dplyr::mutate("T"=T)
results3 <- results3 |> tidyr::drop_na()

# Create a distribution for t (deterministic incubation period)

t_list <- list()

# Save the distribution of deterministic incubation period
for (d in 1:100){
  print(d)
  t_vec <- c()
  for (i in 1:10000){
    row <- sample.int(n = nrow(results3), size = 1)
    t_vec <- append(t_vec, results3[row,7] - log(d) / (0.5*(-results3[row,4]-results3[row,1]-results3[row,2]+sqrt((results3[row,4]-results3[row,1]-results3[row,2])^2+4*results3[row,4]*results3[row,1]*results3[row,3]))))
  }
  t_list[[d]] <- t_vec
}

# Now save the results
saveRDS(t_list, "~/PhD-work/Chapter 2/Data/Parameter distributions/t_T_dist.rds")

# Plot for the paper
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

# Save data in use-able form
saveRDS(data, "~/PhD-work/Chapter 2/Data/Parameter distributions/t_T_TDR.rds")

# Save plot with 600 DPI resolution
tiff("Fig_tT.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

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
    legend.text = element_text(size = 20)      # Legend labels
  )

# Close the dev
dev.off()

# Now get some confidence intervals
c(quantile(t_list[[1]], p = 0.025), quantile(t_list[[1]], p = 0.975))
c(quantile(t_list[[9]], p = 0.025), quantile(t_list[[9]], p = 0.975))
c(quantile(t_list[[20]], p = 0.025), quantile(t_list[[20]], p = 0.975))

# FOR APPENDIX OF THESIS WE CAN DO SOME NICE PLOTS WITH THIS ABOVE

# Now create parameter kernels
ker_alpha <- density(results3[,1],bw=0.007)
ker_beta <- density(results3[,2],bw=0.02)
ker_G <- density(results3[,3],bw=2)
ker_lambda <- density(results3[,4],bw=0.0003)
ker_A <- density(results3[,5],bw=3)
ker_omega <- density(results3[,6],bw=0.003)
ker_C <- density(results3[which(results3[,7]<200),7])
ker_phi <- density(results3[,8])
ker_T <- density(results3[,9])

# Save the kernel densities
save(ker_alpha,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_alpha.RData")
save(ker_beta,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_beta.RData")
save(ker_G,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_G.RData")
save(ker_lambda,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_lambda.RData")
save(ker_A,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_A.RData")
save(ker_omega,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_omega.RData")
save(ker_C,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_C.RData")
save(ker_phi,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_phi.RData")
save(ker_T,file="~/PhD-work/Chapter 2/Data/Parameter distributions/ker_T.RData")

# Plot the kernel densities
plot(ker_alpha,lwd=2,xlim=c(0,0.26),col="black",xlab=expression(alpha),main=expression(paste("Distribution of ", alpha)))
polygon(c(ker_alpha$x,rev(ker_alpha$x)),c(ker_alpha$y,rep(0,length(ker_alpha$y))),col="skyblue",border="black")

plot(ker_beta,lwd=2,xlim=c(0,1),col="black",xlab=expression(beta),main=expression(paste("Distribution of ", beta)))
polygon(c(ker_beta$x,rev(ker_beta$x)),c(ker_beta$y,rep(0,length(ker_beta$y))),col="skyblue",border="black")

plot(ker_G,lwd=2,col="black",xlab=expression("G"),main=expression("Distribution of G"))
polygon(c(ker_G$x,rev(ker_G$x)),c(ker_G$y,rep(0,length(ker_G$y))),col="skyblue",border="black")

plot(ker_lambda,lwd=2,col="black",xlab=expression(lambda),main=expression(paste("Distribution of ", lambda)))
polygon(c(ker_lambda$x,rev(ker_lambda$x)),c(ker_lambda$y,rep(0,length(ker_lambda$y))),col="skyblue",border="black")

plot(ker_A,lwd=2,col="black",xlab=expression("A"),main=expression("Distribution of A"))
polygon(c(ker_A$x,rev(ker_A$x)),c(ker_A$y,rep(0,length(ker_A$y))),col="skyblue",border="black")

plot(ker_w,lwd=2,col="black",xlab=expression(omega),main=expression(paste("Distribution of ", omega)))
polygon(c(ker_w$x,rev(ker_w$x)),c(ker_w$y,rep(0,length(ker_w$y))),col="skyblue",border="black")

plot(ker_C,lwd=2,xlim=c(140,200),col="black",xlab=expression("C"),main=expression("Distribution of C"))
polygon(c(ker_C$x,rev(ker_C$x)),c(ker_C$y,rep(0,length(ker_C$y))),col="skyblue",border="black")

plot(ker_phi,lwd=2,col="black",xlab=expression(phi),main=expression(paste("Distribution of ", phi)))
polygon(c(ker_phi$x,rev(ker_phi$x)),c(ker_phi$y,rep(0,length(ker_phi$y))),col="skyblue",border="black")

plot(ker_T,lwd=2,col="black",xlab=expression("T"),main=expression("Distribution of T"),xaxt="n",xlim=c(0,100000))
axis(1, at = seq(0,100000,by=20000), labels = c(0,20000,40000,60000,80000,"100000"))
polygon(c(ker_T$x,rev(ker_T$x)),c(ker_T$y,rep(0,length(ker_T$y))),col="skyblue",border="black")
