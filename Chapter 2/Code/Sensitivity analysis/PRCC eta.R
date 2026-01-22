load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_eta.RData")

cdf_eta <- cumsum(ker_eta$y)/sum(ker_eta$y)
sample_from_dens_eta <- function(n){
u <- runif(n)
sampled_values <- approx(cdf_eta,ker_eta$x,xout=u)$y
return(sampled_values)
}

results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(results) <- c("eta","time")

xstart <- c(time=0,L=9,M=0)

for (i in 1:100){
print(i)
eta <- sample_from_dens_eta(1)

phi <- 0.249
lambda=0.04050593
A <- 177.95295
w <- 0.192
G <- 62
alpha=0.089
beta=1.088
theta <- ((lambda-alpha-beta)^2+4*lambda*alpha*G)^0.5
x1 <- 0.5*(-lambda-alpha-beta+theta)
params <- c("lambda"=lambda,"alpha"=alpha,"beta"=beta)
threshold <- phi*(x1+lambda)*exp(x1*eta)/theta

MODEL.onestep <- function(x,params){
               L <- x[2]
               M <- x[3]
               alpha <- params["alpha"]
               beta <- params["beta"]
               lambda <- params["lambda"]
               rates <- c(
                              sphago <- alpha*L,
                              fphago <- beta*L,
                              burst <- lambda*M
               )
               total.rates <- sum(rates)
               if (total.rates==0){
                              tau <- -Inf
               }else{
                              tau <- rexp(n=1,rate=total.rates)}
               transitions <- list(
                              successfulPHA <- c(-1,1),
                              failedPHA <- c(-1,0),
                              bursted <- c(G,-1)
               )
               event <- sample.int(n=3,size=1,prob=rates/total.rates)
               x+c(tau,transitions[[event]])
}

MODEL.simul <- function(y,maxstep=10000000){
	        dose <- 9
	        x <- xstart
                output <- array(dim=c(maxstep+1,3))
                colnames(output) <- names(x)
                output[1,] <- x
                k <- 1
                while ((k <= maxstep) && (x["L"] > 0 | x["M"] > 0) && (x["L"] < threshold)){
                               k <- k+1
                               output[k,] <- x <- MODEL.onestep(x,params)}
                t <- output[k,1]
                if (x["L"]!=0){
                   return(c("eta"=eta,x["time"]))
                } else {return(c(NA,NA))}
}

results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,MODEL.simul,mc.cores=5)),na.action="omit")))
}

sens_eta <- results |>
  group_by(eta) |>
  summarise(
    mean_inc = mean(time, na.rm = TRUE),
    occurrences = n()  # Count the number of rows in each group
  )

write.csv(sens_eta, "~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_eta.csv")