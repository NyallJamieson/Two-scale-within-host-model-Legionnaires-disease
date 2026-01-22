results <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(results) <- c("alpha","beta","time")

xstart <- c(time=0,L=9,M=0)

for (i in 1:100){
print(i)

phi <- 0.249
lambda=0.04050593
A <- 177.95295
w <- 0.192
G <- 62
time <- c(2,24,48,72)
logdose <- log(c(13584.15,588828.58,11207102.05,82907236.06))
L0 <- 10^5
yourdata2 <- data.frame(time,logdose)
ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
gamma_pi <- MASS::mvrnorm(n=1,mu=c(coef(ab_model)[1],coef(ab_model)[2]),Sigma=matrix(c(vcov(ab_model)[1],vcov(ab_model)[2],vcov(ab_model)[3],vcov(ab_model)[4]),2,2))
if (min(gamma_pi)<0){
while (min(gamma_pi)<0){
ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
gamma_pi <- MASS::mvrnorm(n=1,mu=c(coef(ab_model)[1],coef(ab_model)[2]),Sigma=matrix(c(vcov(ab_model)[1],vcov(ab_model)[2],vcov(ab_model)[3],vcov(ab_model)[4]),2,2))
}
}
alpha <- unname(gamma_pi[1]*gamma_pi[2])
beta <- unname(gamma_pi[1]*(1-gamma_pi[2]))

theta <- ((lambda-alpha-beta)^2+4*lambda*alpha*G)^0.5
x1 <- 0.5*(-lambda-alpha-beta+theta)
params <- c("lambda"=lambda,"alpha"=alpha,"beta"=beta,"G"=G)

d3 <- c(4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,4*10^5,5*10^4,5*10^4,5*10^4,5*10^4,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,4*10^3,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2,2*10^2)
t3 <- 24*c(1,1,1,1,1,1,1,1,1,2,2,2,1,2,2,2,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4)
t_L <- t3-24
t_R <- t3
rhs <- log(d3)/x1
Data <- data.frame(t_L,t_R,rhs)
vals <- log(d3)/x1
# Create an interval-censored survival object
surv_obj <- survival::Surv(time = t_L+vals, time2 = t_R+vals, type = "interval2")
model <- survival::survreg(surv_obj ~ 1, data = Data, dist = "gaussian")
C <- rnorm(n=1,mean=unname(coef(model)),sd=sqrt((unname(predict(model,se.fit=TRUE)$se.fit[1]))^2+(model$scale)^2))
threshold <- phi*(x1+lambda)*exp(x1*C)/theta

MODEL.onestep <- function(x,params){
               L <- x[2]
               M <- x[3]
               alpha <- params["alpha"]
               beta <- params["beta"]
               lambda <- params["lambda"]
                G <- params["G"]
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
                  return(c("alpha"=alpha,"beta"=beta,x["time"]))
              } else {return(c(NA,NA,NA))}
}

results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,MODEL.simul,mc.cores=2)),na.action="omit")))
}

sens_alpha <- results |>
  dplyr::group_by(alpha) |>
  dplyr::summarise(
    mean_inc = mean(time, na.rm = TRUE),
    occurrences = n()  # Count the number of rows in each group
  )

sens_beta <- results |>
  dplyr::group_by(beta) |>
  dplyr::summarise(
    mean_inc = mean(time, na.rm = TRUE),
    occurrences = n()  # Count the number of rows in each group
  )

write.csv(sens_alpha, "~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_alpha.csv")
write.csv(sens_beta, "~/PhD-work/Chapter 2/Data/Sensitivity analysis/sens_beta.csv")
