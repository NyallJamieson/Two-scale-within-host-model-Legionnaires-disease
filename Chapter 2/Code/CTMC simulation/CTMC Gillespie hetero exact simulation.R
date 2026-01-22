load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_lambda.RData")

cdf_lambda <- cumsum(ker_lambda$y)/sum(ker_lambda$y)
sample_from_dens_lambda <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_lambda,ker_lambda$x,xout=u)$y
  return(sampled_values)
}

lambda_dist <- read.csv("~/PhD-work/Chapter 2/Data/Rupture data/lambda_dist.csv")[,2]

results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(results) <- c("L","time")

for (i in 1:500){

MODEL.onestep <- function(x,params){   
               alpha <- params["alpha"]
               beta <- params["beta"]
	             lambda <- params["lambda"]
	             G <- params["G"]
               L <- x[2]
               M <- x[3]
               rates <- c(
                              sphago <- alpha*L,
                              fphago <- beta*L,
                              burst <- lambda*M
               )
               total.rates <- sum(rates)
               if (total.rates==0)
                              tau <- -Inf
               else
                              tau <- rexp(n=1,rate=total.rates)
               transitions <- list(
                              successfulPHA <- c(-1,1),
                              failedPHA <- c(-1,0),
                              bursted <- c(rpois(lambda=G,n=1),-1)
               )
               event <- sample.int(n=3,size=1,prob=rates/total.rates)
               x+c(tau,transitions[[event]])
}

illcount <- 0
illtime <- c()
varName1 <- "illcount"
varName2 <- "illtime"

MODEL.simul <- function(y,maxstep=10000000){
	     dose <- unname(xstart["L"])
	     x <- xstart
	     lambda <- sample_from_dens_lambda(1)
	     h <- lambda
             A_w <- MASS::mvrnorm(n=1,mu=c(177.95295,0.19242),Sigma=matrix(c(82.49019097,-0.0254690939,-0.0254690939,0.0001933546),2,2))
	     A <- A_w[1]
	     w <- A_w[2]
	     G <- A/(1+(A-1)*exp(-w*(1/lambda-1)))

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

	     p <- c("lambda"=lambda,"alpha"=alpha,"beta"=beta,"G"=G)
	     names(p) <- c("lambda","alpha","beta","G")


	     theta <- sqrt((h-alpha-beta)^2+4*h*alpha*G)
             x1 <- 0.5*(-h-alpha-beta+theta)

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
C <- rnorm(n=1,mean=unname(coef(model)),sd=sqrt((unname(predict(model,se.fit=TRUE)$se.fit[1]))^2+(model$scale)^2))

	     phi <- rbeta(n=1,shape=191.18,shape2=577.26)
	     threshold <- unname(phi*(x1+lambda)*exp(C*x1)/theta)
print(c(alpha,beta,lambda,G,C,phi,threshold))

             illcount.copy <- get(varName1)
             illtime.copy <- get(varName2)
             output <- array(dim=c(maxstep+1,3))
             colnames(output) <- names(x)
             output[1,] <- x
             k <- 1
             while ((k <= maxstep) && (x["L"] > 0 | x["M"] >0) && (x["L"] < threshold)){
                            k <- k+1
                            output[k,] <- x <- MODEL.onestep(x,params=p)}
             t <- unname(output[k,1])
             if (x["L"]!=0){
                  return(c("L"=dose,"time"=t))
             }}

print(i)
xstart <- c(time=0,L=i,M=0)
sim_results <- as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,MODEL.simul,mc.cores=10)),na.action="omit"))
sim_results[,1] <- as.numeric(sim_results[,1])
sim_results[,2] <- as.numeric(sim_results[,2])

results <- dplyr::bind_rows(results,sim_results)}
