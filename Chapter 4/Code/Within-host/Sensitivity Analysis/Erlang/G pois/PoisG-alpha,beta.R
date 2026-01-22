results <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(results) <- c("alpha","beta","time")

xstart <- c(time=0,L=9,M1=0,M2=0,M3=0)

for (i in 1:1000){
  print(i)
  
  phi <- 0.249
  lam=0.04050593
  A <- 177.95295
  w <- 0.192
  G <- 62
  time <- c(2,24,48,72)
  logdose <- log(c(13584.15,588828.58,11207102.05,82907236.06))
  L0 <- 10^5
  yourdata2 <- data.frame(time,logdose)
  ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lam-gamma)^2+4*lam*Pi*gamma*G)+log(0.5*(lam-gamma+((lam-gamma)^2+4*lam*gamma*Pi*G)^0.5)*exp(0.5*time*(-lam-gamma+((lam-gamma)^2+4*lam*gamma*Pi*G)^0.5))-0.5*(lam-gamma-((lam-gamma)^2+4*lam*gamma*Pi*G)^0.5)*exp(0.5*time*(-lam-gamma-((lam-gamma)^2+4*lam*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
  gamma_pi <- MASS::mvrnorm(n=1,mu=c(coef(ab_model)[1],coef(ab_model)[2]),Sigma=matrix(c(vcov(ab_model)[1],vcov(ab_model)[2],vcov(ab_model)[3],vcov(ab_model)[4]),2,2))
  if (min(gamma_pi)<0){
    while (min(gamma_pi)<0){
      ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lam-gamma)^2+4*lam*Pi*gamma*G)+log(0.5*(lam-gamma+((lam-gamma)^2+4*lam*gamma*Pi*G)^0.5)*exp(0.5*time*(-lam-gamma+((lam-gamma)^2+4*lam*gamma*Pi*G)^0.5))-0.5*(lam-gamma-((lam-gamma)^2+4*lam*gamma*Pi*G)^0.5)*exp(0.5*time*(-lam-gamma-((lam-gamma)^2+4*lam*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
      gamma_pi <- MASS::mvrnorm(n=1,mu=c(coef(ab_model)[1],coef(ab_model)[2]),Sigma=matrix(c(vcov(ab_model)[1],vcov(ab_model)[2],vcov(ab_model)[3],vcov(ab_model)[4]),2,2))
    }
  }
  alpha <- unname(gamma_pi[1]*gamma_pi[2])
  beta <- unname(gamma_pi[1]*(1-gamma_pi[2]))
  
  params <- c(alpha=alpha,beta=beta,lambda=0.126,G=59.58144,theta=0.8554645)
  
  time_vec <- seq(0,300,by=0.01)
  time_vec <- round(time_vec, digits=2)
  
  MODEL.onestep <- function(x,params){
    L <- x[2]
    M1 <- x[3]
    M2 <- x[4]
    M3 <- x[5]
    alpha <- params["alpha"]
    beta <- params["beta"]
    lambda <- params["lambda"]
    G <- params["G"]
    theta <- params["theta"]
    rates <- c(
      sphago <- alpha*L,
      fphago <- beta*L,
      trans1 <- lambda*M1,
      trans2 <- lambda*M2,
      burst <- lambda*theta*M3
    )
    total.rates <- sum(rates)
    if (total.rates==0)
      tau <- -Inf
    else
      tau <- rexp(n=1,rate=total.rates)
    event <- sample.int(n=5,size=1,prob=rates/total.rates)
    transitions <- list(
      successfulPHA <- c(-1,1,0,0),
      failedPHA <- c(-1,0,0,0),
      transition1 <- c(0,-1,1,0),
      transition2 <- c(0,0,-1,1),
      rupture <- c(rpois(n=1,lambda=G),0,0,-1)
    )
    x+c(tau,transitions[[event]])}
  
  MODEL.simul <- function(y,maxstep=100000000,threshold=thresh){
    dose <- xstart["L"]
    x <- xstart
    output <- array(dim=c(10000000+1,5))
    colnames(output) <- names(x)
    output[1,] <- x
    k <- 1
    while ((x["L"] > 0 | x["M1"] > 0 | x["M2"] > 0 | x["M3"] >0) && (x["L"] < threshold)){
      k <- k+1
      output[k,] <- x <- MODEL.onestep(x,params)
    }
    t <- output[k,1]
    if (x["L"]!=0){
      return(c("alpha"=alpha,"beta"=beta,x["time"]))
    }
  }
  
  results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,MODEL.simul,mc.cores=128)),na.action="omit")))
}

write.csv(results, "SA_Erlang.Pois.ab.csv")