source("~/PhD-work/Chapter 4/Code/Within-macrophage/Estimate G for Pois in Erlang within-macrophage model.R")

time_vec <- seq(0,150,by=0.01)
time_vec <- round(time_vec, digits=2)

MODEL.onestep <- function(x,params){
  L <- x[2]
  M1 <- x[3]
  M2 <- x[4]
  M3 <- x[5]
  alpha <- params["alpha"]
  beta <- params["beta"]
  lambda <- params["lambda"]
  theta <- params["theta"]
  rates <- c(
    sphago <- alpha*L,
    fphago <- beta*L,
    trans1 <- lambda*M1,
    trans2 <- lambda*M2,
    burst <- lambda*theta*M3
  )
  total.rates <- sum(rates)
  if (total.rates==0){
    tau <- -Inf
  } else {
    tau <- rexp(n=1,rate=total.rates)
  event <- sample.int(n=5,size=1,prob=rates/total.rates)}
  if (event==5){
    time <- round(sum(rexp(n=3,rate=0.126)),digits=2)
    if (time > 300){
      while (time>300 | time<1){
        time <- round(sum(rexp(n=3,rate=0.126)),digits=2)}
    }
    if (time < 1){
      while (time>300 | time<1){
        time <- round(sum(rexp(n=3,rate=0.126)),digits=2)}
    }
    # need to do some tolerance instead of exact point due to R number issues
    size <- round(sample_G_time(sample_size=1,time=time))
    transitions <- list(
      bursted <- c(size,0,0,-1))
    x+c(tau,transitions[[event-4]])}
  else {
    transitions <- list(
      successfulPHA <- c(-1,1,0,0),
      failedPHA <- c(-1,0,0,0),
      transition1 <- c(0,-1,1,0),
      transition2 <- c(0,0,-1,1)
    )
    x+c(tau,transitions[[event]])}
}

MODEL.simul <- function(x,params,maxstep=100000000,threshold=50661){
  output <- array(dim=c(10000000+1,5))
  colnames(output) <- names(x)
  output[1,] <- x
  k <- 1
  while ((x["L"] > 0 | x["M1"] > 0 | x["M2"] > 0 | x["M3"] >0) && (x["L"] < threshold)){
    k <- k+1
    output[k,] <- x <- MODEL.onestep(x,params)
  }
  print("hey")
  t <- output[k,1]
  print(x["L"])
  if (x["L"]!=0){
    print("hi")
    return(t)}
}

params <- c(alpha=0.089,beta=1.088,lambda=0.126,theta=0.8554645)
xstart <- c(time=0,L=1,M1=0,M2=0,M3=0)