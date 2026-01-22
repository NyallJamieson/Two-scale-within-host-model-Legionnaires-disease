params <- c(alpha=0.089,beta=1.088,G=59.58144,theta=0.8554645)

results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(results) <- c("lambda","time")

xstart <- c(time=0,L=9,M1=0,M2=0,M3=0)

for (i in 1:3){
  
  print(i)
  
  lambda <- rnorm(n = 1, mean = 0.1260891, sd = 0.0007486)
  
  time_vec <- seq(0,300,by=0.01)
  time_vec <- round(time_vec, digits=2)
  
  MODEL.onestep <- function(x,params){
    L <- x[2]
    M1 <- x[3]
    M2 <- x[4]
    M3 <- x[5]
    alpha <- params["alpha"]
    beta <- params["beta"]
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
      rupture <- c(60,0,0,-1)
    )
    x+c(tau,transitions[[event]])}
  
  MODEL.simul <- function(y,maxstep=100000000,threshold=50661){
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
      return(c("lambda"=lambda,x["time"]))
    }
  }
  
  results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:5,MODEL.simul,mc.cores=2)),na.action="omit")))
}

write.csv(results, "SA_Erlang.Constant.lambda.csv")
