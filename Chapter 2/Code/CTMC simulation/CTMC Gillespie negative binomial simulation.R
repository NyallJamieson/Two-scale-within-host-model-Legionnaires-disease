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
                              bursted <- c(rnbinom(n=1,mu=62.792,size=32.072),-1)
               )
               event <- sample.int(n=3,size=1,prob=rates/total.rates)
               x+c(tau,transitions[[event]])
}

illcount <- 0
illtime <- c()
varName1 <- "illcount"
varName2 <- "illtime"

MODEL.simul <- function(y,maxstep=10000000,threshold=50661){
	       dose <- xstart["L"]
	       x <- xstart
               illcount.copy <- get(varName1)
               illtime.copy <- get(varName2)
               output <- array(dim=c(maxstep+1,3))
               colnames(output) <- names(x)
               output[1,] <- x
               k <- 1
               while ((k <= maxstep) && (x["L"] > 0 | x["M"] > 0) && (x["L"] < threshold)){
                              k <- k+1
                              output[k,] <- x <- MODEL.onestep(x,params)}
               t <- output[k,1]
               if (x["L"]!=0){
            return(c(dose,x["time"]))
             }
}


params <- c(alpha=0.08944449,beta=1.088102,lambda=0.04050593)
results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(results) <- c("L","time")


for (i in 1:500){
print(i)
xstart <- c(time=0,L=i,M=0)
results <- dplyr::bind_rows(results,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,MODEL.simul,mc.cores=10)),na.action="omit")))
}
