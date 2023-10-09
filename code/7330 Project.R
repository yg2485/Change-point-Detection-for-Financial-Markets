# Create the MCMC function based on t-shirkage prior
# Author: Yanghong Guo


# n -- number of nodes
# m -- defined in the paper
# y -- data

t_shrinkage <- function(y,m){
  n <- length(y)
  # Set prior settings
  a_sigma <- 0.5
  b_sigma <- 0.5
  a_t <- 2
  lambda_0 <- 5
  b_t <- m^2*a_t
  
  # Set initial values
  lambda_initial <- rep(1, n-1)
  sigma_initial <- var(y)
  theta_initial <- y
  
  T <- 5000
  B <- T/2
  
  # Create space to store results
  lambda_store <- matrix(NA, nrow = T, ncol = n-1)
  sigma_store <- rep(NA, T)
  theta_store <- matrix(NA, nrow = T, ncol = n)
  
  lambda <- lambda_initial
  sigma <-  sigma_initial
  theta <-  theta_initial
  
  # Run the Gibbs sampler
  
  for (t in 1:T){
    cat("Running iteration:",t,"\n")
    
    # Update lambda
    for (i in 1:(n-1)){
      lambda[i] <- 1 / rgamma(1, shape = a_t + 1/2 , rate = b_t + (theta[i+1]- theta[i])^2/2/sigma )
    }
    
    # Update sigma
    sse <- sum((y - theta)^2)
    theta_sum <-0
    for (j in 2:n){
      theta_sum <- theta_sum + (theta[j] - theta[j-1])^2/lambda[j-1]
    }
    sigma <- 1 / rgamma(1, shape = a_sigma + n, rate = b_sigma + theta[1]^2/2/lambda_0 + sse/2 + theta_sum/2)
    
    # Update theta
    nv_1 = 1 / (1/sigma + 1/lambda[1]/sigma + 1/lambda_0/sigma)
    mu_1 = nv_1*(y[1]/sigma + theta[2]/lambda[1]/sigma)
    theta[1] <- rnorm(1, mean = mu_1, sd = sqrt(nv_1))
    
    for (k in 2:(n-1)){
      nv = 1 / (1/sigma + 1/lambda[k-1]/sigma + 1/lambda[k]/sigma)
      mu = nv*(y[k]/sigma + theta[k-1]/lambda[k-1]/sigma + theta[k+1]/lambda[k]/sigma)
      theta[k] <- rnorm(1, mean = mu, sd = sqrt(nv))
    }
    
    nv_n = 1 / (1/sigma + 1/lambda[n-1]/sigma)
    mu_n = nv_n*(y[n]/sigma + theta[n-1]/lambda[n-1]/sigma)
    theta[n] <- rnorm(1, mean = mu_n, sd = sqrt(nv_n))
    
    # Store the results in this iteration
    lambda_store[t, ] <- lambda
    sigma_store[t] <- sigma
    theta_store[t,] <- theta
  }
  out <- NULL
  out$lambda <- lambda_store
  out$sigma <- sigma_store
  out$theta <- theta_store
  
  return(out)
}





##################################################################
##### Bitcoin data  2018 - 2020 #####

bitcoin <- read.csv("BTC-USD (1).csv", header = TRUE)
close <- rev(bitcoin$Close)
n_f <- length(close)
y <- rep(NA, n_f-1)
for (i in 1:n_f-1){
  y[i] <- (close[i+1]- close[1])/close[1]
}
n <- n_f - 1
m <- 0.012
#install.packages('metRology')
#library(metRology)
#1-1/2/n
#pt.scaled(sqrt(log(n)/n), 4, mean = 0, sd = 0.012, lower.tail = TRUE, log.p = FALSE)

results <- t_shrinkage(y,m)

plot(results$sigma, type = 'l')
means <- colMeans(results$theta[(B+1):T,])

sigma_post <- sqrt(mean(results$sigma[(B+1):T]))
critical <-qt(1-1/2/n, df = 2*a_t)*m*sigma_post*6

change <- NULL
for (i in 2:n){
  if(abs(means[i]-means[i-1])>critical){
    change <- c(change, i)
  }
}

plot(y, type = 'l', xlab ='Dates' , ylab = 'Cumulative Return for BTC ', xaxt = 'n')
axis(side = 1, at = c(0, 183, 365, 548, 730), labels = c('2018.06','2018.12','2019.06','2019.12', '2020.06'), tcl = -0.2)
title('Estimates by t-shrinkage prior')
lines(means, col = 'red')
abline(v = change, col = 'gray')

##################################################################
##### Dow Jones #####
djia <- read.csv("DJIA.csv", header = TRUE)
close <- rev(djia$Close)
n_f <- length(close)
y <- rep(NA, n_f-1)
for (i in 1:n_f-1){
  y[i] <- (close[i+1]- close[1])/close[1]
}
n <- n_f - 1
m <- 0.009
#install.packages('metRology')
library(metRology)
1-1/2/n
pt.scaled(sqrt(log(n)/n), 4, mean = 0, sd = 0.009, lower.tail = TRUE, log.p = FALSE)

results <- t_shrinkage(y,m)

plot(results$sigma, type = 'l')
plot(results$lambda[, 10], type = 'l')
plot(results$theta[, 15], type = 'l')
means <- colMeans(results$theta[(B+1):T,])

sigma_post <- sqrt(mean(results$sigma[(B+1):T]))
critical <-qt(1-1/2/n, df = 2*a_t)*m*sigma_post

change <- NULL
for (i in 2:n){
  if(abs(means[i]-means[i-1])>critical){
    change <- c(change, i)
  }
}

plot(y, type = 'l', xlab ='Dates' , ylab = 'Cumulative Return for DJI ', xaxt = 'n')
axis(side = 1, at = c(0, 252, 504, 756, 1000), labels = c('2007.01','2008.01','2009.01','2010.01', '2010.12'), tcl = -0.2)
title('Estimates by t-shrinkage prior')
lines(means, col = 'red')
abline(v = change, col = 'gray')


##################################################################
######## Simulated Data ########

# y is simulated data set

results2 <- t_shrinkage(y, m)
plot(y, type = 'l')
means <- colMeans(results2$theta[(B+1):T,])
lines(means, type = 'l', col = 'red')
legend(x = "bottomright", legend=c("Simulated data", "Estimated mean"),
       col=c("black", "red"), lty=1:2, cex=0.8)
title('Estimates by t-shrinkage prior')
