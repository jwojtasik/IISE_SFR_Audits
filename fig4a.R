library(expint)
library(dplyr)

T <- 250
sim <- 10000

k <- 1
theta <- 1

#functions
est_k <- function(x) return(log(x) - digamma(x) - log(mean(yg)) + mean(log(yg)))
f <- function(S, k, theta) return(1-1/gamma(k)*gammainc(k,S/theta)+S/(theta*gamma(k))*gammainc(k-1,S/theta)-tau)
psi1 <- function (S, y=dane) return (1/T * sum(((y-S)^2/y)*1*(y>=S) + 2*(1-tau)*(S-y)))

ev_gamma <- function(S, shape = k, theta = 1) return(1/gamma(shape)*gammainc(shape,S/theta)-S/(theta*gamma(shape))*gammainc(shape-1,S/theta)-1+tau)

est_k <- function(x) return(log(x) - digamma(x) - log(mean(dane)) + mean(log(dane)))

#assumptions
alpha <- .9
taus <- c(.9, .925, .95, .975, .99, .995, .999)

#results 
cis <- matrix(0, ncol = 3, nrow = length(taus))
colnames(cis) <- c('tau', 'M-est 0.9', 'Ml-est 0.9')

tau <- taus[1]


exp_CML <- rep(0, sim)

exp_CI_ml <- matrix(0, nrow = 10, ncol =  length(taus))
exp_CI_m <- matrix(0, nrow = 10, ncol = length(taus))
alphas <- c(0.2, 0.1, 0.05, 0.01)

for(a in 1:length(alphas)){
  alpha <- alphas[a]
  print(a)
  for(row in 1:10){
    print(row)
  for(t in 1:length(taus)){
    tau <- taus[t]
    CML <- rep(0, sim)
    CM <- rep(0, sim)
    ss <- matrix(0, nrow = 5, ncol = 2)
    for(i in 1:5){
      w <- uniroot(ev_gamma, interval = c(0.1, 10*i), shape = 1, theta = 1)
      ss[i,1] <- w$root
      ss[i,2] <- w$f.root
    }
    S <- ss[which.min(ss[,2]), 1]
    for(i in 1:sim){
      dane <- rexp(T)
      te <- mean(dane)
      ss <- matrix(0, nrow = 5, ncol = 2)
      for(j in 1:5){
        w <- uniroot(ev_gamma, interval = c(0.1, 10*j), shape = 1, theta = te)
        ss[j,1] <- w$root
        ss[j,2] <- w$f.root
      }
      S_ml <- ss[which.min(ss[,2]), 1]
      var_ML <- S_ml^2/T  
      if (between(S, S_ml-qnorm(1-alpha/2)*(var_ML^0.5), S_ml+qnorm(1-alpha/2)*(var_ML^0.5)) == TRUE) (CML[i] <- 1)
      
      ss <- matrix(0, nrow = 5, ncol = 2)
      for(j in 1:5){
        w <- optimize(psi1, interval = c(0, 10*j))
        ss[j,1] <- w$minimum
        ss[j,2] <- w$objective
      }
      S_m <- ss[which.min(ss[,2]), 1]
      var_m <- ((1-tau)^2 + 2*tau*exp(-S_m/theta)-exp(-S_m/theta)-2*tau*S_m/theta*gammainc(0, S_m/theta)+S_m*
                  exp(-S_m/theta)-S_m^2*gammainc(0,S_m/theta))/(T*(gammainc(0, S_m/theta))^2)
      
      if (between(S, S_m-qnorm(1-alpha/2)*(var_m^0.5), S_m+qnorm(1-alpha/2)*(var_m^0.5)) == TRUE) (CM[i] <- 1)
      
    }
    exp_CI_ml[row,t] <- mean(CML)
    exp_CI_m[row,t] <- mean(CM)
  }
  }
  
  if(a==1){
    plot(round(colMeans(exp_CI_ml),2)~taus, type= 'l', ylim = c(0.75,1), ylab = 'coverage probability', lwd = 2, xlab = 'fill rate')
    lines(round(colMeans(exp_CI_m),2)~taus, type = 'l', col = 'red', lwd = 2)
  }else{
    lines(round(colMeans(exp_CI_ml),2)~taus, type= 'l', col = 'black', lwd = 2)
    lines(round(colMeans(exp_CI_m),2)~taus, type = 'l', col = 'red', lwd = 2)
  }
  
}
