library(expint)
library(dplyr)

T <- 250
sim <- 10000

k <- 5
theta <- 1

#functions
est_k <- function(x) return(log(x) - digamma(x) - log(mean(yg)) + mean(log(yg)))
f <- function(S, k, theta) return(1-1/gamma(k)*gammainc(k,S/theta)+S/(theta*gamma(k))*gammainc(k-1,S/theta)-tau)
psi1 <- function (S, y=dane) return (1/T * sum(((y-S)^2/y)*1*(y>=S) + 2*(1-tau)*(S-y)))

ev_gamma <- function(S, shape = k, theta = 1) return(1/gamma(shape)*gammainc(shape,S/theta)-S/(theta*gamma(shape))*gammainc(shape-1,S/theta)-1+tau)

est_k <- function(x) return(log(x) - digamma(x) - log(mean(dane)) + mean(log(dane)))

#assumptions

taus <- c(.9, .925, .95, .975, .99, .995, .999)
rows <- 10
CI_ml <- matrix(0, nrow = rows, ncol = length(taus))
CI_m <- matrix(0, nrow = rows, ncol = length(taus))
alphas <- c(0.2, 0.1, 0.05, 0.01)


for(a in 1:length(alphas)){
  alpha <- alphas[a]
  for(row in 1:rows){
  for(t in 1:length(taus)){
    tau <- taus[t]
    Smls <- rep(0, sim)
    Sms <- rep(0, sim)
    CML <- rep(0, sim)
    CM <- rep(0, sim)
    ss <- matrix(0, nrow = 5, ncol = 2)
    for(i in 1:5){
      w <- uniroot(ev_gamma, interval = c(0.1, 20*i), shape = 5, theta = 1)
      ss[i,1] <- w$root
      ss[i,2] <- w$f.root
    }
    S <- ss[which.min(ss[,2]), 1]
  for(i in 1:sim){
    dane <- rgamma(T, shape = k, scale = theta)
    kk <- matrix(0, nrow = 3, ncol = 2)
    for(j in 1:3){
      ww <- uniroot(est_k, interval = c(0.1, quantile(dane, (0.95 + 0.01*j))))
      kk[j,1] <- ww$root
      kk[j,2] <- ww$f.root
    }
    ke <- kk[which.min(kk[,2]), 1]
    te <- mean(dane)/ke
    S_ml <- uniroot(ev_gamma, c(0.1,100), shape = ke, theta = te)$root
    Smls[i] <- S_ml
    h = 1e-10
    gS <- (1/h)*(ev_gamma(S_ml+h, shape = ke, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gk <- (1/h)*(ev_gamma(S_ml, shape = ke+h, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gth <- (1/h)*(ev_gamma(S_ml, shape = ke, theta = te+h) - ev_gamma(S_ml, shape = ke, theta = te))
    var_ML <- (ke*gk^2-2*te*gk*gth+gth^2*te^2*trigamma(ke))/(T*(ke*trigamma(ke)-1)*(gS^2))  
    if (between(S, S_ml-qnorm(1-alpha/2)*(var_ML^0.5), S_ml+qnorm(1-alpha/2)*(var_ML^0.5)) == TRUE) (CML[i] <- 1)
    
    S_m <- optimize(psi1, interval = c(0, 100))$minimum
    Sms[i] <- S_m
    var_m <- (gamma(k)*theta^2*((S_m/theta)^2*gammainc(k-2, S_m/theta)-(2*S_m*tau)/(theta)*gammainc(k-1, S_m/theta)+
                                  gamma(k)*(1-2*tau)*(1-gammainc(k, S_m/theta)/gamma(k))+gamma(k)*tau^2))/(T*(gammainc(k-1, S_m/theta))^2)
    if (between(S, S_m-qnorm(1-alpha/2)*(var_m^0.5), S_m+qnorm(1-alpha/2)*(var_m^0.5)) == TRUE) (CM[i] <- 1)
  
  }
    CI_ml[row, t] <- mean(CML)
    CI_m[row, t] <- mean(CM)
  }
  }
  if(a==1){
    plot(round(colMeans(CI_ml),2)~taus, type= 'l', ylim = c(0.75,1), ylab = 'coverage probability', lwd = 2, xlab = 'fill rate')
    lines(round(colMeans(CI_m),2)~taus, type = 'l', col = 'red', lwd = 2)
  }else{
    lines(round(colMeans(CI_ml),2)~taus, type= 'l', col = 'black', lwd = 2)
    lines(round(colMeans(CI_m),2)~taus, type = 'l', col = 'red', lwd = 2)
  }
  
}

