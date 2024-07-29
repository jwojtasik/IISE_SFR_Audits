library(expint)
library(dplyr)

est_k <- function(x) return(log(x) - digamma(x) - log(mean(yg)) + mean(log(yg)))
f <- function(S, k, theta) return(1-1/gamma(k)*gammainc(k,S/theta)+S/(theta*gamma(k))*gammainc(k-1,S/theta)-tau)
psi1 <- function (S, y=dane) return (1/T * sum(((y-S)^2/y)*1*(y>=S) + 2*(1-tau)*(S-y)))
ev_gamma <- function(S, shape = k, theta = 1) return(1/gamma(shape)*gammainc(shape,S/theta)-S/(theta*gamma(shape))*gammainc(shape-1,S/theta)-1+tau)
est_k <- function(x, dane = est_sample) return(log(x) - digamma(x) - log(mean(dane)) + mean(log(dane)))

### Exponential

k <- 1
theta <- 1

real_tau <- 0.99
real_alpha <- 0.2

T <- 250
N <- 250
sim <- 10000

corrects<- c(0, 0.0005, 0.00075, 0.001, 0.0025, 0.003, 0.004, 0.005, 0.0075, 0.009)
SS_tau_exp <- matrix(0, nrow = sim, ncol = length(corrects))
means_tau_exp <- rep(0, length(corrects))
SS_tau_exp_M <- matrix(0, nrow = sim, ncol = length(corrects))
means_tau_exp_M <- rep(0, length(corrects))

clear_S <- rep(0, length(corrects))
clear_SM <- rep(0, length(corrects))
var_S <- rep(0, length(corrects))
var_SM <- rep(0, length(corrects))

for(i in 1:length(corrects)){
  tau <- real_tau + corrects[i]
  print(tau)
  alpha <- real_alpha
  probM <- rep(0, sim)
  probML <- rep(0, sim)
  SS <- rep(0, sim)
  SSM <- rep(0, sim)
  vSS <- rep(0, sim)
  vSSM <- rep(0, sim)
  for(j in 1:sim){
    est_sample <- rexp(T, rate = 1/theta)
    #ML
    te <- mean(est_sample)
    S_ml <- uniroot(ev_gamma, c(0.1,100), shape = 1, theta = te)$root
    var_ML <- (S_ml^2)/T
    S1 <- S_ml+qnorm(1-alpha/2)*(var_ML^0.5)
    SS_tau_exp[j,i] <- S1
    SS[j] <- S_ml
    vSS[j] <- var_ML
    #M
    S_m <- optimize(psi1, interval = c(0, 100), y = est_sample)$minimum
    var_M <- (gamma(k)*theta^2*((S_m/theta)^2*gammainc(k-2, S_m/theta)-(2*S_m*tau)/(theta)*gammainc(k-1, S_m/theta)+
                                  gamma(k)*(1-2*tau)*(1-gammainc(k, S_m/theta)/gamma(k))+gamma(k)*tau^2))/(T*(gammainc(k-1, S_m/theta))^2)
    S2 <- S_m+qnorm(1-alpha/2)*(var_M^0.5)
    SS_tau_exp_M[j,i] <- S2
    SSM[j] <- S_m
    vSSM[j] <- var_M
     test_sample <- rexp(N, rate = 1/theta)
    if(mean(pmax(test_sample - S1, 0)/test_sample) > 1 - real_tau) (probML[j] <- 1)
    if(mean(pmax(test_sample - S2, 0)/test_sample) > 1 - real_tau) (probM[j] <- 1)
  }
  means_tau_exp[i] <- mean(probML)
  means_tau_exp_M[i] <- mean(probM)
  clear_S[i] <- mean(SS)
  clear_SM[i] <- mean(SSM)
  var_S[i] <- mean(vSS)
  var_SM[i] <- mean(vSSM)
}

corrects_alpha <- c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.19)
SS_alpha_exp <- matrix(0, nrow = sim, ncol = length(corrects_alpha))
means_alpha_exp <- rep(0, length(corrects_alpha))
SS_alpha_exp_M <- matrix(0, nrow = sim, ncol = length(corrects_alpha))
means_alpha_exp_M <- rep(0, length(corrects_alpha))

for(i in 1:length(corrects_alpha)){
  tau <- real_tau
  alpha <- real_alpha - corrects_alpha[i]
  probM <- rep(0, sim)
  probML <- rep(0, sim)

  for(j in 1:sim){
    
    est_sample <- rexp(T, rate = 1/theta)
    #ML
    te <- mean(est_sample)
    S_ml <- uniroot(ev_gamma, c(0.1,100), shape = 1, theta = te)$root
    var_ML <- (S_ml^2)/T
    S1 <- S_ml+qnorm(1-alpha/2)*(var_ML^0.5)

    SS_alpha_exp[j,i] <- S1
    
    #M
    S_m <- optimize(psi1, interval = c(0, 100), y = est_sample)$minimum
    var_M <- (gamma(k)*theta^2*((S_m/theta)^2*gammainc(k-2, S_m/theta)-(2*S_m*tau)/(theta)*gammainc(k-1, S_m/theta)+
                                  gamma(k)*(1-2*tau)*(1-gammainc(k, S_m/theta)/gamma(k))+gamma(k)*tau^2))/(T*(gammainc(k-1, S_m/theta))^2)
    S2 <- S_m+qnorm(1-alpha/2)*(var_M^0.5)
    SS_alpha_exp_M[j,i] <- S2

    test_sample <- rexp(N, rate = 1/theta)
    
    if(mean(pmax(test_sample - S1, 0)/test_sample) > 1 - real_tau) (probML[j] <- 1)
    if(mean(pmax(test_sample - S2, 0)/test_sample) > 1 - real_tau) (probM[j] <- 1)
  }
  means_alpha_exp[i] <- mean(probML)
  means_alpha_exp_M[i] <- mean(probM)

}


### Gamma(5,1)
    
k <- 5
theta <- 1

real_tau <- 0.99
real_alpha <- 0.2

T <- 250
N <- 250
sim <- 10000

corrects <- c(0, 0.0005, 0.00075, 0.001, 0.0025, 0.003, 0.004, 0.005, 0.0075, 0.009)
SS_tau5 <- matrix(0, nrow = sim, ncol = length(corrects))
means_tau5 <- rep(0, length(corrects))
SS_tau5_M <- matrix(0, nrow = sim, ncol = length(corrects))
means_tau5_M <- rep(0, length(corrects))

for(i in 1:length(corrects)){
  tau <- real_tau + corrects[i]
  alpha <- real_alpha
  probML <- rep(0, sim)
  probM <- rep(0, sim)
  for(j in 1:sim){
    est_sample <- rgamma(T, shape = k, scale = theta)
    #ML
    ke <- uniroot(est_k, c(0.001, 100))$root
    te <- mean(est_sample)/ke
    S_ml <- uniroot(ev_gamma, c(0.1,100), shape = ke, theta = te)$root
    h <- 1e-10
    gS <- (1/h)*(ev_gamma(S_ml+h, shape = ke, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gk <- (1/h)*(ev_gamma(S_ml, shape = ke+h, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gth <- (1/h)*(ev_gamma(S_ml, shape = ke, theta = te+h) - ev_gamma(S_ml, shape = ke, theta = te))
    var_ML <- (ke*gk^2-2*te*gk*gth+gth^2*te^2*trigamma(ke))/(T*(ke*trigamma(ke)-1)*(gS^2)) 
    S1 <- S_ml+qnorm(1-alpha/2)*(var_ML^0.5)
    SS_tau5[j,i] <- S1
    #M
    S_m <- optimize(psi1, interval = c(0, 100), y = est_sample)$minimum
    var_M <- (gamma(k)*theta^2*((S_m/theta)^2*gammainc(k-2, S_m/theta)-(2*S_m*tau)/(theta)*gammainc(k-1, S_m/theta)+
                                  gamma(k)*(1-2*tau)*(1-gammainc(k, S_m/theta)/gamma(k))+gamma(k)*tau^2))/(T*(gammainc(k-1, S_m/theta))^2)
    S2 <- S_m+qnorm(1-alpha/2)*(var_M^0.5)
    SS_tau5_M[j,i] <- S2
    
    test_sample <- rgamma(N, shape = k, scale = theta)
    if(mean(pmax(test_sample - S1, 0)/test_sample) > 1 - real_tau) (probML[j] <- 1)
    if(mean(pmax(test_sample - S2, 0)/test_sample) > 1 - real_tau) (probM[j] <- 1)
  }
  means_tau5[i] <- mean(probML)
  means_tau5_M[i] <- mean(probM)
}

corrects_alpha <- c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.19)
SS_alpha5 <- matrix(0, nrow = sim, ncol = length(corrects_alpha))
means_alpha5 <- rep(0, length(corrects_alpha))
SS_alpha5_M <- matrix(0, nrow = sim, ncol = length(corrects_alpha))
means_alpha5_M <- rep(0, length(corrects_alpha))

for(i in 1:length(corrects_alpha)){
  tau <- real_tau
  alpha <- real_alpha - corrects_alpha[i]
  probML <- rep(0, sim)
  probM <- rep(0, sim)
  for(j in 1:sim){
    est_sample <- rgamma(T, shape = k, scale = theta)
    #ML
    ke <- uniroot(est_k, c(0.001, 100))$root
    te <- mean(est_sample)/ke
    S_ml <- uniroot(ev_gamma, c(0.1,100), shape = ke, theta = te)$root
    h <- 1e-10
    gS <- (1/h)*(ev_gamma(S_ml+h, shape = ke, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gk <- (1/h)*(ev_gamma(S_ml, shape = ke+h, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gth <- (1/h)*(ev_gamma(S_ml, shape = ke, theta = te+h) - ev_gamma(S_ml, shape = ke, theta = te))
    var_ML <- (ke*gk^2-2*te*gk*gth+gth^2*te^2*trigamma(ke))/(T*(ke*trigamma(ke)-1)*(gS^2)) 
    S1 <- S_ml+qnorm(1-alpha/2)*(var_ML^0.5)
    SS_alpha5[j,i] <- S1
    
    S_m <- optimize(psi1, interval = c(0, 100), y = est_sample)$minimum
    var_M <- (gamma(k)*theta^2*((S_m/theta)^2*gammainc(k-2, S_m/theta)-(2*S_m*tau)/(theta)*gammainc(k-1, S_m/theta)+
                                  gamma(k)*(1-2*tau)*(1-gammainc(k, S_m/theta)/gamma(k))+gamma(k)*tau^2))/(T*(gammainc(k-1, S_m/theta))^2)
    S2 <- S_m+qnorm(1-alpha/2)*(var_M^0.5)
    SS_alpha5_M[j,i] <- S2
    
    test_sample <- rgamma(N, shape = k, scale = theta)
    if(mean(pmax(test_sample - S1, 0)/test_sample) > 1 - real_tau) (probML[j] <- 1)
    if(mean(pmax(test_sample - S2, 0)/test_sample) > 1 - real_tau) (probM[j] <- 1)
  }
  means_alpha5[i] <- mean(probML)
  means_alpha5_M[i] <- mean(probM)
}
    
### Gamma(30,1)
    
    k <- 30
theta <- 1

real_tau <- 0.99
real_alpha <- 0.2

T <- 250
N <- 250
sim <- 10000

corrects <- c(0, 0.0005, 0.00075, 0.001, 0.0025, 0.003, 0.004, 0.005, 0.0075, 0.009)
SS_tau30 <- matrix(0, nrow = sim, ncol = length(corrects))
means_tau30 <- rep(0, length(corrects))
SS_tau30_M <- matrix(0, nrow = sim, ncol = length(corrects))
means_tau30_M <- rep(0, length(corrects))

for(i in 1:length(corrects)){
  tau <- real_tau + corrects[i]
  alpha <- real_alpha
  probML <- rep(0, sim)
  probM <- rep(0, sim)
  for(j in 1:sim){
    est_sample <- rgamma(T, shape = k, scale = theta)
    #ML
    ke <- uniroot(est_k, c(0.001, 100))$root
    te <- mean(est_sample)/ke
    S_ml <- uniroot(ev_gamma, c(0.1,100), shape = ke, theta = te)$root
    h <- 1e-10
    gS <- (1/h)*(ev_gamma(S_ml+h, shape = ke, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gk <- (1/h)*(ev_gamma(S_ml, shape = ke+h, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gth <- (1/h)*(ev_gamma(S_ml, shape = ke, theta = te+h) - ev_gamma(S_ml, shape = ke, theta = te))
    var_ML <- (ke*gk^2-2*te*gk*gth+gth^2*te^2*trigamma(ke))/(T*(ke*trigamma(ke)-1)*(gS^2)) 
    S1 <- S_ml+qnorm(1-alpha/2)*(var_ML^0.5)
    SS_tau30[j,i] <- S1
    #M
    S_m <- optimize(psi1, interval = c(0, 100), y = est_sample)$minimum
    var_M <- (gamma(k)*theta^2*((S_m/theta)^2*gammainc(k-2, S_m/theta)-(2*S_m*tau)/(theta)*gammainc(k-1, S_m/theta)+
                                  gamma(k)*(1-2*tau)*(1-gammainc(k, S_m/theta)/gamma(k))+gamma(k)*tau^2))/(T*(gammainc(k-1, S_m/theta))^2)
    S2 <- S_m+qnorm(1-alpha/2)*(var_M^0.5)
    SS_tau30_M[j,i] <- S2
    
    test_sample <- rgamma(N, shape = k, scale = theta)
    if(mean(pmax(test_sample - S1, 0)/test_sample) > 1 - real_tau) (probML[j] <- 1)
    if(mean(pmax(test_sample - S2, 0)/test_sample) > 1 - real_tau) (probM[j] <- 1)
  }
  means_tau30[i] <- mean(probML)
  means_tau30_M[i] <- mean(probM)
}

corrects_alpha <- c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.19)
SS_alpha30 <- matrix(0, nrow = sim, ncol = length(corrects_alpha))
means_alpha30 <- rep(0, length(corrects_alpha))
SS_alpha30_M <- matrix(0, nrow = sim, ncol = length(corrects_alpha))
means_alpha30_M <- rep(0, length(corrects_alpha))

for(i in 1:length(corrects_alpha)){
  tau <- real_tau
  alpha <- real_alpha - corrects_alpha[i]
  probML <- rep(0, sim)
  probM <- rep(0, sim)
  for(j in 1:sim){
    est_sample <- rgamma(T, shape = k, scale = theta)
    #ML
    ke <- uniroot(est_k, c(0.001, 100))$root
    te <- mean(est_sample)/ke
    S_ml <- uniroot(ev_gamma, c(0.1,100), shape = ke, theta = te)$root
    h <- 1e-10
    gS <- (1/h)*(ev_gamma(S_ml+h, shape = ke, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gk <- (1/h)*(ev_gamma(S_ml, shape = ke+h, theta = te) - ev_gamma(S_ml, shape = ke, theta = te))
    gth <- (1/h)*(ev_gamma(S_ml, shape = ke, theta = te+h) - ev_gamma(S_ml, shape = ke, theta = te))
    var_ML <- (ke*gk^2-2*te*gk*gth+gth^2*te^2*trigamma(ke))/(T*(ke*trigamma(ke)-1)*(gS^2)) 
    S1 <- S_ml+qnorm(1-alpha/2)*(var_ML^0.5)
    SS_alpha30[j,i] <- S1
    
    S_m <- optimize(psi1, interval = c(0, 100), y = est_sample)$minimum
    var_M <- (gamma(k)*theta^2*((S_m/theta)^2*gammainc(k-2, S_m/theta)-(2*S_m*tau)/(theta)*gammainc(k-1, S_m/theta)+
                                  gamma(k)*(1-2*tau)*(1-gammainc(k, S_m/theta)/gamma(k))+gamma(k)*tau^2))/(T*(gammainc(k-1, S_m/theta))^2)
    S2 <- S_m+qnorm(1-alpha/2)*(var_M^0.5)
    SS_alpha30_M[j,i] <- S2
    
    test_sample <- rgamma(N, shape = k, scale = theta)
    if(mean(pmax(test_sample - S1, 0)/test_sample) > 1 - real_tau) (probML[j] <- 1)
    if(mean(pmax(test_sample - S2, 0)/test_sample) > 1 - real_tau) (probM[j] <- 1)
  }
  means_alpha30[i] <- mean(probML)
  means_alpha30_M[i] <- mean(probM)
}

### Plots    
    
png(filename = paste(paste('alpha', real_tau, sep = '_'), 'png', sep = '.'), width=580, height=470)
plot(means_alpha_exp~corrects_alpha, type = 'l', col = 'red', lwd = 2, 
     main = 'Scenario with modifying the significance level', cex.main = 2,
     xlab = 'decrease in the level of ??',
     ylab = 'risk value', ylim = c(min(means_alpha30, means_alpha30_M,
                                       means_alpha5, means_alpha5_M,
                                       means_alpha_exp, means_alpha_exp_M),
                                  max(means_alpha30, means_alpha30_M,
                                      means_alpha5, means_alpha5_M,
                                      means_alpha_exp, means_alpha_exp_M)))
lines(means_alpha_exp_M~corrects_alpha, type = 'l', col = 'red', lwd = 2, lty = 4)
lines(means_alpha5~corrects_alpha, type = 'l', col = 'blue', lwd = 2)
lines(means_alpha5_M~corrects_alpha, type = 'l', col = 'blue', lwd = 2, lty = 4)
lines(means_alpha30~corrects_alpha, type = 'l', col = 'black', lwd = 2)
lines(means_alpha30_M~corrects_alpha, type = 'l', col = 'black', lwd = 2, lty = 4)
lines(rep(0.1, length(corrects_alpha))~corrects_alpha, type = 'l', col = 'red', lty = 3)
legend("topright", legend=c("Γ(1), ML", "Γ(1), M", 'Γ(5,1), ML', 'Γ(5,1), M', 'Γ(30,1), ML', 'Γ(30,1), M'),
       col=c("red", "red", "blue", 'blue', 'black', 'black'), lwd=rep(2,6), lty = c(rep(c(1,4),3)), cex=0.8)
dev.off()


png(filename = paste(paste('tau', real_tau, sep = '_'), 'png', sep = '.'), width=580, height=470)
plot(means_tau_exp~corrects, type = 'l', col = 'red', lwd = 2, 
     main = 'Scenario with modifying the fill rate', cex.main = 2,
     xlab = 'increase in the level of ??',
     ylab = 'risk value', ylim = c(min(means_tau30, means_tau30_M,
                                       means_tau5, means_tau5_M,
                                       means_tau_exp, means_tau_exp_M),
                                   max(means_tau30, means_tau30_M,
                                       means_tau5, means_tau5_M,
                                       means_tau_exp, means_tau_exp_M)))
lines(means_tau_exp_M~corrects, type = 'l', col = 'red', lwd = 2, lty = 4)
lines(means_tau5~corrects, type = 'l', col = 'blue', lwd = 2)
lines(means_tau5_M~corrects, type = 'l', col = 'blue', lwd = 2, lty = 4)
lines(means_tau30~corrects, type = 'l', col = 'black', lwd = 2)
lines(means_tau30_M~corrects, type = 'l', col = 'black', lwd = 2, lty = 4)
lines(rep(0.1, length(corrects))~corrects, type = 'l', col = 'red', lty = 3)
legend("topright", legend=c("Γ(1), ML", "Γ(1), M", 'Γ(5,1), ML', 'Γ(5,1), M', 'Γ(30,1), ML', 'Γ(30,1), M'),
       col=c("red", "red", "blue", 'blue', 'black', 'black'), lwd=rep(2,6), lty = c(rep(c(1,4),3)), cex=0.8)
dev.off()


supplies_tau <- rbind(colMeans(SS_tau_exp), colMeans(SS_tau_exp_M), colMeans(SS_tau5), 
                      colMeans(SS_tau5_M), colMeans(SS_tau30), colMeans(SS_tau30_M))
colnames(supplies_tau) <- corrects
rownames(supplies_tau) <- c('exp', 'exp_M', 'gamma5', 'gamma5_M', 'gamma30', 'gamma30_M')

supplies_alpha <- rbind(colMeans(SS_alpha_exp), colMeans(SS_alpha_exp_M), colMeans(SS_alpha5), 
                        colMeans(SS_alpha5_M), colMeans(SS_alpha30), colMeans(SS_alpha30_M))

colnames(supplies_alpha) <- corrects_alpha
rownames(supplies_alpha) <- c('exp', 'exp_M', 'gamma5', 'gamma5_M', 'gamma30', 'gamma30_M')
typ <- c('exp', 0, 'gamma5', 0, 'gamma30')
lims <- c(1, 0, 5, 0, 25)
cols_tau <- c(8,0,7,0,6)

for(i in c(1,3,5)){
  
  png(filename = paste(paste('barplot', 'alpha', tau, typ[i], sep = '_'), 'png', sep = '.'), width=680, height=480)
  p <- barplot(supplies_alpha[i:(i+1),], beside = T, ylim = c(lims[i], max(supplies_alpha[i:(i+1),])+1.5), xpd = F,
          xlab = 'decrease in the level of ??')
  text(x = p[,7], y = rep(supplies_alpha[i+1,7] + 0.7, 2), 
       labels = c(round(supplies_alpha[i,7],4), round(supplies_alpha[i+1,7],4)), srt = 90)
  dev.off()
  
  png(filename = paste(paste('barplot', 'tau', tau, typ[i], sep = '_'), 'png', sep = '.'), width=680, height=480)
  p <- barplot(supplies_tau[i:(i+1),], beside = T, ylim = c(lims[i], max(supplies_tau[i:(i+1),])+1.3), xpd = F,
          xlab = 'increase in the level of ??')
  text(x = p[,cols_tau[i]], y = rep(supplies_tau[i+1,cols_tau[i]] + 0.7, 2), 
       labels = c(round(supplies_tau[i,cols_tau[i]],4), round(supplies_tau[i+1,cols_tau[i]],4)), srt = 90)
  dev.off()
  
}