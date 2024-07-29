library(rootSolve)
library(expint)

est_k <- function(x, dane) return(log(x) - digamma(x) - log(mean(dane)) + mean(log(dane)))

#exponential distibution

f_exp <- function(x, lambda = lambda_est) return (exp(-lambda*x) - lambda*x*gammainc(0,lambda*x) - 1 + tau)
var_exp <- function(S, ya) return(S^2/length(ya))

#gamma distribution

f_gamma <- function(S, k = k_est, theta = theta_est) return(1/gamma(k)*gammainc(k,S/theta)-S/(theta*gamma(k))*gammainc(k-1,S/theta)-1+tau)

var_gamma <- function(S, theta, k, ya){
  h <- 0.0000000001
  gS <- (1/h)*(1/gamma(k)*gammainc(k,(S+h)/theta)-(S+h)/(theta*gamma(k))*gammainc(k-1,(S+h)/theta) - 
                 (1/gamma(k)*gammainc(k,S/theta)-S/(theta*gamma(k))*gammainc(k-1,S/theta)))
  gk <- (1/h)*(1/gamma(k+h)*gammainc(k+h,S/theta)-S/(theta*gamma(k+h))*gammainc(k+h-1,S/theta)- 
                 (1/gamma(k)*gammainc(k,S/theta)-S/(theta*gamma(k))*gammainc(k-1,S/theta)))
  gth <- (1/h)*(1/gamma(k)*gammainc(k,S/(theta+h))-S/((theta+h)*gamma(k))*gammainc(k-1,S/(theta+h)) - 
                  (1/gamma(k)*gammainc(k,S/theta)-S/(theta*gamma(k))*gammainc(k-1,S/theta)))
  var <-  (k*gk^2-2*theta*gk*gth+(gth^2)*(theta^2)*trigamma(k))/(length(ya)*(k*trigamma(k)-1)*(gS^2))
  return(var)
}

#M-estimation

loss <- function (S, y) return (((y-S)^2/y)*1*(y>=S) + 2*(1-tau)*(S-y))
psi <- function (S, x=ya) return (1/length(x) * sum(((x-S)^2/x)*1*(x>=S) + 2*(1-tau)*(S-x)))
m <- function (S, y) return (2*(1-tau) - 2*1*(y>S)*((y-S)/y))


#variables and dimension parameters
lambda = 1
theta = 1/lambda
k = 5

#samples parameters
T = 100 #learning sample

#simulation parameters
alpha = 0.2
sym = 10000
tau = 0.9

odp <- matrix(0, ncol = 3, nrow = 10)
wynik_ML <- rep(0, sym)
wynik_M <- rep(0, sym)

#modifying simulation size 

for (i in 1:10) {
  ss = i*T
  odp[i,1] <- ss
  for (j in 1:sym) {
    y <- rgamma(T + ss, shape = k, scale = theta)
    ya <- y[1:T]
    p1 <- T+1; p2 <- ss+T
    yp <- y[p1:p2]
    if(k==1) {
        lambda_est = 1/mean(ya)
        S <- uniroot(f_exp, c(0.01,1000))$root
        var = S^2/length(ya)
        S_m <- optimize(psi, interval = c(0, 100))$minimum
        var_M = ((1-tau)^2 + 2*tau*exp(-S_m/theta)-exp(-S_m/theta)-2*tau*S_m/theta*gammainc(0, S_m/theta)+S_m*
                   exp(-S_m/theta)-S_m^2*gammainc(0,S_m/theta))/(T*(gammainc(0, S_m/theta))^2)
    
    }else {
      k_est <- uniroot(est_k, c(0.001, 100), dane = ya)$root
      theta_est = mean(ya)/k_est
      S <- uniroot(f_gamma, c(0.1,100))$root
      var <- var_gamma(S, theta_est, k_est, ya)
      S_m <- optimize(psi, interval = c(0, 100))$minimum
      var_M <- (gamma(k)*theta^2*((S_m/theta)^2*gammainc(k-2, S_m/theta)-(2*S_m*tau)/(theta)*gammainc(k-1, S_m/theta)+
                gamma(k)*(1-2*tau)*(1-gammainc(k, S_m/theta)/gamma(k))+gamma(k)*tau^2))/(T*(gammainc(k-1, S_m/theta))^2)
    }
    SME <- S+qnorm(1-alpha/2)*(var^0.5)
    SM <- S_m+qnorm(1-alpha/2)*(var_M^0.5)
    wynik_ML[j]<- mean(pmax(yp-SME,0)/yp)
    wynik_M[j]<- mean(pmax(yp-SM,0)/yp)  
  }
  odp[i,2] <- mean(wynik_ML > 1-tau)
  odp[i,3] <- mean(wynik_M > 1-tau)
}

test <- odp
test[,2:3] <- test[,2:3]/100

plot(odp[,1], odp[,3], type = 'l', lwd = 2,
     xlab = 'size of the evaluation sample', 
     ylim = c(.115, max(odp[,2:3])), 
     ylab = 'risk value', col = 'red', 
     main  = paste('P(FR <', paste0(tau,')'), 'based on estimation sample of size T =', T, sep = ' '))
lines(odp[,1], odp[,2], type = 'l', lwd = 2, col = "black")
