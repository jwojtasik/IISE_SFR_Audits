library(expint); library(dplyr); library(latex2exp);library(Brobdingnag)

ev_gamma <- function(S, shape = k, theta = 1, n = N) return((as.numeric(as.brob(S/theta)^(n*shape)*gammainc(1-(n-1)*shape,S/theta))+
                                                               (n*shape-1)*gammainc(shape+1,S/theta)-n*shape*(S/theta)*
                                                               gammainc(shape,S/theta))/(shape*(n*shape-1)*gamma(shape))-1+tau)
ev_gamma_n1 <- function(S, shape = k, theta = 1) return(1/gamma(shape)*gammainc(shape,S/theta)-S/(theta*gamma(shape))*gammainc(shape-1,S/theta)-1+tau)



k = 0.125
taus = c(0.9, 0.95, 0.98, 0.99)
wyniki_ML = matrix(NA, nrow = 25, ncol = length(taus))

### ML-estimation

Ns = c(1,2,5,10,25)
for(N in 1:5){
  for(t in 1:length(taus)){
    tau = taus[t]
    if(Ns[N]*k !=1){
    wyniki_ML[Ns[N], t] <- uniroot(ev_gamma, c(.1,100), n=Ns[N])$root
    }else{
    wyniki_ML[Ns[N], t] <- uniroot(ev_gamma_n1, c(.1,100))$root  
    }
  }
}

#variance
est_k <- function(x) return(log(x) - digamma(x) - log(mean(y_n)) + mean(log(y_n)))

th=1
wyniki_var_ML = matrix(NA, nrow = 25, ncol = length(taus))
Ns = c(1,2,5,10,25)
for(N in 1:5){
  for(t in 1:length(taus)){
    tau = taus[t]
    if(Ns[N]*k !=1){
      S_ml <- uniroot(ev_gamma, c(0.1,100), n=Ns[N])$root
    }else{
      S_ml <- uniroot(ev_gamma_n1, c(0.1,100))$root  
    }
    if(k==1){
      wyniki_var_ML[Ns[N], t] <- S_ml^2/250
      
    }else{
    h = 1e-9
    gS <- (1/h)*(ev_gamma(S_ml+h, shape = k, theta = th, n=Ns[N]) - ev_gamma(S_ml, shape = k, theta = th, n=Ns[N]))
    gk <- (1/h)*(ev_gamma(S_ml, shape = k+h, theta = th, n=Ns[N]) - ev_gamma(S_ml, shape = k, theta = th, n=Ns[N]))
    gth <- (1/h)*(ev_gamma(S_ml, shape = k, theta = th+h, n=Ns[N]) - ev_gamma(S_ml, shape = k, theta = th, n=Ns[N]))
    
    wyniki_var_ML[Ns[N], t] <-  (k*gk^2-2*th*gk*gth+gth^2*th^2*trigamma(k))/(250*(k*trigamma(k)-1)*(gS^2))
    }
  }
}

### M-estimation

result <- function(y, dimension) {
  wynik <- matrix(0, nrow = length(y) - dimension + 1, ncol = dimension)
  for (i in 1:(length(y)-dimension + 1)){
    for (j in 1:dimension){
      t <- i+j-1
      wynik[i,j] <- y[t]}
  }
  return(wynik)
}


psi <- function(S, m_est = m_est_y) {
  
  sum_up <- matrix(0, nrow = nrow(m_est), ncol = ncol(m_est))
  for(i in 1:ncol(m_est)){
    sum_up[,i] <- (m_est[,i] - S)^2*1*(m_est[,i]>=S)
  }
  
  horizon <- rep(NA, nrow(m_est))
  for(i in 1:nrow(m_est)){
  horizon[i] <- sum(sum_up[i,])/sum(m_est[i,]) + 2*(1-tau)*(S-mean(m_est[i,]))  
  }
return(mean(horizon))  
}

wyniki = matrix(NA, nrow = 25, ncol = 4)
wyniki_var = matrix(NA, nrow = 25, ncol = 4)
taus = c(0.9, 0.95, 0.98, 0.99)
Ns = c(1,2,5,10,25)
for(N in 1:5){
  for(t in 1:4){
    tau = taus[t]
    Ss = rep(NA, 2500)
      for(sim in 1:10000){
        y <- rgamma(250, shape = .125, scale = 1)
        m_est_y <- result(y, dimension = Ns[N]) 
        Ss[sim] <- optimize(psi, interval = c(0, 100), m_est = m_est_y)$minimum 
      }
    wyniki[Ns[N], t] <- mean(Ss)
    wyniki_var[Ns[N], t] <- var(Ss)
  }
}

### Plot

par(mar=c(4.2,4.2,2,2), xpd=TRUE)
plot(sqrt(wyniki_var_ML[1,]/250)~wyniki_ML[1,], xlab = 'basestock level', 
     ylab = 'standard error of the estimator', 
     xlim = c(min(wyniki[1,]),1.05*max(wyniki[25,])),
     ylim = c(min(sqrt(wyniki_var[1,]/250)),max(sqrt(wyniki_var[25,]/250))), pch = 1, cex.lab = 1.3)
legend("top", inset=c(0, 0.0025), ncol = 5, legend = c("n=1", "n=2", "n=5", "n=10", "n=25"), 
       pch = c(1, 3, 2, 4, 0), cex = 1.3)
lines(sqrt(wyniki_var_ML[,1]/250)%>% na.omit()~wyniki_ML[,1]%>% na.omit())
points(sqrt(wyniki_var_ML[2,]/250)~wyniki_ML[2,], pch = 3)#, col = 'red')
lines(sqrt(wyniki_var_ML[,2]/250)%>% na.omit()~wyniki_ML[,2]%>% na.omit())
points(sqrt(wyniki_var_ML[5,]/250)~wyniki_ML[5,], pch = 2)#, col = 'blue')
lines(sqrt(wyniki_var_ML[,3]/250)%>% na.omit()~wyniki_ML[,3]%>% na.omit())
points(sqrt(wyniki_var_ML[10,]/250)~wyniki_ML[10,], pch = 4)#, col = 'blue')
lines(sqrt(wyniki_var_ML[,4]/250)%>% na.omit()~wyniki_ML[,4]%>% na.omit())
points(sqrt(wyniki_var_ML[25,]/250)~wyniki_ML[25,], pch = 0)#, col = 'green')

points(sqrt(wyniki_var[1,]/250)~wyniki[1,], col = 'red')
lines(sqrt(wyniki_var[,1]/250)%>% na.omit()~wyniki[,1]%>% na.omit(), col = 'red')
points(sqrt(wyniki_var[2,]/250)~wyniki[2,], pch = 3, col = 'red')
lines(sqrt(wyniki_var[,2]/250)%>% na.omit()~wyniki[,2]%>% na.omit(), col = 'red')
points(sqrt(wyniki_var[5,]/250)~wyniki[5,], pch = 2, col = 'red')
lines(sqrt(wyniki_var[,3]/250)%>% na.omit()~wyniki[,3]%>% na.omit(), col = 'red')
points(sqrt(wyniki_var[10,]/250)~wyniki[10,], pch = 4, col = 'red')
lines(sqrt(wyniki_var[,4]/250)%>% na.omit()~wyniki[,4]%>% na.omit(), col = 'red')
points(sqrt(wyniki_var[25,]/250)~wyniki[25,], pch = 0, col = 'red')

text(0.6, 0.012, TeX(sprintf(paste(r'(\tau)', "= 0.9"))))
text(1.3, 0.025, TeX(sprintf(paste(r'(\tau)', "= 0.95"))))
text(1.8, 0.028, TeX(sprintf(paste(r'(\tau)', "= 0.98"))))
text(2.2, 0.04, TeX(sprintf(paste(r'(\tau)', "= 0.99"))))
