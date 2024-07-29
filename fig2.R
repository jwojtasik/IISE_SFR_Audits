library(expint); library(dplyr); library(latex2exp)

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


par(mar=c(4.2, 4.2, 5, 2), xpd=TRUE)
plot(wyniki[,1], ylim = c(5, 10), ylab = "basestock level", xlab = "review period", cex.lab = 1.3)
lines(wyniki[,1] %>% na.omit()~Ns)
points(wyniki[,2], pch = 3)
lines(wyniki[,2]%>% na.omit()~Ns)
points(wyniki[,3], pch = 2)
lines(wyniki[,3]%>% na.omit()~Ns)
points(wyniki[,4], pch = 0)
lines(wyniki[,4]%>% na.omit()~Ns)
legend("top", inset=c(0, -.3), ncol = 4, legend = c(TeX(sprintf(paste(r'(\tau)', "= 0.9"))), 
                                                    TeX(sprintf(paste(r'(\tau)', "= 0.95 "))), 
                                                    TeX(sprintf(paste(r'(\tau)', "= 0.98 "))), 
                                                    TeX(sprintf(paste(r'(\tau)', "= 0.99 ")))), 
       pch = c(1, 3, 2, 0), cex = 1.3)


#variance

par(mar=c(4.2, 5.5, 5, 2), xpd=TRUE)
plot(wyniki_var[,1], ylim = c(0,.5),  ylab = "variance of the\nbasestock estimator", xlab = "review period", cex.lab = 1.3)
lines(wyniki_var[,1] %>% na.omit()~Ns)
points(wyniki_var[,2], pch = 3)
lines(wyniki_var[,2]%>% na.omit()~Ns)
points(wyniki_var[,3], pch = 2)
lines(wyniki_var[,3]%>% na.omit()~Ns)
points(wyniki_var[,4], pch = 0)
lines(wyniki_var[,4]%>% na.omit()~Ns)
legend("top", inset=c(0, -.3), ncol = 4, legend = c(TeX(sprintf(paste(r'(\tau)', "= 0.9"))), 
                                                    TeX(sprintf(paste(r'(\tau)', "= 0.95 "))), 
                                                    TeX(sprintf(paste(r'(\tau)', "= 0.98 "))), 
                                                    TeX(sprintf(paste(r'(\tau)', "= 0.99 ")))), 
       pch = c(1, 3, 2, 0), cex = 1.3)