library(expint); library(dplyr); library(latex2exp);library(Brobdingnag)

ev_gamma <- function(S, shape = k, theta = 1, n = N) return(as.numeric(((S/theta)^(n*shape)*as.brob(gammainc(1-(n-1)*shape,S/theta))+
                                                               (n*shape-1)*as.brob(gammainc(shape+1,S/theta))-n*shape*(S/theta)*
                                                               as.brob(gammainc(shape,S/theta)))/as.brob(shape*(n*shape-1)*gamma(shape)))-1+tau)



wyniki = matrix(NA, nrow = 25, ncol = 4)

k = 30
th=1
taus = c(0.9, 0.95, 0.98, 0.99)
Ns = c(1,2,5,10,25)
for(N in 1:5){
  for(t in 1:4){
      tau = taus[t]
      wyniki[Ns[N], t] <- uniroot(ev_gamma, c(0.01,36), n=Ns[N])$root
  }
}


png(filename = 'relacja_S_n_duzo.png', width=680, height=480)
par(mar=c(4.2, 4.2, 5, 2), xpd=TRUE)
plot(wyniki[,1], ylim = c(5, 10), ylab = "basestock level", xlab = "review period", cex.lab = 1.3)#, xaxt="n")
#axis(1, at = c(1:5), labels = c(1,2,5,10,25))
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
dev.off()


#variance
est_k <- function(x) return(log(x) - digamma(x) - log(mean(y_n)) + mean(log(y_n)))


wyniki_var = matrix(NA, nrow = 25, ncol = 4)
taus = c(0.9, 0.95, 0.98, 0.99)
Ns = c(1,2,5,10,25)
for(N in 1:5){
  for(t in 1:4){
    tau = taus[t]
    #Ss = rep(NA, 250)
    #for(sim in 1:250){
    #y_n <- rgamma(250, shape = 5, scale = 1)
    #ke <- uniroot(est_k, c(0.01, 100))$root
    #te = mean(y_n)/ke
    S_ml <- uniroot(ev_gamma, c(0.1,12), shape = k, theta = 1, n=Ns[N])$root
    h = 1e-12
    gS <- (1/h)*(ev_gamma(S_ml+h, shape = k, theta = th, n=Ns[N]) - ev_gamma(S_ml, shape = k, theta = th, n=Ns[N]))
    gk <- (1/h)*(ev_gamma(S_ml, shape = k+h, theta = th, n=Ns[N]) - ev_gamma(S_ml, shape = k, theta = th, n=Ns[N]))
    gth <- (1/h)*(ev_gamma(S_ml, shape = k, theta = th+h, n=Ns[N]) - ev_gamma(S_ml, shape = k, theta = th, n=Ns[N]))
    #Ss[sim] <-  (k*gk^2-2*th*gk*gth+gth^2*te^2*trigamma(k))/(length(y_n)*(k*trigamma(k)-1)*(gS^2))
      #}
    wyniki_var[Ns[N], t] <-(k*gk^2-2*th*gk*gth+gth^2*te^2*trigamma(k))/(250*(k*trigamma(k)-1)*(gS^2)) #mean(Ss)
  }
}


par(mar=c(4.2, 5.5, 5, 2), xpd=TRUE)
plot(wyniki_var[,1], ylim = c(0,.18),  ylab = "variance of the\nbasestock estimator", xlab = "review period", cex.lab = 1.3)
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


