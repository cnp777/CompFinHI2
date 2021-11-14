### Finding the optimal epsilon
n = 10000
S0 = K = 1 
T = 1  #Note: use same code but T = 1/2 for Figure 12
sigma = 0.2

eps = seq(0.0002,1,by=0.001)
var = rep(NA, length(eps))

for (i in 1:length(eps)) {
epsilon = eps[i]
set.seed(1)
Z = rnorm(n,0,1)
ST_BS = rep(NA,n)
Digital_delta_LRM = rep(NA,n)
Control1 = rep(NA,n)
Control2 = rep(NA,n)
Digital_delta_mixed = rep(NA,n)

for (j in 1:n)
{
  ST_BS[j] = S0*exp( -0.5*sigma^2*T + sigma*sqrt(T)*Z[j])
  if(ST_BS[j] > K){Digital_delta_LRM[j] = Z[j]/(S0*sigma*sqrt(T)) }else(Digital_delta_LRM[j] = 0) 
  if( abs(ST_BS[j]-K) < epsilon){Control1[j] = (1/(2*epsilon))*ST_BS[j]/S0}else(Control1[j] = 0)
  Control2[j] = min(1,max(0, (ST_BS[j]-K + epsilon)/(2*epsilon)  ))*(Z[j]/(S0*sigma*sqrt(T)))
  Digital_delta_mixed[j] =  Digital_delta_LRM[j] + Control1[j] - Control2[j]
}
var[i] = var(Digital_delta_mixed)/n
}
library(latex2exp)
library(ggplot2)
df_e = data.frame( index = eps, var = var )
ggplot(df_e, aes(x = index, y = var)) + geom_line(color = 'blue') + theme_classic() + 
  theme(panel.grid.major = element_line()) + xlab('Epsilon') + ylab('Variance (per replication)') +
  scale_x_continuous(sec.axis = dup_axis(breaks=c(0.0,0.2,0.4,0.6,0.8,1),labels=NULL,name=NULL), breaks=c(0.0,0.2,0.4,0.6,0.8,1), limits = c(0,1)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks=c(0.0005,0.001,0.0015),labels=NULL, name=NULL), 
                     limits = c(0,0.0015),breaks=c(0.0005,0.001,0.0015), labels = c("0.5","1.0","1.5")) +
  theme(axis.ticks.length=unit(-0.15, "cm"),plot.title = element_text(size = 8)) +
  ggtitle(TeX("$\\times 10^{-3}$"))
ggsave("S5eps.png")

Variance = function(epsilon)
{
  set.seed(1)
  Z = rnorm(n,0,1)
  ST_BS = rep(NA,n)
  Digital_delta_LRM = rep(NA,n)
  Control1 = rep(NA,n)
  Control2 = rep(NA,n)
  Digital_delta_mixed = rep(NA,n)
  
  for (j in 1:n)
  {
    ST_BS[j] = S0*exp( -0.5*sigma^2*T + sigma*sqrt(T)*Z[j])
    if(ST_BS[j] > K){Digital_delta_LRM[j] = Z[j]/(S0*sigma*sqrt(T)) }else(Digital_delta_LRM[j] = 0) 
    if( abs(ST_BS[j]-K) < epsilon){Control1[j] = (1/(2*epsilon))*ST_BS[j]/S0}else(Control1[j] = 0)
    Control2[j] = min(1,max(0, (ST_BS[j]-K + epsilon)/(2*epsilon)  ))*(Z[j]/(S0*sigma*sqrt(T)))
    Digital_delta_mixed[j] =  Digital_delta_LRM[j] + Control1[j] - Control2[j]
  }
  
  return(var(Digital_delta_mixed)/n)
}
# Find minimal value
optimize(Variance, interval = c(0.002,10) , maximum = FALSE) # epsilon = 0.48


### PLOTS 
n = 10000
p = 7

set.seed(1)
samples_BS <- Generate_samples_BS(n,T,K,sigma)

# Beregn theta_hat
X_BS = matrix(NA, nrow = n, ncol = (p+1) );Y_BS = matrix(NA, nrow = n, ncol = (p+1) ) 
for (j in 1:n) {
  for (i in 0:p) {
    X_BS[j,(i+1)] = (samples_BS$S0_BS[j])^(i)
    Y_BS[j,i+1] = i*(samples_BS$S0_BS[j])^(i-1)
  }
}
YY_BS = matrix(NA, nrow = n, ncol = p) # to be used in delta-only regression
for (j in 1:n) {
  for (i in 1:p) {
    YY_BS[j,i] = i*(samples_BS$S0_BS[j])^(i-1)
  }
}

theta_hat_digital_price = solve(t(X_BS)%*%X_BS)%*%t(X_BS)%*%(samples_BS$DigitalT)
theta_hat_digital_half = solve(0.5*t(X_BS)%*%X_BS+0.5*t(Y_BS)%*%Y_BS)%*%(0.5*t(X_BS)%*%(samples_BS$DigitalT)+0.5*t(Y_BS)%*%(samples_BS$Digital_delta_mixed))
theta_hat_digital_delta = solve(t(YY_BS)%*%YY_BS)%*%t(YY_BS)%*%(samples_BS$Digital_delta_mixed)

library(ggplot2)
df_BS = data.frame(samples_BS)

# FIGURE 2.1 
ggplot(df_BS, aes(x=S0_BS,y=DigitalT)) + 
  geom_point(color = 'grey', shape = 1) +
  geom_function(fun=f_theta, colour = "red", args = list(theta_hat = theta_hat_digital_price,p = p)) +
  theme_classic() + xlab('Stock price') + ylab('Digital Call value') +
  geom_function(fun=F_Digital, colour = "black", args = list(t = 0)) +
  scale_x_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL), expand = c(0, 0), limits = c(0.42, 1.8)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL), limits = c(-0.1, 1.2))
ggsave("S5F21.png")

# FIGURE 2.2 
ggplot(df_BS, aes(x=S0_BS))  +
  geom_function(fun=delta_regression, args = list(theta_hat_digital_price, p) , colour = "red") + 
  geom_function(fun=trueDelta_Digital, args = list(t=0) , colour = "black") + 
  theme_classic() + ylab('Digital Call Delta') + xlab("Initial stock price") +
  scale_x_continuous(expand = c(0, 0), limits = c(0.42, 1.8), sec.axis = dup_axis(breaks=NULL, name=NULL)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL), limits = c(-0.2, 2))
ggsave("S5F22.png")

# FIGURE 3 
ggplot(df_BS, aes(x=S0_BS, y=Digital_delta_mixed)) + geom_point(color = 'grey', shape = 1) + 
  theme_classic() + ylab('Digital Call Delta') + xlab("Initial stock price") + 
  ggtitle("Estimated Digital Call Deltas in the Black-Scholes Model") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_function(fun=trueDelta_Digital, args = list(t=0) , colour = "black") + 
  geom_function(fun=delta_regression, args = list(theta_hat_digital_price, p) , colour = "blue", xlim = c(0.65,1.5)) +
  geom_function(fun=delta_regression, args = list(theta_hat=theta_hat_digital_half, p=p), colour = "red", xlim = c(0.65,1.5)) + 
  geom_function(fun=delta_regression_deltaonly, args = list(theta_hat=theta_hat_digital_delta, p=p), colour = "orange",xlim = c(0.65,1.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, 1.6), sec.axis = dup_axis(breaks=NULL, name=NULL)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL), limits = c(-0.2, 4)) +
  annotate(geom="text", x=0.55, y=2.5, label= "Grey o's: Simulated Deltas", color="grey", hjust = 0, size=3) +
  annotate(geom="text", x=0.55, y=2.4, label= "Black: True Delta", color="black", hjust = 0, size=3) +
  annotate(geom="text", x=0.55, y=2.35, label= "Blue: Price only-regression", color="blue", hjust = 0, size=3) +
  annotate(geom="text", x=0.55, y=2.3, label= "Red: Half/half-regression", color="red", hjust = 0, size=3) +
  annotate(geom="text", x=0.55, y=2.25, label= "Orange: Delta only-regression", color="orange", hjust = 0, size=3) 
ggsave("S5F3.png")

# FIGURE 4
rm(list=ls())

simulations = c(500, 1000, 2000, 3000, 4000, 5000, 6000)
polydeg <- c(3,4,5,6,7,8)
HedgeError_half   <- matrix(NA, length(simulations), length(polydeg))
HedgeError_price <- matrix(NA, length(simulations), length(polydeg))

for (q in 1:length(polydeg)){
  for (sim in 1:length(simulations)){
    sigma = 0.2
    dt = 1/52
    d = 1
    w = 1/2
    K = 1
    T = 1
    epsilon = 0.48
    
    p = polydeg[q]
    n = simulations[sim]
    
    S0 = rep(NA, n)
    ST = rep(NA, n)
    
    DigitalT = rep(NA, n)
    Control1 = rep(NA, n)  
    Control2 = rep(NA, n)  
    Digital_delta_LRM   = rep(NA, n)  
    Digital_delta_mixed = rep(NA, n)
    
    set.seed(1)
    N1 <- rnorm(10000, 0, 1)
    set.seed(2)
    N2 <- rnorm(10000,0,1)
    
    thetaHat_half  <- matrix(NA, 52,p+1)
    thetaHat_price <- matrix(NA, 52, p+1)
    
    for (t in 0:51){
      # Idea : Calculate variance minimizing epsilon here using optimize(). 
      
      for (j in 1:n){
        S0[j] <- K*exp(-0.5*sigma^2*(T) + sigma*sqrt(T) * N1[j]) 
        ST[j] <- S0[j]*exp(-0.5*sigma^2*(T-t*dt) + sigma*sqrt(T-t*dt) * N2[j])
        
        if(ST[j] > K){DigitalT[j]=1}else(DigitalT[j]=0)
        if(ST[j] > K){Digital_delta_LRM[j] = N2[j]/(S0[j]*sigma*sqrt(T-t*dt)) }else(Digital_delta_LRM[j] = 0) 
        if( abs(ST[j]-K) < epsilon){Control1[j] = (1/(2*epsilon))*ST[j]/S0[j]}else(Control1[j] = 0)
        Control2[j] = min(1,max(0, (ST[j] - K + epsilon)/(2*epsilon)))*(N2[j]/(S0[j]*sigma*sqrt(T-t*dt)))
        Digital_delta_mixed[j] =  Digital_delta_LRM[j] + Control1[j] - Control2[j]
        
         }
      
      X <- matrix(NA, n, p+1)
      Y <- matrix(NA, n, p+1)
      
      for (i in 0:p){
        for (j in 1:n){
          X[j,i+1] = (S0[j])^(i)
          Y[j,i+1] = i*(S0[j])^(i-1)
        }
      }
      
      thetaHat_half[t+1,]  <- solve(w*t(X)%*%X+(1-w)*t(Y)%*%Y)%*%(w*t(X)%*%DigitalT+(1-w)*t(Y)%*%Digital_delta_mixed)
      thetaHat_price[t+1,] <- solve(t(X)%*%X)%*%t(X)%*%DigitalT
      
    }
    
    Ftheta_half <- function(x,p,r){
      val <- 0
      for (i in 0:p){
        val = val+thetaHat_half[r, i+1]*x^(i)
      }
      return(val)
    }
    
    deltareg_half <- function(x,p,r){
      val <- 0
      for (i in 1:p){
        val = val+thetaHat_half[r, i+1]*i*x^(i-1)
      }
      return(val)
    }
    
    #Price only
    
    Ftheta_price <- function(x,p,r){
      val <- 0
      for (i in 0:p){
        val = val+thetaHat_price[r, i+1]*x^(i)
      }
      return(val)
    }
    
    deltareg_price <- function(x,p,r){
      val <- 0
      for (i in 1:p){
        val = val+thetaHat_price[r, i+1]*i*x^(i-1)
      }
      return(val)
    }
    
    PnL_half  <- rep(NA,1000)
    PnL_price <- rep(NA,1000)
    
    set.seed(1)
    
    for (i in 1:1000){
      
      St <- 1  #S(0)
      Z1 <- rnorm(52,0,1)
      
      Vt_h <- Ftheta_half(1,p,1)
      at_h <- deltareg_half(1,p,1)
      bt_h <- Vt_h-at_h*St
      
      Vt_p <- Ftheta_price(1,p,1)
      at_p <- deltareg_price(1,p,1)
      bt_p <- Vt_p - at_p*St
      
      for (t in 1:51){
        
        St = St*exp(-0.5*sigma^2*dt + sigma * sqrt(dt) * Z1[t])
        
        Vt_h = at_h * St + bt_h
        Vt_p = at_p * St + bt_p
        
        at_h = deltareg_half(St,p,t+1)
        at_p = deltareg_price(St,p,t+1)
        
        bt_h = Vt_h - at_h * St
        bt_p = Vt_p - at_p * St
        
      }
      
      S_T <- St*exp(-0.5*sigma^2*dt + sigma * sqrt(dt) * Z1[52]) 
      
      VT_h = at_h * S_T + bt_h
      VT_p = at_p * S_T + bt_p
      
      PnL_half[i]  <- VT_h - ifelse(S_T > K,1,0) 
      PnL_price[i] <- VT_p - ifelse(S_T > K,1,0) 
      
    }
    
    HedgeError_half[sim, q]  <- sd(PnL_half)/Ftheta_half(1,p,1)   
    HedgeError_price[sim, q] <- sd(PnL_price)/Ftheta_price(1,p,1)  
  }
}

### Black Scholes model: True Values of Digital Call Option
F <- function(x,t){
  
  d = (log(x/K) - 0.5*sigma^2*(T-t))/(sigma*sqrt(T-t))
  
  return( dnorm(d) )
}

TrueDelta <- function(x,t){
  
  d = (log(x/K) - 0.5*sigma^2*(T-t))/(sigma*sqrt(T-t))
  
  return ( dnorm(d)/(x*sigma*sqrt(T-t)) )
}

PnL_true <- rep(NA, 1000)
HedgeError_true <- rep(NA, length(simulations))

for (j in 1:length(simulations)){
  set.seed(1)
  for (q in 1:1000){
    S0 <- 1
    St_TD <- S0
    Vt_TD <- F(1,0)
    at_TD <- TrueDelta(1,0)
    bt_TD <- Vt_TD-at_TD*St_TD
    Z <- rnorm(52,0,1)
    
    for (t in 1:51){
      St_TD = St_TD * exp(-0.5*sigma^2*dt + sigma * sqrt(dt) * Z[t])
      Vt_TD = at_TD * St_TD + bt_TD
      at_TD = TrueDelta(St_TD,t*dt)
      bt_TD = Vt_TD - at_TD * St_TD
    }
    ST_TD <- St_TD * exp(-0.5*sigma^2*dt + sigma * sqrt(dt) * Z[52])
    VT_TD = at_TD * ST_TD + bt_TD
    PnL_true[q] <- VT_TD - ifelse(ST_TD > K,1,0)
  }
  
  HedgeError_true[j] <- sd(PnL_true)/F(1,0)
}

# plot
library(ggplot2)
library(reshape2)

df_half = data.frame(HedgeError_half, row.names = simulations)
colnames(df_half) <- polydeg
df_half$sim <- simulations

df_price = data.frame(HedgeError_price, row.names = simulations)
colnames(df_price) <- polydeg
df_price$sim <- simulations

ggplot(df_half, aes(sim, y=df_half[,1])) + geom_point(color = 'blue', shape = 1) + geom_line(color = 'blue') + 
  geom_point(data=df_half, aes(sim, y=df_half[,2]),color = 'blue', shape = 1) + geom_line(data=df_half, aes(sim, y=df_half[,2]),color = 'blue') +
  geom_point(data=df_half, aes(sim, y=df_half[,3]),color = 'blue', shape = 1) + geom_line(data=df_half, aes(sim, y=df_half[,3]),color = 'blue') +
  geom_point(data=df_half, aes(sim, y=df_half[,4]),color = 'blue', shape = 1) + geom_line(data=df_half, aes(sim, y=df_half[,4]),color = 'blue') +
  geom_point(data=df_half, aes(sim, y=df_half[,5]),color = 'blue', shape = 1) + geom_line(data=df_half, aes(sim, y=df_half[,5]),color = 'blue') +
  geom_point(data=df_half, aes(sim, y=df_half[,6]),color = 'blue', shape = 1) + geom_line(data=df_half, aes(sim, y=df_half[,6]),color = 'blue') +
  geom_point(data=df_price, aes(sim, y=df_price[,1]),color = 'orange', shape = 1) + geom_line(data=df_price, aes(sim, y=df_price[,1]),color = 'orange') +
  geom_point(data=df_price, aes(sim, y=df_price[,2]),color = 'orange', shape = 1) + geom_line(data=df_price, aes(sim, y=df_price[,2]),color = 'orange') +
  geom_point(data=df_price, aes(sim, y=df_price[,3]),color = 'orange', shape = 1) + geom_line(data=df_price, aes(sim, y=df_price[,3]),color = 'orange') +
  geom_point(data=df_price, aes(sim, y=df_price[,4]),color = 'orange', shape = 1) + geom_line(data=df_price, aes(sim, y=df_price[,4]),color = 'orange') +
  geom_point(data=df_price, aes(sim, y=df_price[,5]),color = 'orange', shape = 1) + geom_line(data=df_price, aes(sim, y=df_price[,5]),color = 'orange') +
  geom_point(data=df_price, aes(sim, y=df_price[,6]),color = 'orange', shape = 1) + geom_line(data=df_price, aes(sim, y=df_price[,6]),color = 'orange') +
  geom_hline(aes(yintercept=HedgeError_true[1]), linetype = "dashed", col = "red") + geom_vline(xintercept = 6000, col = "black") +
  xlab('#simulations') + ylab('Hedge error') + theme_classic() + 
  ggtitle("Polynomial digital call hedging") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  annotate(geom="text", x=2000, y=0.8, label= "Orange curves: Price only-regression", color="orange", hjust = 0) +
  annotate(geom="text", x=2000, y=0.75, label= "Blue curves: (Price, Delta)-regression", color="blue", hjust = 0) + 
  annotate(geom="text", x=0, y=0.4, label= "Hedging with true Delta", color="red", hjust = 0, size = 3) +
  annotate(geom="text", x=0, y=0, label= "K = 1, T = 1, vol = 0.2, dt = 1/52, #repetitions = 1000", color="black", hjust = 0, size = 3) +
  annotate(geom="text", x=6100, y=0.7, label= "Poly'deg'", color="blue", hjust = 0, size = 3) + 
  annotate(geom="text", x=6100, y=df_half[7,1], label= "3", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=df_half[7,2], label= "4,5", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=df_half[7,4], label= "6", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=df_half[7,6], label= "7, 8", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,1], label= "3", color="orange", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,3], label= "4, 5", color="orange", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,4], label= "6", color="orange", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,6], label= "7, 8", color="orange", hjust = 0, size = 4) +
  scale_x_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL),breaks=c(0,1000,2000,3000,4000,5000,6000), limits = c(0,6500)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL),breaks=c(0,0.2,0.4,0.6,0.8))

#ggsave("S5F4.png")
