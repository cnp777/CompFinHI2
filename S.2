
# Import functions
source("Functions_HandIn2.R")

# Define parameters to be used in all subproblems 
K = 1
T = 1
sigma = 0.2
n = 10000
p = 7

# Simulation
set.seed(1)
samples_BS <- Generate_samples_BS(n,T,K,sigma)

# REGRESSION
X_BS = matrix(NA, nrow = n, ncol = (p+1) ) # to be used in price-only and half/half regression
for (j in 1:n) {
  for (i in 0:p) {
    X_BS[j,(i+1)] = (samples_BS$S0_BS[j])^(i)
  }
}
Y_BS = matrix(NA, nrow = n, ncol = (p+1) ) # to be used in half/half regression
for (j in 1:n) {
  for (i in 0:p) {
    Y_BS[j,i+1] = i*(samples_BS$S0_BS[j])^(i-1)
  }
}
YY_BS = matrix(NA, nrow = n, ncol = p) # to be used in delta-only regression
for (j in 1:n) {
  for (i in 1:p) {
    YY_BS[j,i] = i*(samples_BS$S0_BS[j])^(i-1)
  }
}

# Calculate theta_hat
anyNA(samples_BS$C_BS);anyNA(X_BS);anyNA(Y_BS);anyNA(YY_BS);anyNA(samples_BS$D_BS_Path)
theta_BS_price = solve(t(X_BS)%*%X_BS)%*%t(X_BS)%*%(samples_BS$C_BS)
theta_BS_half = solve(0.5*t(X_BS)%*%X_BS+0.5*t(Y_BS)%*%Y_BS)%*%(0.5*t(X_BS)%*%(samples_BS$C_BS)+0.5*t(Y_BS)%*%(samples_BS$D_BS_Path))
theta_BS_delta = solve(t(YY_BS)%*%YY_BS)%*%t(YY_BS)%*%(samples_BS$D_BS_Path)

# Plots
library(ggplot2)
df_BS = data.frame(samples_BS)

# FIGURE 2.1
png("S2Figure21.png")
ggplot(df_BS, aes(x=S0_BS,y=C_BS)) + 
  geom_point(color = 'grey', shape = 1) +
  geom_function(fun=F_BS, colour = "black") + 
  geom_function(fun=f_theta, args = list(theta_hat = theta_BS_price, p = p), colour = "red") +
  theme_classic() + scale_x_continuous(expand = c(0, 0), limits = c(0.22, 1.8), sec.axis = dup_axis(breaks=NULL,name=NULL)) + 
  xlab('Stock price') + ylab('Call value') + scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL))
dev.off()

# FIGURE 2.2 
png("S2Figure22.png")
ggplot(df_BS, aes(x=S0_BS)) + geom_function(fun=trueDelta_BS, colour = "black") + 
  geom_function(fun=delta_regression, args = list(theta_hat=theta_BS_price, p = p), colour = "red",) + 
  theme_classic() + scale_x_continuous(expand = c(0, 0), limits = c(0.22, 1.8), sec.axis = dup_axis(breaks=NULL,name=NULL)) + 
  xlab('Stock price') + ylab('Call Delta') + scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL))
dev.off()

# FIGURE 3 
png("S2Figure3.png")
ggplot(df_BS, aes(x=S0_BS, y=D_BS_Path)) + geom_point(color = 'grey', shape = 1) + 
  geom_function(fun=trueDelta_BS, colour = "black",xlim = c(0.5,1.7)) +
  geom_function(fun=delta_regression, args = list(theta_hat=theta_BS_price, p=p), colour = "blue",xlim = c(0.5,1.7)) +
  geom_function(fun=delta_regression, args = list(theta_hat=theta_BS_half, p=p), colour = "red",xlim = c(0.5,1.7)) + 
  geom_function(fun=delta_regression_deltaonly, args = list(theta_hat=theta_BS_delta, p=p), colour = "orange",xlim = c(0.5,1.7)) +
  xlab('Stock price') + ylab('Call Delta') + theme_classic() + 
  ggtitle("Estimated Call Deltas in the Black-Scholes Model") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  annotate(geom="text", x=0.3, y=1.1, label= "Grey o's: Simulated Deltas", color="grey", hjust = 0) +
  annotate(geom="text", x=0.3, y=0.9, label= "Black: True Delta", color="black", hjust = 0) +
  annotate(geom="text", x=0.3, y=0.85, label= "Blue: Price only-regression", color="blue", hjust = 0) +
  annotate(geom="text", x=0.3, y=0.8, label= "Red: Half/half-regression", color="red", hjust = 0) +
  annotate(geom="text", x=0.3, y=0.75, label= "Orange: Delta only-regression", color="orange", hjust = 0) +
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL), limits = c(-0.1,1.2)) +
  scale_x_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL), limits = c(0.3,1.8)) 
dev.off()

# FIGURE 4

rm(list=ls())

simulations = c(500, 1000, 2000, 3000, 4000, 5000, 6000)
polydeg <- c(3,4,5,6,7,8)
HedgeError_half   <- matrix(NA, length(simulations), length(polydeg))
HedgeError_price  <- matrix(NA, length(simulations), length(polydeg))

for (q in 1:length(polydeg)){
  for (sim in 1:length(simulations)){
    sigma = 0.2
    dt = 1/52
    d = 1
    w = 1/2
    K = 1
    T = 1
    
    p = polydeg[q]
    n = simulations[sim]
    
    S0 <- rep(NA, n)
    ST <- rep(NA, n)
    C  <- rep(NA, n)
    D  <- rep(NA, n)
    
    set.seed(1)
    N1 <- rnorm(10000, 0, 1)
    set.seed(2)
    N2 <- rnorm(10000,0,1)
    
    thetaHat_half  <- matrix(NA, 52, p+1)
    thetaHat_price <- matrix(NA, 52, p+1)
    
    for (t in 0:51){
      
      for (j in 1:n){
        S0[j] <- K*exp(-0.5*sigma^2*(T) + sigma*sqrt(T) * N1[j]) 
        ST[j] <- S0[j]*exp(-0.5*sigma^2*(T-t*dt) + sigma*sqrt(T-t*dt) * N2[j])
        C[j] <- max(ST[j]-K,0)
        
        if (ST[j]>= K) {D[j]=ST[j]/S0[j]}else{D[j]=0}
      }
      
      
      X <- matrix(NA, n, p+1)
      Y <- matrix(NA, n, p+1)
      
      for (i in 0:p){
        for (j in 1:n){
          X[j,i+1] = (S0[j])^(i)
          Y[j,i+1] = i*(S0[j])^(i-1)
        }
      }
      
      thetaHat_half[t+1,]  <- solve(w*t(X)%*%X+(1-w)*t(Y)%*%Y)%*%(w*t(X)%*%C+(1-w)*t(Y)%*%D)
      thetaHat_price[t+1,] <- solve(t(X)%*%X)%*%t(X)%*%C
      
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
    
    PnL_half <- rep(NA,1000)
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
      
      PnL_half[i]  <- VT_h - max(S_T-K,0)
      PnL_price[i] <- VT_p - max(S_T-K,0)
      
    }
    
    HedgeError_half[sim, q]  <- sd(PnL_half)/Ftheta_half(1,p,1)   
    HedgeError_price[sim, q] <- sd(PnL_price)/Ftheta_price(1,p,1)  
  }
}

# Hedging with true Delta
F <- function(x,t){
  
  d1 = (log(x/K) + 1/2*sigma^2*(T-t))/(sigma*sqrt(T-t))
  d2 = d1 - sigma*sqrt(T-t)
  
  return( x*pnorm(d1) - K*pnorm(d2) )
}

TrueDelta <- function(x,t){
  d1 = (log(x/K) + 1/2*sigma^2*(T-t))/(sigma*sqrt(T-t))
  return ( pnorm(d1) )
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
    PnL_true[q] <- VT_TD - max(ST_TD - K,0)
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
  ggtitle("Polynomial call hedging") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  annotate(geom="text", x=2000, y=0.4, label= "Orange curves: Price only-regression", color="orange", hjust = 0) +
  annotate(geom="text", x=2000, y=0.38, label= "Blue curves: (Price, Delta)-regression", color="blue", hjust = 0) + 
  annotate(geom="text", x=0, y=0.131, label= "Hedging with true Delta", color="red", hjust = 0, size = 3) +
  annotate(geom="text", x=0, y=0, label= "K = 1, T = 1, vol = 0.2, dt = 1/52, #repetitions = 1000", color="black", hjust = 0, size = 3) +
  annotate(geom="text", x=6100, y=0.3, label= "Poly'deg'", color="blue", hjust = 0, size = 3) + 
  annotate(geom="text", x=6100, y=df_half[7,1], label= "3", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=df_half[7,2], label= "4", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=0.2, label= "5, 6", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=df_half[7,5], label= "7", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=df_half[7,6], label= "8", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,1], label= "3", color="orange", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,2], label= "4", color="orange", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,3], label= "5, 6", color="orange", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,5], label= "7", color="orange", hjust = 0, size = 4) +
  annotate(geom="text", x=6200, y=df_price[7,6], label= "8", color="orange", hjust = 0, size = 4) +
  scale_x_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL),breaks=c(0,1000,2000,3000,4000,5000,6000), limits = c(0,6500)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL),breaks=c(0,0.1,0.2,0.3,0.4), limits = c(0,0.4))

# ggsave("S2Figure4.png")








