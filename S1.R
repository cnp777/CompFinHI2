
# Import functions
source("Functions_HandIn2.R")

# Define parameters to be used in all subproblems 
K = 1
T = 1
sigma = 0.2
n = 10000
p = 7

######################
###      S.1      ####
######################
# --------------------

# Simulation
set.seed(1)
samples_Bachelier = Generate_samples_Bachelier(n, T, K, sigma)
df <- data.frame(S0 = samples_Bachelier$S0, C = samples_Bachelier$C, D = samples_Bachelier$D)

# Regression
X = matrix(NA, nrow = n, ncol = (p+1) ) # to be used in price only and half/half regression
for (j in 1:n) {
  for (i in 0:p) {
    X[j,(i+1)] = (df$S0[j])^(i)
  }
}
Y = matrix(NA, nrow = n, ncol = (p+1) ) # to be used in half/half regression
for (j in 1:n) {
  for (i in 0:p) {
    Y[j,i+1] = i*(df$S0[j])^(i-1)
  }
} 
YY = matrix(NA, nrow = n, ncol = p) # to be used in delta only regression
for (j in 1:n) {
  for (i in 1:p) {
    YY[j,i] = i*(df$S0[j])^(i-1)
  }
}
theta_hat_price = solve(t(X)%*%X)%*%t(X)%*%(df$C)
theta_hat_half = solve(0.5*t(X)%*%X+0.5*t(Y)%*%Y)%*%(0.5*t(X)%*%(df$C)+0.5*t(Y)%*%(df$D))
theta_hat_delta = solve(t(YY)%*%YY)%*%t(YY)%*%(df$D)

# Plots
library(ggplot2)

# FIGURE 2.1
png("S1Figure21.png")
ggplot(df, aes(x=S0,y=C, colour=Functions)) + 
  geom_point(color = 'grey', shape = 1) +
  geom_function(fun=F, colour = "black") + 
  geom_function(fun=f_theta, colour = "red", args = list(theta_hat = theta_hat_price,p = p)) +
  theme_classic() + scale_x_continuous(expand = c(0, 0), limits = c(0.22, 1.8), sec.axis = dup_axis(breaks=NULL,name=NULL)) + 
  xlab('Stock price') + ylab('Call value') + scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL))
dev.off()

# FIGURE 2.2
png("S1Figure22.png")
ggplot(df, aes(x=S0)) + geom_function(fun=trueDelta, colour = "black")  +
  geom_function(fun=delta_regression, args = list(theta_hat_price, p) , colour = "red") + 
  theme_classic() + scale_x_continuous(expand = c(0, 0), limits = c(0.22, 1.8), sec.axis = dup_axis(breaks=NULL,name=NULL)) + 
  ylab('Call Delta') + scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL))
dev.off()

#### FIGURE 3
png("S1Figure3.png")
ggplot(df, aes(x=S0, y=D)) + 
  geom_point(color = 'grey', shape = 1) + 
  geom_function(fun=trueDelta, colour = "black", xlim = c(0.5,2)) +
  geom_function(fun=delta_regression, args = list(theta_hat_price,p), colour = "blue", xlim = c(0.5,2)) +
  geom_function(fun=delta_regression, args = list(theta_hat_half,p), colour = "red", xlim = c(0.5,2)) + 
  geom_function(fun=delta_regression_deltaonly, args = list(theta_hat_delta,p) , colour = "orange", xlim = c(0.5,2)) + 
  xlab('Stock price') + ylab('Call Delta') + theme_classic() + 
  ggtitle("Estimated Call Deltas in the Bachelier Model") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  annotate(geom="text", x=0.3, y=1.1, label= "Grey o's: Simulated Deltas", color="grey", hjust = 0) +
  annotate(geom="text", x=0.3, y=0.9, label= "Black: True Delta", color="black", hjust = 0) +
  annotate(geom="text", x=0.3, y=0.85, label= "Blue: Price only-regression", color="blue", hjust = 0) +
  annotate(geom="text", x=0.3, y=0.8, label= "Red: Half/half-regression", color="red", hjust = 0) +
  annotate(geom="text", x=0.3, y=0.75, label= "Orange: Delta only-regression", color="orange", hjust = 0) +
  scale_x_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL), limits = c(0.3,2)) +
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL), limits = c(-0.1,1.2))
dev.off()

#### FIGURE 4 
rm(list=ls())

simulations = c(500, 1000, 2000, 3000, 4000, 5000, 6000)
polydeg <- c(3,4,5,6,7,8)
HedgeError_half   <- matrix(NA, length(simulations), length(polydeg))
HedgeError_price <- matrix(NA, length(simulations), length(polydeg))

sigma = 0.2
dt = 1/52
d = 1
w = 1/2
K = 1
T = 1

for (q in 1:length(polydeg)){
  for (sim in 1:length(simulations)){

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
    
    thetaHat_half  <- matrix(NA, 52,p+1)
    thetaHat_price <- matrix(NA, 52, p+1)
    
    for (t in 0:51){
      
      for (j in 1:n){
        S0[j] <- K + d*sigma * sqrt(T) * N1[j]
        ST[j] <- S0[j] + sigma * sqrt(T-t*dt) * N2[j]
        C[j] <- max(ST[j]-K,0)
        if (ST[j]>= K) {D[j]=1}else {D[j]=0}
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
        St = St + sigma * sqrt(dt) * Z1[t]
        
        Vt_h = at_h * St + bt_h
        Vt_p = at_p * St + bt_p
        
        at_h = deltareg_half(St,p,t+1)
        at_p = deltareg_price(St,p,t+1)
        
        bt_h = Vt_h - at_h * St
        bt_p = Vt_p - at_p * St
      }
      
      S_T <- St + sigma * sqrt(dt)*Z1[52]
      
      VT_h = at_h * S_T + bt_h
      VT_p = at_p * S_T + bt_p
      
      PnL_half[i]  <- VT_h - max(S_T-K,0)
      PnL_price[i] <- VT_p - max(S_T-K,0)
      
    }
    
    HedgeError_half[sim, q]  <- sd(PnL_half)/Ftheta_half(1,p,1)   
    HedgeError_price[sim, q] <- sd(PnL_price)/Ftheta_price(1,p,1)  
  }
}

## Hedging with true Delta 

F <- function(x,t){
  return((x-K)*pnorm((x-K)/(sigma*sqrt(T-t)))+sigma*sqrt(T-t)*dnorm((x-K)/(sigma*sqrt(T-t))))
}


TrueDelta <- function(x,t){
  return(pnorm((x-K)/(sigma*sqrt(T-t))))
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
      St_TD = St_TD + sigma * sqrt(dt) * Z[t]
      Vt_TD = at_TD * St_TD + bt_TD
      at_TD = TrueDelta(St_TD,t*dt)
      bt_TD = Vt_TD - at_TD * St_TD
    }
    ST_TD <- St_TD + sigma * sqrt(dt)*Z[52]
    VT_TD = at_TD * ST_TD + bt_TD
    PnL_true[q] <- VT_TD - max(ST_TD - K,0)
  }
  
  HedgeError_true[j] <- sd(PnL_true)/F(1,0)
}

# Plot
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
  geom_hline(aes(yintercept=SdErrorTD[1]), linetype = "dashed", col = "red") + geom_vline(xintercept = 6000, col = "black") +
  xlim(0,6500) + ylim(0,0.35) + xlab('#simulations') + ylab('Hedge error') + theme_classic() + 
  ggtitle("Polynomial call hedging") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  annotate(geom="text", x=2000, y=0.35, label= "Orange curves: Price only-regression", color="orange", hjust = 0) +
  annotate(geom="text", x=2000, y=0.33, label= "Blue curves: (Price, Delta)-regression", color="blue", hjust = 0) + 
  annotate(geom="text", x=0, y=0.13, label= "Hedging with true Delta", color="red", hjust = 0, size = 3) +
  annotate(geom="text", x=0, y=0, label= "K = 1, T = 1, vol = 0.2, dt = 1/52, #repetitions = 1000", color="black", hjust = 0, size = 3) +
  annotate(geom="text", x=6000, y=0.3, label= "Poly'deg'", color="blue", hjust = 0, size = 3) + 
  annotate(geom="text", x=6100, y=0.28, label= "3", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=0.2, label= "4, 5", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=0.165, label= "6, 7", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=0.145, label= "8", color="blue", hjust = 0, size = 4) +
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000), sec.axis = dup_axis(breaks=NULL,name=NULL)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL))
ggsave("S1Figure4.png")












