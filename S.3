
rm(list=ls())

###### S3 - Figure 4 

simulations = c(500, 1000, 2000, 3000, 4000, 5000, 6000)
polydeg <- c(3,8)
HedgeError_lrm  <- matrix(NA, length(simulations), length(polydeg))
HedgeError_path <- matrix(NA, length(simulations), length(polydeg))

sigma = 0.2
dt = 1/52
d = 1
w = 1/2
K = 1
T = 1

set.seed(1)
N1 <- rnorm(6000, 0, 1)
set.seed(2)
N2 <- rnorm(6000,0,1)

for (q in 1:length(polydeg)){
  for (sim in 1:length(simulations)){
    
    p = polydeg[q]
    n = simulations[sim]
    
    S0 <- rep(NA, n)
    ST <- rep(NA, n)
    C  <- rep(NA, n)
    D_lrm  <- rep(NA, n)
    D_path <- rep(NA, n)
    
    thetaHat_lrm  <- matrix(NA, 52,p+1)
    thetaHat_path <- matrix(NA, 52, p+1)
    
    for (t in 0:51){
      
      for (j in 1:n){
        S0[j] <- K*exp(-0.5*sigma^2*(T) + sigma*sqrt(T) * N1[j]) 
        ST[j] <- S0[j]*exp(-0.5*sigma^2*(T-t*dt) + sigma*sqrt(T-t*dt) * N2[j])
        C[j] <- max(ST[j]-K,0)
        
        if(ST[j]>= K){D_path[j]=ST[j]/S0[j]}else(D_path[j]=0)
        D_lrm[j] = max(ST[j]-K,0)*(N2[j]/(S0[j]*sigma*sqrt(T)))
      }
      
      
      X <- matrix(NA, n, p+1)
      Y <- matrix(NA, n, p+1)
      
      for (i in 0:p){
        for (j in 1:n){
          X[j,i+1] = (S0[j])^(i)
          Y[j,i+1] = i*(S0[j])^(i-1)
        }
      }
      
      thetaHat_lrm[t+1,]  <- solve(w*t(X)%*%X+(1-w)*t(Y)%*%Y)%*%(w*t(X)%*%C+(1-w)*t(Y)%*%D_lrm)
      thetaHat_path[t+1,] <- solve(w*t(X)%*%X+(1-w)*t(Y)%*%Y)%*%(w*t(X)%*%C+(1-w)*t(Y)%*%D_path)
      
    }
    
    Ftheta_lrm <- function(x,p,r){
      val <- 0
      for (i in 0:p){
        val = val+thetaHat_lrm[r, i+1]*x^(i)
      }
      return(val)
    }
    
    deltareg_lrm <- function(x,p,r){
      val <- 0
      for (i in 1:p){
        val = val+thetaHat_lrm[r, i+1]*i*x^(i-1)
      }
      return(val)
    }
    
    #Price only
    
    Ftheta_path <- function(x,p,r){
      val <- 0
      for (i in 0:p){
        val = val+thetaHat_path[r, i+1]*x^(i)
      }
      return(val)
    }
    
    
    deltareg_path <- function(x,p,r){
      val <- 0
      for (i in 1:p){
        val = val+thetaHat_path[r, i+1]*i*x^(i-1)
      }
      return(val)
    }
    
    PnL_lrm <- rep(NA,1000)
    PnL_path <- rep(NA,1000)
    
    for (i in 1:1000){
      
      St <- 1  #S(0)
      Z1 <- rnorm(52,0,1)
      
      Vt_h <- Ftheta_lrm(1,p,1)
      at_h <- deltareg_lrm(1,p,1)
      bt_h <- Vt_h-at_h*St
      
      Vt_p <- Ftheta_path(1,p,1)
      at_p <- deltareg_path(1,p,1)
      bt_p <- Vt_p - at_p*St
      
      for (t in 1:51){
        
        St = St*exp(-0.5*sigma^2*dt + sigma * sqrt(dt) * Z1[t])
        
        Vt_h = at_h * St + bt_h
        Vt_p = at_p * St + bt_p
        
        at_h = deltareg_lrm(St,p,t+1)
        at_p = deltareg_path(St,p,t+1)
        
        bt_h = Vt_h - at_h * St
        bt_p = Vt_p - at_p * St
        
      }
      
      S_T <- St*exp(-0.5*sigma^2*dt + sigma * sqrt(dt) * Z1[52]) 
      
      VT_h = at_h * S_T + bt_h
      VT_p = at_p * S_T + bt_p
      
      PnL_lrm[i]  <- VT_h - max(S_T-K,0)
      PnL_path[i] <- VT_p - max(S_T-K,0)
      
    }
    
    HedgeError_lrm[sim, q]  <- sd(PnL_lrm)/Ftheta_lrm(1,p,1)   
    HedgeError_path[sim, q] <- sd(PnL_path)/Ftheta_path(1,p,1)  
  }
}

# Hedging with the true delta (Black-Scholes)
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

df_lrm = data.frame(HedgeError_lrm, row.names = simulations)
colnames(df_lrm) <- polydeg
df_lrm$sim <- simulations

df_path = data.frame(HedgeError_path, row.names = simulations)
colnames(df_path) <- polydeg
df_path$sim <- simulations

ggplot(df_lrm, aes(sim, y=df_lrm[,1])) + geom_point(color = 'blue', shape = 1) + geom_line(color = 'blue') + 
  geom_point(data=df_lrm, aes(sim, y=df_lrm[,2]),color = 'blue', shape = 1) + geom_line(data=df_lrm, aes(sim, y=df_lrm[,2]),color = 'blue') +
  geom_point(data=df_path, aes(sim, y=df_path[,1]),color = 'orange', shape = 1) + geom_line(data=df_path, aes(sim, y=df_path[,1]),color = 'orange') +
  geom_point(data=df_path, aes(sim, y=df_path[,2]),color = 'orange', shape = 1) + geom_line(data=df_path, aes(sim, y=df_path[,2]),color = 'orange') +
  geom_hline(aes(yintercept=HedgeError_true[1]), linetype = "dashed", col = "red") + geom_vline(xintercept = 6000, col = "black") +
  xlab('#simulations') + ylab('Hedge error') + theme_classic() + 
  ggtitle("Polynomial call hedging") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  annotate(geom="text", x=0, y=0.35, label= "Orange curves: (Price, Delta)-regression using Pathwise", color="orange", hjust = 0, size = 3) +
  annotate(geom="text", x=0, y=0.4, label= "Blue curves: (Price, Delta)-regression using LRM", color="blue", hjust = 0, size = 3) + 
  annotate(geom="text", x=0, y=0.135, label= "Hedging with true Delta", color="red", hjust = 0, size = 3) +
  annotate(geom="text", x=0, y=0, label= "K = 1, T = 1, vol = 0.2, dt = 1/52, #repetitions = 1000", color="black", hjust = 0, size = 3) +
  annotate(geom="text", x=6000, y=0.7, label= "Poly'deg'", color="blue", hjust = 0, size = 3) + 
  annotate(geom="text", x=6100, y=df_lrm[7,1], label= "3", color="blue", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=df_path[7,1], label= "3", color="orange", hjust = 0, size = 4) +
  annotate(geom="text", x=6100, y=df_lrm[7,2], label= "8", color="blue", hjust = 0, size = 4) + 
  annotate(geom="text", x=6100, y=df_path[7,2], label= "8", color="orange", hjust = 0, size = 4) +
  scale_x_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL),breaks=c(0,1000,2000,3000,4000,5000,6000), limits = c(0,6500)) + 
  scale_y_continuous(sec.axis = dup_axis(breaks=NULL,name=NULL),breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6, 0.7),  limits = c(0,0.7)) 
ggsave("S3Figure4.png")

