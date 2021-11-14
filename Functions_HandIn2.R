########################
# Regression functions #
########################

# Estimated call pricing function
f_theta = function(x, theta_hat, p)
{
  sum = 0
  for (w in 0:p){
    sum = sum + theta_hat[w+1]*x^(w)
  }
  return ( sum )
}
Vectorize(f_theta)

# Estimated Delta / Derivative of the estimated call price function wrt the stock price
delta_regression = function(x, theta_hat, p)
{
  sum = 0
  for (s in 1:p){
    sum = sum + theta_hat[s+1]*s*x^(s-1)
  }
  return ( sum )
}
Vectorize(delta_regression)

# Estimated Delta to be used in delta only regression (as dimensions don't fit in the above function)
delta_regression_deltaonly = function(x, theta_hat, p)
{
  input = 0
  for (j in 1:p){
    input = input + theta_hat[j]*j*x^(j-1)
  }
  return ( input )
}
Vectorize(delta_regression_deltaonly)


#######################
#   BACHELIER MODEL   #
#######################

# True call pricing function 
F = function(x)
{
  return( (x-K)*pnorm((x-K)/(sigma*sqrt(T))) + sigma*sqrt(T)*dnorm((x-K)/(sigma*sqrt(T))) )
} 
Vectorize(F)

# True Delta / Derivative of the true call price function wrt the stock price
trueDelta = function(x)
{
  return ( pnorm((x-K)/(sigma*sqrt(T))) )
}
Vectorize(trueDelta)

# Function to generate n samples 
Generate_samples_Bachelier = function(n, T, K, sigma)
{
  S0 = rep(NA, n)
  ST = rep(NA, n)
  C  = rep(NA, n)
  D  = rep(NA, n)
  
  # Simulation
  N1 = rnorm(n,0,1)
  N2 = rnorm(n,0,1)
  
  for (j in 1:n)
  {
    S0[j] = K + sigma*sqrt(T)*N1[j]
    ST[j] = S0[j] + sigma*sqrt(T)*N2[j]
    C[j]  = max(ST[j]-K,0)
    if(ST[j] >= K){ D[j] = 1 }else( D[j] = 0 )
  }
  output=list(S0 = S0, ST = ST, C = C, D = D)
}
Vectorize(Generate_samples_Bachelier)


#########################
#  BLACK-SCHOLES MODEL  #
#########################

# True call pricing function 
F_BS = function(x)
{
  d1 = (log(x/K) + 1/2*sigma^2*T)/(sigma*sqrt(T))
  d2 = d1 - sigma*sqrt(T)
  
  return( x*pnorm(d1) - K*pnorm(d2) )
} 
Vectorize(F_BS)

# True Delta / Derivative of the true call price function wrt the stock price
trueDelta_BS = function(x)
{
  d1 = (log(x/K) + 1/2*sigma^2*(T))/(sigma*sqrt(T))
  
  return ( pnorm(d1) )
}
Vectorize(trueDelta_BS)

F_Digital = function(x,t)
{
  d1 = (log(x/K) + 1/2*sigma^2*(T-t))/(sigma*sqrt(T-t))
  d2 = d1 - sigma*sqrt(T-t)
  
  return( pnorm(d2) )
} 
Vectorize(F_Digital)

trueDelta_Digital = function(x,t)
{
  d1 = (log(x/K) + 1/2*sigma^2*(T-t))/(sigma*sqrt(T-t))
  d2 = d1 - sigma*sqrt(T-t)
  
  return ( dnorm(d2)/(x*sigma*sqrt(T-t)) )
}
Vectorize(trueDelta_Digital)

# Function to generate n samples 
Generate_samples_BS = function(n, T_, K, sigma)
{
  # Create empty vectors 
  S0_BS = rep(NA,n)
  ST_BS = rep(NA,n)
  
  # Related to the Call
  C_BS      = rep(NA,n)
  D_BS_Path = rep(NA,n)
  D_BS_LRM  = rep(NA,n)
  
  # Related to the Digital Call
  DigitalT = rep(NA, n)
  Control1 = rep(NA,n)
  Control2 = rep(NA,n)
  Digital_delta_LRM   = rep(NA, n)
  Digital_delta_mixed = rep(NA, n)
  epsilon = 0.48 
  
  # Simulation
  N1 = rnorm(n,0,1)
  N2 = rnorm(n,0,1)
  
  for (j in 1:n)
  {
    S0_BS[j] = K*exp( -0.5*sigma^2*T_ + sigma*sqrt(T_)*N1[j])
    ST_BS[j] = S0_BS[j]*exp( -0.5*sigma^2*T_ + sigma*sqrt(T_)*N2[j])
     
    # Related to the Call
    C_BS[j]  = max(ST_BS[j]-K,0)
    if(ST_BS[j]>= K){D_BS_Path[j]=ST_BS[j]/S0_BS[j]}else(D_BS_Path[j]=0)
    D_BS_LRM[j] = max(ST_BS[j]-K,0)*(N2[j]/(S0_BS[j]*sigma*sqrt(T_))) 
    
    # Related to the Digital Call
    if(ST_BS[j] > K){DigitalT[j]=1}else(DigitalT[j]=0)
    if(ST_BS[j] > K){Digital_delta_LRM[j] = N2[j]/(S0_BS[j]*sigma*sqrt(T_)) }else(Digital_delta_LRM[j] = 0) 
    if( abs(ST_BS[j]-K) < epsilon){Control1[j] = (1/(2*epsilon))*ST_BS[j]/S0_BS[j]}else(Control1[j] = 0)
    Control2[j] = min(1,max(0, (ST_BS[j]-K + epsilon)/(2*epsilon)))*(N2[j]/(S0_BS[j]*sigma*sqrt(T)))
    Digital_delta_mixed[j] =  Digital_delta_LRM[j] + Control1[j] - Control2[j]
    
  }
  output = list(S0_BS = S0_BS, ST_BS = ST_BS, C_BS = C_BS, D_BS_Path = D_BS_Path, 
                D_BS_LRM = D_BS_LRM, DigitalT = DigitalT, Digital_delta_LRM = Digital_delta_LRM, 
                Digital_delta_mixed = Digital_delta_mixed)
}
Vectorize(Generate_samples_BS)
