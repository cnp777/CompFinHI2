

### Values from article
# Parameters used in Longstaff & Schwartz article 
S0 <- c(rep(36,4), rep(38,4), rep(40,4), rep(42,4), rep(44,4))
sigma <- rep(c(rep(0.2,2), rep(0.4,2)),5)
mat <- rep(c(1,2),10)
K <- 40
r <- 0.06
# Column 4
FD_AMR <- c(4.478, 4.840, 7.101, 8.508, 3.250, 3.745, 6.148, 7.670, 2.314, 2.885, 5.312, 6.920, 1.617, 2.212, 4.582, 6.248,
            1.110, 1.690, 3.948, 5.647) 
# Column 5
ClosedFormEuropean <- round(EUOption(S0, mat, 0, K, r, sigma, 'Put'),3) 


### Functions 
EUOption <- function(S, mat, t, K, r, sigma, type){
  d1 <- 1 / (sigma * sqrt(mat - t)) * (log(S/K) + (r + sigma^2 / 2) * (mat - t))
  d2 <- d1 - sigma * sqrt(mat - t)
  if (type == 'Call')
  {
    return(pnorm(d1) * S - pnorm(d2) * K * exp(-r * (mat - t)))
  } else
  {
    return(pnorm(-d2) * K * exp(-r * (mat - t)) - pnorm(-d1) * S)
  }
}

####################
####    FD.1    ####
####################

TriDiag_solver <- function(a, b, c, d) { # d = right hand side 
  dim <- length(b) 
  f <- rep(0, dim)
  for (i in 1:(dim - 1)) {
    # Gaussian elimination
    a[i+1] <- a[i+1] / b[i]
    b[i+1] <- b[i+1] - a[i+1]*c[i] 
    # Forward substitution 
    d[i+1] <- d[i+1] - a[i+1]*d[i]
  }
  # Backward substitution
  f[dim] <- d[dim] / b[dim]
  for (i in (dim - 1):1) {
    f[i] <- (d[i] - c[i]*f[i+1]) / b[i]
  }
  return(f)
}

TriDiag_solver_AMR <- function(a, b, c, d, payoff) { # d = right hand side 
  dim <- length(b) 
  f <- rep(0, dim)
  for (i in 1:(dim - 1)) {
    # Gaussian elimination
    a[i+1] <- a[i+1] / b[i]
    b[i+1] <- b[i+1] - a[i+1]*c[i]
    # Forward substitution 
    d[i+1] <- d[i+1] - a[i+1]*d[i]
  }
  # Backward substitution
    f[dim] <- max(d[dim] / b[dim], payoff[dim])
    for (i in (dim - 1):1) {
      f[i] <- max((d[i] - c[i]*f[i+1]) / b[i], payoff[i])
    }
  return(f)
}

CrankNicolson <- function(x_min, x_max, x_steps, t_steps, mat, K, r, sigma, European){
  # Grids
  dt <- mat / (t_steps)
  t_grid <- dt * 0 : t_steps
  dx <- (x_max - x_min) / x_steps
  x_grid <- dx * 0 : x_steps + x_min
  
  # Define elements of the matrices
  mu <- r * x_grid
  sigma_term <- sigma^2 * x_grid^2
  
  A <- 0.25 * dt * (sigma_term / dx^2 - mu / dx)
  B <- -1 - 0.5 * dt * (sigma_term / dx^2 + r)
  C <- 0.25 * dt * (sigma_term / dx^2 + mu / dx)
  B_star <- -1 + 0.5 * dt * (sigma_term / dx^2 + r)
  
  # Fill the (quadratic) matrix on the Right-Hand-Side of Munk eq (19) 
  RHS_m <- diag(B_star[-c(1, (x_steps + 1))])
  RHS_m[row(RHS_m) - col(RHS_m) == 1] <- -A[3:(x_steps)]
  RHS_m[col(RHS_m) - row(RHS_m) == 1] <- -C[2:(x_steps - 1)]
  
  # Define vectors to be used in the Munk Tridiagonal Solver. 
  # These are the inputs in matrix on the LHS of Equation 19
  a_vec <- A[-c(1, x_steps + 1)]
  b_vec <- B[-c(1, x_steps + 1)]
  c_vec <- C[-c(1, x_steps + 1)]
  
  # Define price matrix
  f <- matrix(0, nrow = x_steps + 1, ncol = t_steps + 1)
  
  # Boundary conditions
  f[1,] <- exp(-r*(T - t_grid)) * (K - x_min)
  f[, t_steps + 1] <- pmax(K - x_grid, 0)
  
  # Define empty f_{n+1} vector, i.e. the right-hand-side vector of eq (19)
  RHS_vec <- rep(0, x_steps - 1)
  
  # Do successive backwards iterations, for t = N-1,...,0
  if (European == TRUE){
    for (i in t_steps : 1){ 
  # Remember: the first and last i are not included in the TriDiag solver
      RHS_vec[1] <- -A[2] * f[1, i + 1] - A[2] * f[1, i]
      RHS_vec[x_steps - 1] <- -C[x_steps] * f[x_steps + 1, i + 1] - C[x_steps] * f[x_steps + 1, i]
      
  # Multiply the parts together to obtain the RHS of Munk equation (19).
      RHS <- RHS_m %*% f[-c(1, x_steps + 1), i + 1] + RHS_vec
      
  # Solve to get the "f_n"-vector 
      f[-c(1, x_steps + 1), i] <- TriDiag_solver(a_vec, b_vec, c_vec, RHS)
    }
  } else { # type == American
# Remember: it should in each step be considered whether the option holder wishes to exercise or not
    
    payoff <- K - x_grid[-c(1, (x_steps + 1))]
    
    for (i in t_steps : 1){
      f[1,i] <- max(exp(-r*(T - i * dt)) * f[1, i + 1], K-x_min)
      RHS_vec[1] <- -A[2] * f[1, i + 1] - A[2] * f[1, i]
      RHS_vec[x_steps - 1] <- -C[x_steps] * f[x_steps + 1, i + 1] - C[x_steps] * f[x_steps + 1, i]
      RHS <- RHS_m %*% f[-c(1, x_steps + 1), i + 1] + RHS_vec
      
      f[-c(1, x_steps + 1), i] <- TriDiag_solver_AMR(a_vec, b_vec, c_vec, RHS, payoff)
    }
    
  }
  results <- data.frame(f, row.names = x_grid)
  colnames(results) <- t_grid
  return(results)
}

# Replicate the numbers from Column 5
Sigma02_T2 <- CrankNicolson(x_min=0,x_max=100,x_steps=1000,t_steps = 300, mat=2, K=40, r=0.06, sigma=0.2,European= TRUE)
Sigma04_T2 <- CrankNicolson(x_min=0,x_max=100,x_steps=1000,t_steps = 300, mat=2, K=40, r=0.06, sigma=0.4, European=TRUE)

round(Sigma02_T2[c("36", "38", "40", "42", "44"), c("1", "0")],3) 
round(Sigma04_T2[c("36", "38", "40", "42", "44"), c("1", "0")],3)



####################
####    FD.2    ####
####################

# Test order of convergence for European 
x_steps_error <- 1000
t_steps_error <- 300
x_index <- c("0", "20", "40", "60", "80", "98", "100", "102", "120", "140", "160", "180", "200")
t_index <- c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")

v1EU <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  x_steps_error, 
            t_steps = t_steps_error, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = TRUE)
v2EU <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  x_steps_error/2, 
                       t_steps = t_steps_error, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = TRUE)
v3EU <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  x_steps_error/4, 
                       t_steps = t_steps_error, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = TRUE)
v4EU <- CrankNicolson(x_min = 0, x_max = 200, x_steps =   x_steps_error, 
                       t_steps = t_steps_error/2, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = TRUE)
v5EU <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  x_steps_error, 
                       t_steps = t_steps_error/4, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = TRUE)

v1EU_e <- v1EU[x_index,  t_index]
v2EU_e <- v2EU[x_index,  t_index]
v3EU_e <- v3EU[x_index,  t_index]
v4EU_e <- v4EU[x_index,  t_index]
v5EU_e <- v5EU[x_index,  t_index]

# dx
EU_error_x <- (v2EU_e - v3EU_e) / (v1EU_e - v2EU_e)
library(xtable)
xtable(EU_error_x)

# dt
EU_error_t <- (v4EU_e - v5EU_e) / (v1EU_e - v4EU_e)
xtable(EU_error_t)

# Calculate samlet error and runtime
x_steps_error <- 1000
t_steps_error <- 300
x_index <- c("0", "20", "40", "60", "80", "98", "100", "102", "120", "140", "160", "180", "200")
t_index <- c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")

start.time = proc.time()
v1EU <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  x_steps_error,   t_steps = t_steps_error,   mat = 1, K = 100, r = 0.06, sigma = 0.4, European = TRUE)
v2EU <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  x_steps_error/2, t_steps = t_steps_error,   mat = 1, K = 100, r = 0.06, sigma = 0.4, European = TRUE)
v4EU <- CrankNicolson(x_min = 0, x_max = 200, x_steps =   x_steps_error,  t_steps = t_steps_error/2, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = TRUE)
error_1030 = round(v1EU[x_index, t_index] + (v1EU[x_index, t_index] - v2EU[x_index, t_index])/3 + (v1EU[x_index, t_index] - v4EU[x_index, t_index])/3,5)
runtime.time = proc.time()-start.time 

xtable(error_1030) # with v1

# without v1 (only (10.26) and (10.28) errors)
error_1030n = round((v1EU[x_index, t_index] - v2EU[x_index, t_index])/3 + (v1EU[x_index, t_index] - v4EU[x_index, t_index])/3,5)



####################
####    FD.4    ####
####################

# Reuse functions from FD.1
Sigma02_T2_amr <- CrankNicolson(x_min=0,x_max=100,x_steps=1000,t_steps = 300, mat=2, K=40, r=0.06, sigma=0.2, European = FALSE)
Sigma04_T2_amr <- CrankNicolson(x_min=0,x_max=100,x_steps=1000,t_steps = 300, mat=2, K=40, r=0.06, sigma=0.4, European = FALSE)

round(Sigma02_T2_amr[c("36", "38", "40", "42", "44"), c("1", "0")],3) 
round(Sigma04_T2_amr[c("36", "38", "40", "42", "44"), c("1", "0")],3)

xtable = matrix(ncol = 5, nrow = 20)
xtable[,1] <- S0
xtable[,2] <- sigma
xtable[,3] <- mat
xtable[,4] <- FD_AMR
xtable(xtable)

# Global order of convergence - American
x_steps_error <- 1000
t_steps_error <- 300
x_index <- c("0", "20", "40", "60", "80", "98", "100", "102", "120", "140", "160", "180", "200")
t_index <- c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")


v1AMR <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  4*x_steps_error, t_steps = 4*t_steps_error, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = FALSE)
v2AMR <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  2*x_steps_error, t_steps = 4*t_steps_error, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = FALSE)
v3AMR <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  1*x_steps_error, t_steps = 4*t_steps_error, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = FALSE)
v4AMR <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  4*x_steps_error, t_steps = 2*t_steps_error, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = FALSE)
v5AMR <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  4*x_steps_error, t_steps = 1*t_steps_error, mat = 1, K = 100, r = 0.06, sigma = 0.4, European = FALSE)

v1AMR_e <- v1AMR[x_index,  t_index]
v2AMR_e <- v2AMR[x_index,  t_index]
v3AMR_e <- v3AMR[x_index,  t_index]
v4AMR_e <- v4AMR[x_index,  t_index]
v5AMR_e <- v5AMR[x_index,  t_index]

# dx
AMR_error_x <- (v2AMR_e - v3AMR_e) / (v1AMR_e - v2AMR_e)
library(xtable)
xtable(AMR_error_x)

# dt
AMR_error_t <- (v4AMR_e - v5AMR_e) / (v1AMR_e - v4AMR_e)
xtable(AMR_error_t)

# dt when r = 0
v1AMR0 <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  4*x_steps_error, t_steps = 4*t_steps_error, mat = 1, K = 100, r = 0, sigma = 0.4, European = FALSE)
v2AMR0 <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  2*x_steps_error, t_steps = 4*t_steps_error, mat = 1, K = 100, r = 0, sigma = 0.4, European = FALSE)
v3AMR0 <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  1*x_steps_error, t_steps = 4*t_steps_error, mat = 1, K = 100, r = 0, sigma = 0.4, European = FALSE)
v4AMR0 <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  4*x_steps_error, t_steps = 2*t_steps_error, mat = 1, K = 100, r = 0, sigma = 0.4, European = FALSE)
v5AMR0 <- CrankNicolson(x_min = 0, x_max = 200, x_steps =  4*x_steps_error, t_steps = 1*t_steps_error, mat = 1, K = 100, r = 0, sigma = 0.4, European = FALSE)

v1AMR0_e <- v1AMR0[x_index,  t_index]
v2AMR0_e <- v2AMR0[x_index,  t_index]
v3AMR0_e <- v3AMR0[x_index,  t_index]
v4AMR0_e <- v4AMR0[x_index,  t_index]
v5AMR0_e <- v5AMR0[x_index,  t_index]

AMR0_error_t <- (v4AMR0_e - v5AMR0_e) / (v1AMR0_e - v4AMR0_e)
xtable(AMR0_error_t)


####################
####    FD.5    ####
####################


FDM_Dupire <- function(x_min, x_max, x_steps, mat, spot, r, sigma, t_steps, gamma){
  dt <- mat / (t_steps)
  t_grid <- dt * 0 : t_steps
  dx <- (x_max - x_min) / x_steps
  x_grid <- dx * 0 : x_steps + x_min
  
  # MATRICES
  mu <- r * x_grid
  sigma_CEV <- sigma * x_grid^(gamma - 1)
  sigma_term <- -sigma_CEV^2 * x_grid^2
  
  A <- 0.25 * dt * (sigma_term / dx^2 - mu / dx)
  B <- -1 - 0.5 * dt * (sigma_term / dx^2 + r)
  C <- 0.25 * dt * (sigma_term / dx^2 + mu / dx)
  B_star <- -1 + 0.5 * dt * (sigma_term / dx^2 + r)
  
  # We know LHS now
  LHS_m <- diag(B[-c(1, (x_steps + 1))])
  LHS_m[row(LHS_m) - col(LHS_m) == 1] <- A[3:(x_steps)]
  LHS_m[col(LHS_m) - row(LHS_m) == 1] <- C[2:(x_steps - 1)]
  
  # Values for solver
  a_vec <- -A[-c(1, x_steps + 1)]
  b_vec <- B_star[-c(1, x_steps + 1)]
  c_vec <- -C[-c(1, x_steps + 1)]
  
  # Price matrix
  f <- matrix(0, nrow = x_steps + 1, ncol = t_steps + 1)
  f[1,] <- exp(r*(t_grid)) * spot - x_min
  f[, 1] <- pmax(spot - x_grid, 0)
  
  LHS_vec <- rep(0, x_steps - 1)
  
  # Solve equation
  for (i in 2 : (t_steps + 1)){
    LHS_vec[1] <-  A[2] * f[1, i - 1] + A[2] * f[1, i]
    LHS_vec[x_steps - 1] <- C[x_steps] * f[x_steps + 1, i - 1] + C[x_steps] * f[x_steps + 1, i]
    LHS <- LHS_m %*% f[-c(1, x_steps + 1), i - 1] + LHS_vec
    
    f[-c(1, x_steps + 1), i] <- TriDiag_solver(a_vec, b_vec, c_vec, LHS)
  }
  
  for (i in 2:(t_steps + 1)){
    f[,i] <- exp(-r * t_grid[i]) * f[,i]
  }
  
  result <- data.frame(f, row.names = x_grid)
  colnames(result) <- t_grid
  
  return(result)
}

# Values for plot
dt <- 1 / 300
t_grid <- dt * 0 : 300
dx <- 200 / 1000
x_grid <- dx * 0 : 1000
prices_dupire <- FDM_Dupire(0, 200, 1000, 1, 100, 0.06, 3, 300, 0.5)
put_prices_strike <- x_grid * exp(-0.06) + prices_dupire[, "1"] - 100
put_prices_mat <- 100 * exp(-0.06 * t_grid) + prices_dupire["100",] - 100

# Plots

par(mfrow=c(2,2)) 

plot(x_grid, prices_dupire[, "1"], type = 'l', xlab = 'Strike', ylab = 'Call Price')
lines(x_grid, EUOption(100, 1, 0, x_grid,0.06, 0.2, 'Call'), col = 'red')
legend(125, 100, legend = c('CEV', 'BS'), col = c('black', 'red'), lty = 1)

plot(x_grid, put_prices, type = 'l', xlab = 'Strike', ylab = "Put Price")
lines(x_grid, EUOption(100, 1, 0, x_grid,0.06, 0.2, 'Put'), col = 'red')
legend(5, 85, legend = c('CEV', 'BS'), col = c('black', 'red'), lty = 1)

plot(t_grid, prices_dupire["100",], type = 'l', xlab = 'Time to maturity', ylab = "Call Price")
lines(t_grid, EUOption(100, t_grid, 0, 100 ,0.06, 0.2, 'Call'), col = 'red')
legend(0, 14, legend = c('CEV', 'BS'), col = c('black', 'red'), lty = 1)

plot(t_grid, put_prices_mat, type = 'l', xlab = 'Time to maturity', ylab = "Put Price")
lines(t_grid, EUOption(100, t_grid, 0, 100 ,0.06, 0.2, 'Put'), col = 'red')
legend(0, 8.8, legend = c('CEV', 'BS'), col = c('black', 'red'), lty = 1)










