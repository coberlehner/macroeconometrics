library(vars)
library(expm)
set.seed(1234)

# --------------------------------------------- #
# ---- Write function that simulates data ----- #
# --------------------------------------------- #

simulate_irfs <- function(m, h, n = 5, Ttotal = 750, Tburn = 500) {
  
  # draw coefficient matrices A1 and A2 repeatedly and stop as soon as the eigenvalues
  # lie within the unit circle
  repeat{
  
    # simulate the coefficient matrices
    A1 <- matrix(rnorm(n*n,0,0.2),n,n) # use small sd for fast runtime
    A2 <- matrix(rnorm(n*n,0,0.2),n,n)
    
    # create companion form by stacking the coefficient matrices above an identity and a 0 matrix
    phi <- rbind( 
      cbind(A1,A2),
      cbind(diag(n),matrix(0,n,n))
    )
    
    # check if all eigenvalues lie within the unit circle by checking the maximal eigenvalue
    ev <- eigen(phi)$values
  
    if(max(Mod(ev)) < 1) {
      break
    }
  }
    
  # draw random impact matrix of dimension n*m
  H <- matrix(rnorm(n*m),n,m)
   
  # initialize data for simulation
   Y <- matrix(0,Ttotal,n)
   Y[1,] <- rnorm(n)
   Y[2,] <- rnorm(n)

   # simulate data by making use of the obtained DGP
   for(t in 3:Ttotal) {
     et <- rnorm(m) #draw random errors
     Y[t,] <- A1 %*% Y[t-1,] +
              A2 %*% Y[t-2,] +
              H %*% et
   }
   
   # retain only the last 250 obs
   Y <- Y[(Tburn+1):Ttotal,]
    
   # compute IRFS for horizon h with the formula derived in the previous bullet point
   J <- rbind(diag(n),matrix(0,n,n))
   irf_h <- list()
   
   for (horizon in 0:h) {
     irf_h[[horizon+1]] <- t(J) %*% (phi %^% horizon) %*% J %*% H
   }
   
   # return all simulated matrices, the data and IRFs
   return(list("A1" = A1, "A2" = A2, "H" = H, "Y" = Y, "irf" = irf_h))
   
}

compare_and_plot_irfs <- function(m) {
  
# --------------------------------- #
# --- Estimate Reduced Form VAR --- #
# --------------------------------- #
  
  # simulate data
  data_sim <- simulate_irfs(m=m,h=15)
  y_sim <- data_sim$Y
  irf_true <- data_sim$irf
  
  # estimate a reduced form VAR from simulated data
  rf_var <- VAR(y_sim, p=2)
  sigma <- summary(rf_var)[[3]]
  
  # get coefficient matrices
  n <- 5
  A1 <- matrix(0,n,n)
  A2 <- matrix(0,n,n)
  
  for (row in 1:n) {
    A1[row,] <- summary(rf_var)[[2]][[row]]$coef[1:5]
    A2[row,] <- summary(rf_var)[[2]][[row]]$coef[6:10]
  }
  
  # construct companion form with reduced form coefficient matrices
  phi <- rbind( 
    cbind(A1,A2),
    cbind(diag(n),matrix(0,n,n))
  )
  
  
  # --------------------------------- #
  # ---- zero impact restrictions --- #
  # --------------------------------- #
  
  # perform Cholesky decomposition
  H <- t(chol(sigma))
  
  # use obtained H to compute IRFs
  h <- 15
  J <- rbind(diag(n),matrix(0,n,n))
  irf_zero_rst <- list()
  
  for (horizon in 0:h) {
    irf_zero_rst[[horizon+1]] <- t(J) %*% (phi %^% horizon) %*% J %*% H
  }
  
  
  # -------------------------------- #
  # ------ sign restrictions ------- #
  # -------------------------------- #
  
  # Impose sign restrictions according to true H
  signs <- sign(diag(data_sim$H))
  
  # Draw random orthogonal matrices until one fulfills the restrictions
  accepted <- FALSE
  while(!accepted) {
    
    X <- matrix(rnorm(n^2), n, n)
    Q <- qr.Q(qr(X))
    H_est <- H %*% Q
    
    if(all(sign(diag(H_est)) == signs)) {
      accepted <- TRUE
    }
  }
  
  # use obtained H to compute IRFs
  irf_sign_rst <- list()
  
  for (horizon in 0:h) {
    irf_sign_rst[[horizon+1]] <- t(J) %*% (phi %^% horizon) %*% J %*% H_est
  }
  
  # --------------------------------- #
  # --------- Plot IRFs ------------- #
  # --------------------------------- #
  
  irfs_ts <- function(irfs, i, j) {
    sapply(irfs, function(x) x[i, j])
  }
  
  
  par(mfrow = c(5, 5), mar = c(2, 2, 2, 1))
  
  for (i in 1:5) {
    for (j in 1:5) {
      
      irf_true_ts <- irfs_ts(irf_true, i, j)
      irf_zero_rst_ts <- irfs_ts(irf_zero_rst, i, j)
      irf_sign_rst_ts <- irfs_ts(irf_sign_rst, i, j)
      
      ylim <- range(c(irf_true_ts, irf_zero_rst_ts, irf_sign_rst_ts))
      
      plot(irf_true_ts, type = "l", col = "black", lwd = 2,
           ylim = ylim,
           xlab = "", ylab = "",
           main = bquote( d*y[.(i)] / d*epsilon[.(j)] ))
      
      lines(irf_zero_rst_ts, col = "blue", lwd = 2)
      lines(irf_sign_rst_ts, col = "red", lwd = 2)
      abline(h=0, lty=2)
      
    }
  }
}
