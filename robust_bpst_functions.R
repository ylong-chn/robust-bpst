### Function to compute robust BPST estimators
#' These functions are modified from the the `BPST.est.ho()` function from the `BPST` package  
#' Lily Wang, Ming-Jun Lai, GuanNan Wang, Myungjin Kim, Xinyi Li, Jingru Mu, Jue Wang, Yueying Wang, Shan Yu (2021). 
#' BPST: Bivariate spline over triangulation. R package version 0.1.0.
#' https://github.com/FIRST-Data-Lab/BPST

### 1. Robust BPST S-estimator

## Intermediate functions 
#' update weight 
#' @param r residual vector
#' @param bdp delta parameter in M-scale
#' @param rho rho function for M-scale ("bisquare" is usually used)
#' @param cc consistency constant 
weights_update <- function(r, bdp, rho, cc) {
  ms <- pense::mscale(r, bdp)
  r_tilde <- r / ms
  rw <- robustbase::Mwgt(r_tilde, rho, cc = cc)
  rw / sum(rw * r_tilde^2)
}

#' given weights, compute weighted PLS solution 
#' @param B matrix storing values of basis functions 
#' @param Q2 Q2 matrix obtained from QR decomposition
#' @param K energy function a.k.a. penalty matrix 
#' @param W diagonal weight matrix 
#' @param lambda vector of penalty parameters 
#' @param Y vector of values of response variable  
weighted_PLS <- function(B, Q2, K, W, lambda, Y) {
  solve(t(Q2) %*% (t(B) %*% W %*% B + 2*lambda*K) %*% Q2) %*% t(Q2) %*% t(B) %*% W %*% Y
}

#' calculate the S estimator given lambda and initial solution 
#' @param B matrix storing values of basis functions 
#' @param Q2 Q2 matrix obtained from QR decomposition
#' @param K energy function a.k.a. penalty matrix 
#' @param lambda vector of penalty parameters 
#' @param Y vector of values of response variable 
#' @param rho rho function for M-scale ("bisquare" is usually used)
#' @param bdp delta parameter in M-scale 
#' @param theta_start initial solution for S estimation 
#' @param maxit maximum iteration 
#' @param tol tolerance level for convergence 
robust.BPST.est.s.1 <- function(B, Q2, K, lambda, Y, rho = "bisquare", bdp, theta_start,
                              maxit = 500, tol = 1e-6) {
  
  # preliminary computation 
  cc <- pense::consistency_const(bdp, rho)
  B_til <- as.matrix(B%*%Q2)
  # D <- as.matrix(crossprod(t(crossprod(Q2,as.matrix(K))),Q2))
  
  # initialized iteration counter 
  iter <- 0 
  
  # initial solution 
  theta_update <- as.matrix(theta_start, ncol = length(theta_start))
  
  # initial residual
  r <- drop(Y - B_til %*% theta_update)
  
  while(iter < maxit) {
    
    # update weight 
    w <- weights_update(r, bdp, rho, cc)
    W <- diag(w)/sum(w)
    
    # scale lambda
    # lambda <- lambda / sum(w)
    
    # update parameter 
    par_update <- weighted_PLS(B, Q2, K, W, lambda, Y)
    
    # count iteration 
    iter <- iter + 1
    
    # update residual 
    r <- drop(Y - B_til %*% par_update)
    
    # stopping criterion 
    if(sqrt(sum((par_update - theta_update)^2)) > tol) {
      theta_update <- par_update
    }
    else break 
  }
  
  gamma <- crossprod(t(Q2), par_update)
  
  list(
    iteration = iter, 
    theta = par_update, 
    gamma = crossprod(t(Q2), par_update),
    obj_fun = drop(1/2*(pense::mscale(r, bdp))^2 + lambda * t(gamma) %*% K %*% gamma)
  )
  
}


# robust.BPST.est.s <- function(B, Q2, K, lambda, Y, rho, bdp, theta_start
#                               maxit = 500, tol = 1e-6) {
#   
#   
# }
