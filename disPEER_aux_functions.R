library(assertthat)
library(glmnet)

### functions ###
calc_norm_matrix_d = function(D) {
  D / max(D[lower.tri(D, diag = FALSE)])
}
calc_adj_matrix_d = function(Dnorm, h) {
  A_d = exp(-h * Dnorm ^ 2)
  diag(A_d) = diag(Dnorm)
  A_d
}
construct_adj_matrix_c = function(p, 
                                  centroids_data,
                                  grouping_col) {
  Ac_true = matrix(NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (centroids_data[i, `aparc_aseg_RSN 7`] == 
          centroids_data[j, `aparc_aseg_RSN 7`]) {
        Ac_true[i, j] = 1
      }
      else {
        Ac_true[i, j] = 0
      }
    }
    diag(Ac_true) = 0
  } 
  Ac_true
}
calc_laplacian_c = function(adj) {
  D.diag <- apply(adj, MARGIN = 1, function(row) sum(row))
  D <- diag(D.diag, nrow = nrow(adj), ncol = ncol(adj))
  D - adj
}
calc_norm_laplacian_c = function(L) {
  L.d <- diag(L)
  L.normalized <- L
  p <- dim(L)[2]
  for (u in 1:p) {
    for (v in 1:p) {
      if (L.normalized[u, v] != 0) {
        L.normalized[u, v] <- L.normalized[u,v] / sqrt(L.d[u] * L.d[v])
      }
    }
  }
  L.normalized.d <- rep(1, p)
  L.normalized.d[which(L.d == 0)] <- 0
  diag(L.normalized) <- L.normalized.d
  L.normalized
}
calc_norm_laplacian_d = function(A_D) {
  p = dim(A_D)[1]
  L.norm = matrix(NA, p, p)
  deg = as.vector(rowSums(A_D))
  for (i in 1:p) {
    for (j in 1:p) {
      L.norm[i, j] = -A_D[i, j] / sqrt(deg[i] * deg[j])
    }
  }
  diag(L.norm) = 1
  L.norm
}
calc_eucl_dist = function(regions) {
  n_regions = dim(regions)[1]
  distanceMatrix = matrix(NA, 
                          nrow = n_regions, 
                          ncol = n_regions)
  for (i in 1:n_regions) {
    for (j in 1:n_regions) {
      distanceMatrix[i, j] = sqrt(sum((regions[i, ] - regions[j, ])^2))
    }
  }
  distanceMatrix
}
calc_diss = function(Atrue, Aobs) {
  nom = sum(abs(Aobs - Atrue) > 0)
  denom = 2 * sum(Atrue > 0)
  nom / denom
}
calc_Aobs = function(Atrue, diss) {
  Aobs = Atrue 
  Aobs_graph = graph_from_adjacency_matrix(Aobs, mode = "upper")
  while (round(calc_diss(Atrue, Aobs), 2) != diss) {
    Aobs_graph = rewire(Aobs_graph, 
                        with = keeping_degseq(niter = 1))
    Aobs = as.matrix(as_adjacency_matrix(Aobs_graph))
    diag(Aobs) = 0
    calc_diss(Atrue, Aobs)
  }
  Aobs
}
plot_matrix = function(matrix, which_mat,
                       what_mat, which_hemi,
                       based_on = "") {
  matrix_df = melt(matrix)
  ggplot(data = matrix_df) +
    geom_tile(aes(x = Var2, 
                  y = Var1, 
                  fill = value),
              colour = "grey50") +
    scale_fill_gradient(low = "white", 
                        high = "red",
                        limits = c(0, 1)) +
    labs(x = "Cortical regions", 
         y = "Cortical regions") + 
    # title = paste(which_mat, what_mat, based_on, "matrix", 
    #                "(", which_hemi, "hemisphere )", based_on)) +
    theme_bw() +
    theme(
      axis.text.x = 
        # element_blank(),
        element_text(size = 9,
                     angle = 90,
                     vjust = 0.3),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 11))
}
calc_MSEr = function(b_est, b_real) {
  sum((b_est - b_real) ^ 2) / (sum(b_real ^ 2))
}
calc_MSE = function(b_est, b_real) {
  mean((b_est - b_real) ^ 2)
}
calc_Bias = function(b_est, b_real) {
  mean(abs(b_est - b_real))
}
calc_estimations = function(n, m, 
                            sigma2b, sigma2eps,
                            rC, rD,
                            k, p,
                            Qc_true, Qc_obs,
                            Qd_true, Qd_obs,
                            which_hemi) {
  b = mvrnorm(1, 
              mu = rep(0, p), 
              Sigma = sigma2b * solve(0.1 * diag(p) + 
                                        rC * Qc_true + 
                                        rD * Qd_true))
  Z = matrix(rnorm(n * p), nrow = n, ncol = p)
  X = matrix(rnorm(n * m), nrow = n, ncol = m)
  beta = rep(0, m)
  eta = Z %*% b + X %*% beta
  epsilon = rnorm(n, mean = 0, sd = sqrt(sigma2eps))
  y = eta + epsilon
  
  ### ESTIMATION ###
  ### Ridge estimation: without any additional information ###
  ridge_cv = cv.glmnet(x = cbind(X, Z), y = y, 
                       alpha = 0, intercept = FALSE)
  b_ridge = coef(ridge_cv)[(m + 2):(p + m + 1)] # exclude intercept and betas
  
  Pcx = diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
  y_P = Pcx %*% y
  Z_P = Pcx %*% Z
  eps_tilde = Pcx %*% epsilon
  
  l_tilde_function = function(lambdas) {
    lambdaC = lambdas[1]
    lambdaD = lambdas[2]
    lambdaR = lambdas[3]
    B_lambda =
      lambdaC * Qc_obs +
      lambdaD * Qd_obs +
      lambdaR * diag(p)
    
    n * log(sum(y_P^2) - t(y_P) %*% Z_P %*% solve(B_lambda + t(Z_P) %*% Z_P) %*% t(Z_P) %*% y_P) +
      log(det(B_lambda + t(Z_P) %*% Z_P)) - log(det(B_lambda))
  }
  lambdas_hat = sbplx(x0 = c(1, 1, 1),
                      fn = l_tilde_function,
                      lower = c(10^(-5), 10^(-5), 10^(-5)),
                      upper = c(1e6, 1e6, 1e6))
  lambda_C = lambdas_hat$par[1]
  lambda_D = lambdas_hat$par[2]
  lambda_R = lambdas_hat$par[3]
  
  ### Estimation based on the CONNECTIONS and DISTANCES ###
  b_new = solve(t(Z_P) %*% Z_P + 
                  lambda_C * Qc_obs +
                  lambda_D * Qd_obs +
                  lambda_R * diag(p)) %*% t(Z_P) %*% y_P 
  
  ### Estimation based on the CONNECTIONS ###
  riPEER_1 = riPEER(Q = Qc_obs, 
                    y = y,
                    Z = Z,
                    X = X)
  
  estimations_data_1 = data.table(indexes = 1:p,
                                  b = b,
                                  b_ridge = b_ridge,
                                  b_riPEER = riPEER_1$b.est,
                                  b_disPEER = b_new)
  colnames(estimations_data_1)[5] = "b_disPEER"
  
  estimations_data_2 = data.table(indexes = 1:p,
                                  b_ridge = b_ridge - b,
                                  b_riPEER = riPEER_1$b.est - b,
                                  b_disPEER = b_new - b)
  colnames(estimations_data_2)[4] = "b_disPEER"
  
  
  list(plot_obj_1 = ggplot(data = estimations_data_1, aes(x = indexes)) +
         geom_line(aes(y = b, colour = "True b coefficients"), lty = 2) + 
         geom_point(aes(y = b, colour = "True b coefficients"), lwd = 1.7) +
         geom_line(aes(y = b_ridge, colour = "Ridge estimates"), lty = 2) +
         geom_point(aes(y = b_ridge, colour = "Ridge estimates"), lwd = 1.7) +
         geom_line(aes(y = b_riPEER, colour = "riPEER estimates"), lty = 2) + 
         geom_point(aes(y = b_riPEER, colour = "riPEER estimates"), lwd = 1.7) +
         geom_line(aes(y = b_disPEER, colour = "disPEER estimates"), lty = 2) +
         geom_point(aes(y = b_disPEER, colour = "disPEER estimates"), lwd = 1.7) +
         scale_colour_manual("",
                             breaks = c("True b coefficients", 
                                        "Ridge estimates",
                                        "riPEER estimates", 
                                        "disPEER estimates"),
                             values = c("black", 
                                        "#E7B800",
                                        "#00AFBB", 
                                        "red1")) +
         scale_x_continuous(breaks = 1:p,
                            labels = colnames(Qd_true)) +
         theme(axis.text.x = element_text(size = 9, 
                                          angle = 90,
                                          vjust = 0.3)),
       plot_obj_2 = ggplot(data = estimations_data_2) +
         geom_violin(aes(x = "1", y = b_ridge, fill = "Ridge"), trim = FALSE) +
         geom_violin(aes(x = "2", y = b_riPEER, fill = "riPEER"), trim = FALSE) + 
         geom_violin(aes(x = "3", y = b_disPEER, fill = "disPEER"), trim = FALSE) +
         geom_boxplot(aes(x = "1", y = b_ridge), width = 0.2) +
         geom_boxplot(aes(x = "2", y = b_riPEER), width = 0.2) + 
         geom_boxplot(aes(x = "3", y = b_disPEER), width = 0.2) +
         geom_jitter(aes(x = "1", y = b_ridge)) +
         geom_jitter(aes(x = "2", y = b_riPEER)) + 
         geom_jitter(aes(x = "3", y = b_disPEER)) +
         geom_hline(yintercept = 0, lty = 2) +
         scale_fill_manual(" ",
                           breaks = c("Ridge",
                                      "riPEER", 
                                      "disPEER"),
                           values = c("#E7B800",
                                      "#00AFBB", 
                                      "red1")) +
         ylab("Real and estimate values differences") +
         xlab("Method"),
       MSEr_Ridge = calc_MSEr(b_est = b_ridge,
                              b_real = b),
       MSEr_riPEER = calc_MSEr(b_est = riPEER_1$b.est,
                               b_real = b),
       MSEr_disPEER = calc_MSEr(b_est = b_new,
                                b_real = b),
       lambdas_Ridge = ridge_cv$lambda.min,
       lambdas_riPEER = c(riPEER_1$lambda.Q, 
                          riPEER_1$lambda.R, 
                          riPEER_1$lambda.2),
       lambdas_disPEER = c(lambda_C, lambda_D, lambda_R))
}
calc_multiple_est = function(n, m, 
                             sigma2b, sigma2eps,
                             rC, rD,
                             k, p,
                             Qc_true, Qc_obs,
                             Qd_true, Qd_obs,
                             n_loops) {
  ##### MULTIPLE SIMULATIONS ###
  b_ridge_results = matrix(NA, p, n_loops)
  b_riPEER_results = matrix(NA, p, n_loops)
  b_new_results = matrix(NA, p, n_loops)
  MSEr_results = matrix(NA, n_loops, 3)
  b = mvrnorm(1, 
              mu = rep(0, p), 
              Sigma = sigma2b * solve(0.1 * diag(p) + 
                                        rC * Qc_true + 
                                        rD * Qd_true))
  Z = matrix(rnorm(n * p), nrow = n, ncol = p)
  X = matrix(rnorm(n * m), nrow = n, ncol = m)
  beta = rep(0, m)
  eta = Z %*% b + X %*% beta
  
  for (l in 1:n_loops) {
    epsilon = rnorm(n, mean = 0, sd = sqrt(sigma2eps))
    y = eta + epsilon
    
    ### ESTIMATION ###
    ### Ridge estimation: without any additional information ###
    ridge_cv = cv.glmnet(x = cbind(X, Z), y = y,
                         alpha = 0, intercept = FALSE)
    b_ridge = coef(ridge_cv)[(m + 2):(p + m + 1)] # exclude intercept and betas
    
    Pcx = diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
    y_P = Pcx %*% y
    Z_P = Pcx %*% Z
    eps_tilde = Pcx %*% epsilon
    
    l_tilde_function = function(lambdas) {
      lambdaC = lambdas[1]
      lambdaD = lambdas[2]
      lambdaR = lambdas[3]
      B_lambda =
        lambdaC * Qc_obs +
        lambdaD * Qd_obs +
        lambdaR * diag(p)
      
      n * log(sum(y_P^2) - t(y_P) %*% Z_P %*% solve(B_lambda + t(Z_P) %*% Z_P) %*% t(Z_P) %*% y_P) +
        log(det(B_lambda + t(Z_P) %*% Z_P)) - log(det(B_lambda))
    }
    lambdas_hat = sbplx(x0 = c(1, 1, 1),
                        fn = l_tilde_function,
                        lower = c(10^(-5), 10^(-5), 10^(-5)),
                        upper = c(1e6, 1e6, 1e6))
    lambda_C = lambdas_hat$par[1]
    lambda_D = lambdas_hat$par[2]
    lambda_R = lambdas_hat$par[3]
    
    ### Estimation based on the CONNECTIONS and DISTANCES ###
    b_new = solve(t(Z_P) %*% Z_P + 
                    lambda_C * Qc_obs +
                    lambda_D * Qd_obs +
                    lambda_R * diag(p)) %*% t(Z_P) %*% y_P 
    ### Estimation based on the CONNECTIONS ###
    riPEER_1 = riPEER(Q = Qc_obs, 
                      y = y,
                      Z = Z,
                      X = X)
    
    b_ridge_results[, l] = b_ridge
    b_new_results[, l] = b_new
    b_riPEER_results[, l] = riPEER_1$b.est
    MSEr_results[l, 1] = calc_MSEr(b_est = b_ridge, 
                                   b_real = b)
    MSEr_results[l, 2] = calc_MSEr(b_est = riPEER_1$b.est, 
                                   b_real = b)
    MSEr_results[l, 3] = calc_MSEr(b_est = b_new, 
                                   b_real = b)
  }
  colnames(MSEr_results) = c("Ridge", "riPEER", "disPEER")
  list(median = apply(MSEr_results, 2, median, na.rm = TRUE),
       mean = apply(MSEr_results, 2, mean, na.rm = TRUE))
}
calc_est_MSE = function(n, m, 
                        sigma2b, sigma2eps,
                        rC, rD,
                        k, p,
                        Qc_true, Qc_obs,
                        Qd_true, Qd_obs,
                        n_loops) {
  ##### MULTIPLE SIMULATIONS ###
  b_ridge_results = matrix(NA, n_loops, p)
  b_ols_results = matrix(NA, n_loops, p)
  b_riPEER_results = matrix(NA, n_loops, p)
  b_disPEER_results = matrix(NA, n_loops, p)
  lambdas_dispeer = matrix(NA, n_loops, 3)
  b = mvrnorm(1, 
              mu = rep(0, p), 
              Sigma = sigma2b * solve(0.1 * diag(p) + 
                                        rC * Qc_true + 
                                        rD * Qd_true))
  Z = matrix(rnorm(n * p), nrow = n, ncol = p)
  X = matrix(rnorm(n * m), nrow = n, ncol = m)
  beta = rep(0, m)
  eta = Z %*% b + X %*% beta
  
  for (l in 1:n_loops) {
    epsilon = rnorm(n, mean = 0, sd = sqrt(sigma2eps))
    y = eta + epsilon
    
    ### ESTIMATION ###
    Pcx = diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
    y_P = Pcx %*% y
    Z_P = Pcx %*% Z
    eps_tilde = Pcx %*% epsilon
    
    ### OLS ###
    lm_model = lm(y ~ cbind(X, Z))
    b_ols = as.vector(coef(lm_model)[(m + 2):(p + m + 1)])
    
    ### Ridge estimation: without any additional information ###
    ridge_cv = cv.glmnet(x = cbind(X, Z), y = y,
                         alpha = 0, intercept = FALSE)
    b_ridge = coef(ridge_cv)[(m + 2):(p + m + 1)] # exclude intercept and betas
    
    
    l_tilde_function = function(lambdas) {
      lambdaC = lambdas[1]
      lambdaD = lambdas[2]
      lambdaR = lambdas[3]
      B_lambda =
        lambdaC * Qc_obs +
        lambdaD * Qd_obs +
        lambdaR * diag(p)
      
      n * log(sum(y_P^2) - t(y_P) %*% Z_P %*% solve(B_lambda + t(Z_P) %*% Z_P) %*% t(Z_P) %*% y_P) +
        log(det(B_lambda + t(Z_P) %*% Z_P)) - log(det(B_lambda))
    }
    lambdas_hat = sbplx(x0 = c(1, 1, 1),
                        fn = l_tilde_function,
                        lower = c(10^(-5), 10^(-5), 10^(-5)),
                        upper = c(1e6, 1e6, 1e6))
    lambda_C = lambdas_hat$par[1]
    lambda_D = lambdas_hat$par[2]
    lambda_R = lambdas_hat$par[3]
    
    ### Estimation based on the CONNECTIONS and DISTANCES ###
    b_disPEER = solve(t(Z_P) %*% Z_P + 
                        lambda_C * Qc_obs +
                        lambda_D * Qd_obs +
                        lambda_R * diag(p)) %*% t(Z_P) %*% y_P 
    
    ### Estimation based on the CONNECTIONS ###
    riPEER_1 = riPEER(Q = Qc_obs, 
                      y = y,
                      Z = Z,
                      X = X)
    b_riPEER = riPEER_1$b.est
    
    b_ridge_results[l, ] = b_ridge
    b_ols_results[l, ] = b_ols
    b_riPEER_results[l, ] = b_riPEER
    b_disPEER_results[l, ] = b_disPEER
    lambdas_dispeer[l, ] = lambdas_hat$par
    
  }
  list(t(b_ridge_results),
       t(b_ols_results),
       t(b_riPEER_results),
       t(b_disPEER_results),
       t(b),
       lambdas_dispeer)
}
# split_text = function(text) {
#   splited_text = unlist(strsplit(text, " "))
#   n = length(splited_text)
#   half = round(n / 2)
#   str1 = paste(splited_text[1:half], collapse = " ")
#   str2 = paste(splited_text[(half + 1):n], collapse = " ")
#   cat(str1, str2, sep = "\n")
# }

# h_crossvalidate = function(n, m, 
#                            sigma2b, sigma2eps,
#                            rC, rD, p,
#                            Qc_true, Qc_obs,
#                            k_folds, h_set) {
#   h_n = length(h_set)
#   beta = rep(0, m)
#   Z = matrix(rnorm(n * p), nrow = n, ncol = p)
#   X = matrix(rnorm(n * m), nrow = n, ncol = m)
#   
#   for (i in 1:h_n) {
#     
#     h = h_set[i]
#     Ad_l_true_euc = calc_adj_matrix_d(Dnorm_l_true_euc, h = h)
#     Qd_l_true_euc = calc_norm_laplacian_d(Ad_l_true_euc)
#     Qd_l_obs_euc = Qd_l_true_euc
#     Qd_true = Qd_l_true_euc
#     Qd_obs = Qd_l_obs_euc
#     b = mvrnorm(1, 
#                 mu = rep(0, p), 
#                 Sigma = sigma2b * solve(0.1 * diag(p) + 
#                                           rC * Qc_true + 
#                                           rD * Qd_true))
#     eta = Z %*% b + X %*% beta
#     epsilon = rnorm(n, mean = 0, sd = sqrt(sigma2eps))
#     y = eta + epsilon
#     
#     
#     for (j in 1:k_folds) {
#       ### ESTIMATION ###
#       Pcx = diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
#       y_P = Pcx %*% y
#       Z_P = Pcx %*% Z
#       eps_tilde = Pcx %*% epsilon
#       
#       ### OLS ###
#       lm_model = lm(y ~ cbind(X, Z))
#       b_ols = as.vector(coef(lm_model)[(m + 2):(p + m + 1)])
#       
#       ### Ridge estimation: without any additional information ###
#       ridge_cv = cv.glmnet(x = cbind(X, Z), y = y,
#                            alpha = 0, intercept = FALSE)
#       b_ridge = coef(ridge_cv)[(m + 2):(p + m + 1)] # exclude intercept and betas
#       
#       
#       l_tilde_function = function(lambdas) {
#         lambdaC = lambdas[1]
#         lambdaD = lambdas[2]
#         lambdaR = lambdas[3]
#         B_lambda =
#           lambdaC * Qc_obs +
#           lambdaD * Qd_obs +
#           lambdaR * diag(p)
#         
#         n * log(sum(y_P^2) - t(y_P) %*% Z_P %*% solve(B_lambda + t(Z_P) %*% Z_P) %*% t(Z_P) %*% y_P) +
#           log(det(B_lambda + t(Z_P) %*% Z_P)) - log(det(B_lambda))
#       }
#       lambdas_hat = sbplx(x0 = c(1, 1, 1),
#                           fn = l_tilde_function,
#                           lower = c(10^(-5), 10^(-5), 10^(-5)),
#                           upper = c(1e6, 1e6, 1e6))
#       lambda_C = lambdas_hat$par[1]
#       lambda_D = lambdas_hat$par[2]
#       lambda_R = lambdas_hat$par[3]
#       
#       ### Estimation based on the CONNECTIONS and DISTANCES ###
#       b_disPEER = solve(t(Z_P) %*% Z_P + 
#                           lambda_C * Qc_obs +
#                           lambda_D * Qd_obs +
#                           lambda_R * diag(p)) %*% t(Z_P) %*% y_P 
#       
#       ### Estimation based on the CONNECTIONS ###
#       riPEER_1 = riPEER(Q = Qc_obs, 
#                         y = y,
#                         Z = Z,
#                         X = X)
#       b_riPEER = riPEER_1$b.est
#       
#       b_ridge_results[l, ] = b_ridge
#       b_ols_results[l, ] = b_ols
#       b_riPEER_results[l, ] = b_riPEER
#       b_disPEER_results[l, ] = b_disPEER
#     }
#     
#   }
#   list(t(b_ridge_results),
#        t(b_ols_results),
#        t(b_riPEER_results),
#        t(b_disPEER_results),
#        t(b))
# }



### old code, not needed now
# calc_solver_errors = function(n, m, sigma2b, 
#                               rC, k, p,
#                               Qc_true, Qc_obs,
#                               Qd_true, Qd_obs,
#                               n_loops) {
#   ##### MULTIPLE SIMULATIONS ###
#   errors = 0
#   for (l in 1:n_loops) {
#     Sigma_Z = matrix(NA, p, p)
#     for (i in 1:p) {
#       for (j in 1:p) {
#         Sigma_Z[i, j] = exp(-k * (i - j) ^ 2)
#       }
#     }
#     Z = mvrnorm(n, 
#                 mu = rep(0, p),
#                 Sigma = Sigma_Z)
#     b = mvrnorm(1, 
#                 mu = rep(0, p), 
#                 Sigma = sigma2b * solve(rC * Qc_true + 
#                                           (1 - rC) * Qd_true))
#     X = matrix(rnorm(n * m), nrow = n, ncol = m)
#     beta = runif(m)
#     eta = Z %*% b + X %*% beta
#     epsilon = rnorm(n)
#     y = eta + epsilon
#     
#     ### Estimation ###
#     Pcx = diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
#     y_P = Pcx %*% y
#     Z_P = Pcx %*% Z
#     eps_tilde = Pcx %*% epsilon
#     
#     ### left hemisphere ###
#     l_tilde_function = function(lambdas) {
#       lambdaC = lambdas[1]
#       lambdaD = lambdas[2]
#       lambdaR = lambdas[3]
#       B_lambda =
#         lambdaC * Qc_obs +
#         lambdaD * Qd_obs +
#         lambdaR * diag(p)
#       n * log(t(y_P) %*% y_P -
#                 t(y_P) %*% Z_P %*% solve(B_lambda +
#                                            t(Z_P) %*% Z_P) %*%
#                 t(Z_P) %*% y_P) +
#         log(det(B_lambda + t(Z_P) %*% Z_P)) -
#         log(det(B_lambda))
#     }
#     errors = errors + is.error(sbplx(x0 = c(1, 1, 1),
#                                fn = l_tilde_function,
#                                lower = c(10^(-5), 10^(-5), 10^(-5)),
#                                upper = c(1e6, 1e6, 1e6)))
#   }
#   errors
# }


