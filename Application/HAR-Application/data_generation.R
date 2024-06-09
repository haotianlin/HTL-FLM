#### kernels ####
Matern_kernel <- function(x1,x2, kernel = "Exp", rho = 1){
  Ds <- abs(x1-x2)
  if(kernel=="Exp"){
    Sig_new <- exp(-Ds/rho)
  }
  if(kernel=="M3/2"){ # kernel of Matern Process nu=3/2
    Sig_new <- (1+sqrt(3)*(-Ds/rho))*exp(-sqrt(3)*Ds/rho)
  }
  if(kernel=="M5/2"){ # kernel of Matern Process nu=5/2
    Sig_new <- (1+(sqrt(5)*Ds/rho)+(5*(Ds^2)/(3*rho^2)))*exp(-sqrt(5)*Ds/rho)
  }
  if(kernel=="M7/2"){ # kernel of Matern Process nu=7/2
    Sig_new <- (1+(sqrt(7)*Ds/rho)+(14*(Ds^2)/(5*rho^2))+(7*sqrt(7)*Ds^3)/(15*rho^3))*exp(-sqrt(7)*Ds/rho)
  }
  if(kernel=="Gau"){
    Sig_new <- exp(-(Ds^2)/(2*rho))
  }
  return( Sig_new )
}


self_defined_kernel = function(x1,x2, r1 = r1){
  val = 0
  for(k in 1:100){
    val = val + (1/(k)^{2*r1})*phi_k(x = x1, k = k) *phi_k(x = x2, k = k)
  }
  return(val)
}

sobolev_kernel = function(x1,x2, r1 = 1){
  k2 <- function(x) ((x - 0.5)^2 - 1/12)/2
  k4 <- function(x) ((x - 0.5)^4 - (x - 0.5)^2/2 + 7/240)/24
  k2(x1) * k2(x2) - k4(abs(x1 - x2))
}

sobolev_full_kernel = function(x1,x2){
  
  1 + (x1 - 0.5)*(x2 - 0.5) + sobolev_kernel(x1, x2)
  
}



#### basis ####
phi_k = function(x,k){
  sqrt(2)*cos(k*pi*x)
}


#### generate predictor via Gaussian process ####
X_func_target <- function(x )
{
  val = 0
  
  C = function(s,t){
    Matern_kernel(x1 = s, x2 = t, kernel = "Exp", rho = 1/15)
  }
  
  Sigma = outer(x,x,C)
  
  val <- as.numeric(val + mvrnorm(mu = rep(0,length(x)),Sigma = Sigma))
  return(val)
}

X_func_source <- function(x )
{
  # val <- sin(pi*x)
  val = 0
  
  C = function(s,t){
    Matern_kernel(x1 = s, x2 = t, kernel = "M3/2", rho = 1/15)
  }
  
  
  
  Sigma = outer(x,x,C)
  
  val <- as.numeric(val + mvrnorm(mu = rep(0,length(x)),Sigma = Sigma))
  return(val)
}



#### beta ####
beta_true = function(x){
  
  if(beta_ture_type == 1){
    val <- 0
    f_k <- rep(0,50)
    for(k in 1:50){
      f_k[k] <- 4*sqrt(2)*(-1)^{k-1}*k^{-2}
      val <- val + f_k[k]*phi_k(x,k)
    }
    return(val)
  }else if(beta_ture_type == 2){
    return( 4*cos(3*pi*x) )
  }else if(beta_ture_type == 3){
    return( 4*sin(3*pi*x) + 4*cos(3*pi*x) )
  }else{
    return( 5*(phi_k(x,1) + phi_k(x,2)/2 + phi_k(x,3)/3)  )
  }
  
}

beta_trans = function(x, h){
  val = 0
  f_k = rep(0,100)
  for(k in 1:100){
    R = runif(1, min = -1, 1)
    # R = rsign(1)
    f_k[k] = R*((sqrt(6)*h)/(pi * k^2))
    val = val + f_k[k]*phi_k(x,k)
  }
  return(val + beta_true(x))
}

beta_nontrans = function(x, mu , K){
  Mu = mu(x)
  # Sigma = outer(x, x, K)
  
  Sigma = sapply(x, function(s1) {
    sapply(x, function(s2) {
      K(s1, s2)
    })
  })
  beta = 5*mvrnorm(mu = Mu, Sigma = Sigma) + 5
  return(beta)
}


### generating data ###
gene_transdata = function(N_target, N_source, total_size, trans_index, h, 
                          nontrans_mu = function(x){ return(5*cos(2*pi*x)) },
                          nontran_kernel = function(s,t){ return(exp(-abs(s-t))) },
                          grid, sigma = 0.5){
  
  eta_target = eta_trans = eta_nontrans = NULL
  ##### for target dataset ########
  beta0 = beta_true(grid$pt)
  Xmat_target = Ymat_target = NULL
  for(i in 1:N_target){
    X = X_func_target(x = grid$pt)
    Xmat_target = rbind(Xmat_target, X)
    eta_target = c(eta_target, sum(grid$wt*beta0*X))
    Y = rbinom(1,1,prob = 1/(1+exp(-sum(grid$wt*beta0*X))))
    Ymat_target = c(Ymat_target,Y)
  }
  
  ##### for source dataset ########
  Xmat_auxiliary = Ymat_auxiliary = list()
  Beta_trans = c()
  Beta_nontrans = c()
  for(I in 1:total_size){
    if(I %in% trans_index){
      beta1 = beta_trans(x = grid$pt, h = h)
      Xmat = Ymat = NULL
      for(i in 1:N_source){
        X = X_func_source(x = grid$pt)
        Xmat = rbind(Xmat, X)
        eta_trans = c(eta_trans, sum(grid$wt*beta1*X))
        Y = rbinom(1,1,prob = 1/(1+exp(-sum(grid$wt*beta1*X))))
        Ymat = c(Ymat,Y)
      }
      Xmat_auxiliary[[I]] = Xmat
      Ymat_auxiliary[[I]] = Ymat
      Beta_trans = rbind(Beta_trans,beta1)
      
    }else{
      beta2 = beta_nontrans(x = grid$pt, mu = nontrans_mu, K = nontran_kernel)
      Xmat = Ymat = NULL
      for(i in 1:N_source){
        X = X_func_source(x = grid$pt)
        Xmat = rbind(Xmat, X)
        eta_nontrans = c(eta_nontrans, sum(grid$wt*beta2*X))
        Y = rbinom(1,1,prob = 1/(1+exp(-sum(grid$wt*beta2*X))))
        Ymat = c(Ymat,Y)
      }
      Xmat_auxiliary[[I]] = Xmat
      Ymat_auxiliary[[I]] = Ymat
      Beta_nontrans = rbind(Beta_nontrans,beta2)
    }
    
  }
  names(Xmat_auxiliary) = paste("A", 1:total_size, sep = "")
  names(Ymat_auxiliary) = paste("A", 1:total_size, sep = "")
  
  
  out = list(X_target = Xmat_target,
             Y_target = Ymat_target,
             X_auxiliary = Xmat_auxiliary,
             Y_auxiliary = Ymat_auxiliary,
             trans_index = trans_index, 
             beta_trans = Beta_trans,
             beta_nontrans = Beta_nontrans,
             eta_target = eta_target,
             eta_trans = eta_trans,
             eta_nontrans = eta_nontrans)
  return(out)
  
}


##################### generate only target dataset ###########################
gene_target_data <- function(N_target = 50, grid, sigma = 0.5){
  beta0 = beta_true(grid$pt)
  Xmat_target = Ymat_target = NULL
  eta = NULL
  for(i in 1:N_target){
    X = X_func_target(x = grid$pt)
    Xmat_target = rbind(Xmat_target, X)
    eta = c(eta, sum(grid$wt*beta0*X))
    Y = rbinom(1,1,prob = 1/(1+exp(-sum(grid$wt*beta0*X))))
    Ymat_target = c(Ymat_target,Y)
  }
  out <- list(X = Xmat_target, Y = Ymat_target, eta = eta)
  return(out)
}





#### misclassification ####
GetMisclassification = function(y_true, y_pred){
  sum1 = factor(y_true,levels = c(0,1))
  sum2 = factor(y_pred,levels = c(0,1))
  t = table(sum1, sum2)
  
  accuracy = (t[2,2] + t[1,1]) / (t[2,2] + t[1,2] + t[1,1] + t[2,1])
  misclassfication = 1 - accuracy
  precision = (t[2,2] / (t[2,2] + t[1,2]))
  sensitivity  = t[2,2] / (t[2,2] + t[2,1])
  specificity = t[1,1] / (t[1,1] + t[1,2])
  f1_score = 2*(precision*sensitivity)/(precision + sensitivity)
  
  metrics = data.frame(accuracy = accuracy,
                       misclassfication = misclassfication,
                       precision = precision,
                       sensitivity = sensitivity,
                       specificity = specificity,
                       f1_score = f1_score)
  
  return(misclassfication)
}





#### evaluation metrics ####
Evaluation = function(beta_est, test_target){
  ### calculate prediction error on X_{n+1}
  eta_pred = as.numeric(test_target$X%*%(grid$wt*beta_est))
  # pred_error = mean( (psi_1(eta_pred) - psi_1(test_target$eta))^2 )
  pred_error = mean( (eta_pred - test_target$eta)^2 )
  
  ### calculate the metrics from confusion matrix
  p_pred = 1/(1 + exp(-eta_pred))
  y_pred = ifelse(p_pred > 0.5, 1,0)
  
  sum1 = factor(test_target$Y,levels = c(0,1))
  sum2 = factor(y_pred,levels = c(0,1))
  t = table(sum1, sum2)
  
  accuracy = (t[2,2] + t[1,1]) / (t[2,2] + t[1,2] + t[1,1] + t[2,1])
  misclassfication = 1 - accuracy
  precision = (t[2,2] / (t[2,2] + t[1,2]))
  sensitivity  = t[2,2] / (t[2,2] + t[2,1])
  specificity = t[1,1] / (t[1,1] + t[1,2])
  f1_score = 2*(precision*sensitivity)/(precision + sensitivity)
  
  metrics = data.frame(pred_error = pred_error,
                       accuracy = accuracy,
                       misclassfication = misclassfication,
                       precision = precision,
                       sensitivity = sensitivity,
                       specificity = specificity,
                       f1_score = f1_score)
  
  
  # calculate ROC curve
  df = data.frame(class = test_target$Y, class_pred = y_pred, class_prob = p_pred)
  thresholds <- seq(0,1, by= 0.001)
  TPR <- FPR <- numeric(length(thresholds))
  Positive <- sum(df$class == 1)
  Negative <- sum(df$class == 0)
  for (i in 1:length(thresholds)) {
    # Upper bound probability threshold for equation 10.1
    data_subset <- subset(df, df[ ,'class_prob'] <= thresholds[i])
    # Condition numerator as outlined by indicator function in equation 10.2
    TP <- sum(data_subset[data_subset$class == 1, 'class_prob'] > 0.5) # A
    TN <- sum(data_subset[data_subset$class == 0, 'class_prob'] <= 0.5) # D
    FP <- sum(data_subset[data_subset$class == 0, 'class_prob'] > 0.5) # B
    FN <- sum(data_subset[data_subset$class == 1, 'class_prob'] <= 0.5) # C
    TPR[i] <- 1 - (TP + FN) / Positive # Equation 10.3
    FPR[i] <- 1 - (TN + FP) / Negative # Equation 10.4
  }
  roc.data = data.frame(fpr = FPR, tpr = TPR)
  
  return(list(metrics = metrics,
              roc_data = roc.data))
}











