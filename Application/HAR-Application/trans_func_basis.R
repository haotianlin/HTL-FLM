#### psi functions 
psi_0 = function(x){ log(1+exp(x)) }
psi_1 = function(x){ exp(x)/(1+exp(x)) }
psi_2 = function(x){ exp(x)/((1+exp(x))^2) }


#### obtain basis information
Calculate_basis = function(grid, M){
  
  Sig_RKHS = outer(grid$pt, grid$pt, RKHS_kernel)
  
  n = length(grid$pt)
  gamma_RKHS=eigen(Sig_RKHS)$values
  e_val_RKHS=abs(gamma_RKHS[1:n]/n)
  e_vec_RKHS=eigen(Sig_RKHS)$vectors[,1:n]*sqrt(n)
  m = length(which(e_val_RKHS>0))
  
  
  e_val = e_val_RKHS[1:M]
  e_fun = e_vec_RKHS[,1:M]
  
  return(list(e_val = e_val,
              e_fun = e_fun))
}

#### fit optimal-FGLM algorithm
opt_FGLM = function(xmat, ymat, grid, offset = NULL){
  lambda_candi = 10^(seq(0,-4,length = 100))
  X_comp = xmat%*%(diag(grid$wt))%*%e_fun
  
  # if(dim(X_comp)[1]>=2000){
  #   index = sample(1:dim(X_comp)[1], size = 2000, replace = F)
  #   X_comp = X_comp[index,]
  #   ymat = ymat[index]
  # }
  
  fit = cv.glmnet(x = X_comp, y = ymat, alpha = 0, family = 'binomial', offset = offset, weights = NULL, nfolds = 10, parallel = FALSE,
                  intercept = TRUE, standardize = FALSE, penalty.factor = 1/e_val, lambda=lambda_candi)
  beta_est = e_fun%*%(as.numeric(coef(fit)[-1]))
  alpha_est = as.numeric(coef(fit)[1])
  # print(fit$lambda.min)
  
  # fit = cv.glmnet(x = X_comp, y = ymat, alpha = 0, family = 'binomial', offset = offset, weights = NULL, nfolds = 10, parallel = FALSE,
  #                 intercept = FALSE, standardize = FALSE, penalty.factor = 1/e_val, lambda=lambda_candi)
  # beta_est = e_fun%*%(as.numeric(coef(fit)[-1]))
  # alpha_est = 0
  
  return(list(beta_hat = beta_est,
              alpha_hat = alpha_est))
}


#### fit source-known transfer learning
sktl_glm = function(Data, grid){
  ########################## transfer step ##############################
  X_trans = Data$X_target
  Y_trans = Data$Y_target
  if(length(Data$Y_auxiliary)>=1){
    for(i in 1:length(Data$X_auxiliary)){
      X_trans = rbind(X_trans, Data$X_auxiliary[[i]])
      Y_trans = c(Y_trans, Data$Y_auxiliary[[i]])
    }
  }
  
  # X_trans = Data$X_auxiliary[[1]]
  # Y_trans = Data$Y_auxiliary[[1]]
  # if(length(Data$Y_auxiliary)>=2){
  #   for(i in 1:length(Data$X_auxiliary)){
  #     X_trans = rbind(X_trans, Data$X_auxiliary[[i]])
  #     Y_trans = c(Y_trans, Data$Y_auxiliary[[i]])
  #   }
  # }
  trans_step_res = opt_FGLM(xmat = X_trans, ymat = Y_trans, grid, offset = NULL)
  beta_pooled = trans_step_res$beta_hat
  alpha_pooled = trans_step_res$alpha_hat
  ######################### debiased step ##############################
  offset = as.numeric(Data$X_target%*%(grid$wt*beta_pooled)) + alpha_pooled
  debias_step_res = opt_FGLM(xmat = Data$X_target, ymat = Data$Y_target, grid, offset = offset) 
  beta_delta = debias_step_res$beta_hat
  alpha_delta = debias_step_res$alpha_hat
  
  beta_hat = beta_pooled + beta_delta
  alpha_hat = alpha_pooled + alpha_delta
  return(list(beta_hat = beta_hat,
              beta_pooled = beta_pooled,
              beta_delta = beta_delta,
              alpha_hat = alpha_hat,
              alpha_pooled = alpha_pooled,
              alpha_delta = alpha_delta))
}


#### source aggregation transfer learning 
satl_glm = function(Data, grid, lambda_factor = 1){
  X_target = Data$X_target
  Y_target = Data$Y_target
  
  K = length(Data$Y_auxiliary)
  
  # id1 = sample.split(Y = Y_target, SplitRatio = .5)
  # X_train = subset(X_target, id1 == T); Y_train = subset(Y_target, id1 == T)
  # X_train_center = X_train ;  Y_train_center = Y_train

  # X_test = subset(X_target, id1 == F); Y_test = subset(Y_target, id1 == F)
  # X_test_center = X_test;  Y_test_center = Y_test
  
  X_train_center = X_train = X_target
  Y_train_center = Y_train = Y_target
  
  X_test_center = X_test = xmat_target_test
  Y_test_center = Y_test = ymat_target_test
    
  
  # construct the candidate index sets
  N_T = length(Y_train)
  init_res = opt_FGLM(xmat = X_train, ymat = Y_train, grid, offset = NULL)
  
  # calculate the MSE scores
  Rank_norm = rep(0,K)
  for(k in 1:K){
    X_trans = Data$X_auxiliary[[k]]; Y_trans = Data$Y_auxiliary[[k]]
    N_A = length(Y_trans)
    init_k = opt_FGLM(xmat = X_trans, ymat = Y_trans, grid, offset = NULL)
    # calculate the RKHS norm between target and k-th source 
    Rank_norm[k] = RKHS_norm(beta = (init_res$beta_hat - init_k$beta_hat), grid = grid)
  }
  Tset = list()
  kk_list = unique(rank(Rank_norm))
  for(k in 1:length(kk_list)){
    Tset[[k]] = which(rank(Rank_norm) <= kk_list[k])
  }
  Tset = unique(Tset)
  
  
  # estimate beta for each candidate set
  alpha_T <- rep(0, 1+length(Tset) )
  beta_T = list()
  beta_T[[1]] = as.vector(init_res$beta_hat)
  alpha_T[1] <- init_res$alpha_hat
  for(k in 1:length(Tset)){
    # trans
    DATA1 = list( X_target = X_train,
                  Y_target = Y_train,
                  X_auxiliary = Data$X_auxiliary[Tset[[k]]],
                  Y_auxiliary = Data$Y_auxiliary[Tset[[k]]] )
    trans_res = sktl_glm(Data = DATA1, grid = grid)
    
    beta_T[[k+1]] = as.vector(trans_res$beta_hat)
    alpha_T[k+1] <- trans_res$alpha_hat
  }
  beta_T = beta_T[!duplicated((beta_T))]
  beta_T = as.matrix(as.data.frame(beta_T))
  colnames(beta_T) = seq(1:ncol(beta_T))
  
  # Q-aggregation
  agg_res = Qagg_func(B = beta_T, A = alpha_T, X_test = X_test, Y_test = Y_test, total_step = 100, selection = F, grid = grid)
  
  Staragg_res = Staragg_func(B = beta_T, A = alpha_T, X_test = X_test, Y_test = Y_test, grid = grid, conflevel = 0.001, Data = Data)
  
  return(list(Qagg_res = agg_res,
              Staragg_res = Staragg_res,
              Tset = Tset) )
}



############# Q_aggregation 
Qagg_func = function(B, A, X_test, Y_test, grid , total_step = 100, selection = F){
  if(sum(B==0)==ncol(B)*nrow(B)){
    return(rep(0,nrow(B)))
  }
  K = ncol(B)
  theta_hat = matrix(0,nrow = K, ncol = 1)
  ## low
  Temp = 1/10
  ## mid
  # Temp = 4
  ## high
  # Temp = length(Y_test)
  
  if(selection){
    
  }else{
    for(k in 1:K){
      theta_hat[k,] = exp(-MSE_func(X = X_test, Y = Y_test, beta = B[,k], alpha = A[k], grid = grid)/Temp)
    }
    theta_hat = theta_hat/(sum(theta_hat))
    theta_old = theta_hat
    alpha = as.numeric(A%*%theta_hat)
    beta = B%*%theta_hat
    beta_we = beta
    for(ss in 1:total_step){
      for(k in 1:K){
        theta_hat[k] = exp(-MSE_func(X = X_test, Y = Y_test, beta = B[,k], alpha = A[k], grid = grid)/Temp +
                             1*MSE1_func(X = X_test, Y = Y_test, beta1 = B[,k], beta2 = beta, alpha1 = A[k], alpha2 = alpha, grid = grid)/(Temp) )
      }
      theta_hat = theta_hat/(sum(theta_hat))
      beta = 1/4*B%*%theta_hat + 3/4*beta
      alpha = as.numeric(1/4*A%*%theta_hat + 3/4*alpha)
      if( sum(abs(theta_hat - theta_old))<10^{-4}  ){break}
      theta_old = theta_hat
    }
  }
  
  return(list(beta_hat = beta,
              alpha_hat = alpha,
              theta = theta_hat,
              beta_we = beta_we))
}


############## star-shape aggregation 
Staragg_func = function(B, A, X_test, Y_test, grid, conflevel = 0.01, Data){
  
  # id2 = c(1:length(Y_test))
  # index = sample(x = id2, size = .5*length(id2), replace = F)
  # # index = sample( x = c(1:length(Y_test)), size = .5*length(Y_test), replace = F)
  # X1 = X_test[index,]; Y1 = Y_test[index]
  # X2 = X_test[-index,]; Y2 = Y_test[-index]
  
  
  # id2 = sample.split(Y = Y_test, SplitRatio = .5)
  # X1 = subset(X_test, id2 == T); Y1 = subset(Y_test, id2 == T)
  # X2 = subset(X_test, id2 == F); Y2 = subset(Y_test, id2 == F)
  
  
  X1 = X2 = X_test
  Y1 = Y2 = Y_test
  
  R_func = function(f, a, X, Y, grid){
    eta_pred = X%*%(grid$wt*f) + a
    y_pred = 1/(1+exp(-eta_pred))
    
    RMSE_trans = sum((y_pred - Y)^2)
    return(RMSE_trans)
  }
  
  fnorm_func = function(f, a, X, grid){
    eta_pred = X%*%(grid$wt*f) + a
    y_pred = 1/(1+exp(-eta_pred))
    return( sqrt( sum( y_pred^2 ) ) )
  }
  
  R1 = R2 = rep(0, ncol(B))
  for(k in 1:ncol(B)){
    R1[k] = R_func(f = B[,k], a = A[k], X = X1, Y = Y1, grid = grid)
    R2[k] = R_func(f = B[,k], a = A[k], X = X2, Y = Y2, grid = grid)
  }
  
  f1_index = which.min(R1)
  f1 = B[,f1_index]
  a1 <- A[f1_index]
  
  f1norm = f2norm = rep(0, ncol(B))
  for(k in 1:ncol(B)){
    f1norm[k] = fnorm_func(f = B[,k] - f1, a = A[k], X = X1, grid = grid)
    f2norm[k] = fnorm_func(f = B[,k] - f1, a = A[k], X = X2, grid = grid)
  }
  
  
  # construct preselection elements
  b_func = function(X, Y, B, A){
    pred = matrix(NA, nrow = length(Y), ncol = ncol(B))
    for(k in 1:ncol(B)){
      pred[,k] = X%*%(grid$wt*B[,k]) + A[k]
    }
    b = max( max(Y), max(pred) )
    return( b )
  }
  M = ncol(B);  n0 = length(Y_test)
  b = b_func(X = Data$X_target, Y = Data$Y_target, B = B, A = A)
  Phi = b*sqrt( (log(M) + conflevel)/n0 )
  c = 4*(1+9*b)
  
  F1 = c()
  for(k in 1:M){
    if( R1[k] <= (R1[f1_index] + c*max( Phi*f1norm[k], Phi^2 )) ){
      F1 = c(F1, k)
    }
  }
  
  convex_mse = c()
  for(k in F1){
    lambda_k = min( max(0.5*((R2[k] - R2[f1_index])/(f2norm[k]^2) + 1),0) , 1) 
    convex_mse = c(convex_mse, R_func(f = (1-lambda_k)*B[,k] + lambda_k*f1, a = (1-lambda_k)*A[k] + lambda_k*a1, X = X2, Y = Y2, grid = grid))
  }
  
  
  j_index = F1[which.min(convex_mse)]
  lambda_j = min( max(0.5*((R2[j_index] - R2[f1_index])/(f2norm[j_index]^2) + 1),0) , 1) 
  
  f_titled = lambda_j*f1 + (1 - lambda_j)*B[,j_index]
  a_titled <- lambda_j*a1 + (1 - lambda_j)*A[j_index]
  
  return(list(beta_hat = f_titled,
              alpha_hat = a_titled,
              index = c(j_index, f1_index),
              lambda = c((1-lambda_j), lambda_j)))
}


################ two MSE functions 
MSE_func = function(X, Y, beta, alpha, grid){
  eta_pred = X%*%(grid$wt*beta) + alpha
  y_pred = 1/(1+exp(-eta_pred))

  RMSE_trans = sum((y_pred - Y)^2)
  return(RMSE_trans)
  
  # eta_pred = X%*%(grid$wt*beta) + alpha
  # sum(-Y*eta_pred + psi_0(eta_pred))
  
}

MSE1_func = function(X, Y, beta1, beta2, alpha1, alpha2, grid){
  eta_pred1 = X%*%(grid$wt*beta1) + alpha1
  y_pred1 = 1/(1+exp(-eta_pred1))
  eta_pred2 = X%*%(grid$wt*beta2) + alpha2
  y_pred2 = 1/(1+exp(-eta_pred2))
  
  RMSE_trans = sum( (y_pred1 - y_pred2)^2 )
  return(RMSE_trans)
}
## RKHS related functions ###
RKHS_norm = function(beta, grid){
  KK = outer(grid$pt, grid$pt, RKHS_kernel)
  # KK = outer(grid$pt, grid$pt, sobolev_full_kernel)
  Hnorm = t(beta)%*%solve(KK)%*%beta
  return(Hnorm)
}


############### detection approach
trans_detect = function( Data, grid, Eps = .01,  C0 = 1){
  
  X_target = Data$X_target
  Y_target = Data$Y_target
  
  index = c()
  index2 = c()
  total_MSE = list()
  
  for(i in 1:length(Data$X_auxiliary)){
    X_trans = Data$X_auxiliary[[i]]
    Y_trans = Data$Y_auxiliary[[i]]
    
    DATA = list( X_target = X_target,
                 Y_target = Y_target,
                 X_auxiliary = list(X_trans),
                 Y_auxiliary = list(Y_trans) )
    N_T = length(Y_target)
    N_A = length(Y_trans)
    
    MSE_nontrans = c()
    MSE_trans = c()
    MSE_trans2 = c()
    
    # plot(grid$pt, beta_0, type = 'l', ylim = c(0,15))
    
    for(II in 1:3){
      target_index = sample( x = c(1:N_T), size = .8*N_T, replace = F)
      
      train_X = X_target[target_index, ]
      train_Y = Y_target[target_index]
      test_X = X_target[-target_index, ]
      test_Y = Y_target[-target_index]

      
      
      # non-trans
      nontrans_res = opt_FGLM(xmat = train_X, ymat = train_Y, grid = grid, offset = NULL)
      
      # non-trans MSE
      RMSE_nontrans = mean(  -test_Y*(test_X%*%(grid$wt*nontrans_res$beta_hat) + nontrans_res$alpha_hat) + psi_0(test_X%*%(grid$wt*nontrans_res$beta_hat) + nontrans_res$alpha_hat)   )
      MSE_nontrans = c(MSE_nontrans, RMSE_nontrans)
      
      # lines(grid$pt, nontrans_res$beta_hat, col = 'blue')
      
      # trans
      DATA1 = list( X_target = train_X,
                    Y_target = train_Y,
                    X_auxiliary = list(X_trans),
                    Y_auxiliary = list(Y_trans) )
      trans_res = sktl_glm(Data = DATA1, grid = grid)
      
      # debias1
      RMSE_trans = mean(  -test_Y*(test_X%*%(grid$wt*trans_res$beta_pooled) + trans_res$alpha_pooled) + psi_0(test_X%*%(grid$wt*trans_res$beta_pooled) + trans_res$alpha_pooled)   )
      # RMSE_trans = mean(  -test_Y*(test_X%*%(grid$wt*trans_res$beta_hat) + trans_res$alpha_hat) + psi_0(test_X%*%(grid$wt*trans_res$beta_hat) + trans_res$alpha_hat)   )
      MSE_trans = c(MSE_trans, RMSE_trans)
      # lines(grid$pt, trans_res$beta_pooled, col = 'red') 
      
    }
    
    source_MSE = cbind(MSE_nontrans, MSE_trans)
    total_MSE[[i]] = source_MSE
    source_MSE = Re(source_MSE)
    
    # if( colMeans(source_MSE)[2] - source_MSE[1] <= C0*max(sd(source_MSE[,2]-source_MSE[,1]), Eps)  ){
    #   index = c(index,i)
    # }
    if( colMeans(source_MSE)[1] > (1 + Eps)*colMeans(source_MSE)[2]){
      index = c(index,i)
    }
  }
  
  source_MSE = c()
  for(i in 1:length(total_MSE)){
    source_MSE = rbind(source_MSE,colMeans(total_MSE[[i]]))
  }
  
  # return(index)
  out = list(index = index,
             source_MSE = source_MSE,
             total_MSE = total_MSE)
  return(out)
  
}






