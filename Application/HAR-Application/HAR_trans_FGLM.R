require(extraDistr)
require(gss)
require(pracma)
require(RandomFieldsUtils)
require(MASS)
require(truncnorm)
require(caret)
require(caTools)
require(glmnet)

source("data_generation.R")
source("trans_func_basis.R")

# Cluster specified parameter
# args <-  as.numeric(commandArgs(trailingOnly=TRUE))
# set.seed(2022*args)

load(file = "HAR_processed.RData")

grid.size = dim(All_Data$person1$X_target)[2]; M = grid.size;
grid = list( pt = (1:grid.size)/grid.size, wt = rep(1/grid.size, grid.size) )

Person = names(All_Data)
Kernels = c('Exp', 'M3/2')

eachkernel_res <- list()
for(kernel in Kernels){
  # calculate eigenvalue and eigenfunction based on 
  RKHS_kernel = function(x1, x2){
    Matern_kernel(x1 = x1, x2 = x2, kernel = kernel, rho = 1)
  }
  basis_info = Calculate_basis(grid, M)
  e_val = basis_info$e_val
  e_fun = basis_info$e_fun
  
  eachperson_res = list()
  for(per in Person){
    print(per)
    
    # if(per == 'Person3') break
    
    xmat_target = All_Data[[per]]$X_target
    ymat_target = All_Data[[per]]$Y_target

    train_index = sample.split(Y = ymat_target, SplitRatio = .5)
    xmat_target_train = subset(xmat_target, train_index == T); ymat_target_train = subset(ymat_target, train_index == T)
    xmat_target_test = subset(xmat_target, train_index == F); ymat_target_test = subset(ymat_target, train_index == F)
    
    
    DATA = list(X_target = xmat_target_train, Y_target = ymat_target_train,
                 X_auxiliary = All_Data[[per]]$X_auxiliary, Y_auxiliary = All_Data[[per]]$Y_auxiliary)
    
    labels = matrix(NA, nrow = 7, ncol = length(ymat_target_test))
    rownames(labels) <- c("original", "nontrans", "pooled", "oracle-trans", "Detect Trans" ,  "Q-trans", "star-trans")
    
    etas = matrix(NA, nrow = 6, ncol = length(ymat_target_test))
    rownames(etas) <- c("nontrans", "pooled", "oracle-trans", "Detect Trans" ,  "Q-trans", "star-trans")
    
    # non transfer results
    nontrans_res = opt_FGLM(xmat = DATA$X_target, ymat = DATA$Y_target, grid = grid, offset = NULL)
    
    # pooled and oracle
    sktl_res = sktl_glm(Data = DATA, grid = grid)
    
    #### detection approach
    detect_res = trans_detect(Data = DATA, grid = grid, Eps = 0.01)
    detect_index = detect_res$index
    if(is.null(detect_index)){
      trans_detect_res = nontrans_res
    }else{
      DATA2 = list( X_target = DATA$X_target,
                    Y_target = DATA$Y_target,
                    X_auxiliary = DATA$X_auxiliary[detect_index],
                    Y_auxiliary = DATA$Y_auxiliary[detect_index] )
      trans_detect_res = sktl_glm(Data = DATA2, grid = grid)
    }

    #### QTrans-Osflr
    aggtrans_res = satl_glm(Data = DATA, grid = grid)
    
    etas[1,] = as.numeric(xmat_target_test%*%(grid$wt*nontrans_res$beta_hat)) + nontrans_res$alpha_hat 
    etas[2,] = as.numeric(xmat_target_test%*%(grid$wt*sktl_res$beta_pooled)) + sktl_res$alpha_pooled
    etas[3,] = as.numeric(xmat_target_test%*%(grid$wt*sktl_res$beta_hat)) + sktl_res$alpha_hat
    etas[4,] = as.numeric(xmat_target_test%*%(grid$wt*trans_detect_res$beta_hat)) + trans_detect_res$alpha_hat
    etas[5,] = as.numeric(xmat_target_test%*%(grid$wt*aggtrans_res$Qagg_res$beta_hat)) + aggtrans_res$Qagg_res$alpha_hat
    etas[6,] = as.numeric(xmat_target_test%*%(grid$wt*aggtrans_res$Staragg_res$beta_hat)) + aggtrans_res$Staragg_res$alpha_hat
    
    labels[1,] = as.numeric(ymat_target_test)
    labels[2,] = ifelse( psi_1(etas[1,]) >0.5, 1, 0 )
    labels[3,] = ifelse( psi_1(etas[2,]) >0.5, 1, 0 )
    labels[4,] = ifelse( psi_1(etas[3,]) >0.5, 1, 0 )
    labels[5,] = ifelse( psi_1(etas[4,]) >0.5, 1, 0 )
    labels[6,] = ifelse( psi_1(etas[5,]) >0.5, 1, 0 )
    labels[7,] = ifelse( psi_1(etas[6,]) >0.5, 1, 0 )
    

    
    eachperson_res[[per]] <- list(labels = labels, etas = etas)
    
  }
  eachkernel_res[[kernel]] = eachperson_res
}


# save(eachkernel_res, file = paste("trans_sflr_res",args, ".RData",sep=""))



