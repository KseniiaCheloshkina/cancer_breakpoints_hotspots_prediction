
make_dataset <- function(folder, filename, chr_var, window_var, add_var, new_col_names, data_all){
  
  full_path <- paste(folder, filename, sep="") 
  load(full_path)
  
  if(missing(data_all)){

    all_data <- iter_data
    names(all_data)[names(all_data) == chr_var] <- "chr"
    names(all_data)[names(all_data) == window_var] <- "wind"
    
    for (i in 1:length(add_var)){
      names(all_data)[names(all_data) == add_var[i]] <- new_col_names[i]     
    }
    cols <- c(new_col_names, "chr", "wind", "v_1")
    all_data <- all_data[, cols]
    
  } else {

    names(iter_data)[names(iter_data) == chr_var] <- "chr"
    names(iter_data)[names(iter_data) == window_var] <- "wind"
    
    for (i in 1:length(add_var)){
      names(iter_data)[names(iter_data) == add_var[i]] <- new_col_names[i]     
    }
    cols <- c(new_col_names, "chr", "wind")
    all_data <- merge(data_all, iter_data[, cols], by=c("chr", "wind"))
    
  }
  
  return(all_data)
  
}




lr_train_predict_q <- function(dataset_train, dataset_test, iter_n, fold_n, predictors){
 
  transformers <- lr_train(dataset_train, predictors)
  results <- lr_predict(dataset_train, dataset_test, iter_n, fold_n, 
             predictors, transformers)
  return(results)

}




lr_train <- function(dataset_train, predictors){
  
  # оверзамплим на трейне
  dataset_train_cl1 <- dataset_train[dataset_train$v_1 == 1, ]
  n_rep <- round((nrow(dataset_train) - nrow(dataset_train_cl1)) / nrow(dataset_train_cl1)) - 1
  cl_rep <- dataset_train_cl1[rep(seq_len(nrow(dataset_train_cl1)), n_rep), ]
  
  dataset_train <- rbind(dataset_train, cl_rep)
  dataset_train <- dataset_train[sample(nrow(dataset_train), replace=FALSE), ]
  rownames(dataset_train) <- NULL
  
  fmla <- as.formula(paste0("v_1 ", "~", paste0(predictors, collapse = "+")))
  
  
  trans <- preProcess(dataset_train[, predictors], 
                      method = c("center", "scale"))
  
  dataset_train <- predict(trans, dataset_train) 

  #dataset_test <- predict(trans, dataset_test)
  
  logr <- glm(fmla, data=dataset_train, family='binomial', )

  transformers <- list(trans, logr) 

  return(transformers)
  
}




lr_predict <- function(dataset_train, dataset_test, 
                       iter_n, fold_n, 
                       predictors, transformers){
  
  trans <- transformers[[1]]
  logr <- transformers[[2]]
  
  dataset_test <- predict(trans, dataset_test)
  
  imp <- data.frame(imp=as.numeric(logr$coefficients[predictors]),
                    var=predictors)
  
  imp$iter <- iter_n
  imp$fold <- fold_n
  
  un_dataset <- unique(dataset_train)
  pr_train <- predict(logr, newdata=un_dataset, type="response")
  train_auc <- as.numeric(auc(un_dataset$v_1, pr_train))
  
  # prediction
  pr_test <- predict(logr, newdata=dataset_test, type="response")
  test_auc <- as.numeric(auc(dataset_test$v_1, pr_test))
  
  
  # 5, 10, 15, 20, 25 ,30, 35, 40, 45, 50%
  data_test <- as.data.frame(cbind(as.character(dataset_test$v_1), pr_test))
  data_test$pr_test <- as.numeric(as.character(data_test$pr_test))
  data_test$V1 <- as.numeric(as.character(data_test$V1))
  data_test <- data_test[order(data_test$pr_test, decreasing=T), ]
  test_pos <- sum(data_test$V1)
  
  quant <- 1 - seq(0.5, 0.95, 0.05)
  ind_q <- round(nrow(data_test) * quant)
  
  ratio_data <- as.data.frame(cbind(iter_n, fold_n, quant,
                                    sapply(ind_q, function(x) sum(data_test$V1[1:x]) / test_pos)))
  
  names(ratio_data) <- c("iter", "fold", "quantile", "recall")
  
  
  iter_d <- as.data.frame(cbind(iter_n, fold_n, train_auc, test_auc))
  names(iter_d) <- c("iter", "fold", "train_auc", "test_auc")
  
  res <- list(iter_d, ratio_data, imp)
  
  return(res)
  
}



parse_results <- function(res, n_iters, stat_data, ratio_recall_data, var_imp, file_nm){
  
  stat_data_d <- lapply(res, function(x) x[[1]])
  for (l in 1:n_iters){
    stat_data <- rbind(stat_data, stat_data_d[[l]])
  }  
  
  ratio_recall_data_d <- lapply(res, function(x) x[[2]])
  for (l in 1:n_iters){
    ratio_recall_data <- rbind(ratio_recall_data, ratio_recall_data_d[[l]])
  }  
  
  var_imp_d <- lapply(res, function(x) x[[3]])
  for (l in 1:n_iters){
    var_imp <- rbind(var_imp, var_imp_d[[l]])
  }
  
  
  stat_data <- stat_data[stat_data$iter != 0, ]
  ratio_recall_data <- ratio_recall_data[ratio_recall_data$iter != 0, ]
  var_imp <- var_imp[var_imp$iter != 0, ]
  
  stat_data$dataset <- gsub(".RData","", file_nm)
  var_imp$dataset <- gsub(".RData","", file_nm)
  ratio_recall_data$dataset <- gsub(".RData","", file_nm)
  
  res <- list(stat_data, var_imp, ratio_recall_data)
  return(res)
  
}