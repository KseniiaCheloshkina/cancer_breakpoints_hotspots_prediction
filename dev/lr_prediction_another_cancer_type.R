## Data directories
working_folder <- "..\\"
setwd(working_folder)

q_folder <- ".\\data\\data for model\\quadr_mut\\"
s_folder <- ".\\data\\data for model\\stemloops_mut\\"

res_folder <- ".\\data\\adhoc\\another_cancer\\"

## functions
source_path <- "\\run\\helper_functions.R"
source(source_path)

## libraries
library(caret)
library(pROC)
library(e1071)
library(ggplot2)
library(dplyr)
library(reshape2)
library(foreach)
library(doParallel)


## register parallel processing
n_cores <- 2
registerDoParallel(n_cores)
n_iters <- 5

## Create empty df to add metrics in iterations
stat_data_all <- data.frame(iter=0, fold=0,
                            train_auc=0, test_auc=0, dataset=0)
ratio_recall_data_all <- data.frame(iter=0, fold=0, quantile=0, recall=0, dataset=0)
var_imp_all <- data.frame(iter=0, fold=0, imp=0, var=0, dataset=0)
n_ex <- data.frame(dataset=0, n1=0, n0=0)


all_files <- list.files(s_folder)
all_files <- all_files[grep("RData", all_files)]

train_data_nm <- "blood"
train_data <- all_files[grep(train_data_nm, all_files)]

for (file_n in 1:length(train_data)){
  print(paste("Dataset: №", file_n, "Name: ", gsub(".RData", "", train_data[file_n])), sep="")

  ## create train dataset
  # q & sl
  dataset <- make_dataset(folder=q_folder, filename=train_data[file_n],
                          chr_var="v_3", window_var="v_4", add_var=c("v_2"),
                          new_col_names=c("quadr"), )

  all_data <- make_dataset(folder=s_folder, filename=train_data[file_n],
                           chr_var="v_5", window_var="v_6", add_var=c("v_2", "v_3", "v_4"),
                           new_col_names=c("sl1", "sl2", "sl3"), data_all=dataset)

  all_data$chr <- NULL
  all_data$wind <- NULL
  predictors <- setdiff(names(all_data), "v_1")
  
  
  
  ## create test dataset #1
  # q & sl
  breast_fnm <- gsub(train_data_nm, "breast", train_data[file_n])
  breast_dataset <- make_dataset(folder=q_folder, filename=breast_fnm,
                          chr_var="v_3", window_var="v_4", add_var=c("v_2"),
                          new_col_names=c("quadr"), )
  
  breast_all_data <- make_dataset(folder=s_folder, filename=breast_fnm,
                           chr_var="v_5", window_var="v_6", add_var=c("v_2", "v_3", "v_4"),
                           new_col_names=c("sl1", "sl2", "sl3"), data_all=breast_dataset)
  
  breast_all_data$chr <- NULL
  breast_all_data$wind <- NULL


  ## create test dataset #1
  # q & sl
  pancreatic_fnm <- gsub(train_data_nm, "pancreatic", train_data[file_n])
  pancreatic_dataset <- make_dataset(folder=q_folder, filename=pancreatic_fnm,
                                 chr_var="v_3", window_var="v_4", add_var=c("v_2"),
                                 new_col_names=c("quadr"), )
  
  pancreatic_all_data <- make_dataset(folder=s_folder, filename=pancreatic_fnm,
                                  chr_var="v_5", window_var="v_6", add_var=c("v_2", "v_3", "v_4"),
                                  new_col_names=c("sl1", "sl2", "sl3"), data_all=pancreatic_dataset)
  
  pancreatic_all_data$chr <- NULL
  pancreatic_all_data$wind <- NULL
  
  
  ## model building
  # инициализация : данные для одного датасета
  stat_data_breast <- data.frame(iter=0, fold=0, train_auc=0, test_auc=0)
  ratio_recall_data_breast <- data.frame(iter=0, fold=0, quantile=0, recall=0)
  var_imp_breast <- data.frame(iter=0, fold=0, imp=0, var=0)

  
  stat_data_pancr <- data.frame(iter=0, fold=0, train_auc=0, test_auc=0)
  ratio_recall_data_pancr <- data.frame(iter=0, fold=0, quantile=0, recall=0)
  var_imp_pancr <- data.frame(iter=0, fold=0, imp=0, var=0)
  
  
  
  train_ind_list<-list()
  
  
  # данные о балансе классов
  class_data <- all_data
  class1 <- class_data[class_data$v_1 == 1, ]
  class0 <- class_data[class_data$v_1 == 0, ]
  n_ex <- rbind(n_ex, c(gsub(".RData","", all_files[file_n]), nrow(class1), nrow(class0)))

  
  # сплиты для одного датасета
  res_foreach <- foreach(i=seq(1, n_iters), .packages=c('randomForest','pROC', 'caret')) %dopar% { 
    print(i)
    ## балансировка классов
    # перемешиваем, чтобы разбивать на фолды по порядку
    # class 1

    class1 <- class1[sample(nrow(class1), replace=FALSE), ]
    
    f1_cl1 <- class1[1:(round(nrow(class1) / 3)), ]
    f2_cl1 <- class1[(round(nrow(class1) / 3) + 1):(2 * round(nrow(class1) / 3)), ]
    f3_cl1 <- class1[(2 * round(nrow(class1) / 3) + 1):nrow(class1), ]
    
    train_ind_list[[i]] <- list(rownames(f1_cl1), rownames(f2_cl1), rownames(f3_cl1))

    # class 0
    class0 <- class0[sample(nrow(class0), replace=FALSE), ]
    
    f1_cl0 <- class0[1:(round(nrow(class0) / 3)), ]
    f2_cl0 <- class0[(round(nrow(class0) / 3) + 1):(2 * round(nrow(class0) / 3)), ]
    f3_cl0 <- class0[(2 * round(nrow(class0) / 3) + 1):nrow(class0), ]
    

    # final datasets
    dataset_train_1 <- rbind(f1_cl0, f2_cl0, f1_cl1, f2_cl1)
    dataset_test_1 <- rbind(f3_cl0, f3_cl1)
    
    dataset_train_2 <- rbind(f1_cl0, f3_cl0, f1_cl1, f3_cl1)
    dataset_test_2 <- rbind(f2_cl0, f2_cl1)
    
    dataset_train_3 <- rbind(f3_cl0, f2_cl0, f3_cl1, f2_cl1)
    dataset_test_3 <- rbind(f1_cl0, f1_cl1)
  
    # fold 1
    transformers <- lr_train(dataset_train_1, predictors)
    all_data_breast <- lr_predict(dataset_train_1, breast_all_data, i, 3, 
                          predictors, transformers)
    all_data_pancr <- lr_predict(dataset_train_1, pancreatic_all_data, i, 3, 
                                  predictors, transformers)    
    
    
    stat_data_breast <- rbind(stat_data_breast, all_data_breast[[1]])
    ratio_recall_data_breast <- rbind(ratio_recall_data_breast, all_data_breast[[2]])
    var_imp_breast <- rbind(var_imp_breast, all_data_breast[[3]])
    
    stat_data_pancr <- rbind(stat_data_pancr, all_data_pancr[[1]])
    ratio_recall_data_pancr <- rbind(ratio_recall_data_pancr, all_data_pancr[[2]])
    var_imp_pancr <- rbind(var_imp_pancr, all_data_pancr[[3]])
  
    
    # fold 2
    transformers <- lr_train(dataset_train_2, predictors)
    all_data_breast <- lr_predict(dataset_train_2, breast_all_data, i, 2, 
                                  predictors, transformers)
    all_data_pancr <- lr_predict(dataset_train_2, pancreatic_all_data, i, 2, 
                                 predictors, transformers)    
    
    
    stat_data_breast <- rbind(stat_data_breast, all_data_breast[[1]])
    ratio_recall_data_breast <- rbind(ratio_recall_data_breast, all_data_breast[[2]])
    var_imp_breast <- rbind(var_imp_breast, all_data_breast[[3]])
    
    stat_data_pancr <- rbind(stat_data_pancr, all_data_pancr[[1]])
    ratio_recall_data_pancr <- rbind(ratio_recall_data_pancr, all_data_pancr[[2]])
    var_imp_pancr <- rbind(var_imp_pancr, all_data_pancr[[3]])
    

    # fold 3
    transformers <- lr_train(dataset_train_3, predictors)
    all_data_breast <- lr_predict(dataset_train_3, breast_all_data, i, 1, 
                                  predictors, transformers)
    all_data_pancr <- lr_predict(dataset_train_3, pancreatic_all_data, i, 1, 
                                 predictors, transformers)    
    
    
    stat_data_breast <- rbind(stat_data_breast, all_data_breast[[1]])
    ratio_recall_data_breast <- rbind(ratio_recall_data_breast, all_data_breast[[2]])
    var_imp_breast <- rbind(var_imp_breast, all_data_breast[[3]])
    
    stat_data_pancr <- rbind(stat_data_pancr, all_data_pancr[[1]])
    ratio_recall_data_pancr <- rbind(ratio_recall_data_pancr, all_data_pancr[[2]])
    var_imp_pancr <- rbind(var_imp_pancr, all_data_pancr[[3]])
    
    
    breast_res <- list(stat_data_breast, ratio_recall_data_breast, var_imp_breast)
    pancr_res <- list(stat_data_pancr, ratio_recall_data_pancr, var_imp_pancr)
    out_l <- list(breast_res, pancr_res)
    return(out_l)
  }
  
  res_data_breast <- lapply(res_foreach, function(x) x[[1]])
  res_data_pancr <- lapply(res_foreach, function(x) x[[2]])
  
  # соберем из списка в датафрейм для одного датасета
  # breast
  res_breast <- parse_results(res_data_breast, n_iters, 
                              stat_data_breast, ratio_recall_data_breast, 
                              var_imp_breast, breast_fnm)
  res_pancr <- parse_results(res_data_pancr, n_iters, stat_data_pancr, 
                             ratio_recall_data_pancr, 
                             var_imp_pancr, pancreatic_fnm)

  # add dataset data to general results files
  stat_data_all <- rbind(stat_data_all, res_breast[[1]])
  stat_data_all <- rbind(stat_data_all, res_pancr[[1]])
  var_imp_all <- rbind(var_imp_all, res_breast[[2]])
  var_imp_all <- rbind(var_imp_all, res_pancr[[2]])
  ratio_recall_data_all <- rbind(ratio_recall_data_all, res_breast[[3]])
  ratio_recall_data_all <- rbind(ratio_recall_data_all, res_pancr[[3]])
  
  save(stat_data_all, file=paste(res_folder, "stat_data_all.RData", sep=""))
  save(var_imp_all, file=paste(res_folder, "var_imp_all.RData", sep=""))
  save(ratio_recall_data_all, file=paste(res_folder, "ratio_recall_data_all.RData", sep=""))
  save(n_ex, file=paste(res_folder, "n_ex.RData", sep=""))

} 
