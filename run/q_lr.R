## Data directories
working_folder <- "E:\\Учеба\\Диплом\\Diploma\\scripts\\repo\\"
setwd(working_folder)

q_folder <- ".\\data\\data for model\\quadr_mut\\"
s_folder <- ".\\data\\data for model\\stemloops_mut\\"

res_folder <- ".\\data\\results\\test\\"

## functions
source_path <- "E:\\Учеба\\Диплом\\Diploma\\scripts\\repo\\run\\helper_functions.R"
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
## register parallel processing
n_cores <- 2
registerDoParallel(n_cores)
n_iters <- 15


## Create empty df to add metrics in iterations
stat_data_all <- data.frame(iter=0, fold=0,
                            train_auc=0, test_auc=0, dataset=0)
ratio_recall_data_all <- data.frame(iter=0, fold=0, quantile=0, recall=0, dataset=0)
var_imp_all <- data.frame(iter=0, fold=0, imp=0, var=0, dataset=0)
n_ex <- data.frame(dataset=0, n1=0, n0=0)


all_files <- list.files(s_folder)
all_files <- all_files[grep("RData", all_files)]


for (file_n in 1:length(all_files)){
  print(paste("Dataset: №",file_n,"Name: ",gsub(".RData","",all_files[file_n])),sep="")

  ## create dataset
  # q & sl
  # dataset <- make_dataset(folder=q_folder, filename=all_files[file_n], 
  #                         chr_var="v_3", window_var="v_4", add_var=c("v_2"), 
  #                         new_col_names=c("quadr"), )
  # 
  # all_data <- make_dataset(folder=s_folder, filename=all_files[file_n], 
  #                          chr_var="v_5", window_var="v_6", add_var=c("v_2", "v_3", "v_4"),
  #                          new_col_names=c("sl1", "sl2", "sl3"), data_all=dataset)
  # sl
  all_data <- make_dataset(folder=s_folder, filename=all_files[file_n], 
                           chr_var="v_5", window_var="v_6", add_var=c("v_3", "v_4"),
                           new_col_names=c("coverage_sl_16_50","coverage_sl_6_15"), )
  all_data$chr <- NULL
  all_data$wind <- NULL
  predictors <- setdiff(names(all_data), "v_1")


  ## model building
  # инициализация : данные для одного датасета
  stat_data<-data.frame(iter=0, fold=0,
                        train_auc=0, test_auc=0)
  ratio_recall_data<-data.frame(iter=0, fold=0, quantile=0, recall=0)
  var_imp<-data.frame(iter=0, fold=0, imp=0, var=0)
  
  train_ind_list<-list()
  
  
  # данные о балансе классов
  class_data <- all_data
  class1 <- class_data[class_data$v_1 == 1, ]
  class0 <- class_data[class_data$v_1 == 0, ]
  n_ex <- rbind(n_ex, c(gsub(".RData","", all_files[file_n]), nrow(class1), nrow(class0)))

  
  # сплиты для одного датасета
  res_foreach <- foreach(i=seq(1, n_iters), .packages=c('randomForest','pROC', 'caret')) %dopar% { 

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
    all_data <- lr_train_predict(dataset_train_1, dataset_test_1, i, 3, predictors)

    stat_data <- rbind(stat_data, all_data[[1]])
    ratio_recall_data <- rbind(ratio_recall_data, all_data[[2]])
    var_imp <- rbind(var_imp, all_data[[3]])
  
    
    # fold 2
    all_data <- lr_train_predict(dataset_train_2, dataset_test_2, i, 2, predictors)
    
    stat_data <- rbind(stat_data, all_data[[1]])
    ratio_recall_data <- rbind(ratio_recall_data, all_data[[2]])
    var_imp <- rbind(var_imp, all_data[[3]])
    

    # fold 3
    all_data <- lr_train_predict(dataset_train_3, dataset_test_3, i, 1, predictors)
    
    
    stat_data <- rbind(stat_data, all_data[[1]])
    ratio_recall_data <- rbind(ratio_recall_data, all_data[[2]])
    var_imp <- rbind(var_imp, all_data[[3]])
    
    out_l <- list(stat_data, ratio_recall_data, var_imp)
    return(out_l)
  }  
  
  
  # соберем из списка в датафрейм для одного датасета
  # parse_results <- function(res, n_rsp, stat_data, ratio_recall_data, var_imp, file_n)
  stat_data_d <- lapply(res_foreach, function(x) x[[1]])
  for (l in 1:n_iters){
    stat_data <- rbind(stat_data, stat_data_d[[l]])
  }  
  
  ratio_recall_data_d <- lapply(res_foreach, function(x) x[[2]])
  for (l in 1:n_iters){
    ratio_recall_data <- rbind(ratio_recall_data, ratio_recall_data_d[[l]])
  }  
  
  var_imp_d <- lapply(res_foreach, function(x) x[[3]])
  for (l in 1:n_iters){
    var_imp <- rbind(var_imp, var_imp_d[[l]])
  }
  
  
  stat_data <- stat_data[stat_data$iter != 0, ]
  ratio_recall_data <- ratio_recall_data[ratio_recall_data$iter != 0, ]
  var_imp <- var_imp[var_imp$iter != 0, ]
  
  stat_data$dataset <- gsub(".RData","", all_files[file_n])
  var_imp$dataset <- gsub(".RData","", all_files[file_n])
  ratio_recall_data$dataset <- gsub(".RData","", all_files[file_n])
  
  
  # add dataset data to general results files
  stat_data_all <- rbind(stat_data_all, stat_data)
  var_imp_all <- rbind(var_imp_all, var_imp)
  ratio_recall_data_all <- rbind(ratio_recall_data_all, ratio_recall_data)
  
 save(stat_data_all, file=paste(res_folder, "stat_data_all.RData", sep=""))
 save(var_imp_all, file=paste(res_folder, "var_imp_all.RData", sep=""))
 save(ratio_recall_data_all, file=paste(res_folder, "ratio_recall_data_all.RData", sep=""))
 save(n_ex, file=paste(res_folder, "n_ex.RData", sep=""))

} 
