library(reshape2)
library(dplyr)


#### load data

input_data_path <- "..\\data\\adhoc\\liver\\"
dataset <- read.csv(paste0(input_data_path, "LIRIJP_500kb_coverages_ed4.csv"), sep="\t")
dataset$X <- NULL
target <- "LIRIJP_bp"
predictors <- setdiff(names(dataset), target)


#### make hotspots (from script data_preparation\8_agg.R)


get_hotspots <- function(dataset, target, percentage){
  p0 <- percentage
  percentage <- 1 - p0
  q <- quantile(dataset[, target], percentage)
  new_col_name <- paste0(target, "_", format(p0, scientific = FALSE), "_hotspot")
  dataset[, new_col_name] <- ifelse(dataset[, target] > q, 1, 0)
  
  return(dataset)
}

# 1 option (0.1%)
# 2 option (0.05%)
# 3 option (0.01%)
# 4 option (0.5%)
# 5 option (1%)
dataset <- get_hotspots(dataset, target, 0.001)
dataset <- get_hotspots(dataset, target, 0.0005)
dataset <- get_hotspots(dataset, target, 0.0001)
dataset <- get_hotspots(dataset, target, 0.005)
dataset <- get_hotspots(dataset, target, 0.01)


dataset[, target] <- NULL
htspt <- names(dataset)[grep("hotspot", names(dataset))]
predictors <- setdiff(names(dataset), htspt)


table(dataset[, htspt[3]])



#### remove correlated feats


cormm <- cor(dataset[, predictors], method = "spearman")
cormm <- data.frame(cormm)
cormm$var1 <- row.names(cormm)
corr_matrix <- melt(cormm, id.vars = "var1")
corr_matrix <- corr_matrix[corr_matrix$var1 != corr_matrix$variable, ]
correlated <- corr_matrix[abs(corr_matrix$value) > 0.8, ]
hist(correlated$value)
correlated$var1 <- as.character(correlated$var1)
correlated$variable <- as.character(correlated$variable)


# get wilcoxon test pvalue for each feature
wt <- data.frame()
for (predictor in predictors){
  pval <- wilcox.test(dataset[,predictor] ~ dataset$LIRIJP_bp_0.01_hotspot)$p.value
  wt <- rbind(wt, data.frame("feature" = predictor, "pval" = pval))
}
wt$feature <- as.character(wt$feature)

correlated <- correlated %>% 
  left_join(wt, by=c("var1" = "feature"))
names(correlated)[length(names(correlated))] <- "pval_var1"

correlated <- correlated %>% 
  left_join(wt, by=c("variable" = "feature"))
names(correlated)[length(names(correlated))] <- "pval_variable"


# sequential feature elimination
corr_feat1 <- unique(correlated$var1)

all_corr_feats <- vector()

for (current_feat in corr_feat1){
  
  if (current_feat %in% unique(correlated$var1)){
    
    current_data <- correlated[correlated$var1 == current_feat, ]
    p1 <- unique(current_data$pval_var1)
    idx_min_2 <- which.min(current_data$pval_variable)
    p2 <- current_data$pval_variable[idx_min_2]
    feat2_best <- current_data$variable[idx_min_2]
    
    if (p1 < p2){
      cor_list <- unique(current_data$variable)
    } else {
      cor_list <- unique(current_data$variable)
      cor_list <- setdiff(cor_list, feat2_best)
      cor_list <- c(cor_list, current_feat)
    }
    
    correlated <- correlated[!correlated$var1 %in% cor_list, ]
    correlated <- correlated[!correlated$variable %in% cor_list, ]
    all_corr_feats <- c(all_corr_feats, cor_list)
    
  }  
}

# check
not_corr <- setdiff(predictors, all_corr_feats)
cormm_new <- cor(dataset[, not_corr], method = "spearman")
cormm <- data.frame(cormm_new)
cormm$var1 <- row.names(cormm)
corr_matrix <- melt(cormm, id.vars = "var1")
corr_matrix <- corr_matrix[corr_matrix$var1 != corr_matrix$variable, ]
corr_matrix[abs(corr_matrix$value) > 0.8, ]
hist(corr_matrix$value)

dataset <- dataset[, c(not_corr, htspt)]



#### generate separate files (from script data_preparation\8_agg.R)

f_to_write = "data for pred\\"
level_n <- "500kb_"
htspt <- names(dataset)[grep("hotspot", names(dataset))]
predictors <- setdiff(names(dataset), htspt)

dataset[, htspt] <- as.data.frame(apply(dataset[, htspt], 2, function(x) as.factor(x)))

for (i in 1:length(htspt)){
  y_name <- htspt[i]
  pr <- c(predictors, y_name)
  iter_data <- dataset[, pr]
  names(iter_data)[which(names(iter_data) == y_name)] <- "v_1"
  flph <- paste(input_data_path, f_to_write, level_n, y_name, ".RData", sep="")
  save(iter_data, file=flph)
}







##### train

input_data_path <- "..\\data\\adhoc\\liver\\data for pred\\"
res_folder <- "..\\data\\adhoc\\liver\\res"

## functions
source_path <- "helper_functions.R"
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



n_cores <- 2
registerDoParallel(n_cores)
n_iters <- 15


## Create empty df to add metrics in iterations
stat_data_all <- data.frame(iter=0, fold=0,
                            train_auc=0, test_auc=0, dataset=0)
ratio_recall_data_all <- data.frame(iter=0, fold=0, quantile=0, recall=0, dataset=0)
var_imp_all <- data.frame(iter=0, fold=0, imp=0, var=0, dataset=0)
n_ex <- data.frame(dataset=0, n1=0, n0=0)

all_files <- list.files(input_data_path)
all_files <- all_files[grep("50kb", all_files)]

for (file_n in 1:length(all_files)){
  
  ph <- paste0(input_data_path, all_files[file_n])
  print(paste("Dataset: №",file_n,"Name: ",gsub(".RData","",all_files[file_n])),sep="")
  load(ph)
  all_data <- iter_data
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
    all_data <- lr_train_predict_q(dataset_train_1, dataset_test_1, i, 3, predictors)
    
    stat_data <- rbind(stat_data, all_data[[1]])
    ratio_recall_data <- rbind(ratio_recall_data, all_data[[2]])
    var_imp <- rbind(var_imp, all_data[[3]])
    
    
    # fold 2
    all_data <- lr_train_predict_q(dataset_train_2, dataset_test_2, i, 2, predictors)
    
    stat_data <- rbind(stat_data, all_data[[1]])
    ratio_recall_data <- rbind(ratio_recall_data, all_data[[2]])
    var_imp <- rbind(var_imp, all_data[[3]])
    
    
    # fold 3
    all_data <- lr_train_predict_q(dataset_train_3, dataset_test_3, i, 1, predictors)
    
    
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
