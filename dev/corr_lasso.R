s_folder <- "E:\\Учеба\\Диплом\\Diploma\\data\\stemloops_mut\\data for pred\\"
load(paste0(s_folder, list.files(s_folder))[151])
new_data <- iter_data[, c("v_2", "v_3", "v_4")]
names(new_data) <- c("coverage_15_30","coverage_16_50","coverage_sl_6_15")
corr <- cor(new_data, method="spearman")

predictors <- c("v_2", "v_3", "v_4")
dataset_train_1 <- iter_data
dataset_train_1_cl1 <- dataset_train_1[dataset_train_1$v_1 == 1, ]
n_rep <- round((nrow(dataset_train_1) - nrow(dataset_train_1_cl1)) / nrow(dataset_train_1_cl1)) - 1
cl_rep <- dataset_train_1_cl1[rep(seq_len(nrow(dataset_train_1_cl1)), n_rep), ]

dataset_train_1 <- rbind(dataset_train_1, cl_rep)
dataset_train_1 <- dataset_train_1[sample(nrow(dataset_train_1), replace=FALSE), ]
rownames(dataset_train_1) <- NULL
# library(caret)
trans <- preProcess(dataset_train_1[, predictors], 
                    method = c("center", "scale"))

dataset_train_1 <- predict(trans, dataset_train_1) 
# dataset_test_1 <- predict(trans, dataset_test_1)
table(dataset_train_1$v_1)
# library(glmnet)
logr <- glmnet(x=as.matrix(dataset_train_1[, predictors]), y=dataset_train_1$v_1, 
               family='binomial', alpha=0, lambda = c(1, 0.1, 0.01, 0.001))
max_devratio_ind <- which.max(logr$dev.ratio)
best_lambda <- logr$lambda[max_devratio_ind]
coef(logr, s=best_lambda)