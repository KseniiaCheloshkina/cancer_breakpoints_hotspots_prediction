## libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(xlsx)

## set working folder
current_folder <- "E:\\Учеба\\Диплом\\Diploma\\scripts\\repo\\data\\adhoc\\another_cancer\\"
setwd(current_folder)

## load metrics for blood cancer
data_path <- "../../results/AdditionalFile2.xlsx"
blood_cancer <- read.xlsx(data_path, sheetIndex = 1)
blood_cancer <- blood_cancer %>% filter(Cancer.type == "blood",
                        Model.type == "stemloops and quadruplexes model") %>%
  select(Dataset, Cancer.type, Labeling.type, Aggregation.level, Median.of.test.AUC, 
         Best.probability.quantile.for.threshold.according.to.max.lift.of.recall, 
         Median.recall.for.the.threshokd, Lift.of.recall.for.the.threshold)
names(blood_cancer) <- c("dataset", "cancer_type", "labeling_type", "agg_level",
                         "median_test_auc_blood", "best_prob_q_blood", 
                         "med_recall_blood", "lift_of_recall_blood")
blood_cancer$config <- paste(blood_cancer$agg_level, blood_cancer$labeling_type, sep="_")
blood_cancer$median_test_auc_blood <- as.numeric(as.character(blood_cancer$median_test_auc_blood))

## load testing data

n_ex_files <- "n_ex.RData"
stat_data_all_files <- "stat_data_all.RData" 
ratio_recall_files <- "ratio_recall_data_all.RData"
load(n_ex_files)
load(stat_data_all_files)
load(ratio_recall_files)

## preprocess
n_ex$n1 <- as.numeric(n_ex$n1)
n_ex$n0 <- as.numeric(n_ex$n0)
stat_data_all <- stat_data_all[stat_data_all$dataset != "0", ]
ratio_recall_data_all <- ratio_recall_data_all[ratio_recall_data_all$dataset != "0", ]
n_ex <- n_ex[n_ex$dataset != "0", ]

ratio_recall_data_all <- unique(ratio_recall_data_all)
stat_data_all <- unique(stat_data_all)
stat_data_all$iter <- as.character(stat_data_all$iter)
stat_data_all$fold <- as.character(stat_data_all$fold)

n_ex <- n_ex[, c("dataset", "n1", "n0")]


####### ROC AUC

agg_level_pattern <- "^[0-9]*[a-z]*"
labeling_pattern <- "([0-9]?[.]?[0-9]*)$"

m <- regexpr(agg_level_pattern, stat_data_all$dataset)
stat_data_all$agg_level <- regmatches(stat_data_all$dataset, m)

m <- regexpr(labeling_pattern, stat_data_all$dataset)
stat_data_all$labeling_type <- regmatches(stat_data_all$dataset, m) 

breast_data <- stat_data_all[grep("breast", stat_data_all$dataset), ]
breast_data$cancer_type <- "breast"
names(breast_data)[which(names(breast_data) == "test_auc")] <- "breast"
breast_data$config <- paste(breast_data$agg_level, breast_data$labeling_type, sep="_")
  
pancreatic_data <- stat_data_all[grep("pancr", stat_data_all$dataset), ]
pancreatic_data$cancer_type <- "pancreatic"
names(pancreatic_data)[which(names(pancreatic_data) == "test_auc")] <- "pancreatic"
pancreatic_data$config <- paste(pancreatic_data$agg_level, pancreatic_data$labeling_type, sep="_")

test_cancer <- pancreatic_data %>% 
  select(config, iter, fold, pancreatic) %>%
  inner_join(
    (breast_data %>%
       select(config, iter, fold, breast)
    )
    ) 

full_test_data <- melt(test_cancer, id.vars = c("config", "iter", "fold"), 
                       variable.name = "cancer_type", value.name = "test_auc") %>%
  group_by(config, iter, cancer_type) %>%
  summarize(test_auc = median(test_auc))
# убрали повторяющуюся разметку  
full_test_data <- full_test_data %>%
  filter(!config %in% c("10kb_0.05", "50kb_0.05", "50kb_0.5"))

p <- ggplot(full_test_data, aes(x=test_auc, fill=cancer_type)) +
  geom_density(alpha=0.4) + 
  facet_wrap( ~ config, ncol=4, scales = "free") + 
  theme(legend.position="bottom") + 
  geom_vline(data=blood_cancer, aes(xintercept=median_test_auc_blood),
             linetype="dashed") +
  scale_fill_discrete(name="Cancer type")

p

roc_auc_stat <- full_test_data %>% group_by(cancer_type, config) %>%
  summarize(test_auc = median(test_auc)) %>%
  inner_join(blood_cancer %>% select(config, median_test_auc_blood)) %>%
  mutate(greater_than_blood = as.numeric(test_auc > median_test_auc_blood))
roc_auc_stat$cancer_type <- as.character(roc_auc_stat$cancer_type)

## load metrics for breast and pancreatic cancer
data_path <- "../../results/AdditionalFile2.xlsx"
or_cancer <- read.xlsx(data_path, sheetIndex = 1)
or_cancer <- or_cancer %>% filter((Cancer.type == "breast")|(Cancer.type == "pancreatic"),
                                        Model.type == "stemloops and quadruplexes model") %>%
  select(Dataset, Cancer.type, Labeling.type, Aggregation.level, Median.of.test.AUC)
names(or_cancer) <- c("dataset", "cancer_type", "labeling_type", "agg_level",
                         "median_test_auc_original")
or_cancer$config <- paste(or_cancer$agg_level, or_cancer$labeling_type, sep="_")
or_cancer$median_test_auc_original <- as.numeric(as.character(or_cancer$median_test_auc_original))
or_cancer$cancer_type <- as.character(or_cancer$cancer_type)


roc_auc_stat <- roc_auc_stat %>% left_join(
  (or_cancer %>% select("config", "cancer_type", "median_test_auc_original")), 
  by=c("config", "cancer_type")) %>%
  mutate(greater_than_original = as.numeric(test_auc > median_test_auc_original))


roc_auc_stat %>% group_by(cancer_type) %>%
  summarize(sum(greater_than_blood))

roc_auc_stat %>% group_by(cancer_type) %>%
  summarize(sum(greater_than_original))


write.csv(roc_auc_stat, file = "results.csv", row.names = FALSE)
