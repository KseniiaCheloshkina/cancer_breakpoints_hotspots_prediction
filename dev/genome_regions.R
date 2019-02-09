# require libraries
library(dplyr)
require(data.table)
library(ggplot2)
library(ggpubr)
library(reshape2)

# Change working directory
# current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
current_working_dir <- "E://Учеба//Диплом//Diploma//scripts//repo//dev"
setwd(current_working_dir)

##########  FUNCTIONS
# function to load regions data
load_region <- function(region_path){
  region_df <- read.table(region_path)
  region_df <- region_df[, 1:3]
  names(region_df) <- c("chr", "start", "end")
  region_df['chr'] <- apply(region_df, MARGIN = 1, 
                       FUN = function(x) {gsub(x = x["chr"], 
                                               pattern = "chr", 
                                               replacement = "")})
  region_df <- unique(region_df)
  setDT(region_df)
  
  return(region_df)
}



# function to calculate number of breakpoints (bkpt_df) inside a region (region_df)
count_bkpt <- function(region_df, bkpt_df, region_name){
  
  chr_list <- unique(bkpt_df[, chr])
  merged_data <- data.frame()
  
  for (chr in chr_list){
    print(chr)
    bkpt_df_part <- bkpt_df[chr, on="chr"]
    region_df_part <- region_df[chr, on="chr"]
    
    
    joined_df <- region_df_part[bkpt_df_part, 
                                on=.(start <= chr_bkpt_beg_s, end >= chr_bkpt_beg_e),
                                allow.cartesian=TRUE]
    # unique hits
    joined_df <- unique(joined_df)
    
    final_df <- joined_df %>%
      filter(!is.na(chr)) %>%
      group_by(cancer, chr) %>%
      summarize(n_bkpt_intersected = sum(n_bkpt))
    
    names(final_df)[3] <- paste0("n_bkpt_intersected_", region_name)
    
    merged_data <- rbind(merged_data, data.frame(final_df))
    gc()
  }
  
  return(merged_data)
  
}


############ LOAD BREAKPOINTS DATA

cancer_bkpt_path <- "../data/raw breakpoints/"
all_cancers_data <- list.files(cancer_bkpt_path)
all_cancers_data <- all_cancers_data[
  all_cancers_data != "all_cancer_data_eda.csv"]

bkpt <- data.frame()

for (cancer_type_file in all_cancers_data){
  
  bkpt_part <- read.csv(paste0(cancer_bkpt_path, cancer_type_file))
  cancer_type <- gsub(x = cancer_type_file, 
                      pattern = "_all_data.csv", replacement = "")
  
  bkpt_part["cancer"] <- cancer_type
  bkpt <- rbind(bkpt, bkpt_part)
  
}
# unique donors
bkpt <- bkpt[, c("chr", "chr_bkpt_beg", "chr_bkpt_end", "cancer", "icgc_donor_id")]
bkpt <- unique(bkpt)
bkpt["chr"] <- apply(bkpt, 1, function(x){as.character(x["chr"])})
bkpt <- bkpt %>% group_by(chr, chr_bkpt_beg, chr_bkpt_end, cancer) %>%
  summarize(n_bkpt = n())
bkpt["chr_bkpt_beg_s"] <- bkpt["chr_bkpt_beg"]
bkpt["chr_bkpt_beg_e"] <- bkpt["chr_bkpt_beg"]
bkpt <- bkpt[bkpt["chr"] != "Y", ]
setDT(bkpt)

n_bkpt <- bkpt %>% 
  group_by(cancer, chr) %>%
  summarize(n_bkpt_total = n())




####### READ REGION DATA AND CALCUCATE INTERSECTION STATS
labels_dir <- "../data/adhoc/regions/input/"
all_labels <- list.files(labels_dir)

for (label in all_labels){
  path <- paste0(labels_dir, label)
  print(path)
  region_df <- load_region(path)
  print(head(region_df))
  label <- gsub(pattern = "([0-9]*)$", replacement = "", x = label)
  df_res <- count_bkpt(region_df=region_df, bkpt_df=bkpt, region_name=label)
  n_bkpt <- n_bkpt %>%
    left_join(df_res, by=c("cancer", "chr"))
}

n_bkpt[is.na(n_bkpt)] <- 0

write.csv(n_bkpt, file = "../data/adhoc/regions/results.csv", row.names = FALSE)




######## PLOTS

n_bkpt <- read.csv("../data/adhoc/regions/results.csv")

# by cancer type and region

cancer_region_df <- n_bkpt %>% 
  group_by(cancer) %>%
  summarize(
    n_bkpt_intersected_3UTR = sum(n_bkpt_intersected_3UTR),
    n_bkpt_intersected_5UTR = sum(n_bkpt_intersected_5UTR),
    n_bkpt_intersected_codingExons = sum(n_bkpt_intersected_codingExons),
    n_bkpt_intersected_downstream = sum(n_bkpt_intersected_downstream),
    n_bkpt_intersected_introns = sum(n_bkpt_intersected_introns),
    n_bkpt_intersected_promoters = sum(n_bkpt_intersected_promoters),
    n_bkpt_intersected_WholeGenes = sum(n_bkpt_intersected_WholeGenes)
    ) %>%
  mutate(n_bkpt_total = n_bkpt_intersected_3UTR + 
           n_bkpt_intersected_5UTR + 
           n_bkpt_intersected_codingExons + 
           n_bkpt_intersected_downstream + 
           n_bkpt_intersected_introns + 
           n_bkpt_intersected_promoters + 
           n_bkpt_intersected_WholeGenes
           )

cols <- names(cancer_region_df)
cols <- cols[!cols %in% c("cancer", "n_bkpt_total")]
for (col in cols){
  new_name <- paste0("ratio_", col)
  cancer_region_df[new_name] <- cancer_region_df[col] / cancer_region_df["n_bkpt_total"]
}

ratio_cols <- names(cancer_region_df)[!startsWith(names(cancer_region_df), "n_bkpt")]
gg_cancer_region_df <- melt(cancer_region_df[ratio_cols], id.vars = "cancer", 
     variable.name = "region", value.name = "ratio")
gg_cancer_region_df$region <- gsub(pattern = "ratio_n_bkpt_intersected_",
                                   replacement = "", x = gg_cancer_region_df$region)

g1 <- ggplot(gg_cancer_region_df,
             aes(x=cancer, y=ratio, fill=region)) +
  geom_bar(stat="identity", position="dodge") + 
  theme(legend.position = "bottom") + 
  xlab("Cancer type") +
  ylab("Ratio of breakpoints in the \n region  from total number \n of overlaps")+
  scale_fill_discrete(name="Region")
g1



# by chromosome and region

chr_region_df <- n_bkpt %>% 
  group_by(chr) %>%
  summarize(
    n_bkpt_intersected_3UTR = sum(n_bkpt_intersected_3UTR),
    n_bkpt_intersected_5UTR = sum(n_bkpt_intersected_5UTR),
    n_bkpt_intersected_codingExons = sum(n_bkpt_intersected_codingExons),
    n_bkpt_intersected_downstream = sum(n_bkpt_intersected_downstream),
    n_bkpt_intersected_introns = sum(n_bkpt_intersected_introns),
    n_bkpt_intersected_promoters = sum(n_bkpt_intersected_promoters),
    n_bkpt_intersected_WholeGenes = sum(n_bkpt_intersected_WholeGenes)
  ) %>%
  mutate(n_bkpt_total = n_bkpt_intersected_3UTR + 
           n_bkpt_intersected_5UTR + 
           n_bkpt_intersected_codingExons + 
           n_bkpt_intersected_downstream + 
           n_bkpt_intersected_introns + 
           n_bkpt_intersected_promoters + 
           n_bkpt_intersected_WholeGenes
  )

cols <- names(chr_region_df)
cols <- cols[!cols %in% c("chr", "n_bkpt_total")]
for (col in cols){
  new_name <- paste0("ratio_", col)
  chr_region_df[new_name] <- chr_region_df[col] / chr_region_df["n_bkpt_total"]
}

ratio_cols <- names(chr_region_df)[!startsWith(names(chr_region_df), "n_bkpt")]
gg_chr_region_df <- melt(chr_region_df[ratio_cols], id.vars = "chr", 
                            variable.name = "region", value.name = "ratio")
gg_chr_region_df$region <- gsub(pattern = "ratio_n_bkpt_intersected_",
                                   replacement = "", x = gg_chr_region_df$region)

g2 <- ggplot(gg_chr_region_df,
             aes(x=chr, y=ratio, fill=region)) +
  geom_bar(stat="identity", position="dodge") + 
  theme(legend.position = "bottom") + 
  xlab("Chromosome")  +
  ylab("Ratio of breakpoints in the \n region  from total number \n of overlaps")+
  scale_fill_discrete(name="Region")
g2


# by region

agg_region_df <- n_bkpt %>% ungroup %>%
  select(n_bkpt_intersected_3UTR, n_bkpt_intersected_5UTR, n_bkpt_intersected_codingExons, 
         n_bkpt_intersected_downstream, n_bkpt_intersected_introns, 
         n_bkpt_intersected_promoters, n_bkpt_intersected_WholeGenes) %>%
  summarise_all(sum)

gg_agg_region_df <- melt(agg_region_df, variable.name = "region", value.name = "number")
gg_agg_region_df$region <- gsub(pattern = "n_bkpt_intersected_",
                                replacement = "", x = gg_agg_region_df$region)
sum_bkpt <- sum(gg_agg_region_df$number)
gg_agg_region_df$ratio <- round(gg_agg_region_df$number/sum_bkpt*100, 3)

g3 <- ggplot(gg_agg_region_df,
             aes(x=reorder(region, number, mean), y=number)) +
  geom_bar(stat="identity", position="dodge") + 
  theme(legend.position = "bottom") + 
  xlab("Region")  +
  ylab("Ratio of breakpoints in the region \n  from total number of overlaps")+
  geom_text(aes(y = number + 10000,    
                label = paste0(ratio, '%')),    
            position = position_dodge(width = .9), 
            size = 4)
g3  





ggarrange(g1, g3, g2,  
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)