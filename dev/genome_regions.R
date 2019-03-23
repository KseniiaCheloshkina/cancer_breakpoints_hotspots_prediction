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
  region_df['chr'] <- apply(region_df, 1, function(x){as.character(x["chr"])})
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
    bkpt_df_part <- bkpt_df[chr, on="chr"]
    region_df_part <- region_df[chr, on="chr"]
    
    
    joined_df_s <- region_df_part[bkpt_df_part, 
                                on=.(start <= chr_bkpt_beg_s, end >= chr_bkpt_beg_s),
                                allow.cartesian=TRUE]
    joined_df_s$chr_bkpt_beg_e <- NULL
    
    joined_df_e <- region_df_part[bkpt_df_part, 
                                on=.(start <= chr_bkpt_beg_e, end >= chr_bkpt_beg_e),
                                allow.cartesian=TRUE]
    joined_df_e$chr_bkpt_beg_s <- NULL
    joined_df <- rbind(joined_df_s, joined_df_e)
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

# function to calculate number of genes inside hotspots 
count_hsp <- function(region_df, bkpt_df){
  
  chr_list <- unique(bkpt_df[, chr])
  merged_data <- data.frame()
  
  for (chr in chr_list){
    bkpt_df_part <- bkpt_df[chr, on="chr"]
    region_df_part <- region_df[chr, on="chr"]
    
    
    joined_df_s <- region_df_part[bkpt_df_part, 
                                on=.(start >= chr_bkpt_beg_s, start <= chr_bkpt_beg_e),
                                allow.cartesian=TRUE] %>%
      select(from, chr, to, hsp_nm)

    joined_df_e <- region_df_part[bkpt_df_part, 
                                on=.(end >= chr_bkpt_beg_s, end <= chr_bkpt_beg_e),
                                allow.cartesian=TRUE] %>%
      select(from, chr, to, hsp_nm)
    
    joined_df <- rbind(joined_df_s, joined_df_e)
    # unique hits
    joined_df <- unique(joined_df)
    
    final_df <- joined_df %>%
      filter(!is.na(chr)) %>%
      group_by(chr) %>%
      summarize(n_hsp = n())
    
    null_df <- data.frame("chr" = chr, "n_hsp" = 0, stringsAsFactors = FALSE)
    
    if (nrow(final_df) == 0){
      final_df <- null_df
    } 
    final_df$n_hsp_chr <- nrow(bkpt_df_part)
    merged_data <- rbind(merged_data, data.frame(final_df))
    
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



gg_chr_region_df$chr <- as.character(gg_chr_region_df$chr)
chr_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
               "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
               "21", "22", "X")

g2 <- ggplot(gg_chr_region_df,
             aes(x=chr, y=ratio, fill=region)) +
  geom_bar(stat="identity", position="dodge") + 
  theme(legend.position = "bottom") + 
  xlab("Chromosome")  +
  ylab("Ratio of breakpoints in the \n region  from total number \n of overlaps")+
  scale_fill_discrete(name="Region")+
  scale_x_discrete(limits=chr_order)
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
          ncol = 2, nrow = 2, hjust = -3, vjust = 1)




###### HOTSPOTS

## LOAD 10Kb HOTSPOTS DATA

cancer_bkpt_path <- "../data/preprocessed/breakpoints/structural_mutation_final_10_kb.csv"
bkpt <- read.csv(cancer_bkpt_path)
hsp_names <- names(bkpt)[grep("hotspot", names(bkpt))]
needed_names <- c(c("chr", "window", "from", "to"), hsp_names)
new_bkpt <- bkpt[needed_names]
rm(bkpt)
gc()

new_bkpt["chr"] <- apply(new_bkpt, 1, function(x){as.character(x["chr"])})
new_bkpt["chr_bkpt_beg_s"] <- new_bkpt["from"]
new_bkpt["chr_bkpt_beg_e"] <- new_bkpt["to"]


## READ WHOLE GENES DATA AND CALCUCATE INTERSECTION STATS
labels_path <- "../data/adhoc/regions/input/"
labels <- list.files(labels_path)
labels <- labels[!labels %in% c("promoters020219", "downstream020219")]

all_res <- data.frame()

for (label in labels){
  
  # load region
  path <- paste0(labels_path, label)
  region_df <- load_region(path)
  clean_label <- gsub("020219", "", label)
  print(clean_label)
  
  for (hsp_nm in hsp_names){
    print(hsp_nm)
    # select cols
    needed_names <- c("chr", "from", "to", "chr_bkpt_beg_e", "chr_bkpt_beg_s")
    bkpt_df <- new_bkpt[c(needed_names, hsp_nm)]
    # filter hotspots
    bkpt_df <- bkpt_df[bkpt_df[hsp_nm] == 1, ]
    n_bkpt_df <- nrow(bkpt_df)
    setDT(bkpt_df)
    # count
    df_res <- count_hsp(region_df, bkpt_df)
    df_res$hsp_type <- hsp_nm
    df_res$n_total <- n_bkpt_df
    df_res$label <- clean_label
    
    all_res <- rbind(all_res, data.frame(df_res))
  }
  
}


write.csv(all_res, file = "../data/adhoc/regions/results_regions.csv", row.names = FALSE)






## STATS AND PLOTS


all_res <- read.csv("../data/adhoc/regions/results_regions.csv")

# GENES

genes_results <- all_res[all_res$label == "WholeGenes", ]

cancer_pattern <- "^[a-z]*"
labeling_pattern <- "([0-9]?[.]?[0-9]*)$"

m <- regexpr(cancer_pattern, genes_results$hsp_type)
genes_results$cancer_type <- regmatches(genes_results$hsp_type, m)
m <- regexpr(labeling_pattern, genes_results$hsp_type)
genes_results$labeling_type <- regmatches(genes_results$hsp_type, m) 

genes_results$perc_chr <- genes_results$n_hsp / genes_results$n_hsp_chr
genes_results$perc_total <- genes_results$n_hsp / genes_results$n_total

genes_results$label <- NULL

# by cancer and labeling type
genes_results_c_l <- genes_results %>% 
  group_by(cancer_type, labeling_type, n_total) %>% 
  summarize(n_hsp = sum(n_hsp)) %>%
  mutate(perc = n_hsp/n_total)

g_c_l <- ggplot(genes_results_c_l, aes(cancer_type, labeling_type)) + 
  geom_tile(aes(fill = perc), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue", name="Percentage") + 
  xlab("Cancer type") + 
  ylab("Labeling type")

g_c_l


# by chr and labeling type
genes_results_chr_l <- genes_results %>% 
  group_by(chr, labeling_type, n_total) %>% 
  summarize(n_hsp = sum(n_hsp)) %>%
  mutate(perc = n_hsp/n_total)

genes_results_chr_l$chr <- as.character(genes_results_chr_l$chr)
chr_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
               "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
               "21", "22", "X")

g_chr_l <- ggplot(genes_results_chr_l, aes(chr, labeling_type)) + 
  geom_tile(aes(fill = perc), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue", name="Percentage") + 
  xlab("Chromosome") + 
  ylab("Labeling type")+
  scale_x_discrete(limits=chr_order)

g_chr_l


# by cancer and chromosome for each labeling type
all_g <- ggplot(genes_results, 
                aes(x = labeling_type, y = perc_total, fill = chr)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  facet_grid(~ cancer_type) +
  xlab("Labeling type") + 
  ylab("Percentage") +theme(axis.text.x = element_text(size=10, angle=45)) +
  scale_fill_discrete(name="Chromosome", labels = chr_order)

all_g  

# 
  
ggarrange(all_g, 
          ggarrange(g_c_l, g_chr_l, labels = c("B", "C"), ncol = 2, widths = 4:3),
          nrow = 2, labels = "A" , heights = 4:3)

## other regions

 
m <- regexpr(cancer_pattern, all_res$hsp_type)
all_res$cancer_type <- regmatches(all_res$hsp_type, m)
m <- regexpr(labeling_pattern, all_res$hsp_type)
all_res$labeling_type <- regmatches(all_res$hsp_type, m) 

all_res$perc_chr <- all_res$n_hsp / all_res$n_hsp_chr
all_res$perc_total <- all_res$n_hsp / all_res$n_total

# by cancer_type and region (point - chromosome in dataset)
g1 <- ggplot(all_res, aes(x=cancer_type, y=perc_total, fill=label))+
  geom_boxplot()+
  xlab("Cancer type")+
  scale_fill_discrete(name="Region")+
  ylab("Percentage of \n intersected hotspots from \n total number of hotspots")

all_res$chr <- as.character(all_res$chr)

# by chromosome and region (point -  cancer_type in dataset)
g2 <- ggplot(all_res, aes(x=chr, y=perc_total, fill=label))+
  geom_boxplot()+
  xlab("Chromosome")+
  scale_fill_discrete(name="Region")+
  ylab("Percentage of \n intersected hotspots from \n total number of hotspots")+
  scale_x_discrete(limits=chr_order)

ggarrange(g1, g2, labels = c("A", "B"), nrow = 2)
