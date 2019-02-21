## libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(xlsx)
library(stringr)
require(data.table)

## set working folder
current_folder <- "E:\\Учеба\\Диплом\\Diploma\\scripts\\repo\\data\\adhoc\\cancers_cytoband\\"
setwd(current_folder)



# function to calculate number of genes inside hotspots 
count_hsp_markers <- function(region_df, bkpt_df){
  
  chr_list <- unique(bkpt_df[, chr])
  merged_data <- data.frame()
  
  for (chr in chr_list){
    bkpt_df_part <- bkpt_df[chr, on="chr"]
    region_df_part <- region_df[chr, on="start_chr"]
    
    start_joined_df_1 <- bkpt_df_part[region_df_part, 
                                      on=.( chr_bkpt_beg_s <= begin_start, 
                                            chr_bkpt_beg_e >= begin_start),
                                      allow.cartesian=TRUE] %>%
      select(id, Fusion.gene, start_chr, start_pos, end_chr, 
             end_pos, hsp_nm)
    
    start_joined_df_2 <- bkpt_df_part[region_df_part, 
                                      on=.( chr_bkpt_beg_s <= begin_end, 
                                            chr_bkpt_beg_e >= begin_end),
                                      allow.cartesian=TRUE]  %>%
      select(id, Fusion.gene, start_chr, start_pos, end_chr, 
             end_pos, hsp_nm)
    
    end_joined_df_1 <- bkpt_df_part[region_df_part, 
                                    on=.( chr_bkpt_beg_s <= end_start, 
                                          chr_bkpt_beg_e >= end_start),
                                    allow.cartesian=TRUE]   %>%
      select(id, Fusion.gene, start_chr, start_pos, end_chr, 
             end_pos, hsp_nm)
    
    end_joined_df_2 <- bkpt_df_part[region_df_part, 
                                    on=.( chr_bkpt_beg_s <= end_end, 
                                          chr_bkpt_beg_e >= end_end),
                                    allow.cartesian=TRUE]  %>%
      select(id, Fusion.gene, start_chr, start_pos, end_chr, 
             end_pos, hsp_nm)
    
    final_df <- rbind(start_joined_df_1, start_joined_df_2,
                      end_joined_df_1, end_joined_df_2) 
    final_df <- data.frame(final_df)
    final_df[is.na(final_df[hsp_nm]), hsp_nm] <- 0
    final_df <- final_df %>%
      group_by(id, Fusion.gene, start_chr, start_pos, end_chr, end_pos) %>%
      summarize_all(sum)
    
    final_df$n_hsp_chr <- nrow(bkpt_df_part)
    merged_data <- rbind(merged_data, data.frame(final_df))
    
  }
  
  names(merged_data)[7] <- paste0("n_hsp_", hsp_nm) 
  names(merged_data)[8] <- paste0("n_hsp_chr_", hsp_nm)
  
  return(merged_data)
  
}



### load markers

path <- "Tables S2-S5.xlsx"
data_s2 <- read.xlsx(path, sheetIndex = 1, header = TRUE, startRow = 2)
data_s3 <- read.xlsx(path, sheetIndex = 2, header = TRUE, startRow = 2)
data_s2 <- data_s2[, c("Translocation", "Fusion.gene")]
data_s3 <- data_s3[, c("Translocation", "Fusion.gene")]
data_tr <- rbind(data_s2, data_s3)
loc <- str_split(data_tr$Translocation, "_")
data_tr$start <- unlist(lapply(loc, function(x) x[1]))
data_tr$end <- unlist(lapply(loc, function(x) x[2]))
data_tr <- data_tr[!is.na(data_tr$Translocation), ]

chr_pattern <- "^[0-9]*"
m <- regexpr(chr_pattern, data_tr$start)
data_tr$start_chr <- regmatches(data_tr$start, m)

m <- regexpr(chr_pattern, data_tr$end)
data_tr$end_chr <- regmatches(data_tr$end, m)

pos_pattern <- "[a-z][0-9]*$"
m <- regexpr(pos_pattern, data_tr$start)
data_tr$start_pos <- regmatches(data_tr$start, m)

m <- regexpr(pos_pattern, data_tr$end)
data_tr$end_pos <- regmatches(data_tr$end, m)

data_tr <- data_tr[, c("Fusion.gene", "start_chr", "start_pos", "end_chr", "end_pos")]




### load labeling

path <- "cytoBand.txt"
data_lb <- read.table(path)
data_lb[, 1] <- gsub("chr", "", data_lb[, 1])
names(data_lb) <- c("chr", "start", "end", "subband", "type")

loc <- str_split(data_lb$subband, "[.]")
data_lb$band <- unlist(lapply(loc, function(x) x[1]))
data_lb_fin <- data_lb %>% group_by(chr, band) %>%
  summarize(start = min(start), end = max(end))


#### add cytoband labeling to data

data <- data_tr %>% 
  left_join(data_lb_fin, by = c("start_chr" = "chr", "start_pos" = "band")) %>%
  rename("begin_start" = "start", "begin_end" = "end") %>% 
  left_join(data_lb_fin, by = c("end_chr" = "chr", "end_pos" = "band")) %>%
  rename("end_start" = "start", "end_end" = "end")
data$id <- row.names(data)
setDT(data)


#### load hotspots data

path <- "../../preprocessed/breakpoints/structural_mutation_final_10_kb.csv"
hsp_data <- read.csv(path)
hsp_names <- names(hsp_data)[grep("hotspot_1", names(hsp_data))]
hsp_data <- hsp_data[, c("chr", "from", "to", hsp_names)]
hsp_data$chr <- as.character(hsp_data$chr)
hsp_data["chr_bkpt_beg_s"] <- hsp_data["from"]
hsp_data["chr_bkpt_beg_e"] <- hsp_data["to"]

#### find intersections - for each cancer type hotspots match translocation with hsp

for (i in 1:length(hsp_names)){
  hsp_nm <- hsp_names[i]
  print(hsp_nm)
  # select cols
  needed_names <- c("chr", "from", "to", "chr_bkpt_beg_e", "chr_bkpt_beg_s")
  bkpt_df <- hsp_data[c(needed_names, hsp_nm)]
  # filter hotspots
  bkpt_df <- bkpt_df[bkpt_df[hsp_nm] == 1, ]
  n_bkpt_df <- nrow(bkpt_df)
  setDT(bkpt_df)
  # count
  if (i == 1){
    df_res_all <- count_hsp_markers(data, bkpt_df)
  }
  else{
    df_res <- count_hsp_markers(data, bkpt_df)
    df_res_all <- df_res_all %>% left_join(df_res)
  }
  
}



##### STATS AND PLOTS
hsp_cols <- names(df_res_all)[grep("hotspot", names(df_res_all))]
hsp_cols <- hsp_cols[-grep("chr", hsp_cols)]
# хотспотов
df_res_all$total_hsp <- rowSums(df_res_all[, hsp_cols])
nrow(df_res_all[df_res_all$total_hsp>0,])

# типов рака
for (col in hsp_cols){
  new_col <- paste0(col, "_dummy")
  df_res_all[new_col] <- ifelse(df_res_all[col] > 0, 1, 0)
}

dummy_cols <- names(df_res_all)[grep("dummy", names(df_res_all))]
df_res_all$total_cancers <- rowSums(df_res_all[, dummy_cols])
table(df_res_all$total_cancers)

df_res_all$location <- paste0(df_res_all$start_chr, df_res_all$start_pos, "_", 
                              df_res_all$end_chr, df_res_all$end_pos)
df_res_all$full_location <- ifelse(!is.na(df_res_all$Fusion.gene), 
                                   paste0(
                                     df_res_all$location, " (", 
                                     df_res_all$Fusion.gene, ")"),
                                   df_res_all$location)
int_markers <- df_res_all[df_res_all$total_hsp > 0, ]
write.csv(int_markers,"../cancers_cytoband/results.csv")

int_markers <- melt(int_markers, id.vars="full_location", measure.vars=hsp_cols)
names(int_markers) <- c("translocation", "cancer_type", "n_hsp")
int_markers$cancer_type <- gsub("n_hsp_", "", int_markers$cancer_type)
int_markers$cancer_type <- gsub("_bkpt_hotspot_1", "", int_markers$cancer_type)

int_markers %>% filter(n_hsp > 0) 
# %>% group_by(cancer_type) %>% 
  # summarize(n_markers = n()) %>% arrange(n_markers)
sum(int_markers$n_hsp)

p <- ggplot(int_markers, aes(translocation, cancer_type)) + 
  geom_tile(aes(fill = n_hsp), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue", name="Number of hotspots") + 
  theme(axis.text.x = element_text(size=4, angle=45))+
  xlab("Translocation")+
  ylab("Cancer type")
p
