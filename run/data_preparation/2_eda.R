##### EXPLORATORY DATA ANALYSIS OF BREAKPOINTS


# libraries
library(reshape2)
library(ggplot2)
library(dplyr)

# data
input_folder <- "..//data//raw breakpoints//structural mutation//1_read//"
setwd(input_folder)
all_data <- read.csv("all_cancer_data.csv")
all_data$X <- NULL

# RANGE 
# replace na range with 0
all_data$chr_from_range[is.na(all_data$chr_from_range)] <- 0
all_data$chr_to_range[is.na(all_data$chr_to_range)] <- 0

factor_cols <- c("chr_from","chr_to")
numeric_cols <- c("chr_from_bkpt", "chr_from_range", "chr_to_bkpt","chr_to_range")


sapply(all_data[, factor_cols], levels)
all_data[all_data$chr_from == "MT", ]
# исключим MT хромосому, т.к. 1 разрыв всего


# distributions
dat <- melt(all_data[, numeric_cols])

ggplot(dat[seq(1, 1019636, 100), ], aes(value))+
  geom_histogram()+
  facet_wrap(~variable, scales="free")


ggplot(dat[seq(1, 1019636, 100), ],aes(x=1, y=value))+
  geom_boxplot()+
  facet_wrap(~variable, scales="free")

### RANGE

summary(all_data$chr_to_range)
quantile(all_data$chr_to_range,seq(0.9, 1, 0.01))
quantile(all_data$chr_from_range,seq(0.8, 1, 0.025))

# we will exclude rows with range>10 but later
# let's gather some statistics

nrow(all_data)
# 254908

bad_range <- which(all_data$chr_from_range > 10)
bad_range <- append(bad_range, which(all_data$chr_to_range > 10))
bad_range <- unique(bad_range)
length(bad_range)
# 14469

# bad ranges vs breakpoint type
var_big_range <- as.data.frame(summary(all_data$variant_type[(all_data$chr_to_range > 10)|(all_data$chr_from_range > 10)]))
var_big_range$variant_type <- rownames(var_big_range)
rownames(var_big_range) <- NULL
names(var_big_range)[1] <- 'count'


############################################ transform data for density calculation
all_data<-subset(all_data,select=c("icgc_donor_id","icgc_sample_id","variant_type",
                                   "chr_from","chr_from_bkpt","chr_from_range",
                                   "chr_to","chr_to_bkpt","chr_to_range"))

from_data <- all_data[, c("icgc_donor_id", "icgc_sample_id", "variant_type", "chr_from", "chr_from_bkpt", "chr_from_range")]
to_data <- all_data[, c("icgc_donor_id", "icgc_sample_id", "variant_type", "chr_to", "chr_to_bkpt", "chr_to_range")]
names(from_data) <- c("icgc_donor_id", "icgc_sample_id", "variant_type", "chr", "chr_bkpt", "chr_range")
names(to_data) <- c("icgc_donor_id", "icgc_sample_id", "variant_type", "chr", "chr_bkpt", "chr_range")
data<-rbind(from_data, to_data)
nrow(data)
# 509818



# exclude rows with range > 10
bad_range <- which(data$chr_range > 10)
length(bad_range)
# 22391
data <- data[-bad_range, ]
nrow(data)
# 487427

# exclude MT
mt <- which(data$chr == "MT")
data <- data[-mt, ]
nrow(data)
# 487 425

# expand breakpoint interval with size of range
data$chr_bkpt_beg <- data$chr_bkpt - data$chr_range
data$chr_bkpt_end <- data$chr_bkpt + data$chr_range



################################################ STATS
# gather some stats
length(unique(data$icgc_donor_id))
# 2092
length(unique(data$icgc_sample_id))
# 2234

variants <- data %>% group_by(variant_type) %>% summarize(count=n())

ggplot(variants, aes(x=reorder(variant_type, count),
                     y=count,
                     fill=count,
                     show.legend=FALSE)) +
  geom_bar(stat="identity") +
  coord_flip() +
  geom_text(aes(x=variant_type, y=count, label=count),
          hjust=0.1, 
          size=3,
          color=rgb(100, 100, 100, maxColorValue=255)) +
  theme(legend.position="none") +
  ggtitle("Number of breakpoints by variant type")


chr_s<-data %>% group_by(chr) %>% summarize(count= n())
chr_base<-as.data.frame(cbind(c(as.character(seq(1,22,1)),"X","Y"),
                              c(249250621,243199373,198022430,191154276,
                 180915260,171115067,159138663,146364022,
                 141213431,135534747,135006516,133851895,
                 115169878,107349540,102531392,90354753,
                 81195210,78077248,59128983,63025520,
                 48129895,51304566,155270560,59373566)
               ))
names(chr_base)<-c("chr","length")
chr_s<-merge(chr_s,chr_base,by="chr")
chr_s$length<-as.numeric(as.character(chr_s$length))
chr_s$norm_count<-chr_s$count/chr_s$length
chr_s$norm_count<-format(chr_s$norm_count,digits = 3,scientific = FALSE)
chr_s$norm_count<-as.numeric(chr_s$norm_count)

ggplot(chr_s,aes(x=reorder(chr, count),y=count,
                    fill=count,
                    show.legend=FALSE))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme(legend.position="none")+
  ggtitle("Number of breakpoints by chromosome")

chr_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
               "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
               "21", "22", "X", "Y")
chr_s$chr <- as.character(chr_s$chr)
nrm <- ggplot(chr_s,aes(x=reorder(chr, norm_count),y=norm_count,
                 fill=norm_count,
                 show.legend=FALSE))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme(legend.position="none")+
  xlab("Chromosome")+
  ylab("Normalized number of breakpoints")
nrm
save(nrm, file = "..\\data\\1_D.RData")
chr_s$mb<-floor(chr_s$length/1000000)

# breakpoints location by chromosome
# reorder(chr,chr_bkpt, FUN=median)
data_mb<-data
data_mb$chr_bkpt_mb<-floor(data_mb$chr_bkpt/1000000)

ggplot(data_mb,aes(x=paste("chr",as.factor(chr)),y=chr_bkpt_mb))+
  facet_wrap(~chr,ncol = 6,scales="free")+
  geom_boxplot()+
  ggtitle("Breakpoints locations boxplot in each chromosome (in Mb)")+
  xlab("Chromosome")+ylab("Mb")



data_stat <- data_mb %>% group_by(chr) %>% summarize(med=median(chr_bkpt_mb))
data_stat <-merge(data_stat,chr_s[,c("chr","mb")],by="chr")
data_stat$natural_med<-data_stat$mb/2
data_stat$normalized_diff<-(data_stat$med-data_stat$natural_med)/data_stat$mb
data_stat_p<-data_stat %>% 
  arrange(-normalized_diff) %>%
  mutate(chr = factor(chr, chr)) 

ggplot(data_stat_p,aes(x=chr,y=normalized_diff,fill=chr))+
  geom_bar(stat="identity")+
  geom_text(aes(x=chr,y=normalized_diff, label=scales::percent(normalized_diff)),
            #vjust=0.00003, 
            nudge_x=0.02,nudge_y=0.01,
            size=3)+
  theme(legend.position="none")+
  ggtitle("Deviation from center of chromosome (normalized by length of chromosome) ")+
  ylab("Deviation")


# outlier in Y chromosome...
# distributions do not look like uniform
# check!
ggplot(data,aes(chr_bkpt))+
  facet_wrap(~chr,ncol = 6,scales="free")+
  geom_histogram()+
  ggtitle("Breakpoints locations histograms in each chromosome")




# write new csv
output_folder = "../2_eda_preprocessing/"
setwd(output_folder)
write.csv(data,row.names = FALSE,file="all_cancer_data_eda.csv")







# function for the same transformation OF other files in the future
eda_trans<-function(data){
  
  #replace na range with 0
  a<-c(which(is.na(data$chr_from_range)),which(is.na(data$chr_to_range)))
  if (length(a)!=0){
  data$chr_from_range[is.na(data$chr_from_range)]<-0
  data$chr_to_range[is.na(data$chr_to_range)]<-0
  }
  

  data<-subset(data,select=c("icgc_donor_id","icgc_sample_id","variant_type",
                             "chr_from","chr_from_bkpt","chr_from_range",
                             "chr_to","chr_to_bkpt","chr_to_range"))
  
  a<-data[,c("icgc_donor_id","icgc_sample_id","variant_type","chr_from","chr_from_bkpt","chr_from_range")]
  b<-data[,c("icgc_donor_id","icgc_sample_id","variant_type","chr_to","chr_to_bkpt","chr_to_range")]
  names(a)<-c("icgc_donor_id","icgc_sample_id","variant_type","chr","chr_bkpt","chr_range")
  names(b)<-c("icgc_donor_id","icgc_sample_id","variant_type","chr","chr_bkpt","chr_range")
  data<-rbind(a,b)
  
  # exclude rows with range>10
  bad_range<-which(data$chr_range>10)
  if (length(bad_range)!=0) {data<-data[-bad_range,]}
  
  # REMOVE MT
  mt<-which(data$chr=="MT")
  if (length(mt)!=0) {data<-data[-mt,]}
  
  
  # expand breakpoint interval with size of range
  data$chr_bkpt_beg<-data$chr_bkpt-data$chr_range
  data$chr_bkpt_end<-data$chr_bkpt+data$chr_range
  
  return(data)
}


input_folder <- "..//data//raw breakpoints//structural mutation//1_read//"
setwd(input_folder)
cancer_types<-list.files()
cancer_types<-cancer_types[grep("_all_data.csv",cancer_types)]


for (i in cancer_types){
  data<-read.csv(i)
  data<-eda_trans(data)
  filename<-paste("../2_eda_preprocessing/",i,sep="")
  write.csv(data,row.names = FALSE,file=filename)
}
