########################## EDA

library(ggplot2)
library(dplyr)
library(reshape2)
library(reshape)


## 1

# боксплот по хромосомам положения структур по типам 
# (смотрим, такие же ли отклонения от серединки)
setwd("E:/Учеба/Диплом/Diploma/data/secondary/STEMLOOPS")
s15_30<-read.table("sorted_all_chr_S15-30.bed")
s15_30$str_mb<-floor(s15_30$V2/1000000)
small_dat<-s15_30[seq(1,23757318,10),]

g<-ggplot(small_dat,aes(x=V1,y=str_mb))+
  facet_wrap(~V1,ncol = 6,scales="free")+
  geom_boxplot()+
  ggtitle("Stemloops (S15-30) locations boxplot in each chromosome (in Mb)")+
  xlab("Chromosome")+ylab("Mb")

g

chr_base<-as.data.frame(cbind(c(as.character(seq(1,22,1)),"X","Y"),
                              c(249250621,243199373,198022430,191154276,
                                180915260,171115067,159138663,146364022,
                                141213431,135534747,135006516,133851895,
                                115169878,107349540,102531392,90354753,
                                81195210,78077248,59128983,63025520,
                                48129895,51304566,155270560,59373566)
))
names(chr_base)<-c("chr","length")
chr_base$length<-as.numeric(as.character(chr_base$length))
chr_base$mb<-floor(chr_base$length/1000000)



data_mb <- s15_30
data_stat <- data_mb %>% group_by(V1) %>% summarize(med=median(str_mb))
names(data_stat)[1]<-"chr"
data_stat$chr<-gsub("chr","",data_stat$chr)
data_stat <-merge(data_stat,chr_base[,c("chr","mb")],by="chr")
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
  ggtitle("Deviation from center of chromosome (normalized by length of chromosome).S15-30")+
  ylab("Deviation")



## 2
# отрисовать линейный график плотности
data <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\structural_mutation_final_10_kb.csv")
data <- data[, c("chr", "from", "overall_density", "density_breast", "density_pancreatic")]
data$chr <- as.character(data$chr)
names(data)[3:5]<-c("breakpoints (general)", "breakpoints (breast)", 
                    "breakpoints (pancreatic)")

sec_str <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\STEMLOOPS\\stemloops_density_10k.csv")
sec_str <- sec_str[, c("chr", "from", "coverage_15_30", "coverage_sl_6_15", "coverage_16_50")]
sec_str$chr <- as.character(sec_str$chr)
sec_str$chr <- gsub("chr", "", sec_str$chr)
names(sec_str)[3:5] <- c("medium stemloops", "short stemloops", "long stemloops")



all_data <- merge(sec_str, data, by=c("chr", "from"))
all_data_plot <- melt(all_data, id.vars = c("chr", "from"))

chr_data <- all_data_plot[all_data_plot$chr == "21", ]

ggplot(chr_data,aes(x=from/1000000,y=value,colour=variable))+
  geom_line()+
  facet_grid(variable ~.,scales = "free_y" )+
  ggtitle("Densities of breakpoints and coverage of stemloops across 21 chromosome (in Mb)")+
  xlab("Start position of window")+
  ylab("Density/coverage")



# добавим квадруплексы
sec_str<-read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\QUADRUPLEXES\\quadr_cov_10kb.csv")
sec_str<-sec_str[,c("chr","from","coverage_q")]
sec_str$chr<-as.character(sec_str$chr)
sec_str$chr<-gsub("chr","",sec_str$chr)
names(sec_str)[3]<-c("quadruplexes")



all_data1 <- merge(sec_str, all_data, by=c("chr", "from"))
all_data_plot <- melt(all_data1, id.vars = c("chr", "from"))

chr_data <- all_data_plot[all_data_plot$chr == "21", ]

ggplot(chr_data, aes(x=from/1000000, y=value, colour=variable)) +
  geom_line() +
  facet_grid(variable ~., scales="free_y" ) +
  ggtitle("Densities of breakpoints and coverage of stemloops and quadruplexes across 21 chromosome (in Mb)")+
  xlab("Start position of window")+
  ylab("Density/coverage")






### CORRELATIONS

mutations <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\structural_mutation_final_1_mb.csv")
names(mutations)[which(names(mutations)=="new_wind")]<-"window"

mutations <- mutations[, c("chr", "window", "from", "density_blood", "density_bone",
                           "density_brain", "density_breast", "density_liver",
                           "density_ovary", "density_pancreatic", "density_prostate",
                           "density_skin", "density_uterus", "overall_density" )]

mutations$chr <- as.character(mutations$chr)

# steml <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\STEMLOOPS\\stemloops_density_10k.csv")
# steml <- steml[, c( "chr", "from", "coverage_16_50", "coverage_15_30", "coverage_sl_6_15")]
steml <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\STEMLOOPS\\stemloops_density_1mb.csv")
steml<-steml[,c( "chr" ,"from","coverage_sl_16_50","coverage_sl_15_30","coverage_sl_6_15")]
names(steml)[3:5] <- c("long stem-loops", "medium stem-loops", "short stem-loops")

steml$chr <- as.character(steml$chr)
steml$chr <- gsub("chr", "", steml$chr)

data <- merge(mutations, steml, by=c("chr", "from"))
chr_names <- unique(data$chr)

corr_m <- data.frame(cancer=0, stemloops=0, correlation=0, chr=0)

mut_names <- c("density_blood", "density_bone",
               "density_brain", "density_breast", "density_liver",
               "density_ovary", "density_pancreatic", "density_prostate",
               "density_skin", "density_uterus", "overall_density")
steml_names <- c("long stem-loops", "medium stem-loops", "short stem-loops")

for (i in chr_names){
  chr_data <- data[data$chr == i, ]
  m <- cor(chr_data[, mut_names], chr_data[, steml_names], method="spearman")
  m <- as.data.frame(m)
  m$var1 <- rownames(m)
  rownames(m) <- NULL
  a <- melt(m, id.vars = "var1")
  a$chr <- i
  names(a) <- c("cancer", "stemloops", "correlation", "chr")
  corr_m <- rbind(corr_m, a)
} 

corr_m <- corr_m[-1, ]

corr_m$cancer <- gsub("density_", "", corr_m$cancer)
corr_m$cancer <- gsub("_density", "", corr_m$cancer)
corr_m$stemloops <- gsub("density_", "", corr_m$stemloops)

ggplot(corr_m, aes(correlation))+
  geom_freqpoly(bins = 20)+
  facet_grid(stemloops ~ cancer)
# ggtitle("Histogram of chromosomes correlations between densities of cancer profiles and coverages of stemloops")


ggplot(corr_m, aes(x=gsub("stem-loops", "", stemloops), y=correlation)) +
  geom_boxplot() +
  facet_wrap( ~ cancer) +
  xlab("Stemloops type")