########################## EDA

library(ggplot2)
library(dplyr)
library(reshape2)


## 1

# боксплот по хромосомам положения структур по типам 
# (смотрим, такие же ли отклонения от серединки)
setwd("E://Учеба//Диплом//Diploma//data//secondary//QUADRUPLEXES//")
quadr<-read.table("all_chr_q.bed")
quadr$str_mb<-floor(quadr$V2/1000000)
small_dat<-quadr


g<-ggplot(small_dat,aes(x=V1,y=str_mb))+
  facet_wrap(~V1,ncol = 6,scales="free")+
  geom_boxplot()+
  ggtitle("Quadruplexes locations boxplot in each chromosome (in Mb)")+
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


data_mb <- quadr
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
  ggtitle("Deviation from center of chromosome (normalized by length of chromosome). Quadruplexes")+
  ylab("Deviation")



## 2
# отрисовать линейный график плотности
data<-read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\structural_mutation_final_10_kb.csv")
data<-data[,c("chr","from","overall_density","density_breast","density_pancreatic")]
data$chr<-as.character(data$chr)
names(data)[3:5]<-c("breakpoints (general)","breakpoints (breast)","breakpoints (pancreatic)")

sec_str<-read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\QUADRUPLEXES\\quadr_cov.csv")
sec_str<-sec_str[,c("chr","from","coverage_q")]
sec_str$chr<-as.character(sec_str$chr)
sec_str$chr<-gsub("chr","",sec_str$chr)
names(sec_str)[3]<-c("Quadr")

all_data<-merge(data,sec_str,by=c("chr","from"),all.x = TRUE)
all_data$Quadr[is.na(all_data$Quadr)]<-0

all_data_plot<-melt(all_data,id.vars = c("chr","from"))

chr_data<-all_data_plot[all_data_plot$chr=="21",]

ggplot(chr_data,aes(x=from/1000000,y=value,colour=variable))+
  geom_line()+
  facet_grid(variable ~.,scales = "free_y" )+
  #facet_grid(. ~ variable)+
  ggtitle(" Densities of breakpoints and coverage of quadruplexes across 21 chromosome (in Mb)")+
  xlab("Start position of window")+
  ylab("Density/coverage")




library(ggplot2)
library(reshape2)
library(reshape)

### CORRELATIONS


# mutations <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\structural_mutation_final_10_kb.csv")
mutations <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\structural_mutation_final_1_mb.csv")
names(mutations)[which(names(mutations)=="new_wind")]<-"window"

mutations <- mutations[, c("chr", "window", "from", "density_blood", "density_bone",
                           "density_brain", "density_breast", "density_liver",
                           "density_ovary", "density_pancreatic", "density_prostate",
                           "density_skin", "density_uterus", "overall_density" )]

mutations$chr <- as.character(mutations$chr)


# sec_str <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\QUADRUPLEXES\\quadr_cov_10kb.csv")
sec_str <- read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\QUADRUPLEXES\\quadr_cov_1mb.csv")
sec_str <- sec_str[, c("chr", "from", "coverage_q")]
sec_str$chr <- as.character(sec_str$chr)
sec_str$chr <- gsub("chr", "", sec_str$chr)
names(sec_str)[3]<-c("quadruplexes")


all_data <- merge(mutations, sec_str, by=c("chr", "from"), all.x = TRUE)
all_data$quadruplexes[is.na(all_data$quadruplexes)] <- 0


data <- all_data
chr_names <- unique(data$chr)

corr_m <- data.frame(cancer=0, quadruplexes=0, correlation=0, chr=0)

mut_names <- c("density_blood", "density_bone",
               "density_brain", "density_breast", "density_liver",
               "density_ovary", "density_pancreatic", "density_prostate",
               "density_skin", "density_uterus", "overall_density")
q_names<-c("quadruplexes")

for (i in chr_names){
  chr_data <- data[data$chr == i, ]
  m <- cor(chr_data[, mut_names], chr_data[, q_names], method="spearman")
  m <- as.data.frame(m)
  m$var1 <- rownames(m)
  rownames(m) <- NULL
  a <- melt(m, id.vars = "var1")
  a$chr <- i
  names(a) <- c("cancer", "quadruplexes", "correlation", "chr")
  corr_m <- rbind(corr_m, a)
} 

corr_m <- corr_m[-1, ]

corr_m$cancer <- gsub("density_", "", corr_m$cancer)
corr_m$cancer <- gsub("_density", "", corr_m$cancer)


ggplot(corr_m, aes(correlation)) +
  geom_freqpoly(bins = 20) +
  facet_wrap( ~ cancer, ncol=4)
# ggtitle("Histogram of chromosomes correlations between densities of cancer 
# profiles and coverage of quadruplexes")






ggplot(corr_m, aes(x=cancer, y=correlation)) +
  geom_boxplot()+
  xlab("Quadruplexes")+
  ggtitle("Boxplots of chromosomes correlations between densities of cancer profiles
          and coverages of quadruplexes (1 Mb)")