library(dplyr)
setwd("../data/secondary/QUADRUPLEXES")



# 1 строка на 1 квадруплекс

quadr<-read.table("all_chr_q.bed")
nrow(quadr)
names(quadr)<-c("chr","from","to","seq")
quadr$length_q<-quadr$to - quadr$from
quadr$start_w<-round(quadr$from/10000)*10000
quadr$end_w<-quadr$start_w+10000
quadr$chr<-gsub("chr","",quadr$chr)
quadr<-quadr[quadr$chr!="Y",]


#total q
quadr_n<-quadr %>% group_by(chr,start_w) %>%
  summarize(length_q = sum(length_q),
            count_q = n())
quadr_n$coverage_q<-quadr_n$length_q/10000
names(quadr_n)[2]<-"from"


write.csv(quadr_n,file="..\\data\\preprocessed\\quadruplexes\\quadr_cov.csv")
# не во всех окнах есть квадруплекс, так что потом будет делать left join
