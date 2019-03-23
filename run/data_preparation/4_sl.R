library(dplyr)
setwd("E:/Учеба/Диплом/Diploma/data/secondary/STEMLOOPS")

# 1 строка на 1 стемлуп

str3<-read.table("sorted_all_chr_S6-15.bed")
names(str3)<-c("chr","from","to")
str3$length_sl_6_15<-str3$to - str3$from
str3$start_w<-round(str3$from/10000)*10000
str3$end_w<-str3$start_w+10000

str3_new<-str3 %>% group_by(chr,start_w) %>%
  summarize(length_sl_6_15 = sum(length_sl_6_15),
            count_sl_6_15 = n())
str3_new$coverage_sl_6_15<-str3_new$length_sl_6_15/10000
names(str3_new)[2]<-"from"
# не во всех окнах есть такие стемлупы, так что потом будет делать left join
rm(str3)
gc()


# 1 строка на 1 окно

str1<-read.table("cov_10_S16-50.bed", header = FALSE)
names(str1)<-c("chr","from","to","count_sl_16_50","length_sl_16_50","length_window_16_50","coverage_16_50")

str2<-read.table("cov_10_S15-30.bed", header = FALSE)
names(str2)<-c("chr","from","to","count_sl_15_30","length_sl_15_30","length_window_15_30","coverage_15_30")

sec_str<-merge(str1,str2,by=c("chr","from","to"))
#sec_str<-merge(sec_str,str3,by=c("chr","from","to"))

# уберем лишние хромосомы
chrs<-levels(sec_str$chr)
chrs<-chrs[-grep("chrUn",chrs)]
chrs<-chrs[-grep("_",chrs)]
chrs<-chrs[!chrs %in% c("chrM","chrY")]

sec_str<-sec_str[sec_str$chr %in% chrs,]

sec_str$length_window_16_50<-NULL
sec_str$length_window_15_30<-NULL

#left
sec_str_all<-merge(sec_str,str3_new,by=c("chr","from"),all.x = TRUE)
sec_str_all$count_sl_6_15[is.na(sec_str_all$count_sl_6_15)]<-0
sec_str_all$length_sl_6_15[is.na(sec_str_all$length_sl_6_15)]<-0
sec_str_all$coverage_sl_6_15[is.na(sec_str_all$coverage_sl_6_15)]<-0

rm(str1)
rm(str2)
rm(str3_new)
rm(sec_str)
gc()


sec_str_all$chr<-as.character(sec_str_all$chr)
sec_str_all$chr<-gsub("chr","",sec_str_all$chr)

write.csv(sec_str_all,file="stemloops_density_10k.csv",row.names = FALSE)
