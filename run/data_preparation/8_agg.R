library(dplyr)
library(ggplot2)
library(reshape2)


### This script is used to aggregate 10 kb windows data to lower resolution (100 kb, 500 kb, etc.) 
### Firstly, aggregate breakpoints, then stemloops, then quadruplexes 


# aggregate to window of size
new_wind<-1000000

##  BREAKPOINTS

new_data<-read.csv("..\\data\\preprocessed\\breakpoints\\structural_mutation_final_10_Kb.csv")
new_data$X<-NULL
new_data<-new_data[order(new_data$chr,new_data$window),]
new_data$new_wind<-ceiling(new_data$to/new_wind)


dataset<-new_data %>% group_by(chr,new_wind) %>% 
  summarize(bkpt_in_window_blood = sum(bkpt_in_window_blood),
            bkpt_in_window_bone = sum(bkpt_in_window_bone),
            bkpt_in_window_brain = sum(bkpt_in_window_brain),
            bkpt_in_window_breast = sum(bkpt_in_window_breast),
            bkpt_in_window_liver = sum(bkpt_in_window_liver),
            bkpt_in_window_ovary = sum(bkpt_in_window_ovary),
            bkpt_in_window_pancreatic = sum(bkpt_in_window_pancreatic),
            bkpt_in_window_prostate = sum(bkpt_in_window_prostate),
            bkpt_in_window_skin = sum(bkpt_in_window_skin),
            bkpt_in_window_uterus = sum(bkpt_in_window_uterus),
            from = min(from),
            to = max(to)
  )

## DENSITY

cancer_cols<-grep("bkpt",names(dataset))

for (i in 1:length(cancer_cols)){
  s<-sum(dataset[,cancer_cols[i]])
  n<-ncol(dataset)
  dataset[,n+1]<-dataset[,cancer_cols[i]]/s
  names(dataset)[n+1]<-gsub("bkpt_in_window","density",names(dataset)[cancer_cols[i]])
}
apply(dataset[,15:24],2,sum)




# overall cancer profile
# 1 вероятности появления каждого типа рака
cancer_type_prob_type<-c("brain","blood","bone","uterus","liver","prostate","skin",
                         "ovary","pancreatic","breast")
cancer_type_prob_prob<-c(0.0397145516599441,
                         0.1424139000930810,
                         0.0068259385665529,
                         0.1315544523735650,
                         0.1213155445237360,
                         0.1725100837728820,
                         0.0359913124418244,
                         0.0370772572137760,
                         0.0524356189885200,
                         0.2601613403661190)


cancer_type_prob<-as.data.frame(cbind(cancer_type_prob_type,cancer_type_prob_prob))
names(cancer_type_prob)<-c("type","prob_cancer")
cancer_type_prob$prob_cancer<-as.numeric(as.character(cancer_type_prob$prob_cancer))
cancer_type_prob$type<-as.character(cancer_type_prob$type)


# посчитаем вероятности появления разрывов в определенном окне по формуле полной вероятности
dataset$prod_density_blood<-dataset$density_blood*cancer_type_prob$prob_cancer[cancer_type_prob$type=="blood"]
dataset$prod_density_brain<-dataset$density_brain*cancer_type_prob$prob_cancer[cancer_type_prob$type=="brain"]
dataset$prod_density_bone<-dataset$density_bone*cancer_type_prob$prob_cancer[cancer_type_prob$type=="bone"]
dataset$prod_density_uterus<-dataset$density_uterus*cancer_type_prob$prob_cancer[cancer_type_prob$type=="uterus"]
dataset$prod_density_liver<-dataset$density_liver*cancer_type_prob$prob_cancer[cancer_type_prob$type=="liver"]
dataset$prod_density_prostate<-dataset$density_prostate*cancer_type_prob$prob_cancer[cancer_type_prob$type=="prostate"]
dataset$prod_density_skin<-dataset$density_skin*cancer_type_prob$prob_cancer[cancer_type_prob$type=="skin"]
dataset$prod_density_ovary<-dataset$density_ovary*cancer_type_prob$prob_cancer[cancer_type_prob$type=="ovary"]
dataset$prod_density_pancreatic<-dataset$density_pancreatic*cancer_type_prob$prob_cancer[cancer_type_prob$type=="pancreatic"]
dataset$prod_density_breast<-dataset$density_breast*cancer_type_prob$prob_cancer[cancer_type_prob$type=="breast"]

dataset$overall_density<-dataset$prod_density_blood+dataset$prod_density_brain+dataset$prod_density_bone+
  dataset$prod_density_uterus+dataset$prod_density_liver+dataset$prod_density_prostate+dataset$prod_density_skin+
  dataset$prod_density_ovary+dataset$prod_density_pancreatic+dataset$prod_density_breast


sum(dataset$overall_density)


#write.csv(all_data,"all_cancers.csv")


## BREAKPOINTS HOTSPOTS
#dataset<-new_data
dens<-grep("density",names(dataset))
prods<-grep("prod",names(dataset))
dens<-setdiff(dens,prods)
dens<-setdiff(dens,grep("density_bkpt",names(dataset)))

dataset<-as.data.frame(dataset)

# 1 option (0.1%)

for (i in 1:length(dens)){
  n<-ncol(dataset)
  q<-quantile(dataset[,dens[i]],0.999)
  dataset[,n+1]<-ifelse(dataset[,dens[i]]>q,1,0)
  names(dataset)[n+1]<-paste(gsub("density_","",names(dataset)[dens[i]]),"_bkpt_hotspot_0.1",sep="")
} 



# 2 option (0.05%)

for (i in 1:length(dens)){
  n<-ncol(dataset)
  dataset[,n+1]<-ifelse(dataset[,dens[i]]>quantile(dataset[,dens[i]],0.9995),1,0)
  names(dataset)[n+1]<-paste(gsub("density_","",names(dataset)[dens[i]]),"_bkpt_hotspot_0.05",sep="")
} 



# 3 option (0.01%)

for (i in 1:length(dens)){
  n<-ncol(dataset)
  dataset[,n+1]<-ifelse(dataset[,dens[i]]>quantile(dataset[,dens[i]],0.9999),1,0)
  names(dataset)[n+1]<-paste(gsub("density_","",names(dataset)[dens[i]]),"_bkpt_hotspot_0.01",sep="")
} 

# 4 option (0.5%)
for (i in 1:length(dens)){
  n<-ncol(dataset)
  q<-quantile(dataset[,dens[i]],0.995)
  dataset[,n+1]<-ifelse(dataset[,dens[i]]>q,1,0)
  names(dataset)[n+1]<-paste(gsub("density_","",names(dataset)[dens[i]]),"_bkpt_hotspot_0.5",sep="")
} 



# 5 option (1%)
for (i in 1:length(dens)){
  n<-ncol(dataset)
  q<-quantile(dataset[,dens[i]],0.99)
  dataset[,n+1]<-ifelse(dataset[,dens[i]]>q,1,0)
  names(dataset)[n+1]<-paste(gsub("density_","",names(dataset)[dens[i]]),"_bkpt_hotspot_1",sep="")
} 




htspt<-grep("hotspot",names(dataset))

# общее количество разрывов и по хромосомам в каждом профиле 
# при каждом типе разбивки

for (i in htspt){
  m<-as.data.frame(table(dataset$chr,dataset[,i]))
  m<-m[m$Var2==1,]
  names(m)[3]<-names(dataset)[i]
  
  if (i==htspt[1]){
    m1<-m
  } else {
    m1<-merge(m,m1,by=c("Var1","Var2"))
  }
  
  
}

# всего

total<-as.data.frame(apply(m1[,3:57],2,sum))
total$cancer_type<-rownames(total)
rownames(total)<-NULL
names(total)[1]<-"bkpt_count"
total$labeling[grep("0.1",total$cancer_type)]<-"0.1"
total$labeling[grep("0.01",total$cancer_type)]<-"0.01"
total$labeling[grep("0.05",total$cancer_type)]<-"0.05"
total$labeling[grep("_1",total$cancer_type)]<-"1"
total$labeling[grep("0.5",total$cancer_type)]<-"0.5"

a<-gregexpr("_",total$cancer_type)
b<-lapply(a,function(x) x[[1]])
b<-unlist(b)
total$cancer_type<-substr(total$cancer_type,1,b-1)

library(reshape2)
total_m<-dcast(total,cancer_type ~ labeling, value.var= "bkpt_count" )


# по хромосомам

m1$Var2<-NULL
chr_m<-melt(m1,id.vars = "Var1")

chr_m$labeling[grep("0.1",chr_m$variable)]<-"0.1"
chr_m$labeling[grep("0.01",chr_m$variable)]<-"0.01"
chr_m$labeling[grep("0.05",chr_m$variable)]<-"0.05"
chr_m$labeling[grep("_1",chr_m$variable)]<-"1"
chr_m$labeling[grep("0.5",chr_m$variable)]<-"0.5"


a<-gregexpr("_",chr_m$variable)
b<-lapply(a,function(x) x[[1]])
b<-unlist(b)
chr_m$cancer<-substr(chr_m$variable,1,b-1)

names(chr_m)[1]<-"chr"
chr_m$variable<-NULL

chr_m1<-dcast(chr_m,chr+labeling ~ cancer, value.var = "value")
chr_m1<-chr_m1[order(chr_m1$labeling,chr_m1$chr),]


# не берем 0,01
dataset[,grep("_0.01",names(dataset))]<-NULL
dataset[,grep("_0.05",names(dataset))]<-NULL
dataset[,grep("_0.1",names(dataset))]<-NULL


write.csv(dataset,
          file="..\\data\\preprocessed\\breakpoints\\structural_mutation_final_1_mb.csv")




















## Stemloops




#setwd("~/STEMLOOPS")

setwd("../data/preprocessed/stem-loops")

# read 10kb file
sec_str_all<-read.csv("stemloops_density_10k.csv")


new_wind<-1000000

sec_str<-sec_str_all[order(sec_str_all$chr,sec_str_all$from),]
sec_str$new_wind<-ceiling(sec_str$to/new_wind)

new_sec_str<-sec_str %>% group_by(chr,new_wind) %>% 
  summarize(count_sl_15_30 = sum(count_sl_15_30),
            count_sl_16_50 = sum(count_sl_16_50),
            count_sl_6_15 = sum(count_sl_6_15),
            length_sl_15_30 = sum(length_sl_15_30),
            length_sl_16_50 = sum(length_sl_16_50),
            length_sl_6_15 = sum(length_sl_6_15),
            from = min(from),
            to = max(to)
  )


new_sec_str$coverage_sl_15_30<-new_sec_str$length_sl_15_30/new_wind
new_sec_str$coverage_sl_16_50<-new_sec_str$length_sl_16_50/new_wind
new_sec_str$coverage_sl_6_15<-new_sec_str$length_sl_6_15/new_wind

rm(sec_str)
rm(sec_str_all)
gc()

write.csv(new_sec_str,file="stemloops_density_1mb.csv",row.names = FALSE)
















# stemloops features


library(ggplot2)
library(reshape2)


mutations<-read.csv("..\\data\\preprocessed\\breakpoints\\structural_mutation_final_100kb.csv")
names_mut<-c(c("chr","new_wind","from","to"),names(mutations)[grep("hotspot",names(mutations))])
mutations<-mutations[,names_mut]

mutations$chr<-as.character(mutations$chr)




steml<-read.csv("stemloops_density_100k.csv")
steml<-steml[,c( "chr" ,"from","density_S15.30","density_S16.50","density_S6.15")]
steml$chr<-as.character(steml$chr)
steml$chr<-gsub("chr","",steml$chr)

data<-merge(mutations,steml,by=c("chr","from"))
chr_names<-unique(data$chr)

rm(mutations)
rm(steml)
gc()

str<-grep("density_S",names(data))



for (i in 1:length(chr_names)){
  print(i)
  data_chr<-data[data$chr==chr_names[i],]
  data_chr<-data_chr[order(data_chr$from),]
  n<-nrow(data_chr)
  
  
  for (j in str){
    print(j)
    
    # 3 previous neighbours
    data_chr<-cbind(data_chr,c(0,data_chr[2:n,j]-data_chr[1:(n-1),j]))
    data_chr<-cbind(data_chr,c(0,0,data_chr[3:n,j]-data_chr[1:(n-2),j]))
    data_chr<-cbind(data_chr,c(0,0,0,data_chr[4:n,j]-data_chr[1:(n-3),j]))
    # ranks?
    
    
    # 3 following neighbours
    data_chr<-cbind(data_chr,c(data_chr[1:(n-1),j]-data_chr[2:n,j],0))
    data_chr<-cbind(data_chr,c(data_chr[1:(n-2),j]-data_chr[3:n,j],0,0))
    data_chr<-cbind(data_chr,c(data_chr[1:(n-3),j]-data_chr[4:n,j],0,0,0))
    # ranks?
    
  }
  
  
  names(data_chr)[41:58]<-c(
    paste(names(data)[str[1]],"prev_1",sep="_"),
    paste(names(data)[str[1]],"prev_2",sep="_"),
    paste(names(data)[str[1]],"prev_3",sep="_"),
    
    
    paste(names(data)[str[1]],"fol_1",sep="_"),
    paste(names(data)[str[1]],"fol_2",sep="_"),
    paste(names(data)[str[1]],"fol_3",sep="_"),
    
    paste(names(data)[str[2]],"prev_1",sep="_"),
    paste(names(data)[str[2]],"prev_2",sep="_"),
    paste(names(data)[str[2]],"prev_3",sep="_"),
    
    paste(names(data)[str[2]],"fol_1",sep="_"),
    paste(names(data)[str[2]],"fol_2",sep="_"),
    paste(names(data)[str[2]],"fol_3",sep="_"),
    
    
    paste(names(data)[str[3]],"prev_1",sep="_"),
    paste(names(data)[str[3]],"prev_2",sep="_"),
    paste(names(data)[str[3]],"prev_3",sep="_"),
    
    paste(names(data)[str[3]],"fol_1",sep="_"),
    paste(names(data)[str[3]],"fol_2",sep="_"),
    paste(names(data)[str[3]],"fol_3",sep="_")
    
  )
  
  
  if (i==1){
    new_data_chr<-data_chr
  } else {
    new_data_chr<-rbind(new_data_chr,data_chr)
  }
  
  
}



# 1 2.5 5 quantile

for (j in str){
  
  q<-as.numeric(quantile(new_data_chr[,j],c(0.95,0.975,0.99)))
  new_data_chr<-cbind(new_data_chr,ifelse(new_data_chr[,j]>q[1],1,0))
  new_data_chr<-cbind(new_data_chr,ifelse(new_data_chr[,j]>q[2],1,0))
  new_data_chr<-cbind(new_data_chr,ifelse(new_data_chr[,j]>q[3],1,0))
  
}  


names(new_data_chr)[59:67]<-c(
  paste(names(data)[str[1]],"q5",sep="_"),
  paste(names(data)[str[1]],"q2.5",sep="_"),
  paste(names(data)[str[1]],"q1",sep="_"),
  
  paste(names(data)[str[2]],"q5",sep="_"),
  paste(names(data)[str[2]],"q2.5",sep="_"),
  paste(names(data)[str[2]],"q1",sep="_"),
  
  paste(names(data)[str[3]],"q5",sep="_"),
  paste(names(data)[str[3]],"q2.5",sep="_"),
  paste(names(data)[str[3]],"q1",sep="_")
)



write.csv(new_data_chr,"stemloops_mut_data_100k.csv")










################### QUADRUPLEXES







#### GOOD FORMAT FOR BASE 10 kb 

# template (all windows)
new_data<-read.csv("..\\data\\preprocessed\\breakpoints\\structural_mutation_final_10_Kb.csv")
new_data$X<-NULL

setwd("../data/preprocessed/quadruplexes")

# read 10kb file
sec_str_all<-read.csv("quadr_cov.csv")
sec_str_all$X<-NULL

sec_str_all<-merge(sec_str_all,new_data[,c("chr","from")],
                   by=c("chr","from"),all.y = TRUE)

sec_str_all[,3:5]<-apply(sec_str_all[,3:5],2,function(x) ifelse(is.na(x),0,x))
write.csv(sec_str_all,file="quadr_cov.csv")




### AGGREGATION


library(dplyr)
library(ggplot2)
library(reshape2)


# aggregate
new_wind<-1000000

setwd("../data/preprocessed/quadruplexes")
# read 10kb file
sec_str_all<-read.csv("quadr_cov.csv")

sec_str<-sec_str_all[order(sec_str_all$chr,sec_str_all$from),]
sec_str$X<-NULL
sec_str$to<-sec_str$from+10000
sec_str$new_wind<-ceiling(sec_str$to/new_wind)

new_sec_str<-sec_str %>% group_by(chr,new_wind) %>% 
  summarize(count_q = sum(count_q),
            length_q = sum(length_q),
            from = min(from),
            to = max(to)
  )
new_sec_str$coverage_q<-new_sec_str$length_q/new_wind

rm(sec_str)
rm(sec_str_all)
gc()

write.csv(new_sec_str,file="quadr_cov_1mb.csv",row.names = FALSE)
