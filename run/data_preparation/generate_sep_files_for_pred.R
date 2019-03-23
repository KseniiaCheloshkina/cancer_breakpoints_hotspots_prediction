#### This script is used to generate separate .R files for each configuration and type of predictors 
#### First type - datasets with stemloops coverage and target (breakpoints hotspot)
#### Second type - datasets with quadruplexes coverage and target  (breakpoints hotspot)


############ STEMLOOPS
f_to_write<-"E:\\Учеба\\Диплом\\Diploma\\data\\stemloops_mut\\data for pred\\"
mut_fold<-"E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\"

mut_f<-c("structural_mutation_final_10_kb.csv","structural_mutation_final_100_kb.csv",
         "structural_mutation_final_20_kb.csv","structural_mutation_final_50_kb.csv",
         "structural_mutation_final_500_kb.csv","structural_mutation_final_1_mb.csv")

sec_str_fold<-"E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\STEMLOOPS\\"

sec_str<-c("stemloops_density_100k.csv","stemloops_density_10k.csv" ,
           "stemloops_density_1mb.csv","stemloops_density_20k.csv",
           "stemloops_density_500k.csv","stemloops_density_50k.csv" )


##### 10kb
level_n<-"10kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_10_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"stemloops_density_10k.csv",sep=""))

predictor<-data.frame(
  v1=c("y","coverage_15_30","coverage_16_50","coverage_sl_6_15","chr","from"),
  v2=paste(rep("v",6),seq(1,6),sep="_"),stringsAsFactors = F
)



##### 20kb
level_n<-"20kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_20_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"stemloops_density_20k.csv",sep=""))


predictor<-data.frame(
  v1=c("y","coverage_sl_15_30","coverage_sl_16_50","coverage_sl_6_15","chr","from"),
  v2=paste(rep("v",6),seq(1,6),sep="_"),stringsAsFactors = F
)

########## 50kb
level_n<-"50kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_50_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"stemloops_density_50k.csv",sep=""))


predictor<-data.frame(
  v1=c("y","coverage_sl_15_30","coverage_sl_16_50","coverage_sl_6_15","chr","from"),
  v2=paste(rep("v",6),seq(1,6),sep="_"),stringsAsFactors = F
)
###### 100kb
level_n<-"100kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_100_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"stemloops_density_100k.csv",sep=""))


predictor<-data.frame(
  v1=c("y","coverage_sl_15_30","coverage_sl_16_50","coverage_sl_6_15","chr","from"),
  v2=paste(rep("v",6),seq(1,6),sep="_"),stringsAsFactors = F
)
###### 500kb
level_n<-"500kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_500_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"stemloops_density_500k.csv",sep=""))


predictor<-data.frame(
  v1=c("y","coverage_sl_15_30","coverage_sl_16_50","coverage_sl_6_15","chr","from"),
  v2=paste(rep("v",6),seq(1,6),sep="_"),stringsAsFactors = F
)
###### 1 mb
level_n<-"1mb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_1_mb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"stemloops_density_1mb.csv",sep=""))


predictor<-data.frame(
  v1=c("y","coverage_sl_15_30","coverage_sl_16_50","coverage_sl_6_15","chr","from"),
  v2=paste(rep("v",6),seq(1,6),sep="_"),stringsAsFactors = F
)
###### 




class_data<-merge(mut_file,sec_str_file,by=c("chr","from"))

list_y<-names(class_data)[grep("hotspot",names(class_data))]

class_data<-class_data[,c(list_y,
                          names(class_data)[grep("coverage",names(class_data))],
                          "chr", "from")]
m<-as.data.frame(apply(class_data[,list_y],2,function(x) as.factor(x)))
class_data[,list_y]<-m


for (i in 1:length(list_y)){
  y_name<-list_y[i]
  y<-which(names(class_data)==y_name)
  print(y_name)
  
  predictor_iter<-predictor
  predictor_iter[1,1]<-y_name
  
  
  iter_data<-class_data[,predictor_iter[,1]]
  names(iter_data)<-predictor_iter[,2]
  
  save(iter_data,file=paste(f_to_write,level_n,y_name,".RData",sep=""))

}








######## QUADRUPLEXES

f_to_write<-"E:\\Учеба\\Диплом\\Diploma\\data\\quadr_mut\\data for pred\\"
mut_fold<-"E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\"

mut_f<-c("structural_mutation_final_10_kb.csv","structural_mutation_final_100_kb.csv",
         "structural_mutation_final_20_kb.csv","structural_mutation_final_50_kb.csv",
         "structural_mutation_final_500_kb.csv","structural_mutation_final_1_mb.csv")

sec_str_fold<-"E:\\Учеба\\Диплом\\Diploma\\data\\secondary\\QUADRUPLEXES\\"

sec_str<-c("quadr_cov.csv","quadr_cov_1mb.csv" ,
           "quadr_cov_20kb.csv","quadr_cov_50kb.csv",
           "quadr_cov_100kb.csv","quadr_cov_500kb.csv" )


##### 1mb
level_n<-"1mb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_1_mb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"quadr_cov_1mb.csv",sep=""))


##### 500kb
level_n<-"500kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_500_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"quadr_cov_500kb.csv",sep=""))


##### 100kb
level_n<-"100kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_100_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"quadr_cov_100kb.csv",sep=""))

##### 50kb
level_n<-"50kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_50_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"quadr_cov_50kb.csv",sep=""))


##### 20kb
level_n<-"20kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_20_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"quadr_cov_20kb.csv",sep=""))

##### 10kb
level_n<-"10kb_"
mut_file<-read.csv(paste(mut_fold,"structural_mutation_final_10_kb.csv",sep=""))
mut_file$X<-NULL

sec_str_file<-read.csv(paste(sec_str_fold,"quadr_cov.csv",sep=""))






predictor<-data.frame(
  v1=c("y","coverage_q","chr","from"),
  v2=paste(rep("v",4),seq(1,4),sep="_"),stringsAsFactors = F
)

class_data<-merge(mut_file,sec_str_file,by=c("chr","from"))

list_y<-names(class_data)[grep("hotspot",names(class_data))]

class_data<-class_data[,c(list_y,
                          names(class_data)[grep("coverage",names(class_data))],
                          "chr", "from")]
m<-as.data.frame(apply(class_data[,list_y],2,function(x) as.factor(x)))
class_data[,list_y]<-m



for (i in 1:length(list_y)){
  y_name<-list_y[i]
  y<-which(names(class_data)==y_name)
  print(y_name)
  
  predictor_iter<-predictor
  predictor_iter[1,1]<-y_name
  
  
  iter_data<-class_data[,predictor_iter[,1]]
  names(iter_data)<-predictor_iter[,2]
  
  save(iter_data,file=paste(f_to_write,level_n,y_name,".RData",sep=""))
  
}
