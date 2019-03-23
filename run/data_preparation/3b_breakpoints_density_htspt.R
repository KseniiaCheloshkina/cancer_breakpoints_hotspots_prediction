library(dplyr)
library(ggplot2)
library(reshape2)

## ПЛОТНОСТЬ РАЗРЫВОВ

input_folder<-"E:/Учеба/Диплом/Diploma/data/structural mutation/3_breakpoints"
setwd(input_folder)

cancer_types<-c("blood_all_data.csv","bone_all_data.csv","brain_all_data.csv", "breast_all_data.csv" ,"liver_all_data.csv","ovary_all_data.csv","pancreatic_all_data.csv",
                "prostate_all_data.csv", "skin_all_data.csv","uterus_all_data.csv")

for (j in 1:length(cancer_types)){
  print(j)
  print(cancer_types[j])
  data<-read.csv(cancer_types[j])
  if (j==1){
    data<-subset(data,select = c("chr","window","from","to","bkpt_in_window"))
    names(data)[5]<-paste(names(data)[5],gsub("_all_data.csv","",cancer_types[j]),sep="_")
    all_data<-data
  } else {
    data<-subset(data,select = c("chr","window","bkpt_in_window"))
    
    names(data)[3]<-paste(names(data)[3],gsub("_all_data.csv","",cancer_types[j]),sep="_")
    all_data<-merge(all_data,data,by=c("chr","window"))
  }
  
}


# считаем плотность для каждого типа рака

# не рассматриваем у-хромосому
all_data$chr<-as.character(all_data$chr)
all_data<-all_data[all_data$chr!="Y",]


cancer_cols<-grep("bkpt",names(all_data))

for (i in 1:length(cancer_cols)){
  s<-sum(all_data[,cancer_cols[i]])
  n<-ncol(all_data)
  all_data[,n+1]<-all_data[,cancer_cols[i]]/s
  names(all_data)[n+1]<-gsub("bkpt_in_window","density",names(all_data)[cancer_cols[i]])
}


apply(all_data[,15:24],2,sum)


## overall cancer profile

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
all_data$prod_density_blood<-all_data$density_blood*cancer_type_prob$prob_cancer[cancer_type_prob$type=="blood"]
all_data$prod_density_brain<-all_data$density_brain*cancer_type_prob$prob_cancer[cancer_type_prob$type=="brain"]
all_data$prod_density_bone<-all_data$density_bone*cancer_type_prob$prob_cancer[cancer_type_prob$type=="bone"]
all_data$prod_density_uterus<-all_data$density_uterus*cancer_type_prob$prob_cancer[cancer_type_prob$type=="uterus"]
all_data$prod_density_liver<-all_data$density_liver*cancer_type_prob$prob_cancer[cancer_type_prob$type=="liver"]
all_data$prod_density_prostate<-all_data$density_prostate*cancer_type_prob$prob_cancer[cancer_type_prob$type=="prostate"]
all_data$prod_density_skin<-all_data$density_skin*cancer_type_prob$prob_cancer[cancer_type_prob$type=="skin"]
all_data$prod_density_ovary<-all_data$density_ovary*cancer_type_prob$prob_cancer[cancer_type_prob$type=="ovary"]
all_data$prod_density_pancreatic<-all_data$density_pancreatic*cancer_type_prob$prob_cancer[cancer_type_prob$type=="pancreatic"]
all_data$prod_density_breast<-all_data$density_breast*cancer_type_prob$prob_cancer[cancer_type_prob$type=="breast"]

all_data$overall_density<-all_data$prod_density_blood+all_data$prod_density_brain+all_data$prod_density_bone+
  all_data$prod_density_uterus+all_data$prod_density_liver+all_data$prod_density_prostate+all_data$prod_density_skin+
  all_data$prod_density_ovary+all_data$prod_density_pancreatic+all_data$prod_density_breast


sum(all_data$overall_density)


write.csv(all_data,"all_cancers.csv")





# EDA density

library(dplyr)
library(ggplot2)
library(corrplot)
library(reshape2)

all_data<-read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\3_breakpoints\\all_cancers.csv")
all_data$X<-NULL


#correlations
densities<-subset(all_data,select = c("density_blood","density_bone",
                                      "density_brain","density_breast",
                                      "density_liver","density_ovary",
                                      "density_pancreatic",
                                      "density_prostate", "density_skin" ,
                                      "density_uterus","overall_density"))
a<-cor(densities,method="spearman")
for (i in 1:11){
  a[i,i]<-NA
}
mean(a[1:10,1:10],na.rm=TRUE)
median(a[1:10,1:10],na.rm=TRUE)
max(a[1:10,1:10],na.rm=TRUE)

mean(a[11,],na.rm=TRUE)
median(a[11,],na.rm=TRUE)
max(a[11,],na.rm=TRUE)

new_cols <- vector()
for (col in names(densities)){
  ncol <- gsub("density", "", col)
  ncol <- gsub("_", "", ncol)
  new_cols <- c(new_cols, ncol)
}
names(densities) <- new_cols
corrm <- cor(densities,method="spearman")
g1 <- corrplot(corrm,tl.col ="black")

df <- data.frame(corrm)
df$cancer <- row.names(df)
df_plot <- melt(df, id.vars = "cancer")

g1 <- ggplot(df_plot, aes(x=cancer, y=variable, fill=value))+
  geom_tile() + 
  xlab("Cancer type") + 
  ylab("Cancer type") + 
  scale_fill_continuous(low = "white", high = "steelblue", name="Correlation \n")
g1

save(g1, file = "E:\\Учеба\\Диплом\\Diploma\\scripts\\data for plot\\1_B.RData")



########### Найдем порог для определения breakpoint hotspots
# общий раковый профиль

ggplot(all_data, aes(overall_density))+
  geom_histogram(bins = 500)+
  ggtitle("Histogram of overall density")





















## BREAKPOINTS HOTSPOTS

dens<-grep("density",names(all_data))
prods<-grep("prod",names(all_data))
dens<-setdiff(dens,prods)


# 1 option (0.1%)


for (i in 1:length(dens)){
  n<-ncol(all_data)
  all_data[,n+1]<-ifelse(all_data[,dens[i]]>quantile(all_data[,dens[i]],0.999),1,0)
  names(all_data)[n+1]<-paste(gsub("density_","",names(all_data)[dens[i]]),"_bkpt_hotspot_0.1",sep="")
} 



# 2 option (0.05%)

for (i in 1:length(dens)){
  n<-ncol(all_data)
  all_data[,n+1]<-ifelse(all_data[,dens[i]]>quantile(all_data[,dens[i]],0.9995),1,0)
  names(all_data)[n+1]<-paste(gsub("density_","",names(all_data)[dens[i]]),"_bkpt_hotspot_0.05",sep="")
} 



# 3 option (0.01%)

for (i in 1:length(dens)){
  n<-ncol(all_data)
  all_data[,n+1]<-ifelse(all_data[,dens[i]]>quantile(all_data[,dens[i]],0.9999),1,0)
  names(all_data)[n+1]<-paste(gsub("density_","",names(all_data)[dens[i]]),"_bkpt_hotspot_0.01",sep="")
} 



write.csv(all_data,
          file="E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\structural_mutation_final_2.csv")




## HOTSPOTS EDA

# отрисовать линейный график плотности
data<-read.csv("E:/Учеба/Диплом/Diploma/data/structural mutation/structural_mutation_final_10_kb.csv")
library(ggplot2)

g1 <- ggplot(data[data$chr=="17",],aes(x=from/1000000,y=overall_density))+
  geom_line()+
  xlab("Start position of window")+
  ylab("Density")
save(g1, file="E:\\Учеба\\Диплом\\Diploma\\scripts\\data for plot\\1_C.RData")



new_data<-read.csv("E:\\Учеба\\Диплом\\Diploma\\data\\structural mutation\\structural_mutation_final_10_kb.csv")

htspt<-grep("hotspot",names(new_data))

# общее количество разрывов и по хромосомам в каждом профиле 
# при каждом типе разбивки

for (i in htspt){
  m<-as.data.frame(table(new_data$chr,new_data[,i]))
  m<-m[m$Var2==1,]
  names(m)[3]<-names(new_data)[i]
  
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

##### ЧИСЛО ОКОН!!!!!!!!!!!!!!
# нормализуем на количество окон,а не оснований
chr_base<-as.data.frame(cbind(c(as.character(seq(1,22,1)),"X","Y"),
                              c(24926,24320,19803,19116, 18092, 17112, 15914, 
                                14637, 14122, 13554, 13501, 13386,11517, 10735,
                                10254,  9036,  8120,  7808,  5913,  6303,  4813,  5131, 15528,  5938)
))
names(chr_base)<-c("chr","length")
chr_count_h<-merge(chr_m,chr_base,by="chr")
chr_count_h$length<-as.numeric(as.character(chr_count_h$length))
chr_count_h$norm_value<-chr_count_h$value/chr_count_h$length



chr_count_h$chr <- as.character(chr_count_h$chr)
chr_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
               "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
               "21", "22", "X")
g2 <- ggplot(chr_count_h[chr_count_h$labeling=="0.1",],aes(x=cancer,y=chr,group=1))+
  geom_point(data=chr_count_h[chr_count_h$labeling=="0.1",],
             aes(x=cancer, y=chr,size = norm_value, colour = norm_value))+
  scale_y_discrete(limits=chr_order)+
  xlab("Cancer type")+
  ylab("Chromosome") +  guides(colour=guide_legend(title="Percentage"), 
                               size=guide_legend(title="Percentage"))
  #scale_colour_gradient(low = "white", high = "black")+
  # ggtitle("Number of breakpoint hotspots by chromosome and type of cancer\n1 option (0.1)")
save(g2, file = "E:\\Учеба\\Диплом\\Diploma\\scripts\\data for plot\\1_A.RData")

ggplot(chr_count_h[chr_count_h$labeling=="0.05",],aes(x=cancer,y=chr,group=1))+
  geom_point(data=chr_count_h[chr_count_h$labeling=="0.05",],aes(x=cancer, y=chr,size = norm_value, colour = norm_value))+
  #scale_colour_gradient(low = "white", high = "black")+
  ggtitle("Number of breakpoint hotspots by chromosome and type of cancer\n2 option (0.05)")


ggplot(chr_count_h[chr_count_h$labeling=="0.01",],aes(x=cancer,y=chr,group=1))+
  geom_point(data=chr_count_h[chr_count_h$labeling=="0.01",],aes(x=cancer, y=chr,size = norm_value, colour = norm_value))+
  #scale_colour_gradient(low = "white", high = "black")+
  ggtitle("Number of breakpoint hotspots by chromosome and type of cancer\n3 option (0.01)")









# jaccard similarity of breakpoint hotspots profile
jaccard_sim<-function(x,y){
  cont_table<-as.data.frame(table(x,y))
  #print(cont_table)
  m11<-cont_table$Freq[(cont_table[,1]==1)&(cont_table[,2]==1)]
  m10<-cont_table$Freq[(cont_table[,1]==1)&(cont_table[,2]==0)]
  m01<-cont_table$Freq[(cont_table[,1]==0)&(cont_table[,2]==1)]
  
  jaccard_ind<-m11/(m10+m01+m11)
  return(jaccard_ind)
}


# для каждой строки проверяем долю
col_n<-c("chr","window",names(new_data)[grep("_hotspot_0.01",names(new_data))])
htspt<-new_data[,col_n]

jacc_sim<-data.frame(type1=0,type2=0,jaccard_similarity=0)

for (i in 3:13){
  print(i)
  for (j in 3:13){
    jacc_sim<-rbind(jacc_sim,
                    c(names(htspt)[i],
                      names(htspt)[j],
                      jaccard_sim(htspt[,i],htspt[,j])
                    )
    )
  }
  
}

jacc_sim<-jacc_sim[-1,]
jacc_sim$jaccard_similarity<-as.numeric(jacc_sim$jaccard_similarity)
jacc_sim$type1<-gsub("_bkpt_hotspot_0.01","",jacc_sim$type1)
jacc_sim$type2<-gsub("_bkpt_hotspot_0.01","",jacc_sim$type2)
jacc_sim$type1<-gsub("_density","",jacc_sim$type1)
jacc_sim$type2<-gsub("_density","",jacc_sim$type2)


ggplot(data = jacc_sim, aes(x=type1, y=type2, fill=jaccard_similarity))+
  geom_tile(colour = "grey50")+
  scale_fill_gradient2(low = "white", high = "purple",  
                       space = "Lab", 
                       name="Jaccard\nSimilarity\n")+
  ggtitle("Jaccard similarity of labeled genome (3 option, 0.01)")+
  geom_text(aes(type1, type2, label = format(jaccard_similarity,digits=2)), color = "black", size = 3)+
  xlab("type of cancer / overall")+
  ylab("type of cancer / overall")

