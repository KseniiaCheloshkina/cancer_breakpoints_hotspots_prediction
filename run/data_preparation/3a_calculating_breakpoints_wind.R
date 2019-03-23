library(dplyr)

#### FUNCTION for creating windows of specified length

get_chr_windows<-function(wind=10000){
  # number of bases in chromosomes
  chr_base<-list(c(249250621,243199373,198022430,191154276,
                   180915260,171115067,159138663,146364022,
                   141213431,135534747,135006516,133851895,
                   115169878,107349540,102531392,90354753,
                   81195210,78077248,59128983,63025520,
                   48129895,51304566,155270560,59373566),
                 c(as.character(seq(1,22,1)),"X","Y"))
  
  chr_windows<-data.frame(chr=0,from=0,to=0)
  chr_windows<-chr_windows[-1,]
  # for each chromosome
  for (i in 1:length(chr_base[[2]])){
    print(paste("Chromosome ",chr_base[[2]][i]))
    # generate windows boundaries including right border
    a<-seq(0,chr_base[[1]][i],wind)
    a<-append(a,chr_base[[1]][i])
    a<-unique(a)
    
    from<-a[1:length(a)-1]
    to<-a[2:length(a)]
    
    new_row<-as.data.frame(cbind(rep(chr_base[[2]][i],length(a)-1),from,to))
    names(new_row)<-names(chr_windows)
    chr_windows<-rbind(chr_windows,new_row)
    
  }
  chr_windows$from<-as.numeric(as.character(chr_windows$from))
  chr_windows$to<-as.numeric(as.character(chr_windows$to))
  return(chr_windows)
}



####   FUNCTION to calculate number of breakpoints in window

bkpt <- function(all_data, chr_windows_1){
  chr_windows<-chr_windows_1
  all_data$window_0<-ceiling(all_data$chr_bkpt_beg/window_size)
  all_data$window_1<-ceiling(all_data$chr_bkpt_end/window_size)
  all_data$dif_wind<-all_data$window_1-all_data$window_0
  
  
  a<-all_data %>% group_by(chr,window_0) %>%
    summarize(count_bkpt = n())
  

  b<-all_data[all_data$dif_wind==1,]
  b<-b %>% group_by(chr,window_1) %>%
    summarize(count_bkpt = n())
  

  names(a)[2]<-"window"
  names(b)[2]<-"window"
  quantity<-rbind(a,b)
  quantity<-quantity %>% group_by(chr,window) %>%
    summarize(count_bkpt = sum(count_bkpt))
  

  agg_quantity<-quantity %>% group_by(chr) %>% summarize(count_bkpt = sum(count_bkpt))
  
  

  chr<-unique(chr_windows$chr)
  
  for (i in 1:24){
    chr_windows$window[chr_windows$chr==chr[i]]<-seq(1,nrow(chr_windows[chr_windows$chr==chr[i],]))
  }
  
  chr_windows<-merge(chr_windows,agg_quantity,by="chr",all.x=TRUE)
  names(chr_windows)[5]<-"bkpt_in_chr"
  chr_windows<-merge(chr_windows,quantity,by=c("chr","window"),all.x=TRUE)
  names(chr_windows)[6]<-"bkpt_in_window"
  chr_windows$bkpt_in_window[is.na(chr_windows$bkpt_in_window)]<-0
  chr_windows$density<-chr_windows$bkpt_in_window/chr_windows$bkpt_in_chr
  return(chr_windows)
}





#### CALCULATE BREAKPOINTS DENSITY FOR EACH CANCER TYPE


input_folder<-"../Diploma/data/structural mutation/2_eda_preprocessing"
setwd(input_folder)
output_folder<-"../3_breakpoints/"

# set window size
window_size=10000


chr_windows<-get_chr_windows(window_size)

cancer_types<-list.files()
cancer_types<-cancer_types[-grep("all_cancer_data",cancer_types)]
cancer_types

chr_windows_1<-get_chr_windows(window_size)


for (j in cancer_types){
  data<-read.csv(j)
  data<-bkpt(data,chr_windows_1)
  filename<-paste(output_folder,j,sep="")
  write.csv(data,row.names = FALSE,file=filename)
}



#### CALCULATE STATS

input_folder<-"../Diploma/data/structural mutation/2_eda_preprocessing"
setwd(input_folder)

cancer_types<-list.files()
cancer_types<-cancer_types[-grep("all_cancer_data",cancer_types)]
cancer_types

n_samples<-vector()
n_donors<-vector()
n_bkpt<-vector()

for (j in cancer_types){
  data<-read.csv(j)
  n_samples<-append(n_samples,length(unique(data$icgc_sample_id)))
  n_donors<-append(n_donors,length(unique(data$icgc_donor_id)))
  n_bkpt<-append(n_bkpt,nrow(data))
}

stats<-as.data.frame(cbind(gsub("_all_data.csv","",cancer_types),
             n_samples,
             n_donors,
             n_bkpt))
names(stats)<-c("cancer_type","n_samples","n_donors","n_bkpt")

stats[,2:4]<-apply(stats[,2:4],2,function(x) as.numeric(as.character(x)))





#### REMOVE Y CHROMOSOME FROM DATA

input_folder<-"../Diploma/data/structural mutation/3_breakpoints"
setwd(input_folder)

cancer_types<-c("blood_all_data.csv","bone_all_data.csv","brain_all_data.csv", "breast_all_data.csv" ,"liver_all_data.csv","ovary_all_data.csv","pancreatic_all_data.csv",
                "prostate_all_data.csv", "skin_all_data.csv","uterus_all_data.csv")

for (j in 1:length(cancer_types)){
  print(j)
  print(cancer_types[j])
  data<-read.csv(cancer_types[j])
  data$chr<-as.character(data$chr)
  data<-data[data$chr!="Y",]
  write.csv(data,cancer_types[j])
  
}