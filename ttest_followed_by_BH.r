# Run t-test followed BH
library(argparse)

parser <- ArgumentParser(description='Run t-test followed by BH',prog = "Rscript ttest_followed_by_BH.R")

parser$add_argument('experiment', metavar='experiment.csv', type="character", nargs=1,
                    help='Experiment')

parser$add_argument('counts', metavar='count_file.txt', type="character", nargs=1,
                    help='Count File (PORT output)')

parser$add_argument('result', metavar='result.csv', type="character", nargs=1,
                    help='Out File')


args <- commandArgs(TRUE)

## TESTING
#setwd("~/Documents/DifferentialSplicing/emanuela_refseq_june/")
#args <- c("experiment_IL1B.csv",
#          "annotated_master_list_of_exons_counts_MIN.sense.seq_run_june14.txt",
#          "results_RofevsCtrl_summary2.csv")
###########
if(length(args) < 1) {
  args <- c("-h")
}



args <- parser$parse_args(args)



print(args)

#parser$print_help()

experiment = read.csv(args$experiment,header = TRUE,stringsAsFactors=FALSE)
head(experiment,50)

d = read.csv(args$counts,sep = "\t", header = TRUE, na.strings = "",stringsAsFactors=FALSE)
#CTRL vs. CTRL IL1B
#head(d)
#d = read.csv(pipe("cut -f 1,18-25,42-49 FINAL_master_list_of_exons_counts_MIN.emanuela_NEW.txt"),sep = "\t", header = TRUE, na.strings = "",stringsAsFactors=FALSE)
print("Read data - DONE")
d = d[,c("id",experiment$SampleName)]
#test = d[1,experiment$SampleName[experiment$Group == cond[1]]]
head(d)
cond = unique(experiment$Group)
head(cond)

data=d[rowSums(d[ ,-1])>10,] 
P_Val = c(rep(NA,dim(data)[1]))
#colnames(p_val) = "P_VAL" 
for (i in 1:dim(data)[1]) {
  p_val_cur <- tryCatch( wilcox.test(as.matrix(data[i,experiment$SampleName[experiment$Group == cond[1]]]),as.matrix(data[i,experiment$SampleName[experiment$Group == cond[2]]]),alternative="two.sided",exact=TRUE)$p.value, error=function(x) 1 ) 
  P_Val[i]=p_val_cur 
}

FDR <- p.adjust(P_Val, method="BH")
out_table = cbind(data,P_Val,FDR)

head(out_table[order(out_table$FDR),])

hist(FDR)
write.csv(out_table, args$result)
