# devtools::install_github("trevorld/argparse")
library(argparse)
source("/Users/hayer/Documents/DifferentialSplicing/AnovaSplice/anovasplicecalc.R")

parser <- ArgumentParser(description='Add q-values to output of AnovaSplice',prog = "Rscript add_q_values.R")

parser$add_argument('experiment', metavar='experiment.csv', type="character", nargs=1,
                    help='Experiment')

parser$add_argument('anno', metavar='annotation.gtf', type="character", nargs=1,
                    help='Annotation File')

parser$add_argument('counts', metavar='count_file.txt', type="character", nargs=1,
                    help='Count File (PORT output)')

parser$add_argument('result', metavar='result.csv', type="character", nargs=1,
                    help='Out File')

parser$add_argument('--normalize', dest='normalize', action='store_const',
                    const='normalize', default='not',
                    help='normalize by exon length  (default: raw counts)')

# Normalize option

args <- commandArgs(TRUE)

# TESTING
#setwd("~/Documents/DifferentialSplicing/emanuela_refseq_june/")
#args <- c("experiment_IL1B.csv","fixed_refseq_mm9_small.gtf",
#          "annotated_master_list_of_exons_counts_MIN.sense.seq_run_june14.txt","results_RofevsCtrl_summary2.csv")
if(length(args) < 1) {
  args <- c("-h")
}

args <- parser$parse_args(args)



print(args)

#parser$print_help()

experiment = read.csv(args$experiment,header = TRUE,stringsAsFactors=FALSE)
head(experiment,50)


gtf_file = read.csv(args$anno,sep = "\t", header = FALSE, na.strings = "",stringsAsFactors=FALSE)
cmd = paste("cut -f 3 ", args$anno, " | grep -w exon| wc -l", sep = "")
max_num_exons = as.integer(system(cmd,intern = TRUE))
my_list <- vector("list",max_num_exons)
ensgene <- character(max_num_exons)
exons <- character(max_num_exons)
exons_length <- character(max_num_exons)
num_splice_forms <- integer(max_num_exons)
start_codon_counter <- character(max_num_exons)
#key <- "width"
#value <- 32

#mylist <- list()
#mylist[[ key ]] <- value

k = 1
for (i in seq(dim(gtf_file)[1])) {
  #for (i in seq(100)) {
  if (gtf_file$V3[i] == "start_codon") {
    gene_id = strsplit(gtf_file$V9[i], "[; ]")[[1]][2]
    start_codon_counter[k] = gene_id
  }
  if (gtf_file$V3[i] == "exon") {
    gene_id = strsplit(gtf_file$V9[i], "[; ]")[[1]][2]
    name = paste("exon:", gtf_file$V1[i], ":", gtf_file$V4[i], "-", gtf_file$V5[i], sep = "")
    ensgene[k] = gene_id
    exons[k] = name
    exons_length[k] = gtf_file$V5[i] - gtf_file$V4[i]
    k = k+1
  }
}

mapping = unique(data.frame(exons,ensgene,exons_length,num_splice_forms,stringsAsFactors=FALSE))
head(mapping)
ensgene = unique(ensgene)
head(ensgene)
start_codon_counter = start_codon_counter[start_codon_counter != '']
head(start_codon_counter)
for (ensgene_name in ensgene) {
  mapping$num_splice_forms[mapping$ensgene==ensgene_name] = sum(start_codon_counter==ensgene_name)
}

head(mapping)

rm(exons,exons_length,num_splice_forms,start_codon_counter)

#d = read.csv("/Users/hayer/Documents/DifferentialSplicing/emanuela_refseq/annotated_master_list_of_exons_counts_MIN.study1.txt", 
#             sep = "\t", header = TRUE, na.strings = "",stringsAsFactors=FALSE)
#head(d)

d = read.csv(args$counts,sep = "\t", header = TRUE, na.strings = "",stringsAsFactors=FALSE)
#CTRL vs. CTRL IL1B
head(d)
#d = read.csv(pipe("cut -f 1,18-25,42-49 FINAL_master_list_of_exons_counts_MIN.emanuela_NEW.txt"),sep = "\t", header = TRUE, na.strings = "",stringsAsFactors=FALSE)
print("Read data - DONE")
d = d[,c("id",experiment$SampleName)]

#quit()

# ###### TEST2
###current_exons = mapping$exons[mapping$ensgene == "ENSMUSG00000005763"]
###num_splice_forms = mapping$num_splice_forms[mapping$ensgene == "ENSMUSG00000005763"][1]
###Exon = character(length(current_exons)*dim(experiment)[1])
###Genes = character(length(current_exons)*dim(experiment)[1])
###Group = character(length(current_exons)*dim(experiment)[1])
###Sample = character(length(current_exons)*dim(experiment)[1])
###Counts = integer(length(current_exons)*dim(experiment)[1])
###k = 1
###d2 = d[d$id %in% current_exons,]
###
###for (exon in current_exons) {
###  if (!is.na(as.numeric(d2[d2$id==exon,experiment[,2]])[1])) {
###    Exon[k:(k+dim(experiment)[1]-1)] = exon
###    Genes[k:(k+dim(experiment)[1]-1)] = ensgene_name
###    Group[k:(k+dim(experiment)[1]-1)] = experiment[,1]
###    Sample[k:(k+dim(experiment)[1]-1)] = experiment[,2]
###    Counts[k:(k+dim(experiment)[1]-1)] = as.numeric(d2[d2$id==exon,experiment[,2]]) #log(as.numeric(d2[d2$id==exon,experiment[,2]])/(as.numeric(mapping$exons_length[mapping$exons == exon])/100)+1)
###    k = k+dim(experiment)[1]
###  }
###  
###}
###Exon = Exon[Exon != '']
###Group = Group[Group != '']
###Genes = Genes[1:length(Group)]
###Counts = Counts[1:length(Group)]
###Sample = Sample[1:length(Group)]
####print(Exon)
####print(Group)
####print(Counts)
####allGenes = data.frame(Exon,Genes,Group,Sample,Counts, stringsAsFactors=FALSE)
####allGenes = allGenes[allGenes$Exon != '',]
###
#### TODO
####print("call faith")
###if (length(Exon)== 0 || length(unique(Exon))==1) {
###  pExon[i] = NA
###  pGroup[i] = NA
###  pInteraction[i] = NA
###  NumberOfExons[i] = length(unique(Exon))
###} else {
###  res = twowayanova(data.frame(Exon,Group,Counts, stringsAsFactors=FALSE))
###  #pExon[i] = res[1]
###  pExons = rep(res[1],length(Group))
###  #pGroup[i] = res[2]
###  pGroups = rep(res[2],length(Group))
###  #pInteraction[i] = res[3]
###  pInteractions = rep(res[3],length(Group))
###  #NumberOfExons[i] = length(unique(Exon))
###  data.frame(Exon,Group,Counts,Genes,Sample,pExons,pGroups,pInteractions,num_splice_forms)
###
###  #result = call_FAITHS_function(allGenes)
###  #if (first == TRUE) {
###  #  # TODO: this file should be input too
###  #  write.table(file = "results_RofevsCtrl.csv",data.frame(Exon,Group,Counts,Genes,Sample,pExons,pGroups,pInteractions),row.names = FALSE,sep = ",")
###  #  #print("NINA")
###  #  first = FALSE
###  #} else {
###  #  #print("LALALALALALA")
###  #  write.table(file = "results_RofevsCtrl.csv",data.frame(Exon,Group,Counts,Genes,Sample,pExons,pGroups,pInteractions),append = TRUE,row.names = FALSE,sep = ",",col.names = FALSE)
###  #}
###}#
####if (i == 100) {
####  break
####}
###
# #####

run_everything <- function(ensgene) {
  
  
  #ensname = character(length(ensgene))
  pExon = double(length(ensgene)) 
  pGroup = double(length(ensgene))
  pInteraction =double(length(ensgene))
  NumberOfExons = integer(length(ensgene))
  NumberOfSpliceforms = integer(length(ensgene))
  
  i = 1
  first = TRUE
  for (ensgene_name in ensgene) {
    current_exons = mapping$exons[mapping$ensgene == ensgene_name]
    #print(head(current_exons))
    current_num_splice_forms = mapping$num_splice_forms[mapping$ensgene == ensgene_name][1]
    #current_exons = mapping[mapping$ensgene == ensgene_name,1]
    Exon = character(length(current_exons)*dim(experiment)[1])
    Genes = character(length(current_exons)*dim(experiment)[1])
    Group = character(length(current_exons)*dim(experiment)[1])
    Sample = character(length(current_exons)*dim(experiment)[1])
    Counts = integer(length(current_exons)*dim(experiment)[1])
    k = 1
    d2 = d[d$id %in% current_exons,]
    #print(head(d2))
    for (exon in current_exons) {
      if (!is.na(as.numeric(d2[d2$id==exon,experiment[,2]])[1])) {
        Exon[k:(k+dim(experiment)[1]-1)] = exon
        Genes[k:(k+dim(experiment)[1]-1)] = ensgene_name
        Group[k:(k+dim(experiment)[1]-1)] = experiment[,1]
        Sample[k:(k+dim(experiment)[1]-1)] = experiment[,2]
        #print(exon)
        #print(length(d2[d2$id==exon,experiment[,2]]))
        #print(length(mapping$exons_length[mapping$exons == exon]))
        #if(length(mapping$exons_length[mapping$exons == exon])!= 1) {
        #  print("aaaaaaa")
        #  return()
        #}
        if (args$normalize == "normalize") {
          Counts[k:(k+dim(experiment)[1]-1)] = 
            log(as.numeric(d2[d2$id==exon,experiment[,2]])/
                  (as.numeric(mapping$exons_length[mapping$exons == exon][1])/100)+1) 
        } else {
          Counts[k:(k+dim(experiment)[1]-1)] = as.numeric(d2[d2$id==exon,experiment[,2]])
        }
        k = k+dim(experiment)[1]
      }
      
    }
    Exon = Exon[Exon != '']
    Group = Group[Group != '']
    Genes = Genes[1:length(Group)]
    Counts = Counts[1:length(Group)]
    Sample = Sample[1:length(Group)]
    #print(Exon)
    #print(Group)
    #print(Counts)
    #allGenes = data.frame(Exon,Genes,Group,Sample,Counts, stringsAsFactors=FALSE)
    #allGenes = allGenes[allGenes$Exon != '',]
    
    # TODO
    #print("call faith")
    if (length(Exon)== 0 || length(unique(Exon))==1) {
      pExon[i] = NA
      pGroup[i] = NA
      pInteraction[i] = NA
      NumberOfExons[i] = length(unique(Exon))
      NumberOfSpliceforms[i] = current_num_splice_forms
    } else {
      res = twowayanova(data.frame(Exon,Group,Counts, stringsAsFactors=FALSE))
      #print(head(res))
      pExon[i] = res[1]
      pExons = rep(res[1],length(Group))
      pGroup[i] = res[2]
      pGroups = rep(res[2],length(Group))
      pInteraction[i] = res[3]
      pInteractions = rep(res[3],length(Group))
      NumberOfExons[i] = length(unique(Exon))
      NumberOfSpliceforms[i] = current_num_splice_forms
      #result = call_FAITHS_function(allGenes)
      #if (first == TRUE) {
      #  # TODO: this file should be input too
      #  write.table(file = "results_RofevsCtrl.csv",data.frame(Exon,Group,Counts,Genes,Sample,pExons,pGroups,pInteractions),row.names = FALSE,sep = ",")
      #  #print("NINA")
      #  first = FALSE
      #} else {
      #  #print("LALALALALALA")
      #  write.table(file = "results_RofevsCtrl.csv",data.frame(Exon,Group,Counts,Genes,Sample,pExons,pGroups,pInteractions),append = TRUE,row.names = FALSE,sep = ",",col.names = FALSE)
      #}
    }
    #if (i == 100) {
    #  break
    #}
    if (i %% 1000 == 0 ) {
      print(i)
    }
    i = i +1
    
  }
  
  results = data.frame(
    ensgene,
    pExon, 
    pGroup,
    pInteraction,
    NumberOfExons,
    NumberOfSpliceforms,
    stringsAsFactors=FALSE
  )
  return(results)
}

results = run_everything(ensgene)

write.csv(results, args$result)
