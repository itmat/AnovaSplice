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


args <- commandArgs(TRUE)

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
#key <- "width"
#value <- 32

#mylist <- list()
#mylist[[ key ]] <- value

k = 1
for (i in seq(dim(gtf_file)[1])) {
  #for (i in seq(100)) {
  if (gtf_file$V3[i] == "exon") {
    gene_id = strsplit(gtf_file$V9[i], "[; ]")[[1]][2]
    name = paste("exon:", gtf_file$V1[i], ":", gtf_file$V4[i], "-", gtf_file$V5[i], sep = "")
    ensgene[k] = gene_id
    exons[k] = name
    k = k+1
  }
}

mapping = unique(data.frame(exons,ensgene,stringsAsFactors=FALSE))
ensgene = unique(ensgene)
head(ensgene)
rm(exons)

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

run_everything <- function(ensgene) {
  
  
  #ensname = character(length(ensgene))
  pExon = double(length(ensgene)) 
  pGroup = double(length(ensgene))
  pInteraction =double(length(ensgene))
  NumberOfExons = integer(length(ensgene))
  
  i = 1
  first = TRUE
  for (ensgene_name in ensgene) {
    current_exons = mapping[mapping$ensgene == ensgene_name,1]
    Exon = character(length(current_exons)*dim(experiment)[1])
    Genes = character(length(current_exons)*dim(experiment)[1])
    Group = character(length(current_exons)*dim(experiment)[1])
    Sample = character(length(current_exons)*dim(experiment)[1])
    Counts = integer(length(current_exons)*dim(experiment)[1])
    k = 1
    d2 = d[d$id %in% current_exons,]
    
    for (exon in current_exons) {
      if (!is.na(as.numeric(d2[d2$id==exon,experiment[,2]])[1])) {
        Exon[k:(k+dim(experiment)[1]-1)] = exon
        Genes[k:(k+dim(experiment)[1]-1)] = ensgene_name
        Group[k:(k+dim(experiment)[1]-1)] = experiment[,1]
        Sample[k:(k+dim(experiment)[1]-1)] = experiment[,2]
        Counts[k:(k+dim(experiment)[1]-1)] = as.numeric(d2[d2$id==exon,experiment[,2]])
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
    } else {
      res = twowayanova(data.frame(Exon,Group,Counts, stringsAsFactors=FALSE))
      pExon[i] = res[1]
      pExons = rep(res[1],length(Group))
      pGroup[i] = res[2]
      pGroups = rep(res[2],length(Group))
      pInteraction[i] = res[3]
      pInteractions = rep(res[3],length(Group))
      NumberOfExons[i] = length(unique(Exon))
      #result = call_FAITHS_function(allGenes)
      if (first == TRUE) {
        write.table(file = "results_RofevsCtrl.csv",data.frame(Exon,Group,Counts,Genes,Sample,pExons,pGroups,pInteractions),row.names = FALSE,sep = ",")
        #print("NINA")
        first = FALSE
      } else {
        #print("LALALALALALA")
        write.table(file = "results_RofevsCtrl.csv",data.frame(Exon,Group,Counts,Genes,Sample,pExons,pGroups,pInteractions),append = TRUE,row.names = FALSE,sep = ",",col.names = FALSE)
      }
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
    stringsAsFactors=FALSE
  )
  return(results)
}

results = run_everything(ensgene)

write.csv(results, args$result)


