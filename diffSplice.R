source("/Users/hayer/Documents/DifferentialSplicing/AnovaSplice/anovasplicecalc.R")

setwd("~/Documents/DifferentialSplicing/")
rm(list=ls(all=TRUE))
##### READ GTF FILE ######
gtf_file = read.csv("fixed_mm9_ensembl_genes.gtf",sep = "\t", header = FALSE, na.strings = "",stringsAsFactors=FALSE)
max_num_exons = as.integer(system("cut -f 3 fixed_mm9_ensembl_genes.gtf | grep -w exon| wc -l",intern = TRUE))
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
rm(exons)
# unique(mapping[mapping$exons == "exon:chr1:8435108-8435200",2])
d = read.csv(pipe("cut -f 1-17 FINAL_master_list_of_exons_counts_MIN.emanuela_NEW.txt"),sep = "\t", header = TRUE, na.strings = "",stringsAsFactors=FALSE)

experiment = data.frame(Group = c(rep("Control",8),rep("Treatment",8)), SampleName = colnames(d)[2:17],
                        stringsAsFactors=FALSE)

results = data.frame(
  ensname = character(length(ensgene)),
  pExon = double(length(ensgene)), 
  pGroup = double(length(ensgene)),
  pInteraction =double(length(ensgene)),
  stringsAsFactors=FALSE
)

i = 1
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
  Counts = Counts[1:length(Group)]
  #allGenes = data.frame(Exon,Genes,Group,Sample,Counts, stringsAsFactors=FALSE)
  #allGenes = allGenes[allGenes$Exon != '',]
  
  # TODO
  #print("call faith")
  results[i,] = c(ensgene_name, twowayanova(data.frame(Exon,Group,Counts, stringsAsFactors=FALSE)))
  #result = call_FAITHS_function(allGenes)
  if (i == 1000) {
    break
  }
  i = i +1

}



head(d)

# http://stackoverflow.com/questions/11561856/add-new-row-to-dataframe
insertRow <- function(existingDF, newrow, r) {
  #existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

all_genes = data.frame(
  Exons = character(dim(d)[1]),
  Genes = character(dim(d)[1]),
  Group = character(dim(d)[1]),
  Sample = character(dim(d)[1]),
  Value = integer(dim(d)[1]),
  stringsAsFactors=FALSE
)

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

experiment = data.frame(Group = c(rep("Control",4),rep("Treatment",4)), SampleName = colnames(d)[2:9],
                        stringsAsFactors=FALSE)

k = 1
Exons = character(dim(d)[1])
Genes = character(dim(d)[1])
Group = character(dim(d)[1])
Sample = character(dim(d)[1])
Counts = integer(dim(d)[1])
#for (i in seq(dim(d)[1])) {
for (i in seq(1000)) {
  if (is.wholenumber(i/100)) {print(i) }
  if (dim(mapping[mapping$exons == d$id[i],])[1] == 0) {
    #print("NONE")
  } else {
    #all_genes[k,] = c("1","2","3","4",5)
    #f = strsplit(d$name[i],",")
    for (name in unique(mapping[mapping$exons == d$id[i],2])) {
      #for (row in seq(dim(experiment)[1])) {
      #  Exons[k] = d$id[i]
      #  Genes[k] = name
      #  Group[k] = experiment[row,1]
      #  Sample[k] = experiment[row,2] 
      #  Counts[k] = d[i,experiment[row,2]]
      #  k = k+1
      #}
      Exons[k:(k+dim(experiment)[1]-1)] = d$id[i]
      Genes[k:(k+dim(experiment)[1]-1)] = name
      Group[k:(k+dim(experiment)[1]-1)] = experiment[,1]
      Sample[k:(k+dim(experiment)[1]-1)] = experiment[,2]
      Counts[k:(k+dim(experiment)[1]-1)] = as.numeric(d[i,experiment[,2]])
      k = k+dim(experiment)[1]
      
    }
  }
}

allGenes = data.frame(Exons,Genes,Group,Sample,Counts, stringsAsFactors=FALSE)
allGenes = allGenes[allGenes$Exons != '',]


niner <- function(row_d) {
  
  all_genes <<- insertRow(all_genes,c("1","2","3","4",5),dim(all_genes)[1]+1)  
}


all_genes = rbind(all_genes,c(1,2,3,4,5))
all_genes = insertRow(all_genes,c("1","2","3","4",5),dim(all_genes)[1]+1)


(apply(d, 1, niner))

x <- c(as = "asfef", qu = "qwerty", "yuiop[", "b", "stuff.blah.yech")
# split x on the letter e
strsplit(x, "e")


x <- c(3:5, 11:8, 8 + 0:5)
(ux <- unique(x))


head(c(d$id, d$name))
apply( 2, print)


d_new = apply(d, 2, strsplit, d$name,",")
