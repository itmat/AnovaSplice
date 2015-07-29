setwd("~/Documents/DifferentialSplicing/")

d = read.csv("FINAL_master_list_of_exons_counts_MIN.sense.Kid.txt",sep = "\t", header = TRUE, na.strings = "",stringsAsFactors=FALSE)

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

experiment = data.frame(Group = c(rep("Control",4),rep("Treatment",4)), SampleName = colnames(d)[2:9],
                        stringsAsFactors=FALSE)

k = 1
Exons = character(dim(d)[1]*50)
Genes = character(dim(d)[1]*50)
Group = character(dim(d)[1]*50)
Sample = character(dim(d)[1]*50)
Value = integer(dim(d)[1]*50)
for (i in seq(dim(d)[1])) {
#for (i in seq(10)) {
  if (is.na(d$name[i])) {
    print("NONE")
  } else {
    #all_genes[k,] = c("1","2","3","4",5)
    f = strsplit(d$name[i],",")
    for (name in unique(f[[1]])) {
      for (row in seq(dim(experiment)[1])) {
        Exons[k] = d$id[i]
        Genes[k] = name
        Group[k] = experiment[row,1]
        Sample[k] = experiment[row,2] 
        Value[k] = d[i,experiment[row,2]]
        k = k+1
      }
      
    }
  }
}

allGenes = data.frame(Exons,Genes,Group,Sample,Value, stringsAsFactors=FALSE)



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
