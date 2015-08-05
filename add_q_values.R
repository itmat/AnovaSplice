# Adjust P-values
# devtools::install_github("trevorld/argparse")
library(argparse)

parser <- ArgumentParser(description='Add q-values to output of AnovaSplice',prog = "Rscript add_q_values.R")

parser$add_argument('file', metavar='output.csv', type="character", nargs=1,
                    help='Output of AnovaSplice')

parser$add_argument('file2', metavar='output_fdr.csv', type="character", nargs=1,
                    help='Output with fdr')


args <- commandArgs(TRUE)

if(length(args) < 1) {
  args <- c("-h")
}

args <- parser$parse_args(args)

print(args)

#parser$print_help()

d = read.csv(args$file)

pInteractionFDR = p.adjust(d$pInteraction, method = "BH") 
pGroupFDR = p.adjust(d$pGroup, method = "BH") 
pExonFDR = p.adjust(d$pExon, method = "BH") 

dnew = cbind(d, pExonFDR, pGroupFDR, pInteractionFDR)

write.csv(dnew,args$file2 )