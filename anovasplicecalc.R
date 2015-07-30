# this function should receive an object of type
# list that has at lease three columns:
# "Exon", "Group", and "Counts"

twowayanova <- function(exons){
# this function calculates a two-way ANOVA for a balanced, unpaired design
# check that required column names present  
if(!("Counts" %in% colnames(exons)) | !("Exon" %in% colnames(exons)) | !("Group" %in% colnames(exons)))
{
  stop("Exon, Group, and Counts must be column headers\n", call. = TRUE);
}

# to-do: if not balanced, do not go further
b <- replications(Counts ~ Exon*Group, data = exons)
if(!(typeof(b) == "integer"))
{
  stop("Experiment appears to be unbalanced, ANOVA will not be run\n", call. = TRUE)
}

# run two-way aov
int <- aov(exons$Counts~exons$Exon*exons$Group)

# Order of returned P-values: Exon, Group, Interaction
return(summary(int)[[1]][["Pr(>F)"]][1:3])
}