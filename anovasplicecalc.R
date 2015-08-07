# this function should receive an object of type
# list that has at lease three columns:
# "Exon", "Group", and "Counts"
# for paired test, a fourth column "Sample" is required

twowayanova <- function(exons, aovtype = "unpaired"){
  # this function calculates a two-way ANOVA for a balanced, unpaired design
  # check that required column names present  
  if(!("Counts" %in% colnames(exons)) | !("Exon" %in% colnames(exons)) | !("Group" %in% colnames(exons)))
  {
    stop("Exon, Group, and Counts must be column headers\n", call. = FALSE);
  }
  
  b <- replications(Counts ~ Exon*Group, data = exons)
  if(!(typeof(b) == "integer"))
  {
    stop("Experiment appears to be unbalanced, ANOVA will not be run\n", call. = FALSE)
  }
  
  # select two-way aov test to run
  if(aovtype == "unpaired")
  { int <- aov(exons$Counts~exons$Exon*exons$Group) 
  pvalues <- summary(int)[[1]][["Pr(>F)"]][1:3]
  } 
  else if(aovtype == "paired")
  {
    # two-way aov for paired
    # group variable is within subjects
    if(!("Sample" %in% colnames(exons)))
    {
      stop("Sample must be a column header for the paired ANOVA\n", call. = FALSE);
    }
    #stop("The paired two-way aov is not implemented yet\n", call. = FALSE)
    int <- aov(exons$Counts ~ exons$Exon*exons$Group + Error(exons$Sample/exons$Group))
    pexon <- summary(int)[[3]][[1]][["Pr(>F)"]][[1]]
    pgroup <- summary(int)[[2]][[1]][["Pr(>F)"]][[1]]
    pint <- summary(int)[[3]][[1]][["Pr(>F)"]][[2]]
    pvalues <- c(pexon, pgroup, pint)
    }

  # Order of returned P-values for unpaired: Exon, Group, Interaction
  return(pvalues)
  #return(int)
}