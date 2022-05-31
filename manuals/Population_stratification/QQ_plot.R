library("qqman")

args <- commandArgs(trailingOnly = TRUE)

# check if both input and output file names are provided 
if (length(args)<2)
{
  stop("Please provided input and output file names", call.=FALSE)
}

#input and output file names
input <- args[1]
output <- args[2]

gwasResults <- read.table(input, head=TRUE)
jpeg(output, width=5, height=5, units="in", res=200)

qq(gwasResults$P, main = "Q-Q plot of GWAS p-values", col = "blue4")
dev.off()

