library("qqman")

args <- commandArgs(trailingOnly = TRUE)

assoc <- read.table(args[1], header=T)

png(paste(args[2],".png",sep=""), width=720, height=400)
manhattan(assoc, suggestiveline=F)
dev.off()
