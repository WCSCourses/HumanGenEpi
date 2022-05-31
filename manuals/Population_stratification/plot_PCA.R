# Plot PCA

eigenvec_table=read.table("demo_data.pca.input", header = TRUE)
jpeg("pca_plot.jpeg", width=5, height=5, units="in", res=200)


palette(c("blue4","orange3"))
case=which(eigenvec_table$Pheno=="Case")
control=which(eigenvec_table$Pheno=="Control")
plot(eigenvec_table[3:4],pch=c(4,4),col=1:2,xlab="PC1", ylab="PC2")
points(eigenvec_table$PC1[case],eigenvec_table$PC2[case],pch=4,col=1)
points(eigenvec_table$PC1[control],eigenvec_table$PC2[control],pch=4,col=2)
legend("topleft",pch = c(4,4),col=1:2,legend=c("Case","Control"))
dev.off()
