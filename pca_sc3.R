source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
install.packages("mnormt")
library(pcaMethods)
library(mnormt)
install.packages("~/GITCLONES/GSE72056-melanoma-single-cell-rnaseq/pcaReduce-master/pcaReduce_1.0.tar.gz", repos = NULL, type = "source")
library(pcaReduce)

source("~/GITCLONES/GSE72056-melanoma-single-cell-rnaseq/load.R")

input_data <- t(gse72056_expr)

step1 <- prep(input_data, scale="none", center=TRUE)
#[1] "svd"          "nipals"       "rnipals"      "bpca"         "ppca"         "svdImpute"    "robustPca"   
#[8] "nlpca"        "llsImpute"    "llsImputeAll"
step1_pca <- pca(step1, method="svd", center = FALSE, nPcs=2)
summary(step1_pca)
step1_score <- step1_pca@scores
step1_loading <- step1_pca@loadings


tumor_name <- gse72056_tumor_name # k=19
cell_types <- gse72056_cell_type   # k =7


###plot
# setwd("~/GITCLONES/GSE72056-melanoma-single-cell-rnaseq")
# pdf("test.pdf")
# C <- c("skyblue", "mediumvioletred", "olivedrab3", "aquamarine4", "salmon1",   "orangered3",   "darkseagreen3"  )
plot(step1_score[which(cell_types=="B"),1], step1_score[which(cell_types=="B"),2], type="p",pch=21, bg = "salmon1", col= "salmon2", 
     cex=1.7, xlab="PCA1", ylab="PCA2",ylim=c(-150,200), xlim=c(-220,320),  cex.lab=1.5, cex.axis=1.7, font.lab=2, bty='n')
# points(step1_score[which(cell_types=="CAF"),1], step1_score[which(cell_types=="CAF"),2], type="p",pch=21, bg = "orangered3", col="orangered4",cex=1.7)
# points(step1_score[which(cell_types=="Endothelial"),1], step1_score[which(cell_types=="Endothelial"),2], type="p",pch=21, bg = "plum3", col = "plum4", cex=1.7)
# 
# points(step1_score[which(cell_types=="Macrophage"),1], step1_score[which(cell_types=="Macrophage"),2], type="p",pch=21, bg = "skyblue1", col= "skyblue2", cex=1.7)
# 
# points(step1_score[which(cell_types=="NK"),1], step1_score[which(cell_types=="NK"),2], type="p",pch=21, bg = "olivedrab2", col= "olivedrab3", cex=1.7)
# 
# points(step1_score[which(cell_types=="T"),1], step1_score[which(cell_types=="T"),2], type="p", pch=21, bg = "darkseagreen", col= "darkseagreen4", cex=1.7)
# 
# points(step1_score[which(cell_types=="Tumor"),1], step1_score[which(cell_types=="Tumor"),2], type="p", pch=21, col = "aquamarine3", cex=1.7)
# 
# 
# temp <- legend("topright",  border = NULL,legend = c("B","CAF", "Endothelial", "Macrophage", "NK", "T", "Tumor"), col = c("salmon2", "orangered4", "plum4", "skyblue2", "olivedrab3", "darkseagreen4", "aquamarine3" ),
#                text.width = strwidth("00,000,0"),
#                lty = c(0), pch=c(21,21,21,21,21,21,21), pt.bg=c("salmon2", "orangered4", "plum4", "skyblue2", "olivedrab3", "darkseagreen4", "aquamarine3"), xjust = 0, yjust = 0,bty = "n")
# 




library(pcaMethods)
library(pcaReduce)

biocLite("SC3")
biocLite("scater")
library(pheatmap)
set.seed(1234567)


pollen <- readRDS("~/Desktop/pollen.rds")
pollen



a <- gse72056_expr[apply(gse72056_expr
                         ,1,median) > 0,]
hist(colSums(a),breaks = 500)


# cell annotation
ann <- data.frame(cell_type1 = cell_types)
pd <- new("AnnotatedDataFrame", data = ann)
# cell expression
tmp <- a
colnames(tmp) <- rownames(ann)
# SCESEt object
setwd("~/GITCLONES/GSE72056-melanoma-single-cell-rnaseq")
sceset <- newSCESet(fpkmData = tmp, phenoData = pd, logExprsOffset = 1)

is_exprs(sceset) <- exprs(sceset) > 0.1
sceset <- calculateQCMetrics(sceset)



# pollen <- sc3_prepare(sceset, ks = 2:5)
# pollen <- sc3_estimate_k(pollen)
# pollen@sc3$k_estimation #27



sceset <- sc3(sceset, ks = 27, biology = TRUE, n_cores = 1)
save(sceset,file = "sceset.RData")
plotPCA(sceset, colour_by = "cell_type1")
plotPCA(sceset, colour_by = "sc3_27_clusters")

sc3_plot_consensus(sceset, k = 27, show_pdata = "cell_type1")
p_data <- pData(sceset)

setwd("~/GITCLONES/GSE72056-melanoma-single-cell-rnaseq")
table(sceset@phenoData@data$sc3_27_clusters)
table(sceset@phenoData@data$cell_type1)
pheno_df <- as.data.frame(sceset@phenoData@data)
pheno_df[pheno_df$sc3_27_clusters == 2,1]



#######################################################
setwd("~/GITCLONES/GSE72056-melanoma-single-cell-rnaseq")
sceset_cell <- newSCESet(fpkmData = tmp, phenoData = pd, logExprsOffset = 1)

is_exprs(sceset_cell) <- exprs(sceset_cell) > 0.1
sceset_cell <- calculateQCMetrics(sceset_cell)



# pollen <- sc3_prepare(sceset, ks = 2:5)
# pollen <- sc3_estimate_k(pollen)
# pollen@sc3$k_estimation #27



sceset_cell <- sc3(sceset_cell, ks = 7, biology = TRUE, n_cores = 1)
save(sceset_cell,file = "sceset_cell.RData")
plotPCA(sceset_cell, colour_by = "cell_type1")
plotPCA(sceset_cell, colour_by = "sc3_7_clusters")

table(sceset_cell@phenoData@data$sc3_7_clusters)
table(sceset_cell@phenoData@data$cell_type1)
pheno_df_cell <- as.data.frame(sceset_cell@phenoData@data)
table(pheno_df_cell[pheno_df_cell$sc3_7_clusters == 1,1])
table(pheno_df_cell[pheno_df_cell$sc3_7_clusters == 2,1])
table(pheno_df_cell[pheno_df_cell$sc3_7_clusters == 3,1])
table(pheno_df_cell[pheno_df_cell$sc3_7_clusters == 4,1])
table(pheno_df_cell[pheno_df_cell$sc3_7_clusters == 5,1])
table(pheno_df_cell[pheno_df_cell$sc3_7_clusters == 6,1])
table(pheno_df_cell[pheno_df_cell$sc3_7_clusters == 7,1])



sc3_plot_consensus(sceset_cell, k = 7, show_pdata = "cell_type1")

testde <- get_de_genes(gse72056_expr,gse72056_cell_type)