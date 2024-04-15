# Using SUPPA2 to identify alternative splicing events

# Bolun Li
# Apr 10 2024


.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
library(forcats)
library(Matrix)

# load transcript metadata generated from TALON and integrated with TranDecoder results
TALON_afterqc_orf_secondpass <- readRDS("/data/Choi_lung/scLongreads/Seurat/TALON_afterqc_orf_secondpass.rds")

# load the final version with cell type annotation
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/lr_final.RDS")

# load the original isoform level expression profiles
lr.isoform <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.Isoform.22samples.combined.rds")
ncol(lr)
lr.isoform.sub <- subset(lr.isoform, cell = colnames(lr))
lr.isoform.sub$Sample <- lr$Sample
lr.isoform.sub$Sample_NCI <- lr$Sample_NCI
lr.isoform.sub$Annotation_recluster <- lr$Annotation_recluster
count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass$transcript_name_unique,]
count_mtx@Dimnames[[1]] <- TALON_afterqc_orf_secondpass$annot_transcript_id
dim(count_mtx)
sum(count_mtx)

lr.isoform.sub$Annotation_recluster <- as.character(lr.isoform.sub$Annotation_recluster)
lr.isoform.sub$Celltype <- as.character(lr.isoform.sub$Annotation_recluster)
lr.isoform.sub$res_used <- lr$res_used
lr.isoform.sub$Celltype[which(lr.isoform.sub$Celltype == "Plasma cell/Plasmacytoid DCs")] <- "Plasmacytoid DCs"
lr.isoform.sub$Celltype[which(lr.isoform.sub$Celltype == "Adventitial/peribronchial fibroblasts")] <- "Adventitial_peribronchial fibroblasts"
lr.isoform.sub$Celltype[which(lr.isoform.sub$res_used == 33)] <- "Plasma cells"
table(lr.isoform.sub$Celltype)
lr.isoform.sub$Celltype_Sample <- paste(lr.isoform.sub$Celltype, 
                                        lr.isoform.sub$Sample_NCI, sep = "_")
group <- lr.isoform.sub$Celltype_Sample %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)
x <- as.matrix(count_mtx_sum)
tpm_mtx <- t(t(x)*1e6/colSums(x))
rm(x)
celltype_list <- unique(lr.isoform.sub$Celltype)
head(colnames(tpm_mtx))
for (celltype in celltype_list) {
  idx <- grep(celltype, colnames(tpm_mtx))
  mtx <- tpm_mtx[,idx]
  write.table(mtx,file=paste0("/data/Choi_lung/scLongreads/TALON_workspace/test10s/iso_tpm_formatted_",gsub(" ", "_", celltype),".txt"), sep = "\t", quote = FALSE)
}
