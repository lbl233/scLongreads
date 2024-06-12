# eQTL using long read data

# Bolun Li
# Jun 11 2024


library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
library(forcats)
library(Matrix)


gff = as.data.frame(rtracklayer::import('/data/Choi_lung/lbl/refdata-gex-GRCh38-2020-A/genes/genes.gtf'))
seur <- subset(lr, Celltype == "AT2")
cell_number <- table(seur$Sample)
ind_filtered <- names(cell_number)[which(cell_number < 5)]
table(seur$Sample)
library("tidyverse")
`%nin%` <- Negate(`%in%`)
seur <- subset(seur, Sample %nin% ind_filtered)
count_mtx <- seur@assays$RNA@counts

#extracting the raw counts and metadata
counts = lr@assays$RNA$counts
meta = lr@meta.data[which(lr$Celltype == "AT2"),]
dim(counts)
#only keeping genes that is expressed in at least 10% of all cells
count_mtx = count_mtx[rowSums(count_mtx >0) >= (ncol(count_mtx)/10),]
dim(count_mtx)
#keeping genes with mean counts of greater than 0.1
# counts = counts[rowMeans(counts) > 0.1,]
#log normalization of the counts
# counts = NormalizeData(counts)

group <- as.character(seur$Sample) %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx 
count_mtx_sum <- t(count_mtx_sum)
count_mtx_sum = NormalizeData(count_mtx_sum)
gene_gff <- subset(gff, type == "gene")
rm(gff)
pos <- data.frame(X = gene_gff$gene_id,
                  chr = gene_gff$seqnames, 
                  start = gene_gff$start, 
                  end = gene_gff$start,
                  gene_id = gene_gff$gene_name)
pos_sub <- subset(pos, gene_id %in% rownames(count_mtx_sum))
pos_sub <- pos_sub[-which(duplicated(pos_sub$gene_id)),]
rownames(pos_sub) <- pos_sub$gene_id

count_mtx_sum <- count_mtx_sum[rownames(pos_sub),]
count_mtx_sum@Dimnames[[1]] <- pos_sub$X


centered_mtx <- scale(t(count_mtx_sum))
dim(centered_mtx)
which(centered_mtx=="NaN")
centered_mtx[centered_mtx == "NaN"] <- 0
library(limma)
quantile.norm.mtx <- normalizeQuantiles(t(centered_mtx))

colnames(tss) <- c("#chr", "start", "end", "gene_id")
tss_sub <- tss[which(tss$gene_id %in% rownames(quantile.norm.mtx)),c(1:4)]
bed <- cbind(tss_sub,quantile.norm.mtx[(tss_sub$gene_id),])
write.table(bed, file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/phenotype_gene.bed", sep='\t', quote=F, row.names=FALSE, col.names = T, eol = '\n')

#making expression file for PEER
make_expression_table = function(m){
  ta = m[,-c(1,2,3)]
  t2 = data.frame(t(ta[-1]))
  colnames(t2) = ta[,1]
  return(t2)
}

AT2_expression_peer = make_expression_table(bed)

#saving expression file for PEER
write.table(AT2_expression_peer, file = 'AT2_expression_peer.csv', sep = ',', quote = F, row.names = F, col.names = F)
write.table(AT2_expression_peer, file = '/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/gene_peer.csv', sep = ',', quote = F, row.names = F, col.names = F)

bed <- read.table(gzfile("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/phenotype.bed.gz"),sep="\t", header = FALSE)

AT2_iso20_age_rst <- read.table(gzfile("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_iso20/Standardized_quantiled_age_batch_iso20_5PF.cis_qtl.txt.gz"),sep="\t", header = TRUE)
length(which(AT2_iso20_age_rst$qval < 0.05))
AT2_iso20_age_rst_sub <- subset(AT2_iso20_age_rst, qval < 0.05)

library(ACAT)
rownames(AT2_iso20_age_rst) <- AT2_iso20_age_rst$phenotype_id
rownames(pos) <- pos$X
pos_sub <- subset(pos, X %in% bed$V4)
which(pos_sub$X != bed$V4)
groups <- pos_sub[,c(1,5)]
groups <- groups[order(groups$gene_id),]
write.table(groups,file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/phenotype_groups.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
gene.list <- unique(pos_sub$gene_id)
tmp <- pos_sub %>% group_by(gene_id) %>% summarise(count = length(unique(start)))
table(tmp$count)
tss <- read.table("/data/Choi_lung/scLongreads/TALON_workspace/test10s/tss_gene.bed", sep = "\t", header = FALSE)

file_paths <- list.files("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_gene/", 
                         pattern = "cis_qtl.txt.gz")
rst.list <- list()
for (i in c(1:length(file_paths))) {
  file_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_gene/", file_paths[i])
  covs_used <- str_split_fixed(str_split_fixed(file_paths[i], "__", n = 2)[,2], "\\.", n = 2)[,1]
  rst <- read.table(gzfile(file_path),sep="\t", header = TRUE)
  rst <- subset(rst, qval < 0.05)
  rst.list[[covs_used]] <- rst
}

# short-read
rst <- read.table(gzfile("/data/Choi_lung/TTL/tensor/Output_3PC_PEER/AT2_40PF.cis_qtl.txt.gz"),sep="\t", header = TRUE)
rst_sub <- subset(rst, qval < 0.05)
length(intersect(rst.list$`30PF`$phenotype_id, rst_sub$phenotype_id))
tmp <- left_join(rst.list$`30PF`, rst_sub, by = "phenotype_id")
length(which(tmp$variant_id.x == tmp$variant_id.y))
tmp_sub <- tmp[which(tmp$variant_id.x == tmp$variant_id.y),]
tmp_sub$gene_name <- pos_sub[rownames(tmp_sub),]$gene_id
tmp2 <- Natri_AT2[which(Natri_AT2$feature_id %in% tmp$phenotype_id),]
tmp2 <- subset(tmp2, lfsr < 0.1)
View(tmp2[which(tmp2$snp_rsid %in% tmp$variant_id.x),])

df_sum <- data.frame(covs = names(rst.list),
                     eGenes = sapply(rst.list, nrow))

df_sum$covs <- factor(df_sum$covs, levels = paste0(seq(5,30,5), "PF"))
ggplot(df_sum, aes(x=covs, y=eGenes, fill = covs)) +
  geom_bar(stat="identity") + xlab("covs")+ylab("Number of eGenes") +
  geom_text(aes(label = eGenes), vjust=-1)+
  theme_classic()
file_paths <- list.files("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_iso20/", 
                         pattern = "cis_qtl.txt.gz")
file_paths <- file_paths[1:12]
rst.sig.iso.list <- list()
rst.iso.list <- list()
for (i in c(1:length(file_paths))) {
  file_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_iso20/", file_paths[i])
  covs_used <- str_split_fixed(str_split_fixed(file_paths[i], "_", n = 6)[,6], "\\.", n = 2)[,1]
  rst <- read.table(gzfile(file_path),sep="\t", header = TRUE)
  rst.iso.list[[covs_used]] <- rst
  rst <- subset(rst, qval < 0.05)
  rst.sig.iso.list[[covs_used]] <- rst
}

library(ACAT)
pos_sub <- subset(pos, X %in% rst.iso.list[["5PF"]]$phenotype_id)
phenotype_groups <- pos_sub[,c(1,5)]
colnames(phenotype_groups) <- c("phenotype_id", "gene_id")
QTL_rst_mtx <- rst.iso.list[["5PF"]]
nominate_sGene_by_ACAT <- function(QTL_rst_mtx,     # should contain pval_perm
                                   phenotype_groups # should contain phenotype id and group (gene) id
){
  work_mtx <- left_join(phenotype_groups, QTL_rst_mtx, by = "phenotype_id")
  rst<- work_mtx %>% group_by(gene_id) %>% mutate(pval_aggr_ACAT = ACAT(pval_perm))
  rst <- distinct(rst, gene_id, .keep_all = TRUE)
  rst$FDR <- p.adjust(rst$pval_aggr_ACAT, method = "BH")
  return(rst[,c("gene_id","pval_aggr_ACAT", "FDR")])
}

test <- lapply(rst.iso.list, function(x) nominate_sGene_by_ACAT(x,phenotype_groups))
sig_sGene <- sapply(test, function(x){
  rst <- subset(x, pval_aggr_ACAT < 0.05)
  return(nrow(rst))
})
df_sGene <- data.frame(covs = names(sig_sGene),
                       sGenes = sig_sGene)

df_sGene$covs <- factor(df_sGene$covs, levels = paste0(seq(5,60,5), "PF"))
ggplot(df_sGene, aes(x=covs, y=sGenes, fill = covs)) +
  geom_bar(stat="identity") + xlab("covs")+ylab("Number of sGenes") +
  geom_text(aes(label = sGenes), vjust=-1)+
  theme_classic()
df_sum <- data.frame(covs = names(rst.sig.iso.list),
                     Isoforms = sapply(rst.sig.iso.list, nrow))

df_sum$covs <- factor(df_sum$covs, levels = paste0(seq(5,60,5), "PF"))
ggplot(df_sum, aes(x=covs, y=Isoforms, fill = covs)) +
  geom_bar(stat="identity") + xlab("covs")+ylab("Number of significant Isoforms") +
  geom_text(aes(label = Isoforms), vjust=-1)+
  theme_classic()
