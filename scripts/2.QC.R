# compare long-read and short-read data

# Bolun Li
# Jan 17 2024


.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
load("/data/Choi_lung/scLongreads/Seurat/lr.sc.data.list.21samples.beforeQC.RData")
lr.data.list <- dataset_list
names(lr.data.list)
rm(dataset_list)
load("/data/Choi_lung/scLongreads/Seurat/sr.sc.data.list.samples.beforeQC.RData")
sr.data.list <- illumina.datalist
rm(illumina.datalist)
names(sr.data.list) <- str_split_fixed(illumina_file_names, pattern = "\\.", n = 2)[,1]
names(sr.data.list)

# Uniform the sample ID
sr.data.list <- sr.data.list[-c(3,5,8,21)]
names(sr.data.list) <- c("16_NCI_105_110", "17_NCI_111_116", "18_NCI_117_123_2", "19_NCI_124_129",
                         "1_NCI_1257", "20_NCI_130_135", "21_NCI_136_140", "NCI_17_22",
                         "3_NCI_23_29", "4_NCI_30_35", "5_NCI_36_41", "6_NCI_42_47", "7_NCI_48_54",
                         "8_NCI_56_61", "9_NCI_62_68", "10_NCI_69_74", "11_NCI_75_80", 
                         "12_NCI_81_86", "13_NCI_87_92", "2_NCI_9_16", "14_NCI_93_98", "15_NCI_99_104")
for (i in 1:22) {
  sr.data.list[[i]]$orig.ident <- names(sr.data.list)[i]
}

# merge data list
lr.seur <- merge(lr.data.list[[1]], y = lr.data.list[2:22], add.cell.ids = names(lr.data.list), project = "Long-read")
length(intersect(colnames(sr.data.list$nci117_123_II), colnames(lr.data.list$`18_NCI_117_123_2`)))
ncol(sr.data.list$nci117_123_II)
ncol(lr.data.list$`18_NCI_117_123_2`)
length(intersect(colnames(sr.data.list$nci1257_II), colnames(lr.data.list$`1_NCI_1257`)))
ncol(sr.data.list$nci1257_II)
lr.seur <- merge(lr.data.list[[1]], y = lr.data.list[2:22], add.cell.ids = names(lr.data.list), project = "Long-read")
sr.seur <- merge(sr.data.list[[1]], y = sr.data.list[2:22], add.cell.ids = names(sr.data.list), project = "Short-read")
VlnPlot(lr.seur, features = c("nCount_RNA", "nFeature_RNA", "doublet_RNA",
                           "nCount_Isoforms", "nFeature_Isoforms", "percent.mt"),
        ncol = 3, pt.size = 0) 
VlnPlot(sr.seur, features = c("nCount_RNA", "nFeature_RNA", "Doublet_score"),
        ncol = 3, pt.size = 0, group.by = "orig.ident") 

table(lr.seur$doublet2_RNA, lr.seur$doublet2_isoform)
table(lr.seur$doublet2_isoform)
table(lr.seur$orig.ident)
table(sr.seur$orig.ident)
lr.seur$barcode <- colnames(lr.seur)
sr.seur$barcode <- colnames(sr.seur)
combined.meta <- left_join(lr.seur@meta.data, sr.seur@meta.data,by = "barcode")
saveRDS(combined.meta, file="/data/Choi_lung/scLongreads/Seurat/combined.meta.rds")
table(combined.meta$Sample)
length(which(is.na(combined.meta$Sample)))
pdf("/data/Choi_lung/scLongreads/Seurat/nCount_compare.pdf", width = 8, height = 8)
plot(combined.meta$nCount_RNA.x,combined.meta$nCount_RNA.y)
plot(combined.meta$nFeature_RNA.x,combined.meta$nFeature_RNA.y, col = combined.meta$orig.ident.x, pch=19)
dev.off()
combined.meta$orig.ident <- "1"
combined.meta$orig.ident <- factor(combined.meta$orig.ident.x, levels = unique(combined.meta$orig.ident.x))


combined.meta$Sample_10X <- as.numeric(combined.meta$orig.ident) 
plot(combined.meta$nFeature_RNA.x,combined.meta$nFeature_RNA.y, col = combined.meta$Sample_10X, pch=19)
legend("bottomright", legend = paste("10X-Sample", 1:22), col = 1:22, pch = 19, bty = "n", text.font = 2, cex = 0.5)

df.lr <- as.data.frame(table(combined.meta$orig.ident.x))
df.overlappled <- as.data.frame(table(combined.meta$orig.ident.y))
df.sr <- as.data.frame(table(sr.seur$orig.ident))
df.lr$Overlapped <- df.overlappled$Freq
df.lr$Overlapped.pct <- (100*df.lr$Overlapped)/df.lr$Freq
df.sr$Overlapped <- df.overlappled$Freq
df.sr$Overlapped.pct <- (100*df.sr$Overlapped)/df.sr$Freq
df.sr$Overlapped.pct <- round(df.sr$Overlapped.pct, 2)
hist(df.sr$Overlapped.pct)
ggplot(df.sr, aes(x=Var1, y=Overlapped.pct)) + 
  geom_bar(stat = "identity")+RotatedAxis()
library(ggpubr)
ggbarplot(df.sr, x = "Var1", y = "Overlapped.pct",
          fill = "lightblue", color = "lightblue", 
          palette = "Paired",
          label = TRUE,
          position = position_dodge(0.9))+RotatedAxis()

# Compare features
lr.seur <- NormalizeData(lr.seur, verbose = FALSE)
sr.seur <- NormalizeData(sr.seur, verbose = FALSE)
gene.overlapped <- intersect(rownames(lr.seur), rownames(sr.seur.annot))
barcode.overlapped <- combined.meta$barcode[-which(is.na(combined.meta$Sample))]
nrow(lr.seur)
nrow(sr.seur)
plot(lr.seur@assays$RNA@data[gene.overlapped,"1_NCI_1257_GATGATCCATCATCCC-1"], 
     sr.seur@assays$RNA@data[gene.overlapped,"1_NCI_1257_GATGATCCATCATCCC-1"])
cor(lr.seur@assays$RNA@data[gene.overlapped,"1_NCI_1257_GATGATCCATCATCCC-1"], 
    sr.seur@assays$RNA@data[gene.overlapped,"1_NCI_1257_GATGATCCATCATCCC-1"])

plot(lr.seur@assays$RNA@data["PTPRC",barcode.overlapped], 
     sr.seur@assays$RNA@data["PTPRC",barcode.overlapped])
cor(lr.seur@assays$RNA@data["PTPRC",barcode.overlapped], 
    sr.seur@assays$RNA@data["PTPRC",barcode.overlapped])

Classification.info <- read.table("/data/Choi_lung/scLongreads/B2_percentile/Sample_2_NCI_9_16/2_NCI_9_16_classification.filtered_lite_classification.txt", sep = "\t", header = TRUE)
Classification.info2 <- read.table("/data/Choi_lung/scLongreads/B2_percentile/Sample_10_NCI_69_74/10_NCI_69_74_classification.filtered_lite_classification.txt", sep = "\t", header = TRUE)

Classification.info.filter <- read.table("/data/Choi_lung/scLongreads/B2_percentile/Sample_2_NCI_9_16/2_NCI_9_16_classification_filtered/2_NCI_9_16_classification_filtered.annotated.info.csv", sep = "\t", header = TRUE)
Isoform.annot <- read.table("/data/Choi_lung/scLongreads/B2_percentile/Sample_2_NCI_9_16/2_NCI_9_16_classification_filtered/isoforms_seurat/genes.tsv", sep = "\t", header = FALSE)
Isoform.annot2 <- read.table("/data/Choi_lung/scLongreads/B2_percentile/Sample_10_NCI_69_74/10_NCI_69_74_classification_filtered/isoforms_seurat/genes.tsv", sep = "\t", header = FALSE)

IRF4.iso <- rownames(lr.seur@assays$Isoforms)[grep("IRF4", rownames(lr.seur@assays$Isoforms))]
DefaultAssay(lr.seur) <- "Isoforms"
DotPlot(lr.seur, features = IRF4.iso)+RotatedAxis()

sr.seur.annot <- readRDS("/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/demuxlet_vireo/UMAP.RDS")
sr.annot.data.list <- SplitObject(sr.seur.annot, split.by = "orig.ident")
rm(sr.seur.annot)
sr.annot.data.list <- sr.annot.data.list[-c(1,3)]
names(sr.annot.data.list)
names(sr.annot.data.list) <- c("1_NCI_1257", "2_NCI_9_16",  "NCI_17_22","3_NCI_23_29", "4_NCI_30_35", "5_NCI_36_41", "6_NCI_42_47",
                                "7_NCI_48_54", "8_NCI_56_61", "9_NCI_62_68", "10_NCI_69_74", "11_NCI_75_80", 
                               "12_NCI_81_86", "13_NCI_87_92",  "14_NCI_93_98", "15_NCI_99_104",
                               "16_NCI_105_110", "17_NCI_111_116", "18_NCI_117_123_2", "19_NCI_124_129",
                               "20_NCI_130_135", "21_NCI_136_140")
for (i in 1:22) {
  sr.annot.data.list[[i]]$orig.ident <- names(sr.annot.data.list)[i]
}
gc()
sr.seur.annot <- merge(sr.annot.data.list[[1]], y = sr.annot.data.list[2:22], add.cell.ids = names(sr.annot.data.list), project = "short-read")
table(sr.seur.annot$orig.ident)
sr.seur.annot$barcode <- paste0(str_split_fixed(colnames(sr.seur.annot), "-", n=2)[,1], "-1")
rm(sr.annot.data.list)
gc()

saveRDS(sr.seur.annot@meta.data, file = "/data/Choi_lung/scLongreads/Seurat/annot.meta.rds")

annot.meta <- left_join(lr.seur@meta.data, sr.seur.annot@meta.data,by = "barcode")
ncol(sr.seur.annot)
df <- annot.meta[-which(is.na(annot.meta$seurat_clusters)),]
nrow(df)
df$reads <- df$nCount_RNA.x/1000
df$reads <- df$nCount_RNA.y/1000

pdf("/data/Choi_lung/scLongreads/Seurat/scSaturation_by_cellcluster_sr.pdf", width =25, height = 25)
ggplot(df, aes(x=reads, y=nFeature_RNA.y, colour=orig.ident.x)) + 
  geom_line(linewidth=0.8) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  ylab("Number of genes detected") +
  xlab("Number of reads sampled (k)") +
  theme_bw() +
  theme(text = element_text(size=10), legend.title=element_blank()) +facet_wrap(~seurat_clusters, nrow = 8)
dev.off()
df_sub <- subset(df, seurat_clusters %in% c(1,2,7,15,18,21))
ggplot(df_sub, aes(x=reads, y=nFeature_RNA.x, colour=orig.ident.x)) + 
  geom_line(linewidth=0.8) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(0,150,10)) +
  ylab("Number of genes detected") +
  xlab("Number of reads sampled (k)") +
  theme_bw() +
  theme(text = element_text(size=10), legend.title=element_blank()) +facet_wrap(~seurat_clusters, nrow = 1)

VlnPlot(lr.seur, features = c("percent.rb", "percent.mt"),pt.size = 0) 

tmp <- subset(combined.meta, nFeature_RNA.x > 200 & nFeature_RNA.x < 6000)
tmp <- subset(tmp, percent.mt.x < 25)
filtered.meta <- tmp[-which(is.na(tmp$demux.doublet.call)),]
ncol(sr.seur)
nrow(filtered.meta)

lr.seur.sub <- subset(lr.seur, nFeature_RNA > 200 & nFeature_RNA < 6000)
lr.seur.sub <- subset(lr.seur.sub, percent.mt < 25)
gene_w_zero <- rowSums(lr.seur.sub@assays$RNA@counts == 0)/ncol(lr.seur.sub)
gene_mean_count <- rowSums(lr.seur.sub@assays$RNA@counts)/ncol(lr.seur.sub)
length(which(gene_w_zero > 0.9 ))
nrow(lr.seur.sub)

gene_w_zero <- rowSums(sr.seur.annot@assays$RNA@counts == 0)/ncol(sr.seur.annot)
gene_mean_count <- rowSums(sr.seur.annot@assays$RNA@counts)/ncol(sr.seur.annot)
length(which(gene_w_zero > 0.9 | gene_mean_count < 0.1))

sr.meta <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/annot.meta.rds")
saveRDS(lr.seur, file="/data/Choi_lung/scLongreads/Seurat/lr.sc.data.list.22samples.combined.rds")


#
lr.seur <- readRDS("/data/Choi_lung/scLongreads/Seurat/lr.sc.data.list.22samples.combined.rds")
combine.meta <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/combined.meta.rds")

vireo.rst <- read.table("/data/Choi_lung/scLongreads/CS035134/Sample_2_NCI_9_16/v39/vireo_bcftools_het/donor_ids.tsv", sep = "\t", header = TRUE)
colnames(vireo.rst) <- c("orig.barcode", colnames(vireo.rst)[2:8])
vireo.rst$orig.barcode <- paste0(vireo.rst$cell, "-1")
combine.meta.sub <- subset(combine.meta, orig.ident.x == "2_NCI_9_16")
annot.meta <- left_join(lr.seur@meta.data, sr.meta ,by = "barcode")
annot.meta.sub <- subset(annot.meta, orig.ident.x == "2_NCI_9_16")

NCI_9_16.compare <- left_join(vireo.rst, combine.meta.sub, "orig.barcode")
NCI_9_16.compare <- left_join(annot.meta.sub, vireo.rst, "orig.barcode")
NCI_9_16.compare <- NCI_9_16.compare[-which(is.na(NCI_9_16.compare$seurat_clusters)),]
table(NCI_9_16.compare$donor_id, NCI_9_16.compare$Sample)
df <- as.data.frame(table(combine.meta$orig.ident.x))

DefaultAssay(sr.seur.annot) <- "RNA"
sr.seur.annot <- NormalizeData(sr.seur.annot, verbose = FALSE)

barcode.overlapped <- annot.meta$barcode[-which(is.na(annot.meta$Sample))]
barcode.overlapped <- intersect(lr.seur$barcode, sr.seur.annot$barcode)

sr.seur.annot <- RenameCells(sr.seur.annot, new.names = sr.seur.annot$barcode)
expr.cor <- lapply(barcode.overlapped, function(x)cor(lr.seur@assays$RNA@counts[gene.overlapped,x], 
                                                      sr.seur.annot@assays$RNA@counts[gene.overlapped,x])
)

colnames(sr.seur.annot) <- sr.seur.annot$barcode
sr.seur.annot@assays$RNA@data[10:15,1:5]
plot(lr.seur@assays$RNA@data["PTPRC",barcode.overlapped], 
     sr.seur.annot@assays$RNA@data["PTPRC",barcode.overlapped])
cor(lr.seur@assays$RNA@data["PTPRC",barcode.overlapped], 
    sr.seur.annot@assays$RNA@data["PTPRC",barcode.overlapped])
annot.meta$Sample_individual <- paste(annot.meta$orig.ident.x, annot.meta$Sample, sep = "_")
annot.meta <- annot.meta[-which(is.na(annot.meta$Sample)),]
meta.sum <- annot.meta %>% group_by(Sample_individual) %>% summarise(Gene.median = median(nFeature_RNA.x), 
                                              Gene.min = min(nFeature_RNA.x),
                                              Gene.max = max(nFeature_RNA.x),
                                              Iso.median = median(nFeature_Isoforms),
                                              Iso.min = min(nFeature_Isoforms),
                                              Iso.max = max(nFeature_Isoforms))
df <- as.data.frame(table(annot.meta$Sample_individual))
names(dataset_list) <- name
View(dataset_list$`1_NCI_1257`)
tmp <- dataset_list$`1_NCI_1257`
tmp <- tmp[which(tmp$associated_gene == "SFTPA2"),]
tmp2 <- dataset_list$`3_NCI_23_29`
tmp2 <- tmp2[which(tmp2$associated_gene == "SFTPA2"),]

# Isoform detection
Lung_transcripts <- read.table("/data/lib14/project/scLongread/Lung_specific_trancripts_val_annot.txt", sep = "\t", header = TRUE)
ASTS_genes <- read.table("/data/lib14/project/scLongread/ASTS_genes_lung.txt", sep = "\t", header = TRUE)

file_dir <- "/data/Choi_lung/scLongreads/B2_percentile/"
sample_dir <- list.dirs(file_dir, recursive = FALSE)
name <- list.dirs(file_dir, recursive = FALSE, full.names = FALSE)
name <- str_split_fixed(name, pattern = "_", n = 2)[,2]
name <- name[-23]
dataset_list <- list()
for (i in 1:22) {
  filepath <- paste(file_dir, paste("/Sample_",name[i],"/", name[i], "_classification.filtered_lite_classification.txt", sep = ""), sep = "")
  info <- read.table(filepath, sep = "\t", header = TRUE)
  dataset_list[[i]] <- info
}
ASTS_genes.list <- list()
for (i in 1:22) {
  info <- dataset_list[[i]]
  ASTS_genes.list[[i]] <- info[which(info$associated_gene %in% ASTS_genes$symbol),]
}
genes <- list()
for (i in 1:22) {
  genes[[i]] <- unique(ASTS_genes.list[[i]]$associated_gene)
}
info.list <- list()
for (i in 1:22) {
  info <- dataset_list[[i]]
  info$transcript_id <- str_split_fixed(info$associated_transcript, "_", n = 2)[,1]
  info <- subset(info, transcript_id != "novel")
  info.list[[i]] <- info
}
transcript_ids <- str_split_fixed(Lung_transcripts$transcript_id, "\\.", n = 2)[,1]
for (i in 1:22) {
  info <- info.list[[i]]
  info <- info[which(info$transcript_id %in% transcript_ids),]
  info.list[[i]] <- info
}
names(info.list) <- name
for (i in 1:22) {
  info <- info.list[[i]]
  info$Sample <- name[i]
  info.list[[i]] <- info
}

info.long <- Reduce(bind_rows, info.list)
idx.overlapped <- c()
for (i in 1:233) {
  idx <- which(info.long$transcript_id == transcript_ids[i] & info.long$length == Lung_transcripts$transcript_length[i])
  idx.overlapped <- c(idx.overlapped, idx)
}
tmp <- info.long[idx.overlapped,]
table(tmp$Sample)
tmp1 <-tmp %>% group_by(Sample) %>% count(transcript_id)
table(tmp1$Sample)
