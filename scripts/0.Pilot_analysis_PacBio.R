# Pilot anaylsis of Sample1 by PacBio 11/14/2013

library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)

library(reticulate)
scr <- import("scrublet", convert = FALSE)
py_run_string('import scrublet as scr')

file_dir <- "/data/Choi_lung/scLongreads/B2_percentile/"
sample_dir <- list.dirs(file_dir, recursive = FALSE)
name <- list.dirs(file_dir, recursive = FALSE, full.names = FALSE)
name <- str_split_fixed(name, pattern = "_", n = 2)[,2]
name <- name[-22]
dataset_list <- list()

for (i in 1:length(name)) {
  filepath = paste(file_dir, paste("/Sample_",name[i],"/", name[i], "_classification_filtered", sep = ""), "/genes_seurat", sep = "")
  filepath_isoforms = paste(file_dir, paste("/Sample_",name[i],"/", name[i], "_classification_filtered", sep = ""), "/isoforms_seurat", sep = "")
  # doubletpath = paste(doublet_dir , name[i], sep = "")
  # Seurat
  seur.data = Read10X(data.dir = filepath)
  seur.isoforms.data = Read10X(data.dir = filepath_isoforms)
  dim(seur.data)
  dim(seur.isoforms.data)
  seur = CreateSeuratObject(counts = seur.data, project = name[i], min.cells=5)
  print(ncol(seur))
  head(colnames(seur))
  seur[["Isoforms"]] <- CreateAssayObject(counts = seur.isoforms.data, min.cells=5)
  # save original barcodes from PacBio
  seur$orig.barcode <- colnames(seur)
  # reverse complementary DNA of barcodes
  CB.reverse <- chartr(old="ATGC", new="TACG", colnames(seur))
  CB.reverse <- stri_reverse(CB.reverse)
  CB.lr <- paste(str_split_fixed(CB.reverse, "-", n =2)[,2], "-1", sep = "")
  # Rename the cells using the reversed barcodes
  seur <- RenameCells(seur, new.names = CB.lr)
  # doublet = read.table(paste(doubletpath, '/predicted_doublets.txt', sep = ""))
  # doublet_score = read.table(paste(doubletpath, '/doublet_scores.txt', sep = ""))
  # seur[['doublet_score']] = doublet_score$V1
  # seur[['doublet']] = doublet$V1
  # seur[['Batch']] = Batch[i]
  # detect doublet by Scrublet
  ## Gene level
  scrub <- scr$Scrublet(t(as.matrix(seur@assays$RNA@counts)))
  doublet_scores = scrub$scrub_doublets()
  doublet_scores = py_to_r(doublet_scores)
  seur$doublet_RNA = doublet_scores[[1]]
  seur$doublet2_RNA = doublet_scores[[2]]
  ## Isoform level
  ptm <- proc.time()
  scrub <- scr$Scrublet(t(as.matrix(seur@assays$Isoforms@counts))) # 40GB
  doublet_scores = scrub$scrub_doublets()
  doublet_scores = py_to_r(doublet_scores)
  seur$doublet_isoform = doublet_scores[[1]]
  seur$doublet2_isoform = doublet_scores[[2]]
  print(proc.time() - ptm)
  # seur <- subset(seur, doublet == 0 & nFeature_RNA < 5000)
  seur[["percent.mt"]] = PercentageFeatureSet(seur, pattern = "^MT-")
  seur[['percent.rb']] = PercentageFeatureSet(seur, pattern = '^RP[SL]')
  dataset_list[[i]] <- seur
}
names(dataset_list) <- name
rm("seur.data", "seur.isoforms.data")
save(list="dataset_list", file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.list.21samples.beforeQC.RData")

gc()

print(proc.time()-ptm) # 885.850
table(seur$doublet2_RNA)
table(seur$doublet2)
table(seur$doublet2_RNA, seur$doublet2)
hist(seur$doublet_RNA)
hist(seur$doublet)
VlnPlot(seur, features = c("nCount_RNA", "nFeature_RNA", "doublet_RNA",
                           "nCount_Isoforms", "nFeature_Isoforms", "doublet"),
        ncol = 3, pt.size = 0)
median(seur$nFeature_RNA)
median(seur$nFeature_Isoforms)

seur <- NormalizeData(seur, verbose = FALSE)
seur <- NormalizeData(seur, verbose = FALSE, assay = "Isoforms")
seur <- FindVariableFeatures(seur, verbose = FALSE, nfeatures = 2000)
seur <- ScaleData(seur, features = VariableFeatures(seur), verbose = FALSE)
seur <- RunPCA(seur, features = VariableFeatures(seur), verbose = FALSE)
seur <- FindNeighbors(seur, dims = 1:50)
seur <- RunUMAP(seur, dims = 1:50)
seur <- FindClusters(seur, resolution = .1)
FeaturePlot(seur,features = c("doublet_RNA","doublet"), label = TRUE)
FeaturePlot(seur,features = c("PTPRC","EPCAM", "PECAM1", "SFTPC", 
                              "HOPX", "FOXJ1", "LYVE1", 
                              "CD68", "CD3D", "CD19", 
                              "CD4", "CD79A"),label = TRUE)
DefaultAssay(seur) <- "Isoforms"
FeaturePlot(seur,features = c("PTPRC","PTPRC.4"), label = TRUE)
DotPlot(seur, features = rownames(seur@assays$Isoforms)[grep("EPCAM", rownames(seur@assays$Isoforms))])+RotatedAxis()

# Illumina data
illumina_data_path = "/data/Choi_lung/TTL/Seurat/Final_Batch_Lift"
illumina_file_names <- list.files(illumina_data_path,pattern = ".RDS")
illumina.datalist <- list()

for (i in 1:length(illumina_file_names)) {
  illumina.datalist[[i]] <- readRDS(paste(illumina_data_path,illumina_file_names[i], sep = "/"))
}
save(list="illumina.datalist", file = "/data/Choi_lung/scLongreads/Seurat/sr.sc.data.list.samples.beforeQC.RData")

nci17_22 <- readRDS(paste(illumina_data_path,illumina_file_names[11], sep = "/"))
nci81_86 <- readRDS(paste(illumina_data_path,illumina_file_names[22], sep = "/"))
table(nci17_22$Sample)
summary(nci17_22@meta.data)

save(list = ls(), file = "/data/Choi_lung/scLongreads/pilot_analysis/NCI_17_22.RData")
load("/data/Choi_lung/scLongreads/pilot_analysis/NCI_17_22.RData")
ncol(nci81_86)
ncol(nci17_22)
ncol(seur)

filepath
summary(seur@meta.data)
head(colnames(seur))
head(colnames(nci17_22))
length(intersect(colnames(seur), colnames(nci17_22)))
CB.sr <- colnames(nci17_22)
which(colnames(nci17_22) == "ATCAGGTAGCTGGTGA-1")

CB.lr <- read.table("/data/Choi_lung/scLongreads/CS035134/Sample_NCI_17_22/NCI_17_22_classification_filtered/genes_seurat/barcodes.tsv")

# reverse complementary DNA of barcode
stri_reverse(colnames(seur)[1])
CB.reverse <- chartr(old="ATGC", new="TACG", colnames(seur))
CB.reverse <- stri_reverse(CB.reverse)
CB.lr <- paste(str_split_fixed(CB.reverse, "-", n =2)[,2], "-1", sep = "")
length(intersect(CB.lr, CB.sr))

name <- name[-22]
sample_dir <- sample_dir[-1]
samples_compared <- name[c(2,3,4,12,22)]
samples_compared_dir <- sample_dir[c(2,3,4,12,22)]
names(samples_compared_dir) <- samples_compared
df.list <- list()
for(i in 1:length(samples_compared_dir)){
  df <- read.table(paste0(samples_compared_dir[i],"/v39/saturation.txt"), header = TRUE)
  df$SampleID <- samples_compared[i]
  df$Version <- "v39"
  df.list[[i]] <- df
}

saturation_dir <- list.files("/data/Choi_lung/scLongreads/B2_percentile/saturation/", recursive = FALSE)
df.old.list <- list()

for(i in 1:length(saturation_dir)){
  df.old <- read.table(paste0("/data/Choi_lung/scLongreads/saturation/", saturation_dir[i]), sep = "\t",header = TRUE)
  df.old$SampleID <- name[i]
  df.old$Version <- "v32"
  df.old.list[[i]] <- df.old
}

df_stauration_v39 <- Reduce(rbind, df.list)
df_stauration_v32_1 <- Reduce(rbind, df.old.list)
df_stauration_v32_2 <- Reduce(rbind, df.old.list)
df_stauration_v32_1$Batch <- "B2"
df_stauration_v32_2$Batch <- "B1"
library(reshape2)
df_stauration_v32 <- rbind(df_stauration_v32_1, df_stauration_v32_2)
df_stauration <- melt(df_stauration_v32, id.vars = c("SampleID", "Version", "reads", "Batch"))
df_stauration$reads <- df_stauration$reads/1000000
df_stauration_sub <- subset(df_stauration, SampleID == "10_NCI_69_74")
#linear plot
library(ggplot2)
ggplot(df_stauration, aes(x=reads, y=value, colour=variable)) + 
  geom_line(size=0.8) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(0,60,5)) +
  ylab("Number of genes detected") +
  xlab("Number of reads sampled (millions)") +
  theme_bw() +
  theme(text = element_text(size=10), legend.title=element_blank())+facet_wrap(~SampleID, nrow = 5)

ggplot(df_stauration_sub, aes(x=reads, y=value, colour=variable)) + 
  geom_line(size=0.8) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(0,60,5)) +
  ylab("Number of genes detected") +
  xlab("Number of reads sampled (millions)") +
  theme_bw() +
  theme(text = element_text(size=10), legend.title=element_blank())+facet_wrap(~Batch, nrow = 1)

# Sample 2_NCI_9_16
df <- read.table("/data/Choi_lung/scLongreads/CS035134/Sample_2_NCI_9_16/v39/saturation.txt", header = TRUE)
df$SampleID <- "2_NCI_9_16"
df$Version <- "v39"

df$Batch <- "B1"
df <- melt(df, id.vars = c("SampleID", "Version", "reads", "Batch"))
df$reads <- df$reads/1000000
df_stauration_sub <- subset(df_stauration, SampleID == "2_NCI_9_16")
df_compare <- rbind(df, df_stauration_sub)
df_compare <- subset(df_compare, Batch == "B1")
df_compare$Type <- paste(df_compare$variable, df_compare$Version, sep = "_")
ggplot(df_compare, aes(x=reads, y=value, colour=Type)) + 
  geom_line(size=0.8) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(0,60,5)) +
  ylab("Number of genes detected") +
  xlab("Number of reads sampled (millions)") +
  theme_bw() +
  theme(text = element_text(size=10), legend.title=element_blank())

