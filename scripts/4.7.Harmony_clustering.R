
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(reticulate)
use_condaenv("/data/lib14/conda/envs/scvi-env", required = TRUE)
library(leiden)
library(SeuratWrappers)
library(sceasy)
library(patchwork)


load("/data/Choi_lung/scLongreads/Seurat/lr.sc.data.beforharmony.RData")
dim(lr.seur)
lr.seur@assays$RNA
lr.seur@assays$SCT
DefaultAssay(lr.seur)
lr.seur <- RunPCA(lr.seur, npcs = 30, verbose = F)
lr.seur <- IntegrateLayers(
  object = lr.seur, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = FALSE
)
lr.seur <- FindNeighbors(lr.seur, dims = 1:30, reduction = "harmony")
lr.seur <- FindClusters(lr.seur, resolution = 1.5, graph.name = "SCT_snn", algorithm = 4, method = "igraph")
lr.seur <- FindClusters(lr.seur, resolution = 2, graph.name = "SCT_snn", algorithm = 4, method = "igraph")
lr.seur <- FindClusters(lr.seur, resolution = 2.5, graph.name = "SCT_snn", algorithm = 4, method = "igraph")
lr.seur <- RunUMAP(lr.seur, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(lr.seur, group.by = "predicted.ann_finest_level", split.by = "orig.ident", ncol = 6)+NoLegend()
DimPlot(lr.seur, group.by = "predicted.ann_finest_level", label = TRUE)

saveRDS(lr.seur, file = "/data/Choi_lung/scLongreads/Seurat/final/lr.sc.data.harmony.RDS")

# dimensional reduction and batch effect removal by harmony
lr.seur <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/final/lr.sc.data.harmony.RDS")
lr.final <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/final/lr_final.RDS")
lr.seur$category <- "drop"
lr.seur$category[which(lr.seur$SCT_snn_res.2 %in% c(2,4,6,7,11,12,14,16,17,19,21,25,27,38))] <- "Epithelial"
lr.seur$category[which(lr.seur$SCT_snn_res.2 %in% c(1,5,9,10,31,39))] <- "Endothelial"
lr.seur$category[which(lr.seur$SCT_snn_res.2 %in% c(3,8,13,15,18,20,22,23,24,26,37,40,42,43))] <- "Immune"
lr.seur$category[which(lr.seur$SCT_snn_res.2 %in% c(28,29,33,35))] <- "Stroma"

# resolution 2.5
DotPlot(lr.seur, features = c("EPCAM", "PTPRC", "CLDN5", "COL1A2", "MKI67"))
lr.seur$category <- "drop"
lr.seur$category[which(lr.seur$SCT_snn_res.2.5 %in% c(1,4,6,8,9,15,16,19,20,21,24,25,27,30,31,32,33,40,53))] <- "Epithelial"
lr.seur$category[which(lr.seur$SCT_snn_res.2.5 %in% c(3,7,11,14,18,29,38,41,50))] <- "Endothelial"
lr.seur$category[which(lr.seur$SCT_snn_res.2.5 %in% c(2,5,10,12,13,17,22,23,26,28,34,37,44,49,51,54))] <- "Immune"
lr.seur$category[which(lr.seur$SCT_snn_res.2.5 %in% c(35,36,42,46))] <- "Stroma"
DimPlot(lr.seur, split.by = "category", label = TRUE,ncol = 3)
table(lr.seur$category)
DotPlot(lr.seur, features = c("EPCAM", "PTPRC", "CLDN5", "COL1A2", "MKI67"), group.by = "category")
# dimensional reduction and batch effect removal by scANVI
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/Leiden_clustree.RDS")
lr$category <- "drop"
lr$category[which(lr$originalexp_snn_res.1.5 %in% c(2,3,4,5,11,12,14,18,21,36,42,45))] <- "Epithelial"
lr$category[which(lr$originalexp_snn_res.1.5 %in% c(7,9,10,13,28))] <- "Endothelial"
lr$category[which(lr$originalexp_snn_res.1.5 %in% c(1,6,8,15,16,19,23,30,34,35,37,38,41,44))] <- "Immune"
lr$category[which(lr$originalexp_snn_res.1.5 %in% c(17,20,26))] <- "Stroma"
table(lr$category)
table(lr$originalexp_snn_res.1.5)
View(as.data.frame.matrix(table(lr.seur$SCT_snn_res.2.5, lr.seur$predicted.ann_finest_level)))
which(colnames(lr) != colnames(lr.seur))
View(table(lr$originalexp_snn_res.1.5, lr.seur$SCT_snn_res.2.5))
table(lr$category, lr.seur$category)
lr.final <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final.RDS")
lr_sub <- subset(lr, cells = colnames(lr.final))
lr_sub.invert <- subset(lr, cells = colnames(lr.final), invert = TRUE)
table(lr_sub.invert$category)
lr.seur_sub <- subset(lr.seur, cells = colnames(lr.final))
lr.seur_sub$Annotate_recluster <- lr.final$Annotation_recluster
DimPlot(lr.seur_sub, group.by = "category")
DimPlot(lr.seur_sub.invert, group.by = "category")
DimPlot(lr.seur, label = TRUE)
DimPlot(lr.seur_sub, group.by = "Annotate_recluster", label = TRUE)+NoLegend()
DimPlot(lr.seur_sub, group.by = "Annotate_recluster", split.by = "category", ncol = 3)
DimPlot(lr.seur_sub, split.by = "category", label = TRUE,ncol = 3)

table(lr.seur_sub$category,lr.seur_sub$Annotate_recluster) 
lr.seur_sub.invert <- subset(lr.seur, cells = colnames(lr.final), invert = TRUE)
table(lr.seur_sub.invert$category)
write.table(table(lr$originalexp_snn_res.1.5, lr.seur$SCT_snn_res.2), 
            file = "/data/Choi_lung/scLongreads/Seurat/CellComp_compare_harmony_scANVI.txt", sep = "\t", quote = FALSE)
df_per <- as.data.frame.matrix((table(lr$originalexp_snn_res.1.5, lr.seur$SCT_snn_res.2)))
df_new <- apply(df_per, 1, function(x)x/sum(x))
df_new <- t(df_new)
library(reshape2)
df_new <- melt(df_new)
df_new$scANVI <- df_new$Var1
df_new$category <- "drop"
df_new$category[which(df_new$scANVI %in% c(2,3,4,5,11,12,14,18,21,36,42,45))] <- "Epithelial"
df_new$category[which(df_new$scANVI %in% c(7,9,10,13,28))] <- "Endothelial"
df_new$category[which(df_new$scANVI %in% c(1,6,8,15,16,19,23,30,34,35,37,38,41,44))] <- "Immune"
df_new$category[which(df_new$scANVI %in% c(17,20,26))] <- "Stroma"
df_new$scANVI <- factor(df_new$scANVI, levels = c(2,3,4,5,11,12,14,18,21,36,42,45, 
                                                  7,9,10,13,28,
                                                  1,6,8,15,16,19,23,30,34,35,37,38,41,44,
                                                  17,20,26,
                                                  22,24,25,27,29,31,32,33,39,40,43))

df_new$Harmony <- df_new$Var2
df_new$category2 <- "drop"
df_new$category2[which(df_new$Harmony %in% c(2,4,6,7,11,12,14,16,17,19,21,25,27,38))] <- "Epithelial"
df_new$category2[which(df_new$Harmony %in% c(1,5,9,10,31,39))] <- "Endothelial"
df_new$category2[which(df_new$Harmony %in% c(3,8,13,15,18,20,22,23,24,26,37,40,42,43))] <- "Immune"
df_new$category2[which(df_new$Harmony %in% c(28,29,33,35))] <- "Stroma"
df_new$Percentage <- df_new$value
df_new$Harmony <- factor(df_new$Harmony, levels = c(2,4,6,7,11,12,14,16,17,19,21,25,27,38,
                                                    1,5,9,10,31,39,
                                                    3,8,13,15,18,20,22,23,24,26,40,42,43,
                                                    28,29,33,35,
                                                    30,32,34,36,37,41))
df_new$category2 <- factor(df_new$category2, levels = c("Epithelial", "Endothelial", "Immune", "Stroma", "drop"))
df_new$category <- factor(df_new$category, levels = c("Epithelial", "Endothelial", "Immune", "Stroma", "drop"))
library(ggsci)
library(scales)
ggplot(df_new, aes(scANVI, Harmony, fill= Percentage)) + 
  geom_tile()+ theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))+
  theme_classic() + 
  facet_grid2(.~category, scales = "free_x", space = "free_x", switch = "x",
                                  strip = strip_themed(
                                    background_x = elem_list_rect(
                                      fill = pal_nejm(alpha = .5)(5))))

              