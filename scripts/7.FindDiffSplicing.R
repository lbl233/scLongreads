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
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final.RDS")
load("/data/Choi_lung/scLongreads/Seurat/final/lr.sc.data.afterQC.norm.RData")
# Clean up cell type annotation
lr$Celltype <- as.character(lr$Annotation_recluster)
lr$Celltype[which(lr$Celltype == "Plasma cell/Plasmacytoid DCs")] <- "Plasmacytoid DCs"
lr$Celltype[which(lr$Celltype == "Adventitial/peribronchial fibroblasts")] <- "Adventitial fibroblasts"
lr$Celltype[which(lr$res_used == 30)] <- "Plasma cells"
lr$Celltype[which(lr$Celltype == "Monocyte-derived Mφ")] <- "Interstitial macrophages"
lr$Celltype[which(lr$Celltype == "Alveolar macrophages CCL3+")] <- "Alveolar macrophages CCL3"
lr$Celltype[which(lr$Celltype == "Alveolar macrophages MT+")] <- "Alveolar macrophages MT"
lr$Celltype[which(lr$Celltype == "Non-classical monocytes")] <- "Non_classical monocytes"
lr$Celltype[which(lr$Celltype == "Alveolar Mφ proliferating")] <- "Alveolar Mph proliferating"
table(lr$Celltype)
lr$Celltype <- factor(lr$Celltype, levels = c("AT1", "AT2", "Alveolar transitional cells",
                                              "Club", "Goblet", "Secretory transitional cells",
                                              "Basal", "Multiciliated", "Neuroendocrine", # Epithelial category
                                              "Alveolar macrophages", "Alveolar macrophages CCL3",
                                              "Alveolar macrophages MT", "Alveolar Mph proliferating", 
                                              "Interstitial macrophages", "Classical monocyte", "Non_classical monocytes",
                                              "DC2", "CD4 T cells", "CD8 T cells", "NK T cells", "NK cells","T cell proliferating",
                                              "Mast cells", "B cells", "Plasma cells", "Plasmacytoid DCs", # Immune cells
                                              "EC arterial", "EC venous pulmonary",
                                              "EC venous systemic", "Lymphatic EC", 
                                              "EC aerocyte capillary", "EC general capillary",
                                              "Adventitial fibroblasts",
                                              "Alveolar fibroblasts", "Myofibroblast", 
                                              "SMC", "Mesothelium"))
table(lr$Celltype)
# Epithelial signatures
Epi_sig_list <- c("EPCAM", # Epithelial
                  "HOPX", "SFTA2", #Alveolar
                  "AGER", "RTKN2", "SLC39A8", # AT1
                  "SFTPC", "SFTPD", #AT2
                  "SCGB1A1", "SCGB3A1", #Secretory 
                  "BPIFB1", # club
                  "MUC5AC", "MUC5B", # goblet
                  "SCGB3A2", # transitional cell
                  "KRT5", "KRT17", "S100A2", # Basal
                  "TSPAN1", "FOXJ1", "FAM183A", # Multiciliated
                  "CPE", "CHGA" #Neuroendocrine
)

Imm_sig_list <- c("PTPRC",
                  "FCER1G", "CLEC7A", # Myeloid
                  "CD68", "MARCO", "FABP4",# Macrophage
                  "CCL20", "CCL3", # CCL3+ AMs
                  "MT2A", "MT1E", # MT+ AMs
                  "SPP1", "HAMP", # Monocyte derived Mph (aka interstitial Mph)
                  "FCN1", "S100A12", "CD14","MTSS1", "FCGR3A", # classical and non-classical monocyte
                  "CLEC9A", # DC1
                  "CLEC10A", "CD1E", # DC2
                  "LAD1", "CCL19", #migratory DCs
                  "CCL5", #lymphoid
                  "CD3D", "CD3G", # T cells
                  "CD4", "CD28", "TNFRSF25", # CD4 T cells
                  "CD8A", "CD8B", #CD8
                  "GNLY", "NKG7", #NK 
                  "SLC18A2", "CMA1", #Mast cell
                  'MS4A1','CD79A', # B cell
                  'SMPD3', 'SCT','TNFRSF17', # plasma/ plasma DC
                  "TOP2A", "MKI67" # proliferating
)

# Endothelial signatures
Endo_sig_list <- c("CLDN5", # Endothelial
                   "DKK2", "IGFBP3", #Arterial
                   "ACKR1", "VCAM1", # Venous
                   "CPE", "C7", #Pulmonary
                   "ZNF385D", "OLFM1", #Systemic 
                   "LYVE1", "CCL21", "TFF3", # Lymphatic
                   "CA4", "EDNRB", # Aerocyte
                   "IL7R", "FCN3" # general capilary
)

# Stroma signatures
Stroma_sig_list <- c("COL1A2", "DCN", # Stroma
                     "SCARA5", "PI16", #Adventitial Fibroblasts
                     "COL15A1", "CXCL14",
                     "SPINT2", "LIMCH1", "FGFR4", # Alveolar Fibroblasts
                     "TYRP1", "ITGBL1", #Myofibroblast
                     "MYH11", "ACTA2", #SMC
                     "LAMC3", "PDGFRB",
                     "KLK11", "ITLN1" # Mesothelium
)
lr.seur <- subset(lr.seur, cells = colnames(lr))
lr[["RNA"]] <- CreateAssayObject(counts = as.matrix(lr.seur@assays$RNA@counts))
DefaultAssay(lr) <- "RNA"
lr <- NormalizeData(lr)
p <- DotPlot(lr, features = Reduce(union,c(Epi_sig_list, Imm_sig_list, Endo_sig_list, Stroma_sig_list)), 
             cols = c('lightgrey', 'blue'), scale= TRUE, scale.by = 'radius',
             group.by = "Celltype") + RotatedAxis()
pdf(file = "~/github/scLongread/scLongread/plot/Celltype_marker_expr.pdf", width = 30, height = 10)
p + scale_color_distiller(palette="RdYlBu",type = "div")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major= element_line(color="grey", linetype = "dashed"))+
  theme(axis.title.x=element_blank(),  # X axis title
        axis.title.y=element_blank(),  # Y axis title
        # axis.text.x=element_text(size=10, 
        #                          angle = 45,
        #                          vjust=.6),  # X axis text
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  scale_size_continuous(range = c(1,8))+
  geom_vline(xintercept=c(22.5,63.5,78.5),
             linetype=2,color="darkgreen")
dev.off()
lr$predicted.ann_finest_level
DimPlot(lr, reduction = "scarches",group.by = "Celltype", label = TRUE) + NoLegend()
library(RColorBrewer)
n <- length(unique(lr$Celltype))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
pdf("../plot/UMAP_wCelltypeAnnotation_final.pdf", width = 12)
DimPlot(lr, reduction = "scarches",group.by = "Celltype", label = TRUE,cols = col_vector) + coord_equal() 
DimPlot(lr, reduction = "scarches",group.by = "predicted.ann_finest_level", label = TRUE) + coord_equal() 
dev.off()
DimPlot(lr, reduction = "scarches",group.by = "res_used", label = TRUE) + NoLegend()
saveRDS(lr, "/data/Choi_lung/scLongreads/Seurat/final/lr_final_v2_wRNAexpr.RDS")

# summarize cell proportion across individuals
df_stat <- as.data.frame(table(lr$Celltype, lr$Sample_NCI))
cellnumber_perSample <- as.data.frame(table(lr$Sample_NCI))
colnames(df_stat) <- c("Celltype", "Sample", "Cell number")
df_stat$Cell_num_sum_per_sample <- rep(cellnumber_perSample$Freq, each = 37)

df_stat$`Cell proportion` <- (df_stat$`Cell number`/df_stat$Cell_num_sum_per_sample)*100
sum(df_stat$`Cell proportion`[1:37])
df_stat$Category <- rep(c(rep("Epi", 9), rep("Immune", 17), rep("Endo", 6), rep("Stroma", 5)), 129)
library(ggrepel)
library(ggpubr)
ggboxplot(df_stat, x = "Celltype", y = "Cell proportion", color = "Celltype")+
  facet_wrap(~Category, ncol = 2,scales = "free")+
  RotatedAxis()+NoLegend()
library(tidyverse)

library(gghalves)
categorycolors <- c("coral1","lightslateblue","goldenrod1","lightgray")
pdf("../plot/Cell_prop_perInd.pdf", width = 10)
ggplot(data = df_stat,aes(x=Celltype, y=`Cell proportion`, fill=Category)) +    
  #geom_half_violin(side = "r", color=NA, alpha=0.5,width=0.5) +    
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.5, linewidth=0.5) +    
  geom_half_point_panel(side = "l", shape=21, size=1, color="white") +  
  scale_fill_manual(values = categorycolors) +
  labs(y="Cell Proportion (%)",x=NULL) +   
  theme_classic() + 
  #facet_wrap(~Category, ncol = 2,scales = "free")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "bottom",axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black")) + NoLegend() 
dev.off()
df_stat_sub <- subset(df_stat, `Cell number` >= 5)
table(df_stat_sub$Celltype)
df <- as.data.frame(table(df_stat_sub$Celltype))
ggplot(df, aes(x=Var1, y=Freq,fill=Var1)) +
  geom_bar(stat="identity")+
  scale_fill_manual("Legend", values = col_vector)+
  labs(y="Number of individuals (> 5 cells)",x= "Cell type") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  NoLegend()







# load the original isoform level expression profiles
lr.isoform <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/final/lr.sc.data.Isoform.22samples.combined.rds")
ncol(lr)
lr.isoform.sub <- subset(lr.isoform, cell = colnames(lr))
lr.isoform.sub$Sample <- lr$Sample
lr.isoform.sub$Sample_NCI <- lr$Sample_NCI
lr.isoform.sub$Celltype <- lr$Celltype
count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass$transcript_name_unique,]
count_mtx@Dimnames[[1]] <- TALON_afterqc_orf_secondpass$annot_transcript_id
dim(count_mtx)
sum(count_mtx)

# lr.isoform.sub$Celltype <- paste(lr.isoform.sub$Celltype, 
#                                         lr.isoform.sub$Sample_NCI, sep = "_")
group <- lr.isoform.sub$Celltype %>% fct_inorder()
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
# write separated expression profiles of each cell types at isoform level
for (celltype in celltype_list) {
  idx <- grep(celltype, colnames(tpm_mtx))
  mtx <- tpm_mtx[,idx]
  write.table(mtx,file=paste0("/data/Choi_lung/scLongreads/TALON_workspace/test10s/iso_tpm_formatted_",gsub(" ", "_", celltype),".txt"), sep = "\t", quote = FALSE)
}
