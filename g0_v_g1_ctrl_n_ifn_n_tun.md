## g0_v_g1_ctrl_n_ifn_n_ifn_tun

```{r}
library(Seurat)
library(cowplot)
library(stringr)
library(xlsx)

## Setup Seruat g0_ctrl, g1_ctrl, g0_ifn & g1_ifn
g0.ctrl <- Read10X(data.dir = "g0_ctrl")
g0_ctrl <- CreateSeuratObject(counts = g0.ctrl, project = "g0_ctrl", min.cells = 5)
g0_ctrl$stim <- "G0_CTRL"
g0_ctrl[["percent.mt"]] <- PercentageFeatureSet(g0_ctrl, pattern = "^MT-")
plot1 <- FeatureScatter(g0_ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(g0_ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
g0_ctrl <- subset(g0_ctrl, subset=nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
g0_ctrl <- NormalizeData(object = g0_ctrl, verbose = FALSE)
g0_ctrl <- FindVariableFeatures(object = g0_ctrl, selection.method = "vst", nfeatures = 2000)

g1.ctrl <- Read10X(data.dir = "g1_ctrl")
g1_ctrl <- CreateSeuratObject(counts = g1.ctrl, project = "g1_ctrl", min.cells = 5)
g1_ctrl$stim <- "G1_CTRL"
g1_ctrl[["percent.mt"]] <- PercentageFeatureSet(g1_ctrl, pattern = "^MT-")
plot1 <- FeatureScatter(g1_ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(g1_ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
g1_ctrl <- subset(x = g1_ctrl, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
g1_ctrl <- NormalizeData(object = g1_ctrl, verbose = FALSE)
g1_ctrl <- FindVariableFeatures(object = g1_ctrl, selection.method = "vst", nfeatrues = 2000)

g0.ifn <- Read10X(data.dir = "g0_ifn")
g0_ifn <- CreateSeuratObject(counts = g0.ifn, project = "g0_ifn", min.cells = 5)
g0_ifn$stim <- "G0_IFN"
g0_ifn[["percent.mt"]] <- PercentageFeatureSet(g0_ifn, pattern = "^MT-")
plot1 <- FeatureScatter(g0_ifn, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(g0_ifn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
g0_ifn <- subset(x = g0_ifn, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 6.5)
g0_ifn <- NormalizeData(object = g0_ifn, verbose = FALSE)
g0_ifn <- FindVariableFeatures(object = g0_ifn, selection.method = "vst", nfeatures = 2000)

g1.ifn <- Read10X(data.dir = "g1_ifn")
g1_ifn <- CreateSeuratObject(counts = g1.ifn, project = "g1_ifn", min.cells = 5)
g1_ifn$stim <- "G1_IFN"
g1_ifn[["percent.mt"]] <- PercentageFeatureSet(g1_ifn, pattern = "^MT-")
plot1 <- FeatureScatter(g1_ifn, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(g1_ifn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
g1_ifn <- subset(x = g1_ifn, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 8)
g1_ifn <- NormalizeData(object = g1_ifn, verbose = FALSE)
g1_ifn <- FindVariableFeatures(object = g1_ifn, selection.method = "vst", nfeatrues = 2000)

g0.ifn_tun <- Read10X(data.dir = "g0_ifn_tun")
g0_ifn_tun <- CreateSeuratObject(counts = g0.ifn_tun, project = "g0_ifn_tun", min.cells = 5)
g0_ifn_tun$stim <- "G0_IFN_TUN"
g0_ifn_tun[["percent.mt"]] <- PercentageFeatureSet(g0_ifn_tun, pattern = "^MT-")
plot1 <- FeatureScatter(g0_ifn_tun, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(g0_ifn_tun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
g0_ifn_tun <- subset(x = g0_ifn_tun, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 6)
g0_ifn_tun <- NormalizeData(object = g0_ifn_tun, verbose = FALSE)
g0_ifn_tun <- FindVariableFeatures(object = g0_ifn_tun, selection.method = "vst", nfeatures = 2000)

g1.ifn_tun <- Read10X(data.dir = "g1_ifn_tun")
g1_ifn_tun <- CreateSeuratObject(counts = g1.ifn_tun, project = "g1_ifn_tun", min.cells = 5)
g1_ifn_tun$stim <- "G1_IFN_TUN"
g1_ifn_tun[["percent.mt"]] <- PercentageFeatureSet(g1_ifn_tun, pattern = "^MT-")
plot1 <- FeatureScatter(g1_ifn_tun, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(g1_ifn_tun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
g1_ifn_tun <- subset(x = g1_ifn_tun, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 8)
g1_ifn_tun <- NormalizeData(object = g1_ifn_tun, verbose = FALSE)
g1_ifn_tun <- FindVariableFeatures(object = g1_ifn_tun, selection.method = "vst", nfeatrues = 2000)

## Integrate Conditions
all.anchors <- FindIntegrationAnchors(object.list = list(g0_ctrl, g1_ctrl, g0_ifn, g1_ifn, g0_ifn_tun, g1_ifn_tun), dims = 1:20)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:20)

## Integrate Analysis
DefaultAssay(object = all.combined) <- "integrated"

all.combined <- ScaleData(object = all.combined, verbose = FALSE)
all.combined <- RunPCA(object = all.combined, npcs = 30, verbose = FALSE)

## UMAP and Clustering
all.combined <- RunUMAP(object = all.combined, reduction = "pca", dims = 1:20)
all.combined <- FindNeighbors(object = all.combined, reduction = "pca", dims = 1:20)
all.combined <- FindClusters(all.combined, resolution = 0.2)

## Visualization
p1 <- DimPlot(object = all.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = all.combined, reduction = "umap", label = TRUE)
umap <- plot_grid(p1, p2)
ggsave("umap_combined_plot_res0.2.png", umap, width=53, height=29, units="cm")

umap_c <- DimPlot(object = all.combined, reducton = "umap", split.by = "stim")
ggsave("umap_compare_plot_res0.2.png", umap_c, width=53, height=29, units="cm")

## Identify Sample Specific
samples <- all.combined@graphs$integrated_snn@Dimnames[[1]]
samples <- str_extract(samples, "\\-(.*)")
samples %>% str_replace_all(c("-1_1" = "g0_ctrl_Jan", "-2_1" = "g0_ctrl_Nov", "-3_1" = "g0_dmso_Dec", "-4_1" = "g0_dmso_Nov",
                              "-1_2" = "g1_ctrlA", "-2_2" = "g1_ctrlB", "-3_2" = "g1_ctrlC", "-4_2" = "g1_ctrl_Nov",
                              "-1_3" = "g0_ifn_Dec", "-2_3" = "g0_ifn_Jan", "-1_4" = "g1_ifn_Dec", "-2_4" = "g1_ifn_Jan",
                              "-1_5" = "g0_ifn_tun_Dec", "-2_5" = "g0_ifn_tun_Jan", "-1_6" = "g1_ifn_tun_Dec", "-2_6" = "g1_ifn_tun_Jan")) -> samples
all.combined@meta.data$orig.ident <- samples
p1 <- DimPlot(object = all.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = all.combined, reduction = "umap", label = TRUE)
umap <- plot_grid(p1, p2)
ggsave("umap_combined_sample_plot_res0.2.png", umap, width=53, height=29, units="cm")

umap_s <- DimPlot(object = combined, reduction = "umap", split.by = "orig.ident")
ggsave("umap_split_by_condition.png", umap_s, width=159, height=29, units="cm", limitsize=FALSE)

## Create FindConservedMarkers Function
DefaultAssay(object = all.combined) <- "RNA"
clus = 0
while (clus < length(levels(all.combined@meta.data$seurat_clusters)))
{
    print(paste("Cluster", clus, sep=" "))
    markers <- FindConservedMarkers(object = all.combined, ident.1=clus, grouping.var="stim", verbose=FALSE)
    if (clus == 0) {
        write.xlsx(markers, file="g0_v_g1_ctrl_n_ifn_n_tun_markers_res0.2.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=FALSE)
    }
    else {
        write.xlsx(markers, file="g0_v_g1_ctrl_n_ifn_n_tun_markers_res0.2.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=TRUE)
    }
    clus = clus + 1
}

save.image(file="res0.2.RData")

combined <- RenameIdents(all.combined, '0' = 'Mesenchyme', '1' = 'Mesenchyme', '2' = 'Mesenchyme', '3' = 'Mesenchyme', '4' = 'Cycling Cells',
                         '5' = 'Mesenchyme', '6' = 'Mesenchyme', '7' = 'Tubule', '8' = 'Neuron', '9' = 'Mesenchyme', '10' = 'Artifact',
                         '11' = 'Glomerular Epithelial Cells', '12' = 'Endothelial Cells')

umap <- DimPlot(object = combined, label = FALSE, cols = c('#A6E9A9', '#49C64E', '#0000FF', '#FFEE20', '#D3D3D3',  '#FF0000', '#A600FF'), pt.size=1)
ggsave('umap_cluster_colored_w_label_res0.2.png', umap, width=40, height=40, units='cm')

DefaultAssay(combined) <- "integrated"

###############################################################
## Recluster Proximal and Distal Tubule
cluster7 <- SubsetData(combined, ident.use='Tubule', do.clean=T)
cluster7 <- NormalizeData(cluster7)
cluster7 <- FindVariableFeatures(cluster7, selection.method = 'vst')
cluster7 <- ScaleData(cluster7)
cluster7 <- RunPCA(cluster7, npcs=30)
cluster7 <- RunUMAP(cluster7, reduction='pca', dims = 1:20)
cluster7 <- FindNeighbors(cluster7, reduction='pca', dims=1:20)
cluster7 <- FindClusters(cluster7, resolution=0.1)

## Visualization
p1 <- DimPlot(object = cluster7, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = cluster7, reduction = "umap", label = TRUE)
umap <- plot_grid(p1, p2)
ggsave("umap_proximal_and_distal_tubule_res0.1.png", umap, width=53, height=29, units="cm")

umap_c <- DimPlot(object = cluster7, reducton = "umap", split.by = "stim")
ggsave("umap_sample_compare_proximal_and_distal_tubule_res0.1.png", umap_c, width=53, height=29, units="cm")

p1 <- DimPlot(object = cluster7, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = cluster7, reduction = "umap", label = TRUE)
umap <- plot_grid(p1, p2)
ggsave("umap_combined_sample_plot_proximal_and_distal_tubule_res0.1.png", umap, width=53, height=29, units="cm")

all.clus7.markers <- FindAllMarkers(object = cluster7, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25, grouping.var = "stim")

## Create XSLS File
clus = 0
while (clus < length(levels(cluster7@meta.data$seurat_clusters)))
{
    all.clus7.markers %>% dplyr::filter(cluster==clus & p_val_adj <= 0.05) -> val
    if (clus == 0) {
        rownames(val) <- val$gene
        write.xlsx(val, file="cluster7_markers_res0.1.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=FALSE)
    }
    else {
        rownames(val) <- val$gene
        write.xlsx(val, file="cluster7_markers_res0.1.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=TRUE)
    }
    clus = clus + 1
}

## Create and Use FindConservedMarkers Function
## DefaultAssay(object = cluster7) <- "RNA"
## clus = 0
## while (clus < length(levels(cluster7@meta.data$seurat_clusters)))
## {
##     print(paste("Cluster", clus, sep=" "))
##     markers <- FindConservedMarkers(object = cluster7, ident.1=clus, grouping.var="stim", verbose=FALSE)
##     if (clus == 0) {
##         write.xlsx(markers, file="proximal_and_distal_tubule_markers_res0.1.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=FALSE)
##     }
##     else {
##         write.xlsx(markers, file="proximal_and_distal_tubule_markers_res0.1.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=TRUE)
##     }
##     clus = clus + 1
## }

umap <- DimPlot(object = cluster7, label = FALSE, cols = c('#0000FF', '#b1b1ff', '#808080', '#808080','#808080','#808080','#808080','#808080','#808080'),
                pt.size=1)
ggsave('umap_cluster_7_colored_res0.1.png', umap, width=40, height=40, units='cm')

##########################################################################################
##########################################################################################

## Perform Diff Expressed Analysis (NOTE: Cluster names weren't changed the are NUMBERS)
clus_7 <- subset(x = all.combined, idents = "7")
Idents(object = clus_7) <- "stim"
avg.clus7 <- log1p(x = AverageExpression(object = clus_7, verbose = FALSE)$RNA)
avg.clus7$gene <- rownames(x = avg.clus7)
avg.clus7 %>% dplyr::mutate(diff = G1_IFN - G0_IFN) %>% dplyr::arrange(desc(diff)) -> clus7.list

clus_11 <- subset(x = all.combined, idents = "11")
Idents(object = clus_11) <- "stim"
avg.clus11 <- log1p(x = AverageExpression(object = clus_11, verbose = FALSE)$RNA)
avg.clus11$gene <- rownames(x = avg.clus11)
avg.clus11 %>% dplyr::mutate(diff = G1_IFN - G0_IFN) %>% dplyr::arrange(desc(diff)) -> clus11.list # Visualize the difference

clus_12 <- subset(x = all.combined, idents = "12")
Idents(object = clus_12) <- "stim"
avg.clus12 <- log1p(x = AverageExpression(object = clus_12, verbose = FALSE)$RNA)
avg.clus12$gene <- rownames(x = avg.clus12)
avg.clus12 %>% dplyr::mutate(diff = G1_IFN - G0_IFN) %>% dplyr::arrange(desc(diff)) -> clus12.list

## General Relabel and Set celltype.stim w/ celltype
all.combined$celltype.stim <- paste(Idents(object = all.combined), all.combined$stim, sep = "_")
all.combined$celltype <- Idents(object = all.combined)
Idents(object = all.combined) <- "celltype.stim"

## Subcluster Tubule and 2 GEC Clusters
Idents(object = all.combined) <- "seurat_clusters"
clus7_11_12 <- SubsetData(all.combined, ident.use = c("7", "11", "12"))
Idents(object = clus7_11_12) <- "stim"
avg.clus7_11_12 <- log1p(x = AverageExpression(object = clus7_11_12, verbose = FALSE)$RNA)
avg.clus7_11_12$gene <- rownames(x = avg.clus7_11_12)

## Create Marker List for Each Cluster
all.markers <- FindAllMarkers(object = all.combined)
Idents(object = all.combined) <- "stim"
avg.all <- log1p(x = AverageExpression(object = all.combined, verbose = FALSE)$RNA)
avg.all$gene <- rownames(x = avg.all)
avg.all %>% dplyr::mutate(diff = G0_CTRL - G1_CTRL) %>% dplyr::arrange(desc(diff)) -> all.list

## To View Genes w/ Largest abs(avg_logFC)
p1 <- ggplot(avg.all, aes(G0_CTRL, G1_CTRL)) + geom_point(size=3) + ggtitle("G0_CTRL vs. G1_CTRL") +
    geom_point(data=avg.all[c("RPS4Y1", "DLK1", "LY6E", 'S100A6', 'NPNT', 'COL1A2', 'ITM2B', 'CHCHD2', 'HLA-C', 'LETMD1', 'PLPP3', 'UCHL1', 'DES',
                                 'OLFML3', 'HLA-A', 'LDHA', 'CRELD2', 'MYL9', 'TPM2', 'EPSTI1', 'ASS1', 'LMCD1', 'NDFIP1', 'FOS', 'TMEM98', 'TSPO'),],
               color="red")
p1 <- LabelPoints(plot = p1, points = c('RPS4Y1', 'DLK1', 'LY6E', 'S100A6', 'NPNT', 'COL1A2', 'ITM2B', 'CHCHD2', 'HLA-C', 'LETMD1'))
ggsave("all_clusters_g0_ctrl_v_g1_ctrl_scatter_plot.png", p1, width=20, height=20, units="cm")
write.csv(c12, file="cluster12_g0_v_g1_ifn_tun_res0.5_markers.csv", row.names = c12$genes)

p1 <- ggplot(avg.all, aes(G0_CTRL, G1_CTRL)) + geom_point(size=3) + ggtitle("G0_CTRL vs. G1_CTRL")
ggsave("all_clusters_g0_ctrl_v_g1_ctrl_scatter_plot.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.all, aes(G0_IFN, G1_IFN)) + geom_point(size=3) + ggtitle("G0_IFN vs. G1_IFN")
ggsave("all_clusters_g0_ifn_v_g1_ifn_scatter_plot.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.all, aes(G0_IFN_TUN, G1_IFN_TUN)) + geom_point(size=3) + ggtitle("G0_IFN_TUN vs. G1_IFN_TUN")
ggsave("all_clusters_g0_ifn_tun_v_g1_ifn_tun_scatter_plot.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.all, aes(G0_CTRL, G0_IFN)) + geom_point(size=3) + ggtitle("G0_CTRL vs. G0_IFN")
ggsave("all_clusters_g0_ctrl_v_g0_ifn_scatter_plot.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.all, aes(G0_CTRL, G0_IFN_TUN)) + geom_point(size=3) + ggtitle("G0_CTRL vs. G0_IFN_TUN")
ggsave("all_clusters_g0_ctrl_v_g0_ifn_tun_scatter_plot.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.all, aes(G1_CTRL, G1_IFN)) + geom_point(size=3) + ggtitle("G1_CTRL vs. G1_IFN")
ggsave("all_clusters_g1_ctrl_v_g1_ifn_scatter_plot.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.all, aes(G1_CTRL, G1_IFN_TUN)) + geom_point(size=3) + ggtitle("G1_CTRL vs. G1_IFN_TUN")
ggsave("all_clusters_g1_ctrl_v_g1_ifn_tun_scatter_plot.png", p1, width=20, height=20, units="cm")

###########################################################################################################
##############  CLUS_5_7_10  ##############################################################################
###########################################################################################################
## Subcluster Tubule and 2 GEC Clusters
clus5_7_10 <- SubsetData(all.combined, ident.use = c("5", "7", "10"))
Idents(object = clus5_7_10) <- "stim"
avg.clus5_7_10 <- log1p(x = AverageExpression(object = clus5_7_10, verbose = FALSE)$RNA)
avg.clus5_7_10$gene <- rownames(x = avg.clus5_7_10)

p1 <- ggplot(avg.clus5_7_10, aes(G0_CTRL, G1_CTRL)) + geom_point() + ggtitle("G0_CTRL vs. G1_CTRL")
ggsave("clus5_7_10_g0_ctrl_v_g1_ctrl_scatter_plot.png", p1, width=53, height=29, units="cm")

p1 <- ggplot(avg.clus5_7_10, aes(G0_IFN, G1_IFN)) + geom_point() + ggtitle("G0_IFN vs. G1_IFN")
ggsave("clus5_7_10_g0_ifn_v_g1_ifn_scatter_plot.png", p1, width=53, height=29, units="cm")

p1 <- ggplot(avg.clus5_7_10, aes(G0_IFN_TUN, G1_IFN_TUN)) + geom_point() + ggtitle("G0_IFN_TUN vs. G1_IFN_TUN")
ggsave("clus5_7_10_g0_ifn_tun_v_g1_ifn_tun_scatter_plot.png", p1, width=53, height=29, units="cm")

p1 <- ggplot(avg.clus5_7_10, aes(G0_CTRL, G0_IFN)) + geom_point() + ggtitle("G0_CTRL vs. G0_IFN")
ggsave("clus5_7_10_g0_ctrl_v_g0_ifn_scatter_plot.png", p1, width=53, height=29, units="cm")

p1 <- ggplot(avg.clus5_7_10, aes(G0_CTRL, G0_IFN_TUN)) + geom_point() + ggtitle("G0_CTRL vs. G0_IFN_TUN")
ggsave("clus5_7_10_g0_ctrl_tun_v_g0_ifn_tun_scatter_plot.png", p1, width=53, height=29, units="cm")

p1 <- ggplot(avg.clus5_7_10, aes(G1_CTRL, G1_IFN)) + geom_point() + ggtitle("G1_CTRL vs. G1_IFN")
ggsave("clus5_7_10_g1_ctrl_v_g1_ifn_scatter_plot.png", p1, width=53, height=29, units="cm")

p1 <- ggplot(avg.clus5_7_10, aes(G1_CTRL, G1_IFN_TUN)) + geom_point() + ggtitle("G1_CTRL vs. G1_IFN_TUN")
ggsave("clus5_7_10_g1_ctrl_v_g1_ifn_tun_scatter_plot.png", p1, width=53, height=29, units="cm")

###########################################################################################################
##############  CLUS_7  ###################################################################################
###########################################################################################################
Idents(all.combined) <- all.combined@meta.data$seurat_clusters

## Subcluster Cluster 7
clus7 <- SubsetData(all.combined, ident.use = c("7"))
Idents(object = clus7) <- "stim"
avg.clus7 <- log1p(x = AverageExpression(object = clus7, verbose = FALSE)$RNA)
avg.clus7$gene <- rownames(x = avg.clus7)

p1 <- ggplot(avg.clus7, aes(G0_IFN, G1_IFN)) + geom_point(size=3) + ggtitle("G0_IFN vs. G1_IFN") +
    geom_point(data=avg.clus7[c('MT1H', 'SLC2A3', 'LGALS3', 'RPS4Y1', 'DDIT4', 'LY6E', 'IGF2', 'CLDN1', 'WFDC2', 'CD24', 'PSMA2', 'ANXA4',
                                'JAG1', 'MAL', 'MPC2', 'CLU', 'FXYD2'),], color="red")
p1 <- LabelPoints(plot = p1, points = c('MT1H', 'SLC2A3', 'LGALS3', 'RPS4Y1', 'DDIT4', 'LY6E', 'IGF2', 'CLDN1', 'WFDC2', 'CD24', 'PSMA2', 'ANXA4',
                                        'JAG1', 'MAL', 'MPC2', 'CLU', 'FXYD2'))
ggsave("clus7_g0_ifn_v_g1_ifn_scatter_plot_w_labels.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.clus7, aes(G0_IFN, G1_IFN)) + geom_point(size=3) + ggtitle("G0_IFN vs. G1_IFN") +
    geom_point(data=avg.clus7[c('MT1H', 'SLC2A3', 'LGALS3', 'RPS4Y1', 'DDIT4', 'LY6E', 'IGF2', 'CLDN1', 'WFDC2', 'CD24', 'PSMA2', 'ANXA4',
                                'JAG1', 'MAL', 'MPC2', 'CLU', 'FXYD2'),], color="red")
ggsave("clus7_g0_ifn_v_g1_ifn_scatter_plot_wout_labels.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.clus7, aes(G0_IFN_TUN, G1_IFN_TUN)) + geom_point(size=3) + ggtitle("G0_IFN_TUN vs. G1_IFN_TUN") +
    geom_point(data=avg.clus7[c('MYL9', 'RPS4Y1', 'KRT17', 'ANXA1', 'COMT', 'EMP3', 'CLU', 'CLDN1', 'RARRES3', 'JAG1', 'FXYD2', 'MT1F'),], color="red")
p1 <- LabelPoints(plot = p1, points = c('MYL9', 'RPS4Y1', 'KRT17', 'ANXA1', 'COMT', 'EMP3', 'CLU', 'CLDN1', 'RARRES3', 'JAG1', 'FXYD2', 'MT1F'))
ggsave("clus7_g0_ifn_tun_v_g1_ifn_tun_scatter_plot_w_labels.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.clus7, aes(G0_IFN_TUN, G1_IFN_TUN)) + geom_point(size=3) + ggtitle("G0_IFN_TUN vs. G1_IFN_TUN") +
    geom_point(data=avg.clus7[c('MYL9', 'RPS4Y1', 'KRT17', 'ANXA1', 'COMT', 'EMP3', 'CLU', 'CLDN1', 'RARRES3', 'JAG1', 'FXYD2', 'MT1F'),], color="red")
ggsave("clus7_g0_ifn_tun_v_g1_ifn_tun_scatter_plot_wout_labels.png", p1, width=20, height=20, units="cm")

###########################################################################################################
##############  CLUS_10  ###################################################################################
###########################################################################################################
## Subcluster Cluster 10
clus10 <- SubsetData(all.combined, ident.use = c("10"))
Idents(object = clus10) <- "stim"
avg.clus10 <- log1p(x = AverageExpression(object = clus10, verbose = FALSE)$RNA)
avg.clus10$gene <- rownames(x = avg.clus10)

p1 <- ggplot(avg.clus10, aes(G0_IFN, G1_IFN)) + geom_point(size=3) + ggtitle("G0_IFN vs. G1_IFN") +
    geom_point(data=avg.clus10[c('RPS4Y1', 'LY6E', 'LDHA', 'MYL9', 'DLK1', 'S100A6', 'POSTN', 'RGS5', 'SERPINE3', 'ASS1', 'SOX4', 'WT1', 'EPHB2'),],
               color="red")
p1 <- LabelPoints(plot = p1, points = c('RPS4Y1', 'LY6E', 'LDHA', 'MYL9', 'DLK1', 'S100A6', 'POSTN', 'RGS5', 'SERPINE3', 'ASS1', 'SOX4', 'WT1', 'EPHB2'))
ggsave("clus10_g0_ifn_v_g1_ifn_scatter_plot_w_labels.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.clus10, aes(G0_IFN, G1_IFN)) + geom_point(size=3) + ggtitle("G0_IFN vs. G1_IFN") +
    geom_point(data=avg.clus10[c('RPS4Y1', 'LY6E', 'LDHA', 'MYL9', 'DLK1', 'S100A6', 'POSTN', 'RGS5', 'SERPINE3', 'ASS1', 'SOX4', 'WT1', 'EPHB2'),],
               color="red")
ggsave("clus10_g0_ifn_v_g1_ifn_scatter_plot_wout_labels.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.clus10, aes(G0_IFN_TUN, G1_IFN_TUN)) + geom_point(size=3) + ggtitle("G0_IFN_TUN vs. G1_IFN_TUN") +
    geom_point(data=avg.clus10[c('DLK1', 'LY6E', 'COL1A2', 'NPNT', 'S100A6', 'MYL9', 'ITM2B', 'HLA-A', 'NFIB', 'NUPR1', 'GABRA2', 'ASS1', 'CHCHD2',
                                 'UCHL1'),], color="red")
p1 <- LabelPoints(plot = p1, points = c('DLK1', 'LY6E', 'COL1A2', 'NPNT', 'S100A6', 'MYL9', 'ITM2B', 'HLA-A', 'NFIB', 'NUPR1', 'GABRA2', 'ASS1', 'CHCHD2',
                                        'UCHL1'))
ggsave("clus10_g0_ifn_tun_v_g1_ifn_tun_scatter_plot_w_labels.png", p1, width=20, height=20, units="cm")

p1 <- ggplot(avg.clus10, aes(G0_IFN_TUN, G1_IFN_TUN)) + geom_point(size=3) + ggtitle("G0_IFN_TUN vs. G1_IFN_TUN") +
    geom_point(data=avg.clus10[c('DLK1', 'LY6E', 'COL1A2', 'NPNT', 'S100A6', 'MYL9', 'ITM2B', 'HLA-A', 'NFIB', 'NUPR1', 'GABRA2', 'ASS1', 'CHCHD2',
                                 'UCHL1'),], color="red")
ggsave("clus10_g0_ifn_tun_v_g1_ifn_tun_scatter_plot_wout_labels.png", p1, width=20, height=20, units="cm")


all.combined$celltype.stim <- paste(Idents(object = all.combined), all.combined$stim, sep = "_")
all.combined$celltype <- Idents(object = all.combined)
Idents(object = all.combined) <- "celltype.stim"

## Create Marker List for Each Cluster (7, 11, 12)
for (clus in c('7', '11', '12'))
{
    for (sample in c('1', '2', '3', '4'))
    {
        if (sample == '1') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = all.combined, ident.1 = paste(clus, 'G0_IFN', sep='_'), ident.2 = paste(clus, 'G1_IFN', sep='_'), verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "g0_v_g1_ifn", append=FALSE)
        }
        else if (sample == '2') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = all.combined, ident.1 = paste(clus, 'G0_CTRL', sep='_'), ident.2 = paste(clus, 'G0_IFN', sep='_'), verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "g0_ctrl_v_g0_ifn", append=TRUE)
        }
        else if (sample == '3') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = all.combined, ident.1 = paste(clus, 'G1_CTRL', sep='_'), ident.2 = paste(clus, 'G1_IFN', sep='_'), verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "g1_ctrl_v_g1_ifn", append=TRUE)
        }
        else if (sample == '4') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = all.combined, ident.1 = paste(clus, 'G0_IFN_TUN', sep='_'), ident.2 = paste(clus, 'G1_IFN_TUN', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "g0_ifn_tun_v_g1_ifn_tun", append=TRUE)
        }
    }
}

## Cluster 14 isn't present in all samples so write.xlsx fails thus write.csv needed to be used
clus = '14'
for (sample in c('1', '2', '3', '4'))
    {
        if (sample == '1') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = all.combined, ident.1 = paste(clus, 'G0_IFN', sep='_'), ident.2 = paste(clus, 'G1_IFN', sep='_'), verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.csv(c, file=paste(cluster, 'g0_v_g1_ifn_res0.2_markers.csv', sep='_'), row.names=TRUE)
        }
        else if (sample == '2') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = all.combined, ident.1 = paste(clus, 'G0_CTRL', sep='_'), ident.2 = paste(clus, 'G0_IFN', sep='_'), verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.csv(c, file=paste(cluster, 'g0_ctrl_v_g0_ifn_res0.2_markers.csv', sep='_'), row.names=TRUE)
        }
        else if (sample == '3') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = all.combined, ident.1 = paste(clus, 'G1_CTRL', sep='_'), ident.2 = paste(clus, 'G1_IFN', sep='_'), verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.csv(c, file=paste(cluster, 'g1_ctrl_v_g1_ifn_res0.2_markers.csv', sep='_'), row.names=TRUE)
        }
        else if (sample == '4') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = all.combined, ident.1 = paste(clus, 'G0_IFN_TUN', sep='_'), ident.2 = paste(clus, 'G1_IFN_TUN', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.csv(c, file=paste(cluster, 'g0_v_g1_ifn_tun_res0.2_markers.csv', sep='_'), row.names=TRUE)
        }
    }

## ## An example of what type of analysis was performed
## cluster14 <- FindMarkers(object = all.combined, ident.1 = "14_G0_IFN", ident.2 = "14_G1_IFN", verbose = FALSE)
## cluster14 <- FindMarkers(object = all.combined, ident.1 = "14_G0_CTRL", ident.2 = "14_G0_IFN", verbose = FALSE)
## cluster14 <- FindMarkers(object = all.combined, ident.1 = "14_G1_CTRL", ident.2 = "14_G1_IFN", verbose = FALSE)
## cluster14 <- FindMarkers(object = all.combined, ident.1 = "14_G0_IFN_TUN", ident.2 = "14_G1_IFN_TUN", verbose = FALSE)

plots <- VlnPlot(object = gec_g0_g1_ctrl, features = c('WFDC2', 'CLDN4'), group.by='cluster.id', pt.size=0)
plots <- VlnPlot(object = gec_g0_g1_ctrl, features = c('SLC3A1', 'HNF1B'), group.by='cluster.id', pt.size=0)
plots <- VlnPlot(object = gec_g0_g1_ctrl, features = c('HNF4A', 'CUBN'),  group.by='cluster.id', pt.size=0)
plots <- VlnPlot(object = gec_g0_g1_ctrl, features = c('LRP2', 'DPP4'), group.by='cluster.id', pt.size=0)
plots <- VlnPlot(object = gec_g0_g1_ctrl, features = c('OSR1', 'WT1'), group.by='cluster.id', pt.size=0)
plots <- VlnPlot(object = gec_g0_g1_ctrl, features = c('PODXL', 'PTPRO'), group.by='cluster.id', pt.size=0)
plots <- VlnPlot(object = gec_g0_g1_ctrl, features = c('NPHS1', 'MAFB'), group.by='cluster.id', pt.size=0)
plots <- VlnPlot(object = gec_g0_g1_ctrl, features = c('CENPF', 'STMN2', 'GAP43'), group.by='cluster.id', pt.size=0)


## General Violin Plot for Certain Marker Genes
cluster.id <- unname(combined@active.ident)
combined@meta.data$cluster.id <- cluster.id

plots <- VlnPlot(object = combined, features = c('PODXL', 'PTPRO'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_podxl&ptpro_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('NPHS1', 'MAFB'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_nphs1&mafb_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('SLC3A1', 'TMEM176A'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_slc3a1&tmem176a_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('DPP4', 'CLDN4'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_dpp4&cldn4_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('WFDC2', 'IGFBP3'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_wfdc2&igfbp3_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('POSTN', 'GAP43'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_postn&gap43_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('STMN2', 'TOP2A'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_stmn2&top2a_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('CENPF', 'CD34'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_cenpf&cd34_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('TIE1'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_tie1_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('CUBN', 'LRP2'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_cubn_lrp2_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c('CDH2', 'CDH1'), group.by='cluster.id', pt.size=0,
                 cols=c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF'))
ggsave("vlnplot_cdh2_cdh1_res0.2.png", width=53, height=29, units="cm")


## Gene Specific FeaturePlots for Diff Exp
featureplot <- FeaturePlot(object = combined, features = c("PODXL"), cols = c("grey", "red"), split.by="stim", pt.size=0.2)
ggsave("featureplot_podxl_byStim_res0.2.png", featureplot, width=53, height=29, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("PODXL"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_podxl_res0.2.png", featureplot, width=53, height=29, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("APOL1"), cols = c("grey", "red"), split.by="stim", pt.size=0.2)
ggsave("featureplot_apol1_res0.2.png", featureplot, width=53, height=29, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("AHDC1"), cols = c("grey", "red"), split.by="stim", pt.size=0.2)
ggsave("featureplot_ahdc1_res0.2.png", featureplot, width=53, height=29, units="cm")

## Set Idents Prior to Create Diff Exp Plots
Idents(combined) <- combined@meta.data$stim

## Diff Exp Violin Plots for AHDC1
plots <- VlnPlot(object = combined, features = c("AHDC1"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_ahdc1_g0ifn_v_g1ifn_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("AHDC1"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_ahdc1_g0ifn_v_g1ifn_res0.2.png", width=53, height=29, units="cm")

#############################################################
## Diff Exp Violin Plots for PHGDH
plots <- VlnPlot(object = combined, features = c("PHGDH"), group.by="cluster.id", idents = c("G0_CTRL", "G0_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_ahdc1_g0ifn_v_g1ifn_res0.2.png", width=53, height=29, units="cm")

######################################################################################################################################
## Diff Exp Violin Plots for APOL1
plots <- VlnPlot(object = combined, features = c("APOL1"), , group.by="cluster.id", idents = c("G0_IFN", "G0_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_apol1_g0ifn_v_g0ifn_res0.2.png", width=40, height=40, units="cm")

###
plots <- VlnPlot(object = combined, features = c("APOL1"), , group.by="cluster.id", idents = c("G0_IFN", "G1_CTRL"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_apol1_g0ifn_v_g1ctrl_res0.2.png", width=40, height=40, units="cm")

###
plots <- VlnPlot(object = combined, features = c("APOL1"), , group.by="cluster.id", idents = c("G0_IFN", "G0_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_apol1_g0ifn_v_g0ifn_res0.2.png", width=40, height=40, units="cm")

###
plots <- VlnPlot(object = combined, features = c("APOL1"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_apol1_g0ifn_v_g1ifn_res0.2.png", width=40, height=40, units="cm")

###
plots <- VlnPlot(object = combined, features = c("APOL1"), , group.by="cluster.id", idents = c("G1_IFN", "G0_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_apol1_g1ifn_v_g0ifntun_res0.2.png", width=40, height=40, units="cm")

###
plots <- VlnPlot(object = combined, features = c("APOL1"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_apol1_g0ifntun_v_g1ifntun_res0.2.png", width=40, height=40, units="cm")

#####################################################################################################################################
## Diff Exp Violin Plots for CXCL9, CXCL11, MUG13, SNORA14B, PPP1R18
plots <- VlnPlot(object = combined, features = c("CXCL9", "CXCL11"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_cxcl9&cxcl11_g0ifn_v_g1ifn_res0.2.png", width=53, height=29, units="cm")

###
plots <- VlnPlot(object = combined, features = c("MUG13", "SNORA14B"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_mug13&snora14b_g0ifn_v_g1ifn_res0.2.png", width=53, height=29, units="cm")

###
plots <- VlnPlot(object = combined, features = c("PPP1R18"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_ppp1r18_g0ifn_v_g1ifn_res0.2.png", width=53, height=29, units="cm")

###
plots <- VlnPlot(object = combined, features = c("CXCL9", "CXCL11"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_cxcl9&cxcl11_g0ifntun_v_g1ifntun_res0.2.png", width=53, height=29, units="cm")

###
plots <- VlnPlot(object = combined, features = c("MUG13", "SNORA14B"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_mug13&snora14b_g0ifntun_v_g1ifntun_res0.2.png", width=53, height=29, units="cm")

###
plots <- VlnPlot(object = combined, features = c("PPP1R18"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_ppp1r18_g0ifntun_v_g1ifntun_res0.2.png", width=53, height=29, units="cm")



## Diff Exp Violin Plots G0_v_G1_IFN
plots <- VlnPlot(object = combined, features = c("UCHL1", "CHCHD2"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_uchl1&chchd2_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("ASS1", "UBD"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_ass1&ubd_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("SPARC", "PODXL"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_sparc&podxl_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("WT1", "PTPRO"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_wt1&ptpro_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("HLA-A", "HLA-C"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_hla-a&hla-c_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("S100A6"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_s100a6_res0.2.png", width=53, height=29, units="cm")

## Diff Exp VlnPlot EDEM1, IL6, HSPA5
plots <- VlnPlot(object = combined, features = c("EDEM1", "HSPA5"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"),
                 split.by="stim", pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_edem1_hspa5_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("EDEM1", "HSPA5"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"),
                 split.by="stim", pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_edem1_hspa5_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("IL6"), , group.by="cluster.id", idents = c("G0_IFN", "G1_IFN"),
                 split.by="stim", pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_il6_res0.2.png", width=53, height=29, units="cm")


## UMAP FeaturePlot By Markers
featureplot <- FeaturePlot(object = combined, features = c("PODXL", "NPHS1", "SPARC", "NPSH2"), cols = c("grey", "red"), split.by=NULL, pt.size=0.2)
ggsave("featureplot_podxl&nphs1&sparc&npsh2_res0.2.png", featureplot, width=53, height=29, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("SLC3A1", "EPCAM", "CD34", "TIE1"), cols = c("grey", "red"), split.by=NULL, pt.size=0.2)
ggsave("featureplot_slc3a1&epcam&cd34&tie1_res0.2.png", featureplot, width=53, height=29, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("GAP43"), cols = c("grey", "red"), split.by=NULL, pt.size=0.2)
ggsave("featureplot_gap43_res0.2.png", featureplot, width=53, height=29, units="cm")

## Diff Exp Violing Plots G0_v_G1_IFN_TUN
plots <- VlnPlot(object = combined, features = c("UCHL1", "CHCHD2"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_uchl1&chchd2_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("ASS1", "UBD"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_ass1&ubd_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("SPARC", "PODXL"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_sparc&podxl_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("WT1", "PTPRO"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_wt1&ptpro_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("HLA-A", "HLA-C"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_hla-a&hla-c_res0.2.png", width=53, height=29, units="cm")

plots <- VlnPlot(object = combined, features = c("S100A6"), , group.by="cluster.id", idents = c("G0_IFN_TUN", "G1_IFN_TUN"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_s100a6_res0.2.png", width=53, height=29, units="cm")

## Diff Exp Violing Plots G0_v_G1_CTRL
plots <- VlnPlot(object = combined, features = c("PODXL", "EPCAM"), , group.by="cluster.id", idents = c("G0_CTRL", "G1_CTRL"), split.by="stim",
                 pt.size=0, combine=FALSE)
plots <- lapply(
    X = plots,
    FUN = function(p) p + ggplot2::scale_fill_manual(values = c('grey', 'red'))
)
CombinePlots(plots = plots, legend = 'right')
ggsave("diffexp_vlnplot_podxl&epcam_res0.2.png", width=53, height=29, units="cm")

## DotPlot
genes <- c('ACTN4', 'ANLN', 'ARHGAP24', 'CD2AP', 'CFI', 'COL4A3', 'COL4A4', 'E2F3', 'INF2', 'LMX1B', 'PAX2', 'TRPC6', 'WT1', 'COQ6', 'CRB2',
           'ITGB4', 'LAMA5', 'MYO1E', 'NPHP4', 'NUP133', 'NUP160', 'NUP85', 'TTC21B', 'ADCK4', 'ALG1', 'ANKFY1', 'MAFB', 'ARHGDIA',
           'CUBN', 'DGKE', 'DHTKD1', 'DLC1', 'EMP2', 'FAT1', 'GAPVD1', 'LAMB2', 'NPHS1', 'NPHS2', 'NUP107', 'NUP205', 'NUP93', 'PLCE1',
           'PDSS2', 'PMM2', 'PTPRO', 'SCARB2', 'XPO5', 'ZMPSTE24')

DotPlot(combined, features=genes, dot.scale = 8, split.by = 'stim', cols = c('blue', 'red', 'green', 'navy', 'orange', 'purple', 'limegreen')) + RotatedAxis()

DotPlot(combined, features=genes, dot.scale = 8, split.by = c('cluster.id'), cols = c('blue', 'red', 'green', 'navy', 'orange', 'purple')) + RotatedAxis()

group.by=c('Early GEC', 'Proximal and Distal Tubule', 'Glomerular Epithelial Cells', 'Endothelial Cells')


########################
### Monocle Analysis ###
########################
library(monocle)
library(Seurat)
library(reticulate)
library(hdf5r)
library(dplyr)

##########################################################
## Subset Clusters
sub_g0_ifn_tun <- SubsetData(all.combined, ident.use = c("0", "4", "5", "7", "8", "10", "14"))
new.cluster.ids <- c( "Mesenchyme", "Proliferating Cells", "Early GEC", "Proximal and Distal Tubule", "Neuron", "Glomerular Epithelial Cells",
                     "Endothelial Cells")
names(x = new.cluster.ids) <- c( "0", "4", "5", "7", "8", "10", "14")
sub_g0_ifn_tun <- RenameIdents(object = sub_g0_ifn_tun, new.cluster.ids)

## Normalize Data
sub_g0_ifn_tun <- NormalizeData(sub_g0_ifn_tun)

## Find Variable Features
sub_g0_ifn_tun <- FindVariableFeatures(sub_g0_ifn_tun, selection.method = "vst", nfeatures=2000)

## ScaleData
all.genes <- rownames(sub_g0_ifn_tun)
sub_g0_ifn_tun <- ScaleData(object = sub_g0_ifn_tun)

## Perfrom Linear Dimensional Reduction
sub_g0_ifn_tun <- RunPCA(object = sub_g0_ifn_tun, features = VariableFeatures(object = sub_g0_ifn_tun), pcs.compute=30)

## Cluster Cells
sub_g0_ifn_tun <- FindNeighbors(object = sub_g0_ifn_tun, dims = 1:5)
sub_g0_ifn_tun <- FindClusters(object = sub_g0_ifn_tun, resolution = 0.6)

## Create UMAP / tSNE Plots
sub_g0_ifn_tun <- RunUMAP(object = sub_g0_ifn_tun, dims = 1:30)

## Extract data, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(sub_g0_ifn_tun@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sub_g0_ifn_tun@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

## Construct Monocle CDS
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

## Estimate Size Factors and Dispersion
HSMM <- monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

## Filter low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

## Identify CellTypes (Gating Functions)
podxl <- row.names(subset(fData(HSMM), gene_short_name == "PODXL"))  ## Glomerular
tmem176a  <- row.names(subset(fData(HSMM), gene_short_name == "TMEM176A")) ## Proximal Tubule
tmem176b  <- row.names(subset(fData(HSMM), gene_short_name == "TMEM176B")) ## Proximal Tubule
epcam <- row.names(subset(fData(HSMM), gene_short_name == "EPCAM")) ## Proximal Tubule
igfbp3 <- row.names(subset(fData(HSMM), gene_short_name == "IGFBP3")) ## Mesenchyme
postn <- row.names(subset(fData(HSMM), gene_short_name == "POSTN")) ## Mesenchyme
fap <- row.names(subset(fData(HSMM), gene_short_name == "FAP")) ##
dnabj9 <- row.names(subset(fData(HSMM), gene_short_name == "DNABJ9")) ## ER STRESS
hspa5 <- row.names(subset(fData(HSMM), gene_short_name == "HSPA5")) ## ER STRESS

## Initiate Cell Type Classification
cth <- newCellTypeHierarchy()

## Create Gating Functions
cth <- addCellType(cth, "Glomerular", classify_func = function(x)
{ x[podxl,] >= 1 })
cth <- addCellType(cth, "Proximal Tubule", classify_func = function(x)
{ x[epcam,] >= 1 | x[tmem176a,] >= 1 | x[tmem176b,] >= 1 })
cth <- addCellType(cth, "Mesenchyme", classify_func = function(x)
{  x[fap, ] >= 1 })
cth <- addCellType(cth, "ER STRESS", classify_func = function(x)
{ x[hspa5, ] >= 1 })

## ClassifyCells via Gating Function
HSMM <- classifyCells(HSMM, cth, 0.1)

## Perform Supervised Clustering & Reduce Dimension of Cells
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 30,
                        norm_method = 'log', reduction_method = 'tSNE',
                        residualModelFormulaStr = "~num_genes_expressed",
                        verbose = T)

HSMM1 <- clusterCells(HSMM, num_clusters = 4)

## Seurat Clusters into Monocle
sub_g0_ifn_tun <- SubsetData(all.combined, ident.use = c("0", "4", "5", "7", "8", "10", "14"))
new.cluster.ids <- c( "Mesenchyme", "Proliferating Cells", "Early GEC", "Proximal and Distal Tubule", "Neuron", "Glomerular Epithelial Cells",
                     "Endothelial Cells")
names(x = new.cluster.ids) <- c("0", "4", "5", "7", "8", "10", "14")
sub_g0_ifn_tun <- RenameIdents(object = sub_g0_ifn_tun, new.cluster.ids)
DimA <- t(sub_g0_ifn_tun@reductions$umap@cell.embeddings)
rownames(DimA) <- 1:2
HSMM1@reducedDimA <- DimA
metadata <- pData(HSMM1)
metadata$Cluster <- sub_g0_ifn_tun@active.ident
pData(HSMM1) <- metadata

## Trajectory Analysis Unsupervised Clustering
HSMM3 <- detectGenes(HSMM1, min_expr = 0.1)
fData(HSMM3)$use_for_ordering <- fData(HSMM3)$num_cells_expressed > 0.05 * ncol(HSMM3)

## Set Components (num_dim) to Top 30
HSMM4 <- reduceDimension(HSMM3,
                         max_components = 2,
                         norm_method = 'log',
                         num_dim = 3,
                         reduction_method = 'tSNE',
                         verbose = T)

HSMM4 <- clusterCells(HSMM4)

## Extract Cluster Distinguishing Clusters
clustering_DEG_genes <-
    differentialGeneTest(HSMM4[expressed_genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores=24)

HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

HSMM5 <-
    setOrderingFilter(HSMM4,
                      ordering_genes = HSMM_ordering_genes)

HSMM6 <- reduceDimension(HSMM5, method = 'DDRTree')

HSMM6 <- orderCells(HSMM6)

## Attach Correct IDs
cluster.id <- unname(sub_g0_ifn_tun@active.ident)
HSMM6@phenoData@data$cluster.id <- cluster.id

colors <- c('#A6E9A9', '#49C64E', '#FCB0B0', '#0000FF', '#FFEE20', '#FF0000', '#A600FF')

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "cluster.id", cell_size = .5) +
    scale_color_manual(values = colors)

ggsave(filename="g0_v_g1_ctrl_ifn_n_tun_trajectory_plot.png", trajectory_plot, width=53, height=29, units="cm")

quit("yes")

#############################################################################
########## Monocle Analysis All Glomerular & Epithelial Cells ################
#############################################################################

## Subset Clusters
sub_gEpi <- SubsetData(all.combined, ident.use = c("5", "7", "10", "14"))
new.cluster.ids <- c("Early GEC", "Proximal and Distal Tubule", "Glomerular Epithelial Cells", "Endothelial Cells")
names(x = new.cluster.ids) <- c("5", "7", "10", "14")
sub_gEpi <- RenameIdents(object = sub_gEpi, new.cluster.ids)

## Normalize Data
sub_gEpi <- NormalizeData(sub_gEpi)

## Find Variable Features
sub_gEpi <- FindVariableFeatures(sub_gEpi, selection.method = "vst", nfeatures=2000)

## ScaleData
all.genes <- rownames(sub_gEpi)
sub_gEpi <- ScaleData(object = sub_gEpi)

## Perfrom Linear Dimensional Reduction
sub_gEpi <- RunPCA(object = sub_gEpi, features = VariableFeatures(object = sub_gEpi), pcs.compute=30)

## Cluster Cells
sub_gEpi <- FindNeighbors(object = sub_gEpi, dims = 1:5)
sub_gEpi <- FindClusters(object = sub_gEpi, resolution = 0.6)

## Create UMAP / tSNE Plots
sub_gEpi <- RunUMAP(object = sub_gEpi, dims = 1:30)

## Extract data, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(sub_gEpi@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sub_gEpi@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

## Construct Monocle CDS
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

## Estimate Size Factors and Dispersion
HSMM <- monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

## Filter low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

## Identify CellTypes (Gating Functions)
podxl <- row.names(subset(fData(HSMM), gene_short_name == "PODXL"))  ## Glomerular
tmem176a  <- row.names(subset(fData(HSMM), gene_short_name == "TMEM176A")) ## Proximal Tubule
tmem176b  <- row.names(subset(fData(HSMM), gene_short_name == "TMEM176B")) ## Proximal Tubule
epcam <- row.names(subset(fData(HSMM), gene_short_name == "EPCAM")) ## Proximal Tubule
igfbp3 <- row.names(subset(fData(HSMM), gene_short_name == "IGFBP3")) ## Mesenchyme
postn <- row.names(subset(fData(HSMM), gene_short_name == "POSTN")) ## Mesenchyme
fap <- row.names(subset(fData(HSMM), gene_short_name == "FAP")) ##
dnabj9 <- row.names(subset(fData(HSMM), gene_short_name == "DNABJ9")) ## ER STRESS
hspa5 <- row.names(subset(fData(HSMM), gene_short_name == "HSPA5")) ## ER STRESS

## Initiate Cell Type Classification
cth <- newCellTypeHierarchy()

## Create Gating Functions
cth <- addCellType(cth, "Glomerular", classify_func = function(x)
{ x[podxl,] >= 1 })
cth <- addCellType(cth, "Proximal Tubule", classify_func = function(x)
{ x[epcam,] >= 1 | x[tmem176a,] >= 1 | x[tmem176b,] >= 1 })
cth <- addCellType(cth, "Mesenchyme", classify_func = function(x)
{  x[fap, ] >= 1 })
cth <- addCellType(cth, "ER STRESS", classify_func = function(x)
{ x[hspa5, ] >= 1 })

## ClassifyCells via Gating Function
HSMM <- classifyCells(HSMM, cth, 0.1)

## Perform Supervised Clustering & Reduce Dimension of Cells
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 30,
                        norm_method = 'log', reduction_method = 'tSNE',
                        residualModelFormulaStr = "~num_genes_expressed",
                        verbose = T)

HSMM1 <- clusterCells(HSMM, num_clusters = 4)

## Seurat Clusters into Monocle
sub_gEpi <- SubsetData(all.combined, ident.use = c("5", "7", "10", "14"))
new.cluster.ids <- c("Early GEC", "Proximal and Distal Tubule", "Glomerular Epithelial Cells", "Endothelial Cells")
names(x = new.cluster.ids) <- c("5", "7", "10", "14")
sub_gEpi <- RenameIdents(object = sub_gEpi, new.cluster.ids)
DimA <- t(sub_gEpi@reductions$umap@cell.embeddings)
rownames(DimA) <- 1:2
HSMM1@reducedDimA <- DimA
metadata <- pData(HSMM1)
metadata$Cluster <- sub_gEpi@active.ident
pData(HSMM1) <- metadata

## Trajectory Analysis Unsupervised Clustering
HSMM3 <- detectGenes(HSMM1, min_expr = 0.1)
fData(HSMM3)$use_for_ordering <- fData(HSMM3)$num_cells_expressed > 0.05 * ncol(HSMM3)

## Set Components (num_dim) to Top 30
HSMM4 <- reduceDimension(HSMM3,
                         max_components = 2,
                         norm_method = 'log',
                         num_dim = 3,
                         reduction_method = 'tSNE',
                         verbose = T)

HSMM4 <- clusterCells(HSMM4)

## Extract Cluster Distinguishing Clusters
clustering_DEG_genes <-
    differentialGeneTest(HSMM4[expressed_genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores=24)

HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

HSMM5 <-
    setOrderingFilter(HSMM4,
                      ordering_genes = HSMM_ordering_genes)

HSMM6 <- reduceDimension(HSMM5, method = 'DDRTree')

HSMM6 <- orderCells(HSMM6)

## Attach Correct IDs
cluster.id <- unname(sub_gEpi@active.ident)
HSMM6@phenoData@data$cluster.id <- cluster.id

colors <- c('#FCB0B0', '#0000FF', '#FF0000', '#A600FF')

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "cluster.id", cell_size = .5) +
    scale_color_manual(values = colors)

ggsave(filename="all_glomerular_n_epithelial_cells_cluster_plot.png", trajectory_plot, width=53, height=29, units="cm")

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "orig.ident", cell_size = .5)

ggsave(filename="all_glomerular_n_epithelial_cells_sample_date_specific_plot.png", trajectory_plot, width=53, height=29, units="cm")

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "stim", cell_size = .5)

ggsave(filename="all_glomerular_n_epithelial_cells_sample_specific_plot.png", trajectory_plot, width=53, height=29, units="cm")

quit("yes")


#############################################################################
########## Monocle Analysis Only Glomerular Epithelial Cells ################
#############################################################################

## Subset Clusters
sub_gEpi <- SubsetData(all.combined, ident.use = c("10"))
new.cluster.ids <- c("Glomerular Epithelial Cells")
names(x = new.cluster.ids) <- c("10")
sub_gEpi <- RenameIdents(object = sub_gEpi, new.cluster.ids)

## Normalize Data
sub_gEpi <- NormalizeData(sub_gEpi)

## Find Variable Features
sub_gEpi <- FindVariableFeatures(sub_gEpi, selection.method = "vst", nfeatures=2000)

## ScaleData
all.genes <- rownames(sub_gEpi)
sub_gEpi <- ScaleData(object = sub_gEpi)

## Perfrom Linear Dimensional Reduction
sub_gEpi <- RunPCA(object = sub_gEpi, features = VariableFeatures(object = sub_gEpi), pcs.compute=30)

## Cluster Cells
sub_gEpi <- FindNeighbors(object = sub_gEpi, dims = 1:5)
sub_gEpi <- FindClusters(object = sub_gEpi, resolution = 0.6)

## Create UMAP / tSNE Plots
sub_gEpi <- RunUMAP(object = sub_gEpi, dims = 1:30)

## Extract data, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(sub_gEpi@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sub_gEpi@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

## Construct Monocle CDS
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

## Estimate Size Factors and Dispersion
HSMM <- monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

## Filter low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

## Identify CellTypes (Gating Functions)
podxl <- row.names(subset(fData(HSMM), gene_short_name == "PODXL"))  ## Glomerular
tmem176a  <- row.names(subset(fData(HSMM), gene_short_name == "TMEM176A")) ## Proximal Tubule
tmem176b  <- row.names(subset(fData(HSMM), gene_short_name == "TMEM176B")) ## Proximal Tubule
epcam <- row.names(subset(fData(HSMM), gene_short_name == "EPCAM")) ## Proximal Tubule
igfbp3 <- row.names(subset(fData(HSMM), gene_short_name == "IGFBP3")) ## Mesenchyme
postn <- row.names(subset(fData(HSMM), gene_short_name == "POSTN")) ## Mesenchyme
fap <- row.names(subset(fData(HSMM), gene_short_name == "FAP")) ##
dnabj9 <- row.names(subset(fData(HSMM), gene_short_name == "DNABJ9")) ## ER STRESS
hspa5 <- row.names(subset(fData(HSMM), gene_short_name == "HSPA5")) ## ER STRESS

## Initiate Cell Type Classification
cth <- newCellTypeHierarchy()

## Create Gating Functions
cth <- addCellType(cth, "Glomerular", classify_func = function(x)
{ x[podxl,] >= 1 })
cth <- addCellType(cth, "Proximal Tubule", classify_func = function(x)
{ x[epcam,] >= 1 | x[tmem176a,] >= 1 | x[tmem176b,] >= 1 })
cth <- addCellType(cth, "Mesenchyme", classify_func = function(x)
{  x[fap, ] >= 1 })
cth <- addCellType(cth, "ER STRESS", classify_func = function(x)
{ x[hspa5, ] >= 1 })

## ClassifyCells via Gating Function
HSMM <- classifyCells(HSMM, cth, 0.1)

## Perform Supervised Clustering & Reduce Dimension of Cells
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 30,
                        norm_method = 'log', reduction_method = 'tSNE',
                        residualModelFormulaStr = "~num_genes_expressed",
                        verbose = T)

HSMM1 <- clusterCells(HSMM, num_clusters = 4)

## Seurat Clusters into Monocle
sub_gEpi <- SubsetData(all.combined, ident.use = c("10"))
new.cluster.ids <- c("Glomerular Epithelial Cells")
names(x = new.cluster.ids) <- c("10")
sub_gEpi <- RenameIdents(object = sub_gEpi, new.cluster.ids)
DimA <- t(sub_gEpi@reductions$umap@cell.embeddings)
rownames(DimA) <- 1:2
HSMM1@reducedDimA <- DimA
metadata <- pData(HSMM1)
metadata$Cluster <- sub_gEpi@active.ident
pData(HSMM1) <- metadata

## Trajectory Analysis Unsupervised Clustering
HSMM3 <- detectGenes(HSMM1, min_expr = 0.1)
fData(HSMM3)$use_for_ordering <- fData(HSMM3)$num_cells_expressed > 0.05 * ncol(HSMM3)

## Set Components (num_dim) to Top 30
HSMM4 <- reduceDimension(HSMM3,
                         max_components = 2,
                         norm_method = 'log',
                         num_dim = 3,
                         reduction_method = 'tSNE',
                         verbose = T)

HSMM4 <- clusterCells(HSMM4)

## Extract Cluster Distinguishing Clusters
clustering_DEG_genes <-
    differentialGeneTest(HSMM4[expressed_genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores=24)

HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

HSMM5 <-
    setOrderingFilter(HSMM4,
                      ordering_genes = HSMM_ordering_genes)

HSMM6 <- reduceDimension(HSMM5, method = 'DDRTree')

HSMM6 <- orderCells(HSMM6)

## Attach Correct IDs
cluster.id <- unname(sub_gEpi@active.ident)
HSMM6@phenoData@data$cluster.id <- cluster.id

colors <- c('#FF0000')

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "cluster.id", cell_size = .5) +
    scale_color_manual(values = colors)

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "orig.ident", cell_size = 2.5)

ggsave(filename="glomerular_epithelial_cells_sample_date_specific_plot.png", trajectory_plot, width=53, height=29, units="cm")

quit("yes")

#########################################################################################
########## Monocle Analysis All Glomerular Epithelial Cells & Mesenchyme ################
#########################################################################################
library(monocle)
library(Seurat)
library(reticulate)
library(hdf5r)
library(dplyr)

## Seurat Clusters
mes <- SubsetData(all.combined, ident.use = c("0", "1"))
new.cluster.ids <- c("Mesenchyme 1", "Mesenchyme 2")
names(x = new.cluster.ids) <- c("0", "1")
mes <- RenameIdents(object = mes, new.cluster.ids)
Idents(mes) <- mes@meta.data$stim
mes <- SubsetData(mes, ident.use = c("G0_CTRL"))
Idents(mes) <- mes@meta.data$seurat_clusters
new.cluster.ids <- c("Mesenchyme 1", "Mesenchyme 2")
names(x = new.cluster.ids) <- c("0", "1")
mes <- RenameIdents(object = mes, new.cluster.ids)

sub_gEpi <- SubsetData(all.combined, ident.use = c("5", "7", "10", "14"))
new.cluster.ids <- c("Early GEC", "Proximal and Distal Tubule", "Glomerular Epithelial Cells", "Endothelial Cells")
names(x = new.cluster.ids) <- c("5", "7", "10", "14")
sub_gEpi <- RenameIdents(object = sub_gEpi, new.cluster.ids)
Idents(sub_gEpi) <- sub_gEpi@meta.data$stim
sub_gEpi <- SubsetData(sub_gEpi, ident.use = c("G0_CTRL", "G0_IFN", "G1_IFN", "G0_IFN_TUN", "G1_IFN_TUN"))
Idents(sub_gEpi) <- sub_gEpi@meta.data$seurat_clusters
new.cluster.ids <- c("Early GEC", "Proximal and Distal Tubule", "Glomerular Epithelial Cells", "Endothelial Cells")
names(x = new.cluster.ids) <- c("5", "7", "10", "14")
sub_gEpi <- RenameIdents(object = sub_gEpi, new.cluster.ids)

sub_gEpi_n_Mes <- merge(mes, y = sub_gEpi)

## Normalize Data
sub_gEpi_n_Mes <- NormalizeData(sub_gEpi_n_Mes)

## Find Variable Features
sub_gEpi_n_Mes <- FindVariableFeatures(sub_gEpi_n_Mes, selection.method = "vst", nfeatures=2000)

## ScaleData
all.genes <- rownames(sub_gEpi_n_Mes)
sub_gEpi_n_Mes <- ScaleData(object = sub_gEpi_n_Mes)

## Perfrom Linear Dimensional Reduction
sub_gEpi_n_Mes <- RunPCA(object = sub_gEpi_n_Mes, features = VariableFeatures(object = sub_gEpi_n_Mes), pcs.compute=30)

## Cluster Cells
sub_gEpi_n_Mes <- FindNeighbors(object = sub_gEpi_n_Mes, dims = 1:5)
sub_gEpi_n_Mes <- FindClusters(object = sub_gEpi_n_Mes, resolution = 0.6)

## Create UMAP / tSNE Plots
sub_gEpi_n_Mes <- RunUMAP(object = sub_gEpi_n_Mes, dims = 1:30)

## Extract data, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(sub_gEpi_n_Mes@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sub_gEpi_n_Mes@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

## Construct Monocle CDS
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

## Estimate Size Factors and Dispersion
HSMM <- monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

## Filter low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

## Identify CellTypes (Gating Functions)
podxl <- row.names(subset(fData(HSMM), gene_short_name == "PODXL"))  ## Glomerular
tmem176a  <- row.names(subset(fData(HSMM), gene_short_name == "TMEM176A")) ## Proximal Tubule
tmem176b  <- row.names(subset(fData(HSMM), gene_short_name == "TMEM176B")) ## Proximal Tubule
epcam <- row.names(subset(fData(HSMM), gene_short_name == "EPCAM")) ## Proximal Tubule
igfbp3 <- row.names(subset(fData(HSMM), gene_short_name == "IGFBP3")) ## Mesenchyme
postn <- row.names(subset(fData(HSMM), gene_short_name == "POSTN")) ## Mesenchyme
fap <- row.names(subset(fData(HSMM), gene_short_name == "FAP")) ##
dnabj9 <- row.names(subset(fData(HSMM), gene_short_name == "DNABJ9")) ## ER STRESS
hspa5 <- row.names(subset(fData(HSMM), gene_short_name == "HSPA5")) ## ER STRESS

## Initiate Cell Type Classification
cth <- newCellTypeHierarchy()

## Create Gating Functions
cth <- addCellType(cth, "Glomerular", classify_func = function(x)
{ x[podxl,] >= 1 })
cth <- addCellType(cth, "Proximal Tubule", classify_func = function(x)
{ x[epcam,] >= 1 | x[tmem176a,] >= 1 | x[tmem176b,] >= 1 })
cth <- addCellType(cth, "Mesenchyme", classify_func = function(x)
{  x[fap, ] >= 1 })
cth <- addCellType(cth, "ER STRESS", classify_func = function(x)
{ x[hspa5, ] >= 1 })

## ClassifyCells via Gating Function
HSMM <- classifyCells(HSMM, cth, 0.1)

## Perform Supervised Clustering & Reduce Dimension of Cells
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                        norm_method = 'log', reduction_method = 'tSNE',
                        residualModelFormulaStr = "~num_genes_expressed",
                        verbose = T)

HSMM1 <- clusterCells(HSMM, num_clusters = 2)

## Seurat Clusters into Monocle
mes <- SubsetData(all.combined, ident.use = c("0", "1"))
new.cluster.ids <- c("Mesenchyme 1", "Mesenchyme 2")
names(x = new.cluster.ids) <- c("0", "1")
mes <- RenameIdents(object = mes, new.cluster.ids)
Idents(mes) <- mes@meta.data$stim
mes <- SubsetData(mes, ident.use = c("G0_CTRL"))
Idents(mes) <- mes@meta.data$seurat_clusters
new.cluster.ids <- c("Mesenchyme 1", "Mesenchyme 2")
names(x = new.cluster.ids) <- c("0", "1")
mes <- RenameIdents(object = mes, new.cluster.ids)

sub_gEpi <- SubsetData(all.combined, ident.use = c("5", "7", "10", "14"))
new.cluster.ids <- c("Early GEC", "Proximal and Distal Tubule", "Glomerular Epithelial Cells", "Endothelial Cells")
names(x = new.cluster.ids) <- c("5", "7", "10", "14")
sub_gEpi <- RenameIdents(object = sub_gEpi, new.cluster.ids)
Idents(sub_gEpi) <- sub_gEpi@meta.data$stim
sub_gEpi <- SubsetData(sub_gEpi, ident.use = c("G0_CTRL", "G0_IFN", "G1_IFN", "G0_IFN_TUN", "G1_IFN_TUN"))
Idents(sub_gEpi) <- sub_gEpi@meta.data$seurat_clusters
new.cluster.ids <- c("Early GEC", "Proximal and Distal Tubule", "Glomerular Epithelial Cells", "Endothelial Cells")
names(x = new.cluster.ids) <- c("5", "7", "10", "14")
sub_gEpi <- RenameIdents(object = sub_gEpi, new.cluster.ids)

sub_gEpi_n_Mes <- merge(mes, y = sub_gEpi)

## Normalize Data
sub_gEpi_n_Mes <- NormalizeData(sub_gEpi_n_Mes)

## Find Variable Features
sub_gEpi_n_Mes <- FindVariableFeatures(sub_gEpi_n_Mes, selection.method = "vst", nfeatures=2000)

## ScaleData
all.genes <- rownames(sub_gEpi_n_Mes)
sub_gEpi_n_Mes <- ScaleData(object = sub_gEpi_n_Mes)

## Perfrom Linear Dimensional Reduction
sub_gEpi_n_Mes <- RunPCA(object = sub_gEpi_n_Mes, features = VariableFeatures(object = sub_gEpi_n_Mes), pcs.compute=30)

## Cluster Cells
sub_gEpi_n_Mes <- FindNeighbors(object = sub_gEpi_n_Mes, dims = 1:5)
sub_gEpi_n_Mes <- FindClusters(object = sub_gEpi_n_Mes, resolution = 0.6)

## Create UMAP / tSNE Plots
sub_gEpi_n_Mes <- RunUMAP(object = sub_gEpi_n_Mes, dims = 1:30)

DimA <- t(sub_gEpi_n_Mes@reductions$umap@cell.embeddings)
rownames(DimA) <- 1:2
HSMM1@reducedDimA <- DimA
metadata <- pData(HSMM1)
metadata$Cluster <- sub_gEpi_n_Mes@active.ident
pData(HSMM1) <- metadata

## Trajectory Analysis Unsupervised Clustering
HSMM3 <- detectGenes(HSMM1, min_expr = 0.1)
fData(HSMM3)$use_for_ordering <- fData(HSMM3)$num_cells_expressed > 0.05 * ncol(HSMM3)

## Set Components (num_dim) to Top 30
HSMM4 <- reduceDimension(HSMM3,
                         max_components = 2,
                         norm_method = 'log',
                         num_dim = 3,
                         reduction_method = 'tSNE',
                         verbose = T)

HSMM4 <- clusterCells(HSMM4)

## Extract Cluster Distinguishing Clusters clustering_DEG_genes <-
clustering_DEG_genes <-  differentialGeneTest(HSMM4[expressed_genes,], fullModelFormulaStr = '~Cluster', cores=24)

HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

HSMM5 <-
    setOrderingFilter(HSMM4,
                      ordering_genes = HSMM_ordering_genes)

HSMM6 <- reduceDimension(HSMM5, method = 'DDRTree')

HSMM6 <- orderCells(HSMM6)

## Attach Correct IDs
cluster.id <- unname(sub_gEpi_n_Mes@active.ident)
HSMM6@phenoData@data$cluster.id <- cluster.id

colors <- c('#FCB0B0', '#0000FF', '#FF0000', '#A600FF')

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "cluster.id", cell_size = .5) +
    scale_color_manual(values = colors)

ggsave(filename="all_glomerular_n_epithelial_cells_cluster_plot.png", trajectory_plot, width=53, height=29, units="cm")

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "orig.ident", cell_size = .5)

ggsave(filename="all_glomerular_n_epithelial_cells_sample_date_specific_plot.png", trajectory_plot, width=53, height=29, units="cm")

trajectory_plot <-
    plot_cell_trajectory(HSMM6, color_by = "stim", cell_size = .5)

ggsave(filename="all_glomerular_n_epithelial_cells_sample_specific_plot.png", trajectory_plot, width=53, height=29, units="cm")

quit("yes")

#############################################################################
## Monocle 3 Analysis
#############################################################################
library(Seurat)
library(monocle3)
library(reticulate)
library(hdf5r)
library(dplyr)

## Set Identities to cluster.id
Idents(combined) <- combined@meta.data$cluster.id

## Extract counts, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(combined@assays$integrated@data), 'sparseMatrix')
pd <- combined@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

## Construct Monocle CDS
cds <- new_cell_data_set(data,
                         cell_metadata = pd,
                         gene_metadata = fData)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)

## Use cell.cluster Labels w/ Monocle3 UMAP
plot_cells(cds, color_cells_by="cluster.id")

plot_cells(cds, color_cells_by="seurat_clusters", group_label_size=4)

cds <- cluster_cells(cds)

cds <- learn_graph(cds)

plot_cells(cds, color_cells_by = "cluster.id",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           graph_label_size=1.5)

plot_cells(cds,
           color_cells_by = "stim",
           graph_label_size=1.5,
           label_groups_by_cluster=FALSE,
           label_cell_groups=FALSE)

g1ifntun <- cds[,colData(cds)$stim == "G1_IFN_TUN"]
plot_cells(g1ifntun, color_cells_by="cluster.id",
           label_groups_by_cluster=FALSE, cell_size=1,
           group_label_size=4)

g0ifntun <- cds[,colData(cds)$stim == "G0_IFN_TUN"]
plot_cells(g0ifntun, color_cells_by="cluster.id",
           label_groups_by_cluster=FALSE, cell_size=1,
           group_label_size=4)

g0ctrl <- cds[,colData(cds)$stim == "G0_CTRL"]
plot_cells(g0ctrl, color_cells_by="cluster.id",
           label_groups_by_cluster=FALSE, cell_size=1,
           group_label_size=4)

markers = c("PODXL", "NPHS1", "EPCAM", "PTPRO", "LRP2")
plot_cells(cds, genes=markers, label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

g0ctrl <- cds[,colData(cds)$stim == "G0_CTRL"]
plot_cells(g0ctrl, genes=markers, color_cells_by="cluster.id",
           label_groups_by_cluster=FALSE, label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

g1ctrl <- cds[,colData(cds)$stim == "G1_CTRL"]
plot_cells(g1ctrl, genes=markers, color_cells_by="cluster.id",
           label_groups_by_cluster=FALSE, label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

g0ifntun <- cds[,colData(cds)$stim == "G0_IFN_TUN"]
plot_cells(g0ifntun, genes=markers, color_cells_by="cluster.id",
           label_groups_by_cluster=FALSE, label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

g1ifntun <- cds[,colData(cds)$stim == "G1_IFN_TUN"]
plot_cells(g1ifntun, genes=markers, color_cells_by="cluster.id",
           label_groups_by_cluster=FALSE, label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


#############################################################################
## SWNE Plot w/ Seurat v3.0
library(Seurat)
library(swne)

## Extract integrated expression matrix
aligned.counts <- as.matrix(GetAssayData(all.combined, assay = "integrated"))

dim(aligned.counts)

hist(as.matrix(aligned.counts))

## Reconstruct matrix to remove negative values
aligned.counts <- t(apply(aligned.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))

### Run NMF
k <- 20
n.cores <- 16
nmf.res <- RunNMF(aligned.counts, k = k, alpha = 0, init = "ica", n.cores = n.cores)

pc.emb <- Embeddings(se.obj, "pca")
snn <- CalcSNN(t(pc.emb), k = 20, prune.SNN = 0)

swne.embedding <- EmbedSWNE(nmf.res$H, snn, alpha.exp = 1.5, snn.exp = 1, n_pull = 3)

################################################################################
library(Seurat)
library(pheatmap)
library(gplots)
library(dplyr)
library(RColorBrewer)

## Heatmap w/ Integrated Average Expressions
Idents(combined) <- "stim"

combined.averages <- AverageExpression(combined, assays='integrated', return.seurat=TRUE, add.ident="cluster.id")
cad = c("PODXL", "NPHS1", "NPHS2", "MAFB", "PTPRO", "POSTN", "CXCL12", "IGFBP3", "COL3A1")
genes = c("PODXL", "NPHS1", "NPHS2", "MAFB", "PTPRO")
## DoHeatmap(combined.averages, assay="integrated", features=cad,
##           cells=c("G0_IFN_TUN_Glomerular Epithelial Cells", "G1_IFN_TUN_Glomerular Epithelial Cells"))

## G0_CTRL vs G1_CTRL
avgs <- combined.averages@assays$integrated@scale.data
avgs.df <- as.data.frame(avgs)
test <- avgs.df[c('PODXL', 'NPHS1', 'NPHS2', 'MAFB', 'PTPRO', 'CUBN', 'SLC3A1', 'LRP2', 'DPP4', 'TMEM176A', 'WFDC2',
                  'CLDN4', 'TIE1', 'CD31', 'TOP2A', 'CENPF', 'POSTN', 'OSR1', 'SIX2', 'PAX2', 'WT1', 'GAP43', 'STMN2'),
                c('G0_CTRL_Mesenchyme', 'G0_CTRL_Glomerular Epithelial Cells', 'G0_CTRL_Cycling Cells', 'G0_CTRL_Early GEC',
                  'G0_CTRL_Proximal and Distal Tubule', 'G0_CTRL_Endothelial Cells', 'G0_CTRL_Neuron',
                  'G1_CTRL_Mesenchyme', 'G1_CTRL_Early GEC', 'G1_CTRL_Cycling Cells', 'G1_CTRL_Glomerular Epithelial Cells',
                  'G1_CTRL_Proximal and Distal Tubule', 'G1_CTRL_Neuron', 'G1_CTRL_Endothelial Cells')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)

## G1_CTRL vs G1_IFN
test <- avgs.df[c('PODXL', 'NPHS1', 'NPHS2', 'MAFB', 'PTPRO', 'CUBN', 'SLC3A1', 'LRP2', 'DPP4', 'TMEM176A', 'WFDC2',
                  'CLDN4', 'TIE1', 'CD31', 'TOP2A', 'CENPF', 'POSTN', 'OSR1', 'SIX2', 'PAX2', 'WT1', 'GAP43', 'STMN2'),
                c('G1_CTRL_Mesenchyme', 'G1_CTRL_Glomerular Epithelial Cells', 'G1_CTRL_Cycling Cells', 'G1_CTRL_Early GEC',
                  'G1_CTRL_Proximal and Distal Tubule', 'G1_CTRL_Endothelial Cells', 'G1_CTRL_Neuron',
                  'G1_IFN_Mesenchyme', 'G1_IFN_Early GEC', 'G1_IFN_Cycling Cells', 'G1_IFN_Glomerular Epithelial Cells',
                  'G1_IFN_Proximal and Distal Tubule', 'G1_IFN_Neuron', 'G1_IFN_Endothelial Cells')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)

## Heatmap w/ RNA Average Expressions
Idents(combined) <- "stim"

combined.averages <- AverageExpression(combined, assays='RNA', return.seurat=TRUE, add.ident="cluster.id")
avgs <- combined.averages@assays$RNA@scale.data
avgs.df <- as.data.frame(avgs)

## G0_CTRL vs G1_CTRL
test <- avgs.df[c('PODXL', 'NPHS1', 'NPHS2', 'MAFB', 'PTPRO', 'CUBN', 'SLC3A1', 'LRP2', 'DPP4', 'TMEM176A', 'WFDC2',
                  'CLDN4', 'TIE1', 'CD31', 'TOP2A', 'CENPF', 'POSTN', 'OSR1', 'SIX2', 'PAX2', 'WT1', 'GAP43', 'STMN2'),
                c('G0_CTRL_Mesenchyme', 'G0_CTRL_Glomerular Epithelial Cells', 'G0_CTRL_Cycling Cells', 'G0_CTRL_Early GEC',
                  'G0_CTRL_Proximal and Distal Tubule', 'G0_CTRL_Endothelial Cells', 'G0_CTRL_Neuron',
                  'G1_CTRL_Mesenchyme', 'G1_CTRL_Early GEC', 'G1_CTRL_Cycling Cells', 'G1_CTRL_Glomerular Epithelial Cells',
                  'G1_CTRL_Proximal and Distal Tubule', 'G1_CTRL_Neuron', 'G1_CTRL_Endothelial Cells')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)

## G1_CTRL vs G1_IFN
test <- avgs.df[c('PODXL', 'NPHS1', 'NPHS2', 'MAFB', 'PTPRO', 'CUBN', 'SLC3A1', 'LRP2', 'DPP4', 'TMEM176A', 'WFDC2',
                  'CLDN4', 'TIE1', 'CD31', 'TOP2A', 'CENPF', 'POSTN', 'OSR1', 'SIX2', 'PAX2', 'WT1', 'GAP43', 'STMN2'),
                c('G1_CTRL_Mesenchyme', 'G1_CTRL_Glomerular Epithelial Cells', 'G1_CTRL_Cycling Cells', 'G1_CTRL_Early GEC',
                  'G1_CTRL_Proximal and Distal Tubule', 'G1_CTRL_Endothelial Cells', 'G1_CTRL_Neuron',
                  'G1_IFN_Mesenchyme', 'G1_IFN_Early GEC', 'G1_IFN_Cycling Cells', 'G1_IFN_Glomerular Epithelial Cells',
                  'G1_IFN_Proximal and Distal Tubule', 'G1_IFN_Neuron', 'G1_IFN_Endothelial Cells')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)


## DotPlot Function
Idents(combined) <- "stim"
sample <- c("G0 CTRL", "G1 CTRL", "G0 IFN", "G1 IFN", "G0 IFN TUN", "G1 IFN TUN")
names(sample) <- levels(combined)
combined <- RenameIdents(combined, sample)
combined@meta.data$sample <- Idents(combined)

## APOL1 DotPlot
DotPlot(g0_g1_ifn, features='APOL1', dot.scale = 20, cols = c("white", "red"))

## APOL Genes
genes <- c('APOL1', 'APOL2', 'APOL3', 'APOL4', 'APOL5', 'APOL6')
Idents(combined) <- "sample"
g0_g1_ctrl <- SubsetData(combined, ident.use = c("G0 CTRL", "G1 CTRL"))
Idents(g0_g1_ctrl) <- "cluster.id"
DotPlot(g0_g1_ctrl, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

g0_g1_ifn <- SubsetData(combined, ident.use = c("G0 IFN", "G1 IFN"))
Idents(g0_g1_ifn) <- "cluster.id"
DotPlot(g0_g1_ifn, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## UPR Genes
genes <- c('SYVN1', 'ERP44', 'DNAJB9', 'PDIA4', 'PDIA5', 'PDIA6', 'HERPUD2', 'HSPA5', 'DDIT3', 'ERN1', 'PPP1R15A', 'EDEM1',
           'GOSR2', 'GSK3A', 'IGFBP1', 'HSP90B1', 'SEC31A', 'SERP1')
g0_g1_ifn_tun <- SubsetData(combined, ident.use = c("G0 IFN TUN", "G1 IFN TUN"))
Idents(g0_g1_ifn_tun) <- "cluster.id"
upr <- SubsetData(g0_g1_ifn_tun, ident.use = c("Glomerular Epithelial Cells", "Proximal and Distal Tubule", "Endothelial Cells"))
DotPlot(upr, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## Apoptosis Genes
genes <- c('ASM1', 'BAD', 'BAK1', 'BAX', 'BCL2', 'BCL10', 'BIK', 'BIRC8', 'CARD8', 'CARD19', 'CASP1', 'CASP3',
           'CASP9', 'CRADD', 'DIABLO', 'FADD', 'LCN2', 'PUMA', 'SARP2', 'SERPINB9', 'TGFB1', 'TGFB2', 'TNFAIP8', 'TRADD', 'TRAIL')
g0_g1_ifn_tun <- SubsetData(combined, ident.use = c("G0 IFN TUN", "G1 IFN TUN"))
Idents(g0_g1_ifn_tun) <- "cluster.id"
apoptosis <- SubsetData(g0_g1_ifn_tun, ident.use = c("Glomerular Epithelial Cells", "Proximal and Distal Tubule", "Endothelial Cells"))
DotPlot(apoptosis, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## Dediff Genes
genes <- c('PODXL', 'NPHS1', 'NPHS2', 'PTPRO', 'MFAB')
g0_g1_ifn_tun <- SubsetData(combined, ident.use = c("G0 IFN TUN", "G1 IFN TUN"))
Idents(g0_g1_ifn_tun) <- "cluster.id"
apoptosis <- SubsetData(g0_g1_ifn_tun, ident.use = c("Glomerular Epithelial Cells", "Early GEC"))
DotPlot(apoptosis, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## Dediff Genes
genes <- c('PODXL', 'NPHS1', 'NPHS2', 'PTPRO', 'MFAB')
g0_g1_ifn <- SubsetData(combined, ident.use = c("G0 IFN", "G1 IFN"))
Idents(g0_g1_ifn) <- "cluster.id"
apoptosis <- SubsetData(g0_g1_ifn, ident.use = c("Glomerular Epithelial Cells", "Early GEC"))
DotPlot(apoptosis, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## Dediff Genes
genes <- c('SLC3A1', 'EPCAM', 'CUBN', 'LRP2', 'WFDC2', 'CLDN4', 'DPP4', 'TMEM176A')
g0_g1_ifn_tun <- SubsetData(combined, ident.use = c("G0 IFN TUN", "G1 IFN TUN"))
Idents(g0_g1_ifn_tun) <- "cluster.id"
apoptosis <- SubsetData(g0_g1_ifn_tun, ident.use = c("Proximal and Distal Tubule"))
DotPlot(apoptosis, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## Dediff Genes
genes <- c('SLC3A1', 'EPCAM', 'CUBN', 'LRP2', 'WFDC2', 'CLDN4', 'DPP4', 'TMEM176A')
g0_g1_ifn <- SubsetData(combined, ident.use = c("G0 IFN", "G1 IFN"))
Idents(g0_g1_ifn) <- "cluster.id"
apoptosis <- SubsetData(g0_g1_ifn, ident.use = c("Proximal and Distal Tubule"))
DotPlot(apoptosis, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## Dediff Genes
genes <- c('TIE1', 'CD31')
g0_g1_ifn_tun <- SubsetData(combined, ident.use = c("G0 IFN TUN", "G1 IFN TUN"))
Idents(g0_g1_ifn_tun) <- "cluster.id"
apoptosis <- SubsetData(g0_g1_ifn_tun, ident.use = c("Endothelial Cells"))
DotPlot(apoptosis, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## Dediff Genes
genes <- c('TIE1', 'CD31')
g0_g1_ifn <- SubsetData(combined, ident.use = c("G0 IFN", "G1 IFN"))
Idents(g0_g1_ifn) <- "cluster.id"
apoptosis <- SubsetData(g0_g1_ifn, ident.use = c("Endothelial Cells"))
DotPlot(apoptosis, features=genes, dot.scale = 20, split.by = 'sample', cols = c('red', 'red')) + RotatedAxis()

## Color Coded VlnPlot w/ Max Y = 3
colors <- c('#A6E9A9', '#49C64E', '#FFA8AA', '#0000FF', '#FFEE20', '#FF0000', '#A600FF')
subset <- SubsetData(combined, ident.use = c("G0 CTRL", "G1 CTRL", "G0 IFN", "G1 IFN"))
VlnPlot(object = subset, features = c('PODXL', 'NPHS1'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('NPHS2', 'MAFB'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('PTPRO', 'CUBN'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('SLC3A1', 'LRP2'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('DPP4', 'TMEM176A'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('WFDC2', 'CLDN4'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('TIE1', 'CD31'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('TOP2A', 'CENPF'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('POSTN', 'OSR1'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('SIX2', 'PAX2'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('WT1', 'GAP43'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)
VlnPlot(object = subset, features = c('STMN2'), group.by='cluster.id', pt.size=0, y.max=3, cols=colors)

## Heatmap w/ Average Expressions
Idents(combined) <- "stim"

combined.averages <- AverageExpression(combined, assays='integrated', return.seurat=TRUE, add.ident="cluster.id")

## UPR Genes
avgs <- combined.averages@assays$integrated@scale.data
avgs.df <- as.data.frame(avgs)
test <- avgs.df[c('SYVN1', 'ERP44', 'DNAJB9', 'PDIA4', 'PDIA5', 'PDIA6', 'HERPUD2', 'HSPA5', 'DDIT3', 'ERN1', 'PPP1R15A', 'EDEM1',
                  'GOSR2', 'GSK3A', 'IGFBP1', 'HSP90B1', 'SEC31A', 'SERP1'),
                c('G0_IFN_Glomerular Epithelial Cells', 'G1_IFN_Glomerular Epithelial Cells', 'G0_IFN_TUN_Glomerular Epithelial Cells',
                  'G1_IFN_TUN_Glomerular Epithelial Cells')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)

## Apoptosis Genes
avgs <- combined.averages@assays$integrated@scale.data
avgs.df <- as.data.frame(avgs)
test <- avgs.df[c('BAD', 'BAK1', 'BAX', 'BCL2', 'BCL10', 'BIK', 'CASP1', 'CASP3', 'CASP9', 'CRADD', 'DIABLO', 'FADD', 'LCN2',
                  'PUMA', 'SARP2', 'SERPINB9', 'TGFB1', 'TGFB2', 'TNFAIP8', 'TRADD'),
                c('G0_IFN_TUN_Glomerular Epithelial Cells', 'G1_IFN_TUN_Glomerular Epithelial Cells')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)

## Dediff Heatmap
avgs <- combined.averages@assays$integrated@scale.data
avgs.df <- as.data.frame(avgs)
test <- avgs.df[c('PODXL', 'NPHS1', 'MFAB'),
                c('G0_IFN_Glomerular Epithelial Cells', 'G1_IFN_Glomerular Epithelial Cells',
                  'G0_IFN_Early GEC', 'G1_IFN_Early GEC', 'G0_IFN_TUN_Glomerular Epithelial Cells',
                  'G1_IFN_TUN_Glomerular Epithelial Cells', 'G0_IFN_TUN_Early GEC',
                  'G1_IFN_TUN_Early GEC')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)

## FeatureScatter
g0_ifn_tun <- SubsetData(combined, ident.use = c("G0_IFN_TUN"))
Idents(g0_ifn_tun) <- "cluster.id"
cells <- SubsetData(g0_ifn_tun, ident.use = c("Glomerular Epithelial Cells", "Early GEC", "Proximal and Distal Tubule",
                                              "Endothelial Cells"))
DefaultAssay(cells) <- "integrated"
FeatureScatter(cells, feature1="PODXL", feature2="NPHS1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="PODXL", feature2="COL3A1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="NPHS1", feature2="COL3A1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="PODXL", feature2="DDIT3", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="NPHS1", feature2="DDIT3", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="UBD", feature2="APOL1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="PPP1R18", feature2="APOL1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)

g1_ifn_tun <- SubsetData(combined, ident.use = c("G1_IFN_TUN"))
Idents(g1_ifn_tun) <- "cluster.id"
cells <- SubsetData(g1_ifn_tun, ident.use = c("Glomerular Epithelial Cells", "Early GEC", "Proximal and Distal Tubule",
                                              "Endothelial Cells"))
DefaultAssay(cells) <- "integrated"
FeatureScatter(cells, feature1="PODXL", feature2="NPHS1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="PODXL", feature2="COL3A1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="NPHS1", feature2="COL3A1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="PODXL", feature2="DDIT3", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="NPHS1", feature2="DDIT3", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="UBD", feature2="APOL1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)
FeatureScatter(cells, feature1="PPP1R18", feature2="APOL1", group.by="cluster.id", cols=c('#FFA8AA', '#0000FF', '#FC1212', '#A600FF'), pt.size=3)


## Convert To Loop
library(loomR)

DefaultAssay(combined) <- "integrated"
features.combined <- FindVariableFeatures(combined)
features.combined@graphs <- list()
combined.loom <- as.loom(features.combined, assay='integrated', filename='combined.loom', verbose=TRUE)
combined.loom$close_all()

Idents(combined) <- "sample"
subset <- SubsetData(combined, ident.use = c("G0 CTRL", "G1 CTRL", "G0 IFN", "G1 IFN"))
features.subset <- FindVariableFeatures(subset)
subset.loom <- as.loom(features.subset, filename='subset.loom', verbose=TRUE)
subset.loom$close_all()

```

```{python}
## Import Relevant Packages
import scanpy
import scanpy.api as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
from scipy.stats import mode
from collections import Counter

sc.settings.verbosity = 3
sc.set_figure_params(color_map='viridis')
sc.logging.print_versions()

color = ['#A6E9A9', '#FFA5A5', '#FF0000', '#A600FF', '#FFEE20', '#0000FF']

## Read In Data
combined = sc.read_loom('./combined.loom')

## Subset Out Artifact
subset = combined[combined.obs['cluster_id'].isin(['Mesenchyme', 'Cycling Cells', 'Early GEC', 'Proximal and Distal Tubule', 'Neuron',
                                                       'Glomerular Epithelial Cells', 'Endothelial Cells'])]
g0_g1_ctrl_ifn = subset[subset.obs['stim'].isin(['G0_CTRL', 'G1_CTRL', 'G0_IFN', 'G1_IFN'])]

## Reorder Categorical Variables
g0_g1_ctrl_ifn.obs['cluster_id'].cat.reorder_categories(['Glomerular Epithelial Cells', 'Early GEC', 'Proximal and Distal Tubule',
                                                         'Endothelial Cells', 'Cycling Cells', 'Mesenchyme', 'Neuron'], inplace=True)

axes = sc.pl.stacked_violin(g0_g1_ctrl_ifn, ['PODXL', 'NPHS1', 'SLC3A1', 'LRP2', 'DPP4', 'TMEM176A', 'WFDC2', 'CLDN4', 'TIE1',
                                             'TOP2A', 'CENPF', 'POSTN', 'OSR1', 'PAX2', 'GAP43', 'STMN2'],
                            groupby='cluster_id', show=False, row_palette=['#FF0000', '#FFA8AA', '#0000FF', '#A600FF', '#49C64E',
                                                                           '#A6E9A9', '#FFEE20'])
for ax in axes:
    ax.set_ylim(0,4)







```
