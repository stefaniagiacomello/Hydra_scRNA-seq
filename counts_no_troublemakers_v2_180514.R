setwd("Documents/WABI/i_adameyko_1606/tpm/tpm_no_troublemakers_v2/")

library(Seurat)
library(dplyr)

cells_to_filter <- read.delim("../../pca/qc_summary.filtered_cells.txt", sep="\t")  
SS38 <- read.delim("SS_038.rsem_counts.txt", sep = ",")  
SS39 <- read.delim("SS_039.rsem_counts.txt", sep = ",")  
SS40 <- read.delim("SS_040.rsem_counts.txt", sep = ",")

#Merge the count tables 
count <- merge(SS38, SS39)
count_total <- merge(count, SS40)

#Keep only good cells
cells_to_remove <- row.names(cells_to_filter)
count_total_filtered <- count_total[ , !(names(count_total) %in% row.names(cells_to_filter))]
rownames(count_total_filtered) = count_total_filtered[ , 1]
count_total_filtered = count_total_filtered[ , -1]

#Clean-up the Global Environment
remove(count, count_total, SS38, SS39, SS40, cells_to_filter)

#Remove zero read genes
z<-which(rowSums(count_total_filtered)==0)
if (length(z>0)){
  count_total_filtered<-count_total_filtered[-z,]
}

write.table(count_total_filtered, file="counts_total_filtered_180514.txt", col.names = NA, sep="\t", quote = F)
count_total_filtered <- read.table("counts_total_filtered_180514.txt", header = T, row.names = 1, sep="\t")

metafile <- matrix(nrow = 1016, ncol=4)
row.names(metafile) = colnames(count_total_filtered)
metafile[ , 1] = row.names(metafile)
metafile = as.data.frame(metafile)
colnames(metafile) = c("Cell_ID", "Plate", "Legend", "Color")
metafile$Plate <- ifelse(grepl("SS_038", metafile$Cell_ID, ignore.case = T), "SS_038",
                         ifelse(grepl("SS_039", metafile$Cell_ID, ignore.case = T), "SS_039", "SS_040"))
metafile$Legend <- ifelse(grepl("SS_038", metafile$Cell_ID, ignore.case = T), "Neuron - DsRED+",
                          ifelse(grepl("SS_039", metafile$Cell_ID, ignore.case = T), "Stem Cell - GFP+", "DsREDweak_GFPweak"))
metafile$Color <- ifelse(grepl("SS_038", metafile$Cell_ID, ignore.case = T), "#cc2a36",
                         ifelse(grepl("SS_039", metafile$Cell_ID, ignore.case = T), "#98DF8A", "#e39e54"))
write.table(metafile, file="metafile_count_180514.txt", sep="\t", quote = F, row.names = T, col.names = NA)

R <- read.table("counts_total_filtered_180514.txt", sep="\t", row.names=1, header=T)

M <- read.table("metafile_count_180514.txt",sep="\t", header=T, row.names=1)

col_x <- M$Plate
l2cols <- c("red", "green", "orange")[as.integer(factor(col_x, levels = c("SS_038", "SS_039", "SS_040")))]
l2syms <- c(17, 16, 15)[as.integer(factor(col_x, levels = c("SS_038", "SS_039", "SS_040")))]


ercc <- grep("ERCC",rownames(R))
R <- R[-ercc , ]

data <- CreateSeuratObject(raw.data=R, 
                           min.cells = 3, min.genes = 5000, 
                           project = "hydra", meta.data = M)

VlnPlot(object = data, features.plot = c("nGene"), nCol = 1)
VlnPlot(object = data, features.plot = c("nGene", "nUMI"), nCol = 2, group.by = "Plate")
GenePlot(object = data, gene1 = "nUMI", gene2 = "nGene")
data <- SetAllIdent(object = data, id = "Legend")


scale.factor <- mean(colSums(R))
data <- NormalizeData(object = data, normalization.method = "LogNormalize", 
                      scale.factor = scale.factor)

data <- FindVariableGenes(object = data, mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.5, x.high.cutoff = 50, y.cutoff = 0.5, do.plot=T)

data <- ScaleData(object = data, vars.to.regress = c("nGene"), display.progress=T)
#dataB <- ScaleData(object = data, vars.to.regress = c("nGene","Plate"), display.progress=T)

data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:5,  genes.print = 5)
VizPCA(object = data, pcs.use = 1:4)

p1 <- PCAPlot(object = data, dim.1 = 1, dim.2 = 2, do.return=T,group.by="Plate")
#p2 <- PCAPlot(object = dataB, dim.1 = 1, dim.2 = 2, do.return=T,group.by="Plate")

library(ggplot2)
grid.arrange(p1, ncol=2)
plot(p1)


PCAPlot(object = data, dim.1 = 1, dim.2 = 2)
data <- ProjectPCA(object = data, do.print = FALSE)
PCHeatmap(object = data, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

data <- JackStraw(object = data, num.replicate = 100)
pdf(file="JackStraw.pdf", height = 20)
JackStrawPlot(object = data, PCs = 1:20)
dev.off()

PCElbowPlot(object = data)

set.seed(1)
data <- FindClusters(object = data, reduction.type = "pca", dims.use = 1:20, 
                      resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = data)

data <- RunTSNE(object = data, dims.use = 1:20, do.fast = TRUE)
pdf(file="tSNE_Seraut_cl_180514.pdf")
TSNEPlot(object = data,do.label = T, pt.size = 1)
dev.off()

pdf(file="tSNE_Seraut_plate_180514.pdf")
TSNEPlot(object = data, group.by = "Plate")
dev.off()

saveRDS(data, file = "count_seraut_180514.rds")
data <- readRDS("count_seraut_180514.rds")

cluster5.markers <- FindMarkers(object = data, ident.1 = 5, min.pct = 0.25)
print(head(cluster5.markers),5)
VlnPlot(object = data, features.plot = rownames(cluster5.markers)[1:6],nCol=3,size.title.use=10)
FeaturePlot(object = data, features.plot = rownames(cluster5.markers)[1:6], cols.use = c("yellow", "red","black"), reduction.use = "tsne")

cluster11.markers <- FindMarkers(object = data, ident.1 = 11, min.pct = 0.25)
print(head(cluster11.markers),5)
VlnPlot(object = data, features.plot = rownames(cluster11.markers)[1:6],nCol=3,size.title.use=10)
FeaturePlot(object = data, features.plot = rownames(cluster11.markers)[1:6], cols.use = c("yellow", "red","black"), reduction.use = "tsne")

FeaturePlot(object = data, features.plot = c("cluster102707", "cluster43207"), cols.use = c("green", "blue")) 
                                       

cluster0.markers <- FindMarkers(object = data, ident.1 = 0, min.pct = 0.25)
print(x= head(x=cluster0.markers,n=10))
saveRDS(cluster0.markers, file="cluster0.markers_180521.rds")
cluster0.markers <- readRDS(file="cluster0.markers_180521.rds")
write.table(data.frame(cluster0.markers), file="cluster0_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster1.markers <- FindMarkers(object = data, ident.1 = 1, min.pct = 0.25)
print(x= head(x=cluster1.markers,n=10))   
saveRDS(cluster1.markers, file="cluster1.markers_180521.rds")
cluster1.markers <- readRDS(file="cluster1.markers_180521.rds")
write.table(data.frame(cluster1.markers), file="cluster1_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster2.markers <- FindMarkers(object = data, ident.1 = 2, min.pct = 0.25)
print(x= head(x=cluster2.markers,n=10))  
saveRDS(cluster2.markers, file="cluster2.markers_180514.rds")
cluster2.markers <- readRDS(file="cluster2.markers_180514.rds")
write.table(data.frame(cluster2.markers), file="cluster2_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster3.markers <- FindMarkers(object = data, ident.1 = 3, min.pct = 0.25)
print(x= head(x=cluster3.markers,n=10)) 
saveRDS(cluster3.markers, file="cluster3.markers_180514.rds")
cluster3.markers <- readRDS(file="cluster3.markers_180514.rds")
write.table(data.frame(cluster3.markers), file="cluster3_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster4.markers <- FindMarkers(object = data, ident.1 = 4, min.pct = 0.25)
print(x= head(x=cluster4.markers,n=10))  
saveRDS(cluster4.markers, file="cluster4.markers_180514.rds")
cluster4.markers <- readRDS(file="cluster4.markers_180514.rds")
write.table(data.frame(cluster4.markers), file="cluster4_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster5.markers <- FindMarkers(object = data, ident.1 = 5, min.pct = 0.25)
print(x= head(x=cluster5.markers,n=10))
saveRDS(cluster5.markers, file="cluster5.markers_180514.rds")
cluster5.markers <- readRDS(file="cluster5.markers_180514.rds")
write.table(data.frame(cluster5.markers), file="cluster5_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster6.markers <- FindMarkers(object = data, ident.1 = 6, min.pct = 0.25)
print(x= head(x=cluster6.markers,n=10))  
saveRDS(cluster6.markers, file="cluster6.markers_180514.rds")
cluster6.markers <- readRDS(file="cluster6.markers_180514.rds")
write.table(data.frame(cluster6.markers), file="cluster6_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster7.markers <- FindMarkers(object = data, ident.1 = 7, min.pct = 0.25)
print(x= head(x=cluster7.markers,n=10))  
saveRDS(cluster7.markers, file="cluster7.markers_180514.rds")
cluster7.markers <- readRDS(file="cluster7.markers_180514.rds")
write.table(data.frame(cluster7.markers), file="cluster7_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster8.markers <- FindMarkers(object = data, ident.1 = 8, min.pct = 0.25)
print(x= head(x=cluster8.markers,n=100))  
saveRDS(cluster8.markers, file="cluster8.markers_180514.rds")
cluster8.markers <- readRDS(file="cluster8.markers_180514.rds")
write.table(data.frame(cluster8.markers), file="cluster8_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster9.markers <- FindMarkers(object = data, ident.1 = 9, min.pct = 0.25)
print(x= head(x=cluster9.markers,n=10))  
saveRDS(cluster9.markers, file="cluster9.markers_180514.rds")
cluster9.markers <- readRDS(file="cluster9.markers_180514.rds")
write.table(data.frame(cluster9.markers), file="cluster9_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster10.markers <- FindMarkers(object = data, ident.1 = 10, min.pct = 0.25)
print(x= head(x=cluster10.markers,n=10))  
saveRDS(cluster10.markers, file="cluster10.markers_180514.rds")
cluster10.markers <- readRDS(file="cluster10.markers_180514.rds")
write.table(data.frame(cluster10.markers), file="cluster10_markers_180521.txt", sep="\t", quote = F, col.names = NA)

cluster11.markers <- FindMarkers(object = data, ident.1 = 11, min.pct = 0.25)
print(x= head(x=cluster11.markers,n=10))  
saveRDS(cluster11.markers, file="cluster11.markers_180514.rds")
cluster11.markers <- readRDS(file="cluster11.markers_180514.rds")
write.table(data.frame(cluster11.markers), file="cluster11_markers_180521.txt", sep="\t", quote = F, col.names = NA)

pdf(file="top_gene_per_cl_180514.pdf")
FeaturePlot(object = data, features.plot = c("cluster74791", "cluster59896", "cluster40223", "cluster219771", 
                                             "cluster61768", "cluster52149", "cluster103279", "cluster85839", "cluster172891", "cluster116435", "cluster103058", "cluster180452"), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()

#Make tSNE plots with genes requested by Alex: these genes are a subgroup of the ones analyzed in the plots below. 
pdf(file="tSNE_genes_specified_by_Alex_180913.pdf")
FeaturePlot(object = data, features.plot = c("cluster115188", "cluster103279", "cluster185625", "cluster150163", 
                                             "cluster131995", "cluster60981", "cluster62692", "cluster209661"), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()

pdf(file="tSNE_genes_specified_by_Alex_181030_1.pdf")
FeaturePlot(object = data, features.plot = c("cluster115188", "cluster103279", "cluster185625", "cluster150163", "cluster131995", "cluster60981"),
                                             cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()

pdf(file="tSNE_genes_specified_by_Alex_181030_2.pdf")
FeaturePlot(object = data, features.plot = c("cluster62692", "cluster209661", "cluster21977", "cluster8358", "cluster35264", "cluster20828"),
                                            cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()

pdf(file="tSNE_genes_specified_by_Alex_181030_3.pdf")
FeaturePlot(object = data, features.plot = c("cluster40223", "cluster46777", "cluster95001", "cluster195321", "cluster74791", "cluster20477"), 
            cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()

pdf(file="tSNE_genes_specified_by_Alex_181030_4.pdf")
FeaturePlot(object = data, features.plot = c("cluster65816", "cluster59896", "cluster103058", "cluster172891", "cluster180452", "cluster43015"), 
            cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()


pdf(file="tSNE_mouse_gene_190513.pdf")
FeaturePlot(object = data, features.plot = c("cluster8724"), 
            cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()



cluster_neuron.markers <- FindMarkers(object = data, ident.1 = c(5,6,7,8,9,10,11), ident.2 = c(0, 1, 2, 3, 4), 
                                min.pct = 0.25)
print(x= head(x=cluster_neuron.markers,n=10))  
saveRDS(cluster_neuron.markers, file="clusters_neurons.markers_180611.rds")
cluster_neuron.markers <- readRDS(file="clusters_neurons.markers_180611.rds")
write.table(data.frame(cluster_neuron.markers), file="cluster_neuron.markers_vs_rest_180611.txt", sep="\t", quote = F, col.names = NA)

cluster_stem.markers <- FindMarkers(object = data, ident.1 = c(1,4), ident.2 = c(0, 2, 3, 5,6,7,8,9,10,11), 
                                      min.pct = 0.25)
print(x= head(x=cluster_stem.markers,n=10))  
saveRDS(cluster_stem.markers, file="clusters_stem.markers_180611.rds")
cluster_stem.markers <- readRDS(file="clusters_stem.markers_180611.rds")
write.table(data.frame(cluster_stem.markers), file="cluster_stem.markers_vs_rest_180611.txt", sep="\t", quote = F, col.names = NA)

#Save what cells belong to which tSNE cluster
seurat_cl <- as.data.frame(data@ident)
colnames(seurat_cl)=c("Cluster")
write.table(seurat_cl, file="seurat_clusters_identity_180514.txt", col.names = NA, sep = "\t")

data.markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)

FeaturePlot(object = data, features.plot = c("cluster105605"), cols.use = c("grey", "blue"), reduction.use = "tsne")

top20 <- data.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

#Save some data
data_matrix <- data@scale.data
data_matrix=as.data.frame(data_matrix)
write.table(data_matrix, file="../../tpm/tpm_no_troublemakers_v2/seurat_scaled_matrix.txt", col.names = NA, quote = F, sep="\t")

variable_genes_whole_transcriptome <- data@var.genes
write.table(variable_genes_whole_transcriptome, file="list_var_genes_Seurat_whole_trans_180514.txt")

cells_used_for_tSNE_whole_transcriptome <- data@data@Dimnames

pca_seurat <- data@dr$pca@cell.embeddings
saveRDS(pca_seurat, file="pca_seurat_180517.rds")

tSNE_coordinates_whole_transcriptome <- data@dr$tsne@cell.embeddings
write.table(tSNE_coordinates_whole_transcriptome, file="list_tSNE_coords_Seurat_whole_trans_180514.txt")

seurat_cells_genes <- as.data.frame(data@data@Dimnames)

#Import the data file created
data <- readRDS("count_seraut_180514.rds")

metafile <- read.table("metafile_count_180514.txt", sep = "\t", row.names = 1, header = T)
metafile=metafile[rownames(pca_seurat), ]
write.table(metafile, file="metafile_seurat_180517.txt", quote = F, col.names = NA, sep="\t")
metafile <- read.table("metafile_seurat_180517.txt", sep = "\t", row.names = 1, header = T)

col_x <- metafile$Plate
l2cols <- c("#cc2a36", "#98DF8A", "#e39e54")[as.integer(factor(col_x, levels = c("SS_038", "SS_039", "SS_040")))]
l2syms <- c(17, 16, 15)[as.integer(factor(col_x, levels = c("SS_038", "SS_039", "SS_040")))]
pdf(file="seurat_PCA_tpm_no_troublemakers_filtered_180517.pdf")
plot(pca_seurat, type="p", pch=l2syms, col=l2cols, xlab='Comp.1', ylab='Comp.2')
legend("topleft", legend=c("Neuron - DsRED+", "Stem Cell - GFP+", "DsREDweak_GFPweak"), col=c("#cc2a36", "#98DF8A", "#e39e54"), pch=c(17,16,15))
dev.off()

#Calculate genes per cell in Seurat filtered object
#Create filtered Seurat matrix:
R <- read.table("counts_total_filtered_180514.txt", sep="\t", row.names=1, header=T)
seurat_cells_genes <- as.data.frame(data@data@Dimnames)
seurat_genes <- as.data.frame(seurat_cells_genes[1])
rownames(seurat_genes)=seurat_genes$c..cluster0....cluster1....cluster10....cluster100....cluster1000...
seurat_cells <- seurat_cells_genes[2]
R_filtered <- R[rownames(seurat_genes) , rownames(pca_seurat) ]

#Number genes expressed on Seurat filtered matrix 
R_filtered_genes_per_cell <- as.data.frame(apply(R_filtered, 2, function(R_filtered) sum(R_filtered >0))) 
colnames(R_filtered_genes_per_cell) = "genes_per_cell"
max(R_filtered_genes_per_cell)
# 37287
min(R_filtered_genes_per_cell)
# 4967
maxColorValue <- 37287
palette <- colorRampPalette(c("yellow","red", "black"))(maxColorValue)
#Color legend
pdf(file="color_legend_genes_180517.pdf")
col_genes = palette[cut(R_filtered_genes_per_cell$genes_per_cell, maxColorValue)]
z=matrix(1:100,nrow=1)
x=1
y=seq(0,37287,len=100)
image(x,y,z,col=palette,axes=FALSE,xlab="",ylab="", aspect="iso")
axis(2)
dev.off()

pdf(file="seurat_PCA_genespercell_180517.pdf")
plot(pca_seurat, col=palette[cut(R_filtered_genes_per_cell$genes_per_cell, maxColorValue)], cex=1, xlab="", pch=l2syms, ylab="", main ="Number of genes per cell")
legend("topleft", legend=c("Neuron - DsRED+", "Stem Cell - GFP+", "DsREDweak_GFPweak"), bty = 'n', col=c("grey", "grey", "grey"), pch=c(17,16,15))
dev.off()

#tSNE
pdf(file="seurat_tSNE_180517.pdf")
plot(tSNE_coordinates_whole_transcriptome,pch=l2syms,col=l2cols,main="tSNE colored by plate - all genes", xlab = "tSNE1", ylab="tSNE2")
legend("topright", legend=c("Neuron - DsRED+", "Stem Cell - GFP+", "DsREDweak_GFPweak"), bty = 'n', col=c("#cc2a36", "#98DF8A", "#e39e54"), pch=c(17,16,15))
dev.off()

pdf(file="seurat_tSNE_genespercell_180517.pdf")
plot(tSNE_coordinates_whole_transcriptome,pch=l2syms,col=palette[cut(R_filtered_genes_per_cell$genes_per_cell, maxColorValue)],main="tSNE colored by number of genes per cell", xlab = "tSNE1", ylab="tSNE2")
legend("topright", legend=c("Neuron - DsRED+", "Stem Cell - GFP+", "DsREDweak_GFPweak"), bty = 'n', col=c("#cc2a36", "#98DF8A", "#e39e54"), pch=c(17,16,15))
dev.off()

l2cols_seurat <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F")[as.integer(factor(seurat_cl$Cluster, levels = c("0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")))]

pdf(file="seurat_tSNE_cluster_180517.pdf")
plot(tSNE_coordinates_whole_transcriptome,pch=l2syms,col=l2cols_seurat,main="tSNE colored by cluster", xlab = "tSNE1", ylab="tSNE2")
legend("topright", legend=c("Neuron - DsRED+", "Stem Cell - GFP+", "DsREDweak_GFPweak"), bty = 'n', col=c("black", "black", "black"), pch=c(17,16,15))
dev.off()

#Violin plots of selected genes by Alex
pdf(file="violin_plots_1.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster115188", "cluster103279", "cluster185625", "cluster150163"))
dev.off()

pdf(file="violin_plots_2.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster131995", "cluster60981", "cluster62692", "cluster209661"))
dev.off()

pdf(file="violin_plots_3.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster21977", "cluster8358", "cluster35264", "cluster20828"))
dev.off()

pdf(file="violin_plots_4.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster40223", "cluster46777", "cluster95001", "cluster195321"))
dev.off()

pdf(file="violin_plots_5.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster74791", "cluster20477", "cluster65816", "cluster59896"))
dev.off()

pdf(file="violin_plots_6.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster103058", "cluster172891", "cluster180452", "cluster43015"))
dev.off()

pdf(file="violin_plots_7.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster41630", "cluster68091", "cluster40725", "cluster76720"))
dev.off()

pdf(file="violin_plots_8.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster176765", "cluster119726", "cluster183995"))
dev.off()

pdf(file="violin_plots_9.pdf")
VlnPlot(object = data, nCol = 2, cols.use=c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#98DF8A", "#FF9896", "#9467BD", "#C5B0D5", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F"), features.plot = c("cluster8724", "cluster145885", "cluster150163", "cluster202122"))
dev.off()


###MAke count table for Alex with cell names with cluster identity accordingly to his nomenclature######
clusters <- data.frame(data@ident)
library(plyr)
clusters$cluster_name <- revalue(clusters$data.ident, c("0"="Nb", "1"="SC1", "2"="GC", "3"="Nc", "4"="SC2", "5"="N1", "6"="N2", "7"="N3", "8"="N4", "9"="N5", "10"="N6", "11"="N7"))
clusters$cell_cluster <- paste(rownames(clusters), clusters$cluster_name, sep="_")
raw_counts <- data@raw.data
colnames(raw_counts_2) <- clusters$cell_cluster
write.table(raw_counts_2, file="raw_count_table_seurat_filter_cluster_ID_191009_SG.txt", quote = F, col.names = NA, sep="\t")

