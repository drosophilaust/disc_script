### This script provide the workflow that how we read 10x output data and perform substream analysis.
### 

### This workflow used Seurat 2.3.4. You can download this version by the following commands:
# install.packages('devtools')
# devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
# You might need to restart R after this step.

library(Seurat)
library(dplyr)

# load 10x data and create seurat project
# frt82b is the wild type sample
frt82b.data<-Read10X(data.dir="~/path/to/10x/output/outs/filtered_gene_bc_matrices/drosophilagenome/")
frt82b <- CreateSeuratObject(raw.data = frt82b.data, min.cells = 10,project = "MAPS")

# quality check
GenePlot(object = frt82b, gene1 = "nUMI", gene2 = "nGene",pch.use = 16,cex.use = 0.5,col.use = "red")

#
frt82b <- RenameCells(frt82b,add.cell.id = "FRT82B")

# filter cells based on numbers of gene expression
frt82b <- FilterCells(object = frt82b, subset.names = c("nGene"), 
                      low.thresholds = 1000, high.thresholds = 4000)

# Normalizing the data
frt82b <- NormalizeData(object = frt82b)

frt82b <- FindVariableGenes(object = frt82b, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.05, x.high.cutoff = 3, y.cutoff = 0.5)

frt82b <- ScaleData(object = frt82b, vars.to.regress = c("nUMI"),do.par = T,num.cores = 8)


#Perform linear dimensional reduction
frt82b <- RunPCA(object = frt82b, pc.genes = frt82b@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5)

# heatmap for each Principle Component
PCHeatmap(object = frt82b, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)


# Run Non-linear dimensional reduction (tSNE & UMAP)
frt82b <- RunTSNE(object = frt82b, dims.use = 1:20, do.fast = TRUE)

frt82b <- RunUMAP(object = frt82b,genes.use = frt82b@var.genes)

# Find sub clusters
frt82b <- FindClusters(object = frt82b, reduction.type = "pca", dims.use = 1:20, 
                       resolution = 0.3, print.output = 0, save.SNN = TRUE,force.recalc = T)

# plot tSNE
TSNEPlot(object = frt82b, do.label = T)

# plot UMAP
DimPlot(object = frt82b,reduction.use = "umap")


# violin plot of expression in different sub groups
VlnPlot(object = frt82b, features.plot = c("nub","twi","Sox15","hth","tsh","Ubx","odd","dpp","wg","ci","ap","en","hh"))


# plot gene expression based on umap
FeaturePlot(object = frt82b, features.plot = c("nub","twi","Sox15","hth","tsh","Ubx","odd","dpp","wg","ci","ap","en","hh"),
            cols.use = c("#CCCCFF","red"),reduction.use = "umap")



