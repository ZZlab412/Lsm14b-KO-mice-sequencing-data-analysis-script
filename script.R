###LSM14B controls oocyte mRNA storage and stability to ensure female fertility

########  Single-Library Analysis with CellRanger in server
########  CellRanger version=3.1.0

nohup /home/.../10xgenomics/cellranger-3.1.0/cellranger count --id=SampleName \
--localcores 40 \
--transcriptome= .../refdata-cellranger-mm10 \
--fastqs=/home/.../Singlecell_Rawdata \
--sample=SampleName &

#############################################################################################
         #########///////=========Seurat version=4.0.2========\\\\\\\#############
#############################################################################################
### Detail information see online vignettes of Seurat
### https://satijalab.org/seurat/vignettes.html

library(dplyr)
library(Seurat)
library(patchwork)

###Load Data
Seurat_Object <- Read10X(data.dir = ".../filtered_feature_bc_matrix")
Seurat_Object <- CreateSeuratObject(counts = Seurat_Object, project = "Seurat_Object", min.cells = 3, min.features = 200)

###Value setting of filteration
Seurat_Object <- subset(Seurat_Object, subset =  nFeature_RNA > 1500 & nFeature_RNA < 7500 & percent.mt < 25)

###Perform integration
Object.anchors <- FindIntegrationAnchors(object.list = list(Object1,Object2), dims = 1:15)
Object.integrated <- IntegrateData(anchorset = Object.anchors, dims = 1:15)
DefaultAssay(object = Object.integrated) <- "integrated"
Object.integrated <- ScaleData(Object.integrated, verbose = FALSE)
Object.integrated <- RunPCA(Object.integrated, npcs = 30, verbose = FALSE)
Object.integrated <- RunUMAP(Object.integrated, reduction = "pca", dims = 1:15)

#### UMAP and Clustering
Object.integrated <- FindNeighbors(Object.integrated, reduction = "pca", dims = 1:20)
Object.integrated <- FindClusters(Object.integrated, resolution = 0.5)
p1 <- DimPlot(Object.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Object.integrated, reduction = "umap", group.by = "seurat_clusters")

###find markers
DefaultAssay(Object.integrated) <- "RNA"
Object.integrated.All.Markers <- FindAllMarkers(Object.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Object.integrated<-BuildClusterTree(Object.integrated)
Tool(object = Object.integrated, slot = 'BuildClusterTree')
PlotClusterTree(Object.integrated)

###subset single cell type
cell_type <- subset(Object.integrated, idents = c("cell_type"))

###Analysis of differentially expressed genes in cell clusters
Differentially.expressed.genes <- FindMarkers(cell_type, ident.1 = "feature_1", ident.2 = "feature_2", verbose = FALSE,logfc.threshold = 0.25)


#############################################################################################
         #########///////=========Monocle version=2.18.0========\\\\\\\#############
#############################################################################################
### Detail information see online vignettes of Monocle
### http://cole-trapnell-lab.github.io/monocle-release/docs/
library(Seurat)
library(monocle)

Subset_Seurat <- readRDS(file = ".../Seurat.rds")

###Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(Subset_Seurat@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Subset_Seurat@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle<- newCellDataSet(data,
                           phenoData = pd,
                           featureData = fd,
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
						   
Subset_Seurat <- FindVariableFeatures(Subset_Seurat, selection.method = "mean.var.plot", nfeatures = 2000)
VariableFeatures(Subset_Seurat)
seurat_var_genes = VariableFeatures(Subset_Seurat)
head(seurat_var_genes)
monocle <- estimateSizeFactors(monocle)
monocle <- estimateDispersions(monocle)
seurat_var = setOrderingFilter(monocle, seurat_var_genes)
plot_ordering_genes(seurat_var)
cds <- reduceDimension(seurat_var, method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds)
plot_cell_trajectory(cds, color_by = "seurat_clusters")


#############################################################################################
         #########///////=========DESeq2 version=1.30.1========\\\\\\\#############
#############################################################################################
library(DESeq2)
Data <- read.table("./Data.txt")

DES_data<- DESeqDataSetFromMatrix(countData = Data,
                                 colData = sampleTable,
                                 design = ~Group)               
DES_data <- DESeq(DES_data)
results <- results(DES_data, pAdjustMethod = "fdr", alpha = 0.05)


#############################################################################################
         #########///////=========clusterProfiler version=3.18.1========\\\\\\\#############
#############################################################################################
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)

data$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(data),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")

kegg_enrich <- enrichKEGG(gene = names(data),
                          organism = 'mouse',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)
				 
go_enrich <- enrichGO(gene = names(data),
                      OrgDb = 'org.Mm.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)		
