# Alpha_cell_SMOC1
github code repository for Human Pancreatic α-Cell Heterogeneity and Trajectory Inference Analyses Reveal SMOC1 as a β-Cell Dedifferentiation Gene. I will be as precise as possible, but if you still have questions, please reach out to me through issue or message. 

# Preliminary Works
We have applied several QC parameters, and loop them in this code. We began with SoupX"ed" file and DoubletFinder with exactly same parameters. Please refer to the respository under scRNA-snRNA Marker for more detailed code for SuopX. 

```R
library(ROCR)
library(Matrix)
library(KernSmooth)
library(ROCR)
library(parallel)
library(AnnotationHub)
library(ensembldb)
library(DoubletFinder)


in_data_dir <- ("/your/directory/here")
samples <- dir(in_data_dir)

#################################################################################
#################################################################################
## Seurat Object Creation Again



seurat_list <- lapply(samples, function(sample){
  cur_data <- Read10X(paste0(in_data_dir,"/",sample))
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
  )
  cur_seurat$SampleID <- sample
  
  cur_seurat$log10GenesPerUMI <- log10(cur_seurat$nFeature_RNA) / log10(cur_seurat$nCount_RNA)
  cur_seurat$mitoRatio <- PercentageFeatureSet(object = cur_seurat, pattern = "^MT-")
  cur_seurat$mitoRatio <- cur_seurat@meta.data$mitoRatio / 100
  cur_seurat@meta.data -> metadata
  
  filtered_seurat <- subset(x = cur_seurat, 
                            subset= (nCount_RNA >= 500) & 
                              (nFeature_RNA >= 250) & 
                              (log10GenesPerUMI > 0.80) & 
                              (mitoRatio < 0.20))
  
  
  counts <- GetAssayData(object = filtered_seurat, slot = "counts")
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  filtered_counts <- counts[keep_genes, ]
  
  # FIX: Align metadata to filtered cells
  filtered_meta <- filtered_seurat@meta.data[colnames(filtered_counts), , drop = FALSE]
  filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_meta)
  
  # Continue downstream
  filtered_seurat <- SCTransform(filtered_seurat)
  filtered_seurat <- FindVariableFeatures(filtered_seurat)
  filtered_seurat <- ScaleData(filtered_seurat)
  filtered_seurat <- RunPCA(filtered_seurat)
  filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:10)
  filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.5)
  filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:10)
  
  # DoubletFinder pipeline
  sweep.res.list_kidney <- paramSweep(filtered_seurat, PCs = 1:10, sct = TRUE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  
  homotypic.prop <- modelHomotypic(filtered_seurat$seurat_clusters)
  nExp_poi <- round(0.10 * ncol(filtered_seurat))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  filtered_seurat <- doubletFinder(
    filtered_seurat,
    PCs = 1:10,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp_poi.adj,
    reuse.pANN = NULL,
    sct = TRUE
  )
  
  
  filtered_seurat@meta.data -> metaclean_maxst
  
  
  names(metaclean_maxst)[grepl("DF.classifications", names(metaclean_maxst))] <- "Doubletfinder"
  print(paste("Finished DoubletFinder for ",sample))
  filtered_seurat@meta.data <- metaclean_maxst
  
  return(filtered_seurat)
})
```


# Merging and Gene & Cell Level Macro Filtering
```
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
merged_seurat@meta.data -> metadata





filtered_maxst <- subset(x = merged_seurat, 
                          subset= (nCount_RNA >= 500) & 
                            (nFeature_RNA >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))


filtered_maxst <- JoinLayers(filtered_maxst)
counts <- GetAssayData(object = filtered_maxst, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_maxst <- CreateSeuratObject(filtered_counts, meta.data = filtered_maxst@meta.data)
filtered_maxst@meta.data -> metaclean


```

# Cell Phase scoring and Integration 
```

#cell phase scoring
seurat_phase <- NormalizeData(filtered_maxst)
# Load cell cycle markers

library(AnnotationHub)
library(ensembldb)


cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)


# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)


# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")



# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
# Identify the most variable genes
DefaultAssay(seurat_phase) <-"RNA"
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
# Scale the counts
seurat_phase <- NormalizeData(seurat_phase)
seurat_phase <- ScaleData(seurat_phase)
# Perform PCA
seurat_phase <- RunPCA(seurat_phase)
# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca", group.by="Phase")
options(future.globals.maxSize = 256000 * 1024^2)
split_seurat <- SplitObject(seurat_phase, split.by = "SampleID")


for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

#Select Anchors & Integration Features 
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 2000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
```
# Defining Cell type 
Cell type identification, DEGs, pathway enridchment were done exactly same method from our previous publication and share the code repository from previous publication (scRNA-snRNA-Marker) 

# Extract the Seurat based Coordinate and Cell Types for Embedding 
For the consistency, we have extracted the Seurat defined UMAP coordinate and cell type identity. I have utilized code from Sam Morabito. [Visit Original Code](https://smorabit.github.io/blog/2021/velocyto/)

```

# save metadata table:
seurat_integrated$barcode <- colnames(seurat_integrated)
seurat_integrated$UMAP_1 <- seurat_integrated@reductions$umap@cell.embeddings[,1]
seurat_integrated$UMAP_2 <- seurat_integrated@reductions$umap@cell.embeddings[,2]
write.csv(seurat_integrated@meta.data, file='./metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_integrated, assay='RNA', slot='counts')
writeMM(counts_matrix, file='./counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_integrated@reductions$pca@cell.embeddings, file='./pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='./gene_names.csv',
  quote=F,row.names=F,col.names=F)

```

# RNA Velocity & Trajectory Analysis 
I simply utilized the Document frmo scVelo. You can access it from this link [scVelo Docs](https://scvelo.readthedocs.io/en/stable/getting_started.html)
```python
# import libraries
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

#Load & Save

X = io.mmread("/media/randy/playground/RNA/alpha_ADA/NC_review_2nd_round/Cell Only/counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("/media/randy/playground/RNA/alpha_ADA/NC_review_2nd_round/Cell Only/metadata.csv")

# load gene names:
with open('/media/randy/playground/RNA/alpha_ADA/NC_review_2nd_round/Cell Only/gene_names.csv','r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("/media/randy/playground/RNA/alpha_ADA/NC_review_2nd_round/Cell Only/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['celltype'], frameon=False, save=True)

# save dataset as anndata format
adata.write('seurat_integrated.h5ad')

```
# Then Process the Loom Data 
```python

import anndata as ad
import os

# Set verbosity and plot style
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

# Define base directory and sample IDs
base_dir = "/your/folder/containing/samples"
sample_ids = ["B2C", "B3C", "B5C", "B2N", "B3N", "B5N"]

# Load, clean barcodes, and make var names unique
ldata_list = []
for idx, sample in enumerate(sample_ids, start=1):
    loom_path = os.path.join(base_dir, sample, "velocyto", f"{sample}.loom")
    ldata = ad.read_loom(loom_path)

    # Clean barcodes (remove prefix, replace -1 with _N)
    ldata.obs.index = [
        bc.split(':')[1].replace('-1', f'_{idx}') for bc in ldata.obs.index
    ]
    
    # Make var names unique to avoid column collision
    ldata.var_names_make_unique()
    
    ldata_list.append(ldata)

# Concatenate all loom files
ldata = ldata_list[0].concatenate(ldata_list[1:], batch_key=None)

for idx in range(1, len(sample_ids)+1):
    ldata.obs.index = ldata.obs.index.str.replace(f"_{idx}-{idx-1}", f"-1_{idx}", regex=False)


```
# Cluster Extraction and Velocity Calculation 
For more focused analysis between alpha and beta cell subclusters, we have extracted clusters of interest and performed the RNA velocity and PAGA analysis accordingly. 

```python

# merge matrices into the original adata object. UnboundLocalError may happnes, but you can move on
adata = scv.utils.merge(adata, ldata)
adata.write('seurat_loom_added.h5ad')


adata_albe = adata[adata.obs["celltype"].isin(["alpha_1","alpha_2","alpha_3","alpha_4","AB","beta_1","beta_2","beta_3"])]


adata_cell = adata_albe[adata_albe.obs["seq_type"].isin(["scRNA"])]
adata_nuc = adata_albe[adata_albe.obs["seq_type"].isin(["snRNA"])]

# pre-process (neighbor recalculated, maybe cause problem) 
sc.pp.neighbors(adata_albe, n_neighbors=30, n_pcs=50)  # Compute neighbors without the unsupported argument
scv.pp.filter_and_normalize(adata_albe)
scv.pp.moments(adata_albe, n_neighbors=30, n_pcs=50)
scv.tl.velocity(adata_albe, mode="stochastic")
scv.tl.velocity_graph(adata_albe)
scv.tl.rank_velocity_genes(adata_albe, groupby='celltype', min_corr=.3)

df = pd.DataFrame(adata_albe.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv("velocity_top_alpha_integ_new.csv")
pd.options.display.html.use_mathjax = False 

## Velocity Confidence
scv.tl.velocity_confidence(adata_albe)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata_albe, c=keys, cmap='coolwarm', perc=[5, 95])

scv.pl.velocity_graph(adata_albe, threshold=.1, color='celltype')

x, y = scv.utils.get_cell_transitions(adata_albe, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata_albe, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata_albe, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

scv.tl.velocity_pseudotime(adata_albe)
scv.pl.scatter(adata_albe, color='velocity_pseudotime', cmap='gnuplot')

########################## PAGA
# this is needed due to a current bug - bugfix is coming soon.
adata_albe.uns['neighbors']['distances'] = adata_albe.obsp['distances']
adata_albe.uns['neighbors']['connectivities'] = adata_albe.obsp['connectivities']
scv.tl.paga(adata_albe, groups='celltype')
df = scv.get_df(adata_albe, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata_albe, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, title="scRNA/snRNA Integrated Stochastic")


```
# Trajectory-Pseudotime Gene Expression 
To get the genes of the high correlation in the trajectory, we refer to the PAGA from above and obtained the genes with high trajectory correlation 
```

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)


ab_beta3 <- subset(albe_v4, ident=c("AB","beta_3","beta_2","beta_1"))



gene_annotation <- as.data.frame(rownames(ab_beta3@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(ab_beta3@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(ab_beta3@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = ab_beta3@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"


New_matrix <- ab_beta3@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(ab_beta3@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

library(monocle3)

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition



list_cluster <- ab_beta3@active.ident
names(list_cluster) <- ab_beta3@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-ab_beta3@reductions[["umap"]]@cell.embeddings
cds_from_seurat@reduce_dim_aux$gene_loadings <- ab_beta3@reductions[["pca"]]@feature.loadings

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4)

# Choose the starting point of the cell. If this doesn't work set options(shiny.useragg = FALSE)
cds_from_seurat <- order_cells(cds_from_seurat, reduction_method="UMAP")

# Plot
plot_cells( cds = cds_from_seurat,color_cells_by = "pseudotime",show_trajectory_graph = FALSE) + 
  ggtitle("AB-beta_1 Transition", subtitle="Base = PAGA Origin")




# only choose cell works 
pr_graph_test_res <- graph_test(cds_from_seurat, neighbor_graph="principal_graph", cores=32, verbose = FALSE) # Group relative to pseudotime

hits      <- as.data.frame(pr_graph_test_res[pr_graph_test_res$gene_short_name,])
hits$pass <- hits$morans_I > 0.1 & hits$q_value < 0.001 



hits_most <- hits

hits_most <- hits_most %>% dplyr::filter(pass=="TRUE")
hits_most <- hits_most %>% top_n(n=100, morans_I)



hit30 <- hits_most %>% top_n(n=30, morans_I)
hit50 <- hits_most %>% top_n(n=50, morans_I)



write.table(hits,file="AB_beta_1_transition.txt",col.names=T,row.names=F,sep="\t",quote=F)




pt.matrix <- as.matrix(cds_from_seurat@assays@data$counts[match(hit30$gene_short_name,rowData(cds_from_seurat)[,1]),order(pseudotime(cds_from_seurat))])


#pt.matrix <- exprs(cds_from_seurat)[match(genes,rownames(rowData(cds_from_seurat))),order(pseudotime(cds_from_seurat))]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- rownames(hit30);
ht_ins_vitro <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 10),
  km = 1,
  row_title_rot                = 5,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  column_title = "Expression Patterns Clustering from AB to beta_1 Transition")
print(ht_ins_vitro)





pt.matrix <- as.matrix(cds_from_seurat@assays@data$counts[match(hits_most$gene_short_name,rowData(cds_from_seurat)[,1]),order(pseudotime(cds_from_seurat))])


#pt.matrix <- exprs(cds_from_seurat)[match(genes,rownames(rowData(cds_from_seurat))),order(pseudotime(cds_from_seurat))]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- rownames(hits_most);
ht_ins_vitro <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 7),
  km = 1,
  row_title_rot                = 5,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  column_title = "Expression Patterns Clustering from AB to beta_1 Transition")
print(ht_ins_vitro)
```

# Gene Commonality Analysis 
Commonality genes are simple intersections in each trajectory. We have compiled the genes in the each forward and reverse trajectories and obtained the common genes either alpha or beta cell trajectory. 

```

devtools::install_github("gaospecial/ggVennDiagram")

common_tobeta <- list("TA"=t1_betaup,
                      "TB"=t2_alpha_down,
                      "TC"=t3_alpha_down,
                      "TD"=t4_alpha_down)

common_toalpha <- list("TA"=t1_betadown,
                      "TB"=t2_alpha_up,
                      "TC"=t3_alpha_up,
                      "TD"=t4_alpha_up)

ggVennDiagram(common_tobeta, label_alpha = 0) + ggtitle("Commonality to Beta Cell Trajectory") + scale_fill_gradient(low="white",high = "red")
ggVennDiagram(common_toalpha, label_alpha = 0) + ggtitle("Commonality to Alpha Cell Trajectory") + scale_fill_gradient(low="white",high = "red")


alpha_up_common <- list("T12"=intersect(t1_betadown, t2_alpha_up),
                        "T23"=intersect(t2_alpha_up, t3_alpha_up),
                        "T34"=intersect(t3_alpha_up, t4_alpha_up),
                        "T13"=intersect(t1_betadown, t3_alpha_up),
                        "T14"=intersect(t1_betadown, t4_alpha_up),
                        "T24"=intersect(t2_alpha_up, t4_alpha_up),
                        "T123"=intersect(intersect(t1_betadown, t2_alpha_up), t3_alpha_up),
                        "T234"=intersect(intersect(t2_alpha_up, t3_alpha_up), t4_alpha_up),
                        "T134"=intersect(intersect(t1_betadown, t3_alpha_up), t4_alpha_up),
                        "All"=intersect(intersect(intersect(t1_betadown, t3_alpha_up), t2_alpha_up), t4_alpha_up))

alpha_up_common$T12 -> C_AB
write.csv(as.data.frame(C_AB), file="./paper_table/C_AB.csv")




beta_up_common <- list("T12"=intersect(t1_betaup, t2_alpha_down),
                        "T23"=intersect(t2_alpha_down, t2_alpha_down),
                        "T34"=intersect(t3_alpha_down, t4_alpha_down),
                        "T13"=intersect(t1_betaup, t3_alpha_down),
                        "T14"=intersect(t1_betaup, t4_alpha_down),
                        "T24"=intersect(t2_alpha_down, t4_alpha_down),
                        "T123"=intersect(intersect(t1_betaup, t2_alpha_down), t3_alpha_down),
                        "T234"=intersect(intersect(t2_alpha_down, t3_alpha_down), t4_alpha_down),
                        "T134"=intersect(intersect(t1_betaup, t3_alpha_down), t4_alpha_down),
                        "All"=intersect(intersect(intersect(t1_betaup, t3_alpha_down), t2_alpha_down), t4_alpha_down))
```




