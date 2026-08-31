# ########################################################
# This script aims to analyse the heterogeneity of cells
#  using dimension reduction techniques
# ########################################################


# DIMENSIONAL REDUCTION (TSNE/UMAP)
###################################

# ..........................................................................................................
## @knitr heterogeneity_dimReduc
# ..........................................................................................................

nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
if(DIMREDUC_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_dimreduc = nbPC
}

sc10x = RunUMAP( sc10x, dims = 1:nbPC_dimreduc);
sc10x = RunTSNE( sc10x, dims = 1:nbPC_dimreduc);

# Save resulting coordinates for all cells as 'tsv' files
# write.table( Embeddings(sc10x, reduction = "umap"),
#              file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_umap.tsv")),
#              quote = FALSE,
#              row.names = TRUE,
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");
# 
# write.table( Embeddings(sc10x, reduction = "tsne"),
#              file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_tsne.tsv")),
#              quote = FALSE,
#              row.names = TRUE,
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");
# 

# ..........................................................................................................
## @knitr dimReduc_ggplot_covariables
# ..........................................................................................................

## Choose a safe reduction to use (robust, no error possible)
## Priority: useReduction > umap > pca

reduc <- "umap"
if (exists("useReduction")) {
  reduc <- useReduction
}
if (!reduc %in% Reductions(sc10x)) {
  if ("umap" %in% Reductions(sc10x)) {
    reduc <- "umap"
  } else if ("pca" %in% Reductions(sc10x)) {
    reduc <- "pca"
  } else {
    stop("No valid dimensional reduction found in sc10x")
  }
}

# Use the package HemaScribe from https://github.com/RabadanLab/HemaScribe 

sc10x <- HemaScribe(
  input     = sc10x,
  prefilter = 0,     
  reference = "WT"    
);

sc10x <- HemaScape(
  input = sc10x,
  vars  = "orig.ident",  
  sigma = 0.1            
);

sc10x$combined.annot <- as.character(sc10x$combined.annot)
sc10x$combined.annot[sc10x$combined.annot == "Immature_neutrophil"] <- "Immature_neutro"
sc10x$combined.annot <- factor(sc10x$combined.annot)

# Mettre à jour les Idents
Idents(sc10x) <- "combined.annot"

# Sélectionner les cellules dans la zone haute de l'EE (branche EryP)
ee_coords <- as.data.frame(Embeddings(sc10x, reduction = "EE"))
colnames(ee_coords) <- c("EE_1", "EE_2")

# Cellules dans la branche haute (ajuste les coordonnées selon ton plot)
cells_eryp_zone <- rownames(ee_coords[ee_coords$EE_1 > 0.35 & 
                                        ee_coords$EE_1 < 0.65 & 
                                        ee_coords$EE_2 > 0.75, ])

# Voir de quels clusters elles viennent
table(sc10x$combined.annot[cells_eryp_zone])
table(sc10x$orig.ident[cells_eryp_zone])

# Identifier les cellules contaminantes dans la zone EryP
cells_to_remove <- rownames(ee_coords[ee_coords$EE_1 > 0.35 & 
                                        ee_coords$EE_1 < 0.65 & 
                                        ee_coords$EE_2 > 0.75 &
                                        !sc10x$combined.annot %in% c("EryP", "MPP2"), ])

cat("Cellules contaminantes à supprimer :", length(cells_to_remove), "\n")
table(sc10x$combined.annot[cells_to_remove])

# Supprimer ces cellules
cells_to_keep <- Cells(sc10x)[!Cells(sc10x) %in% cells_to_remove]
sc10x <- subset(sc10x, cells = cells_to_keep)

cat("Cellules après nettoyage :", ncol(sc10x), "\n")

## Plot cells by sample of origin

print(
  DimPlot(sc10x, reduction = "EE", group.by = "orig.ident") +
    ggtitle("Map of cells by sample of origin")
)


## QC exploration plots

if (QC_EXPLORATION_MODE == TRUE) {
  print(
    DimPlot(sc10x, reduction = reduc, group.by = "outlier") +
      ggtitle("Map of cells by \nQC filtering status")
  )
}


## % Ribosomal genes

print(
  FeaturePlot(sc10x, reduction = "EE", features = "percent.ribo") +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    ggtitle("Map of cells with level of\npercentage of ribosomal genes")
)

if (QC_EXPLORATION_MODE == TRUE) {
  print(
    DimPlot(sc10x, reduction = reduc, group.by = "outlier.percent.ribo") +
      ggtitle("Map of cells by filter status\non % ribosomal genes")
  )
}


## % Mitochondrial genes

print(
  FeaturePlot(sc10x, reduction = "umap", features = "percent.mito") +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    ggtitle("Map of cells with level of\npercentage of mitochondrial genes")
)

if (QC_EXPLORATION_MODE == TRUE) {
  print(
    DimPlot(sc10x, reduction = reduc, group.by = "outlier.percent.mito") +
      ggtitle("Map of cells by filter status\non % mitochondrial genes")
  )
}


## RNA counts (UMIs)

print(
  FeaturePlot(sc10x, reduction = "umap", features = "nCount_RNA") +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    ggtitle("Map of cells with RNA counts")
)

if (QC_EXPLORATION_MODE == TRUE) {
  print(
    DimPlot(sc10x, reduction = reduc, group.by = "outlier.nCount_RNA") +
      ggtitle("Map of cells by filter status\non nCount_RNA value")
  )
}


## Number of detected genes

print(
  FeaturePlot(sc10x, reduction = "EE", features = "nFeature_RNA") +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    ggtitle("Map of cells with\nnumber of detected genes")
)

if (QC_EXPLORATION_MODE == TRUE) {
  print(
    DimPlot(sc10x, reduction = reduc, group.by = "outlier.nFeature_RNA") +
      ggtitle("Map of cells by filter status\non nFeature_RNA value")
  )
}


