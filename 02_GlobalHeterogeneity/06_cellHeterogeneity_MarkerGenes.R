# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

# MARKER GENES
##############

## @knitr heterogeneity_markerGenes_old

# # # Identify marker genes
# markers = FindAllMarkers( object          = sc10x,
#                           test.use        = FINDMARKERS_METHOD,
#                           only.pos        = FINDMARKERS_ONLYPOS,
#                           min.pct         = FINDMARKERS_MINPCT,
#                           logfc.threshold = FINDMARKERS_LOGFC_THR,
#                           return.thresh   = FINDMARKERS_PVAL_THR,
#                           random.seed     = SEED,
#                           verbose         = .VERBOSE);



# sc10x$cluster_cond <- paste0(sc10x$seurat_clusters, "_", sc10x$orig.ident)
# Idents(sc10x) <- "cluster_cond"
# # 
# # Test cluster par cluster
# # markers <- FindMarkers(
# #    object          = sc10x,
# #    ident.1         = "12_WT_NT", 
# #    ident.2         = "12_WT_MCSF",   
# #    test.use        = FINDMARKERS_METHOD,
# #    only.pos        = FALSE,       
# #    min.pct         = FINDMARKERS_MINPCT,
# #    logfc.threshold = FINDMARKERS_LOGFC_THR,
# #    random.seed     = SEED,
# #    verbose         = .VERBOSE
# #  )
# # 
# # 
# # 
# # 

## @knitr heterogeneity_markerGenes

# Convertir en character AVANT ifelse pour éviter la conversion en numéros
sc10x$combined.annot <- as.character(sc10x$combined.annot)

# 1. Fusion MPP3
sc10x$combined.annot <- ifelse(
  sc10x$combined.annot %in% c("FcG_neg_MPP3", "FcG_pos_MPP3"), 
  "MPP3", 
  sc10x$combined.annot
)

# 2. Filtrage < 100 cellules
annot_counts_init <- table(sc10x$combined.annot)
clusters_to_keep  <- names(annot_counts_init[annot_counts_init >= 100])

# 3. Subset et refactorisation
sc10x <- subset(sc10x, subset = combined.annot %in% clusters_to_keep)
sc10x$combined.annot <- factor(sc10x$combined.annot, levels = clusters_to_keep)

# ==============================================================================
# ANALYSE DIFFÉRENTIELLE (DÉG) PAR POPULATION
# ==============================================================================
# Annotation and Condition dynamique
sc10x$celltype_condition <- paste(as.character(sc10x$combined.annot), 
                                  sc10x$orig.ident, sep = "_")
Idents(sc10x) <- "celltype_condition"

# Vérifier
head(levels(Idents(sc10x)))

# On garde TOUTES les populations restantes (déjà filtrées à >= 100 cellules)
cell_types_keep <- levels(sc10x$combined.annot)

message("--- Statut des populations ---")
message("Populations conservées (>= 100 cellules) : ", paste(cell_types_keep, collapse = ", "))
message("Populations exclues (< 100 cellules) : ", paste(names(annot_counts_init[annot_counts_init < 100]), collapse = ", "))
message("-------------------------------------------------")

list_results <- list()

# Comparaison dynamique basée sur tes variables SAMPLE_id1 et SAMPLE_id2
for (ct in cell_types_keep) {
  
  # Génération dynamique des identifiants Seurat
  id_1 <- paste0(ct, "_", SAMPLE_id1)
  id_2 <- paste0(ct, "_", SAMPLE_id2)
  
  message("Analyse différentielle pour : ", ct, " (", SAMPLE_id2, " vs ", SAMPLE_id1, ")...")
  
  if (id_2 %in% levels(Idents(sc10x)) && id_1 %in% levels(Idents(sc10x))) {
    
    n_2 <- sum(Idents(sc10x) == id_2)
    n_1 <- sum(Idents(sc10x) == id_1)
    
    # Sécurité interne à Seurat : au moins 10 cellules par bras pour le test statistique
    if (n_2 >= 10 && n_1 >= 10) {
      
      res <- FindMarkers(
        object          = sc10x,
        ident.1         = id_2,   # Groupe d'intérêt (KO)
        ident.2         = id_1,   # Groupe contrôle (WT)
        test.use        = FINDMARKERS_METHOD,
        only.pos        = FALSE, 
        min.pct         = FINDMARKERS_MINPCT,
        logfc.threshold = FINDMARKERS_LOGFC_THR,
        verbose         = .VERBOSE
      )
      
      if (nrow(res) > 0) {
        message("  -> ", nrow(res), " DEG trouvés (", SAMPLE_id2, ": ", n_2, " cells, ", SAMPLE_id1, ": ", n_1, " cells)")
        res$gene    <- rownames(res)
        res$cluster <- ct 
        list_results[[ct]] <- res
      } else {
        message("  -> Aucun gène ne passe les filtres pour ", ct)
      }
      
    } else {
      message("  -> Ignoré : sous-groupe trop petit (", SAMPLE_id2, ": ", n_2, ", ", SAMPLE_id1, ": ", n_1, ")")
    }
    
  } else {
    message("  -> Manque de cellules dans l'une des conditions pour le type : ", ct)
  }
}

markers <- do.call(rbind, list_results)

# Tri des marqueurs pour la table
topMarkers_table = by( markers, markers[["cluster"]], function(x) {
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  return( if(is.null( FINDMARKERS_SHOWTOP_TABLE)) x else head( x, n = FINDMARKERS_SHOWTOP_TABLE));
});

# Tri des marqueurs pour le heatmap
topMarkers_heatmap = by( markers, markers[["cluster"]], function(x) {
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  return( if(is.null( FINDMARKERS_SHOWTOP_HEATMAP)) x else head( x, n = FINDMARKERS_SHOWTOP_HEATMAP));
});

topMarkers_tableDF = do.call( rbind, topMarkers_table);
topMarkers_heatmapDF = do.call( rbind, topMarkers_heatmap);

topMarkers_tableDT = topMarkers_tableDF[c("gene", "cluster", "avg_log2FC", "p_val_adj")]


## @knitr heterogeneity_markerGenes_table
datatable(topMarkers_tableDT,
          class = "compact",
          filter = "top",
          rownames = FALSE,
          selection = 'none',
          colnames = c("Gene", "Cell Type", paste0("Log2FC (", SAMPLE_id2, " vs ", SAMPLE_id1, ")"), "Adj. Pvalue"),
          extensions = c('Buttons', 'Select'),
          options = list(dom = "Blfrtip",
                         autoWidth = FALSE,
                         buttons = list(
                           list(extend = 'copy', text = '📋 Copy'),
                           list(extend = 'print', text = '🖨️ Print'),
                           list(extend = 'excel', text = '📊 Excel'),
                           list(extend = 'csv', text = '📄 CSV'),
                           list(extend = 'pdf', text = '📕 PDF')
                         ), fixedHeader = TRUE)) %>%
  formatStyle(columns = "avg_log2FC", 
              target = "cell",
              background = styleColorBar(data = range(topMarkers_tableDT[["avg_log2FC"]]), 'lightblue', angle = -90),
              backgroundSize = '95% 50%',
              backgroundRepeat = 'no-repeat',
              backgroundPosition = 'center') %>%
  formatStyle(columns = "cluster",
              backgroundColor = styleEqual(names(combined_cols),
                                           scales::alpha(combined_cols, 0.3)))


## @knitr heterogeneity_markerGenes_heatmap_mean

# # 1. Configuration des identités sur les annotations combinées directement
Idents(sc10x) <- "combined.annot"
clusterID <- Idents(sc10x)

# Définition directe de cluster_set à partir des niveaux de combined.annot effectivement présents
cluster_set <- intersect(ordre_hematopoietique, levels(factor(sc10x$combined.annot)))

# Extraire les gènes uniques requis pour éviter les duplications lors de l'extraction
topMarkersGenes <- unique(topMarkers_heatmapDF[["gene"]])

# OPTIMISATION RAM : On extrait uniquement la matrice pour les gènes cibles !
expMat <- as.matrix(GetAssayData(sc10x, slot = "data")[topMarkersGenes, , drop = FALSE])

# Réalignement et tri des cellules selon l'ordre des clusters
clusterOrdering <- order(clusterID)
expMat          <- expMat[, clusterOrdering, drop = FALSE]
clusterID       <- clusterID[clusterOrdering]

# 2. Calcul robuste des moyennes par cluster
# ........................................................................
mean_list <- list()

for(cluster_id in cluster_set) {
  # Cellules appartenant au cluster courant
  clusters_cells <- names(clusterID)[which(clusterID == cluster_id)]
  
  if(length(clusters_cells) > 0) {
    # Calcul de la moyenne pour ce cluster (uniquement pour les gènes associés à ce cluster)
    genes_this_cluster <- topMarkers_heatmapDF[topMarkers_heatmapDF$cluster == cluster_id, "gene"]
    
    # Sécurité au cas où aucun gène n'est trouvé pour un cluster donné
    if(length(genes_this_cluster) > 0) {
      
      all_means_for_these_genes <- sapply(cluster_set, function(cl) {
        cells_cl <- names(clusterID)[which(cl == clusterID)]
        if(length(cells_cl) > 0) {
          rowMeans(expMat[genes_this_cluster, cells_cl, drop = FALSE])
        } else {
          rep(0, length(genes_this_cluster))
        }
      })
      
      # Si un seul gène est présent, on force le format matrice
      if(is.vector(all_means_for_these_genes)) {
        all_means_for_these_genes <- matrix(all_means_for_these_genes, nrow = 1, dimnames = list(genes_this_cluster, cluster_set))
      }
      
      # Ajout d'un préfixe pour éviter les doublons
      rownames(all_means_for_these_genes) <- paste0(cluster_id, ".", genes_this_cluster)
      mean_list[[cluster_id]] <- all_means_for_these_genes
    }
  }
}

# Fusion verticale de tous les blocs de moyennes calculés
meanExpMat <- do.call(rbind, mean_list)

# 3. Préparation rigoureuse des couleurs et annotations
# ........................................................................
row_clusters <- sub("\\..*$", "", rownames(meanExpMat))

annotation_row_df <- data.frame(
  Markers = factor(row_clusters, levels = cluster_set), 
  stringsAsFactors = FALSE, 
  row.names = rownames(meanExpMat)
)

annotation_col_df <- data.frame(
  Cluster = factor(cluster_set, levels = cluster_set), 
  stringsAsFactors = FALSE, 
  row.names = cluster_set
)

# /!\ ALIGNEMENT STRICT ET SÉCURISÉ DES COULEURS /!\
clean_cols <- combined_cols[cluster_set]

# Sécurité si une couleur manque
clean_cols[is.na(clean_cols)] <- "#D3D3D3" 
names(clean_cols) <- cluster_set

# Construction de la liste finale que pheatmap va valider
annotation_colors_list <- list(
  Markers = clean_cols,
  Cluster = clean_cols
)

# 4. Génération des Heatmaps et du DotPlot
# ............................................................................

cat("\n \n")
print(pheatmap(expMat,
               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(80),
               breaks = seq(-4.2, 4.2, 0.1),
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "row",
               annotation_row = annotation_row_df,
               annotation_col = data.frame(Cluster = factor(as.character(clusterID), levels = cluster_set), stringsAsFactors = FALSE, row.names = colnames(expMat)),
               annotation_colors = annotation_colors_list,
               show_colnames = FALSE,
               fontsize_row = 8,
               main = "Scaled normalized expression by cell\nof top markers"));

cat("\n \n")
print(DotPlot(sc10x, 
              group.by = "combined.annot",
              features = topMarkersGenes,
              cols = "RdBu") +
        scale_y_discrete(limits = rev(intersect(ordre_hematopoietique, unique(sc10x$combined.annot)))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Scaled normalized mean expression of topmarker genes in clusters"));

cat("\n \n")
genes_clean_labels <- sub("^.*\\.", "", rownames(meanExpMat))
print(pheatmap(meanExpMat,
               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "row",
               annotation_row = annotation_row_df,
               annotation_col = annotation_col_df,
               annotation_colors = annotation_colors_list,
               show_colnames = TRUE,
               labels_row = genes_clean_labels,
               fontsize_row = 8,
               main = "Scaled mean normalized expression\nby cluster of top markers"));


## @knitr heterogeneity_markerGenes_expression_projection

# Déterminer la réduction à utiliser
reduction_to_use <- if(exists("useReduction")) useReduction else "EE"

invisible(lapply(names(topMarkers_heatmap), function(clusterName) {
  if(!(clusterName %in% cluster_set)) return(NULL)
  
  cat("\n \n")
  cat("#### Cl. <span style='border-radius: 3px; border: 3px solid ", combined_cols[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n")
  cat("\n \n")
  
  # Highlight du cluster
  highlightClusterPlot(clusterName, seuratObject = sc10x, reduction = reduction_to_use)
  
  # FeaturePlot pour chaque gène du cluster
  invisible(lapply(topMarkers_heatmap[[clusterName]][["gene"]], function(featureName) {
    
    gene_fc <- topMarkers_heatmap[[clusterName]][
      topMarkers_heatmap[[clusterName]]$gene == featureName, "avg_log2FC"]
    
    cols <- if(gene_fc > 0) c("lightgrey", "red") else c("lightgrey", "blue")
    
    print(FeaturePlot(sc10x,
                      features  = featureName,
                      reduction = reduction_to_use,
                      cols      = cols) +
            theme(axis.title.x    = element_blank(),
                  axis.title.y    = element_blank(),
                  legend.position = "none") +
            ggtitle(paste0(featureName, if(gene_fc > 0) " up" else " down")))
  }))
  
  cat(" \n \n")
}));


## @knitr heterogeneity_markerGenes_expression_violin

# Appliquer l'ordre hématopoïétique
sc10x$combined.annot <- factor(sc10x$combined.annot, 
                               levels = intersect(ordre_hematopoietique, 
                                                  unique(sc10x$combined.annot)))
Idents(sc10x) <- "combined.annot"

# Plot expression values of marker genes as violinplot for each cluster
invisible(lapply(names(topMarkers_heatmap), function(clusterName) {
  if(!(clusterName %in% cluster_set)) return(NULL)
  
  cat("\n \n")
  cat("#### Cl. <span style='border-radius: 3px; border: 3px solid ", combined_cols[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");
  cat("\n \n")
  
  # Remind cluster name in an empty figure
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
  text(x = 0.5, y = 0.5, paste("Cluster", clusterName), cex = 2, col = combined_cols[clusterName]);
  
  # Violinplot en utilisant la palette combinée
  invisible(lapply(topMarkers_heatmap[[clusterName]][["gene"]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = combined_cols, slot = "data"));
  
  cat(" \n \n"); 
}));

## @knitr celltype_per_clusters 

# Annotation des clusters (progéniteurs)

# 
# sc10x <- RenameIdents(sc10x, CLUSTER_ANNOTATION_KO_NT_KO_MCSF)
# 
# sc10x$celltype <- Idents(sc10x)
# 
# DimPlot(sc10x, label = TRUE, group.by = "celltype")


# Couleurs pour les deux conditions
# cols <- c("KO_NT" = "#66C2A5", "KO_MCSF" = "#FC8D62")
# 
# 
# composition_df <- sc10x[[]] %>%
#   
#   group_by(seurat_clusters, orig.ident) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   
#   
#   group_by(orig.ident) %>%
#   mutate(freq_per_cond = n / sum(n)) %>% 
#   ungroup() %>%
#   
#   
#   group_by(seurat_clusters) %>%
#   mutate(pct_normalized = freq_per_cond / sum(freq_per_cond) * 100) %>%
#   ungroup() %>%
#   
#   
#   mutate(celltype = CLUSTER_ANNOTATION_KO_NT_KO_MCSF[as.character(seurat_clusters)])
# 
# 
# ggplot(composition_df, aes(x = pct_normalized, y = celltype, fill = orig.ident)) +
#   geom_bar(stat = "identity") +
#   geom_vline(xintercept = 50, linetype = "solid", color = "black", linewidth = 0.8) +
#   scale_fill_manual(values = cols) +
#   labs(title = "Proportions normalisées (Biais de taille d'échantillon retiré)",
#        x = "% de contribution normalisée", 
#        y = "Cell type") +
#   theme_minimal() +
#   theme(legend.position = "bottom")

# Hemascribe github

# sc10x <- HemaScribe(
#   input     = sc10x,
#   prefilter = 0,     
#   reference = "WT"    
# )
# 
# sc10x <- HemaScape(
#   input = sc10x,
#   vars  = "orig.ident",  
#   sigma = 0.1            
# )

# Variable loading 

broad_pops <- unique(sc10x$broad.annot)
# broad_cols <- setNames(as.vector(pals::alphabet(length(broad_pops))), broad_pops)

fine_pops <- sort(unique(sc10x$fine.annot))
# fine_cols <- setNames(as.vector(pals::polychrome(length(fine_pops))), fine_pops)

combined_pops <- sort(unique(sc10x$combined.annot))
# combined_cols <- setNames(as.vector(pals::alphabet(length(combined_pops))), combined_pops)

# Broad

DimPlot(sc10x, group.by = "broad.annot", cols = broad_cols, label = TRUE, repel = TRUE) + 
  ggtitle("Broad Annotation")

DimPlot(sc10x, 
        reduction = "EE",     
        group.by = "broad.annot", 
        label = FALSE, 
        cols = broad_cols,
        pt.size = 0.8) +            
  ggtitle("EE Broad Annotation") +
  theme_minimal()


# Fine 

DimPlot(sc10x, group.by = "fine.annot", cols = fine_cols, label = TRUE, repel = TRUE) + 
  ggtitle("Fine Annotation (HSPC Subtypes)")

DimPlot(sc10x, 
        reduction = "EE",     
        group.by = "fine.annot", 
        label = FALSE, 
        cols = fine_cols,
        pt.size = 0.8) +            
  ggtitle("EE Fine Annotation (HSPC Subtypes)") +
  theme_minimal()

# Combined

# Supprimer les clusters avec moins de 50 cellules

# Supprimer les cellules du cluster 6
#Idents(sc10x) <- "seurat_clusters"
# sc10x <- subset(sc10x, idents = 7, invert = TRUE)

# Vérifier
table(sc10x$combined.annot)


#sc10x <- subset(sc10x$combined.annot, drop(sc10x$combined.annot))

annot_counts <- table(sc10x$combined.annot)
clusters_to_remove <- names(annot_counts[annot_counts < 100])

print(clusters_to_remove)

# # Fusionner FcG_neg_MPP3 et FcG_pos_MPP3 en MPP3
# sc10x$combined.annot <- as.character(sc10x$combined.annot)
# sc10x$combined.annot[sc10x$combined.annot %in% c("FcG_neg_MPP3", "FcG_pos_MPP3")] <- "MPP3"
# sc10x$combined.annot <- factor(sc10x$combined.annot)

sc10x <- subset(sc10x, combined.annot %in% clusters_to_remove, invert = TRUE)

clustersCount <- as.data.frame(table(Cluster = sc10x[["combined.annot"]]), responseName = "CellCount")

# Garder seulement les clusters avec >= 70 cellules
clustersCount <- clustersCount[clustersCount$CellCount >= 70, ]

datatable(clustersCount,
          class = "compact",
          rownames = FALSE,
          colnames = c("Cluster", "Nb. Cells"),
          options = list(dom = "<'row'rt>",
                         autoWidth = FALSE,
                         columnDefs = list(list(targets = 0:(ncol(clustersCount)-1),
                                                className = 'dt-center')),
                         orderClasses = FALSE,
                         paging = FALSE,
                         processing = TRUE,
                         scrollCollapse = TRUE,
                         scroller = TRUE,
                         scrollX = TRUE,
                         scrollY = "525px",
                         stateSave = TRUE))

DimPlot(sc10x, group.by = "combined.annot", cols = combined_cols, label = TRUE, repel = TRUE) + 
  ggtitle("Combinded Annotation")

DimPlot(sc10x, 
        reduction = "EE",     
        group.by = "combined.annot", 
        label = FALSE, 
        cols = combined_cols,
        pt.size = 0.4) +            
  ggtitle("EE Combined Annotation") +
  theme_minimal()

p <- DimPlot(sc10x, 
             reduction = "EE",     
             group.by = "combined.annot", 
             label = FALSE, 
             cols = combined_cols,
             pt.size = 0.7) +            
  ggtitle("EE Combined Annotation") +
  theme_minimal() +
  theme(legend.position = "none")

# Extraire l'ordre exact des couleurs depuis le plot lui-même
plot_data <- p$data
cluster_color_map <- combined_cols[as.character(unique(plot_data$combined.annot))]

LabelClusters(plot = p, 
              id = "combined.annot",
              repel = TRUE,
              size = 4,
              box = TRUE,
              fill = cluster_color_map,
              color = "black",
              alpha = 0.8);

## @knitr scoring

# On ouvre l'onglet principal "Scoring" avec l'option tabset pour ses sous-onglets
cat("\n# Scoring {.tabset .tabset-fade}\n")

apoptosis_genes <- c("Bax", "Bak1", "Bcl2", "Bcl2l1", "Mcl1", "Bbc3", "Pmaip1", "Trp53")
inflammation_genes <- c("Stat1", "Irf1", "Irf7", "Isg15", "Oas1a", "Nfkb1", "Rela", "Tnfaip3", "Cxcl2")
aging_genes <- c("Cdc42", "Clu", "Apoe", "Selp", "Nupr1", "Cdkn2a")
quiescence_genes <- c("Cdkn1a", "Cdkn1b", "G0s2", "Fgd5", "Mecom", "Mitf")  

# Apoptose
sc10x <- AddModuleScore(sc10x, features = list(apoptosis_genes), name = "Apoptose_score")

# Inflammation
sc10x <- AddModuleScore(sc10x, features = list(inflammation_genes), name = "Inflammation_score")

# Aging (Vieillissement)
sc10x <- AddModuleScore(sc10x, features = list(aging_genes), name = "Aging_score")

# Quiescence
sc10x <- AddModuleScore(sc10x, features = list(quiescence_genes), name = "Quiescence_score")

# 1. On extrait les métadonnées dans un tableau propre pour ggplot
df_meta <- sc10x@meta.data

# 2. On vire les cellules qui n'ont pas de nom (les NA) pour nettoyer le plot
df_meta <- df_meta[!is.na(df_meta$combined.annot), ]

# 4. On applique cet ordre (inversé avec rev() pour que HSC soit en haut sur l'axe Y)
df_meta$combined.annot <- factor(df_meta$combined.annot, levels = rev(ordre_hematopoietique))

# Sous-onglet 1 : Quiescence
cat("\n## Quiescence \n")
p1 <- ggplot(df_meta, aes(x = Quiescence_score1, y = combined.annot, fill = orig.ident)) +
  geom_density_ridges(scale = 2.5, rel_min_height = 0.01, alpha = 0.5, color = "black") + 
  scale_fill_manual(values = c("KO_NT" = "#FC8D62", "KO_MCSF" = "#2ca25f")) + 
  ggtitle("Quiescence Score") +
  theme_minimal() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())
print(p1)
cat("\n")

# Sous-onglet 2 : Inflammation
cat("\n## Inflammation \n")
p2 <- ggplot(df_meta, aes(x = Inflammation_score1, y = combined.annot, fill = orig.ident)) +
  geom_density_ridges(scale = 2.5, rel_min_height = 0.01, alpha = 0.5, color = "black") + 
  scale_fill_manual(values = c("KO_NT" = "#FC8D62", "KO_MCSF" = "#2ca25f")) + 
  ggtitle("Inflammation Score") +
  theme_minimal() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())
print(p2)
cat("\n")

# Sous-onglet 3 : Aging
cat("\n## Aging \n")
p3 <- ggplot(df_meta, aes(x = Aging_score1, y = combined.annot, fill = orig.ident)) +
  geom_density_ridges(scale = 2.5, rel_min_height = 0.01, alpha = 0.5, color = "black") + 
  scale_fill_manual(values = c("KO_NT" = "#FC8D62", "KO_MCSF" = "#2ca25f")) + 
  ggtitle("Aging Score") +
  theme_minimal() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())
print(p3)
cat("\n")

# Sous-onglet 4 : Apoptosis
cat("\n## Apoptosis \n")
p4 <- ggplot(df_meta, aes(x = Apoptose_score1, y = combined.annot, fill = orig.ident)) +
  geom_density_ridges(scale = 2.5, rel_min_height = 0.01, alpha = 0.5, color = "black") + 
  scale_fill_manual(values = c("KO_NT" = "#FC8D62", "KO_MCSF" = "#2ca25f")) + 
  ggtitle("Apoptosis Score") +
  theme_minimal() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())
print(p4)
cat("\n")

# # 1 & 2 : Tes étapes de calcul restent STRICTEMENT identiques
# global_abundance <- sc10x[[]] %>%
#   group_by(combined.annot) %>%
#   summarise(n_total = n(), .groups = "drop") %>%
#   mutate(real_pct_of_dataset = round(n_total / sum(n_total) * 100, 1))
# 
# composition_df <- sc10x[[]] %>%
#   group_by(combined.annot, orig.ident) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(orig.ident) %>%
#   mutate(freq_in_sample = n / sum(n)) %>%
#   ungroup() %>%
#   group_by(combined.annot) %>%
#   mutate(pct_normalized = freq_in_sample / sum(freq_in_sample) * 100) %>%
#   left_join(global_abundance, by = "combined.annot") %>%
#   ungroup()
# 
# # 3. Graphique avec l'ajout des barres d'abondance à droite
# cols <- c("WT_NT" = "#66C2A5", "WT_MCSF" = "#377EB8")
# 
# ggplot(composition_df, aes(x = pct_normalized, y = combined.annot, fill = orig.ident)) +
#   # Tes barres de composition originales
#   geom_bar(stat = "identity", width = 0.7) +
#   
#   # AJOUT : Les barres d'abondance globale (en gris à droite)
#   # Elles commencent à 102 et s'arrêtent à 102 + le pourcentage
#   geom_rect(aes(xmin = 102, xmax = 102 + real_pct_of_dataset, 
#                 ymin = as.numeric(as.factor(combined.annot)) - 0.3, 
#                 ymax = as.numeric(as.factor(combined.annot)) + 0.3),
#             fill = "grey80", inherit.aes = FALSE,
#             data = composition_df[!duplicated(composition_df$combined.annot), ]) +
#   
#   # On garde le texte au bout des nouvelles barres grises
#   geom_text(aes(x = 103 + real_pct_of_dataset, y = combined.annot, label = paste0(real_pct_of_dataset, "%")), 
#             hjust = 0, size = 3, color = "black", inherit.aes = FALSE,
#             data = composition_df[!duplicated(composition_df$combined.annot), ]) + 
#   
#   geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
#   scale_fill_manual(values = cols) +
#   # On élargit l'axe X pour accommoder les barres grises (jusqu'à ~150 selon ton plus gros cluster)
#   scale_x_continuous(limits = c(0, 150), breaks = c(0, 25, 50, 75, 100)) +
#   labs(title = "Analyse de la composition par Cluster",
#        subtitle = "Gauche : Balance normalisée WT MCSF | Droite : Abondance globale (% du dataset)",
#        x = "Contribution relative (normalisée)",
#        y = "Cluster / Cell Type",
#        fill = "Condition") +
#   theme_minimal()



# # 3. Graphique avec l'ajout des barres d'abondance à droite
# cols <- c("WT_NT" = "#66C2A5", "KO_NT" = "#FC8D62")   # KO orange
# 
# ggplot(composition_df, aes(x = pct_normalized, y = combined.annot, fill = orig.ident)) +
#   # Tes barres de composition originales
#   geom_bar(stat = "identity", width = 0.7) +
#   
#   
#   
#   geom_rect(aes(xmin = 102, xmax = 102 + real_pct_of_dataset, 
#                 ymin = as.numeric(as.factor(combined.annot)) - 0.3, 
#                 ymax = as.numeric(as.factor(combined.annot)) + 0.3),
#             fill = "grey80", inherit.aes = FALSE,
#             data = composition_df[!duplicated(composition_df$combined.annot), ]) +
#   
#   
#   geom_text(aes(x = 103 + real_pct_of_dataset, y = combined.annot, label = paste0(real_pct_of_dataset, "%")), 
#             hjust = 0, size = 3, color = "black", inherit.aes = FALSE,
#             data = composition_df[!duplicated(composition_df$combined.annot), ]) + 
#   
#   geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
#   scale_fill_manual(values = cols) +
#   
#   scale_x_continuous(limits = c(0, 150), breaks = c(0, 25, 50, 75, 100)) +
#   labs(title = "Cluster Composition Analysis",
#        subtitle = "Left: Normalized NT/MCSF balance | Right: Overall abundance (% of dataset)",
#        x = "Relative contribution",
#        y = "Cluster / Cell Type",
#        fill = "Condition") +
#   theme_minimal()

# DotPlot(
#   sc10x,
#   features = CLUSTERS_MARKERS
# ) + RotatedAxis()

# Remettre les numéros de clusters
# Idents(sc10x) <- sc10x$seurat_clusters

## @knitr pseudotime_hemascape

# Variable loading 

branch_cols <- setNames(as.vector(pals::polychrome(length(unique(sc10x$branch_pred)))), 
                        unique(sc10x$branch_pred))
density_cols <- setNames(as.vector(pals::glasbey(length(unique(sc10x$density_cluster_pred)))), 
                         unique(sc10x$density_cluster_pred))

n_segments <- length(unique(sc10x$branch_segment_clusters_pred))
calc_cols <- colorRampPalette(pals::alphabet())
segment_cols <- setNames(calc_cols(n_segments), 
                         unique(sc10x$branch_segment_clusters_pred))


# Branch prediction
DimPlot(sc10x, reduction = "umap", group.by = "branch_pred", cols = branch_cols, label = TRUE, repel = TRUE) + 
  ggtitle("HemaScape: Branch Prediction")

DimPlot(sc10x, 
        reduction = "EE", 
        group.by = "branch_pred", 
        label = FALSE,
        cols = branch_cols,
        pt.size = 0.4) + 
  ggtitle("EE HemaScape: Branch Prediction") +
  theme_minimal()

# Density
DimPlot(sc10x, reduction = "umap", group.by = "density_cluster_pred", label = TRUE)

DimPlot(sc10x, 
        reduction = "EE", 
        group.by = "density_cluster_pred", 
        label = FALSE,
        cols = density_cols,
        pt.size = 0.4) + 
  ggtitle("EE HemaScape: Density Prediction") +
  theme_minimal()


# Prediction
DimPlot(sc10x, reduction = "umap", group.by = "branch_segment_clusters_pred", label = TRUE)

DimPlot(sc10x, 
        reduction = "EE", 
        group.by = "branch_segment_clusters_pred", 
        label = FALSE,
        cols = segment_cols,
        pt.size = 0.4) + 
  ggtitle("EE HemaScape: Segment Clusters") +
  theme_minimal()


FeaturePlot(sc10x, 
            features = "pseudotime_pred", 
            reduction = "umap") +
  scale_color_gradientn(
    colours = c("blue", "cyan", "yellow", "red", "darkred"),
    na.value = "grey90"
  ) +
  labs(title = "Pseudotime", color = "pseudotime") +
  theme_minimal()

FeaturePlot(sc10x, 
            features = "pseudotime_pred", 
            reduction = "EE") +
  scale_color_gradientn(
    colours = c("blue", "cyan", "yellow", "red", "darkred"),
    na.value = "grey90"
  ) +
  labs(title = "Pseudotime", color = "pseudotime") +
  theme_minimal()

# G2M 
g2m_genes <- unlist(CELL_CYCLE_G2MPHASE_GENELIST)
sc10x <- AddModuleScore(sc10x, features = list(g2m_genes), name = "G2M_score")

# S phase
s_genes <- unlist(CELL_CYCLE_SPHASE_GENELIST)
sc10x <- AddModuleScore(sc10x, features = list(s_genes), name = "S_score")

# Feature plots avec taille des points à 0.8
FeaturePlot(sc10x, features = "G2M_score1", reduction = "EE", pt.size = 0.8) +
  ggtitle("G2M Phase Score") +
  scale_color_gradientn(colours = c("blue", "white", "red"))

FeaturePlot(sc10x, features = "S_score1", reduction = "EE", pt.size = 0.8) +
  ggtitle("S Phase Score") +
  scale_color_gradientn(colours = c("blue", "white", "red"))


## @knitr genes_of_interest

# Appliquer l'ordre hématopoïétique aux Idents
sc10x$combined.annot <- factor(sc10x$combined.annot, 
                               levels = intersect(ordre_hematopoietique, 
                                                  unique(sc10x$combined.annot)))
Idents(sc10x) <- "combined.annot"

genes_of_interest <- c("Pu.1", "Cebpb", "Sirt3")  # noms en format souris

# Vérifier que les gènes existent dans l'objet
genes_found <- genes_of_interest[genes_of_interest %in% rownames(sc10x)]
genes_not_found <- genes_of_interest[!genes_of_interest %in% rownames(sc10x)]

# Puis lancer les violins
invisible(lapply(genes_found, function(gene) {
  
  cat("\n \n")
  cat("#### <span style='border-radius: 3px; border: 3px solid black; padding:0px 2px'>", gene, "</span>\n")
  cat("\n \n")
  
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste("Gene:", gene), cex = 2, col = "black")
  
  invisible(violinFeatureByCluster(gene, seuratObject = sc10x, clustersColor = combined_cols, slot = "data"))
  
  cat(" \n \n")
}))