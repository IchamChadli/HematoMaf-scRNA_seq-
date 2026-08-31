# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# # CLUSTER FUNCTIONAL ENRICHMENT
# ################################
 
## @knitr heterogeneity_markerGenes_functional_enrichment

# for( clusterName in levels( markers$cluster)){
# # Vérifier que markers$cluster contient bien les noms combined.annot
# levels(markers$cluster)
#   # Add a sub-tab with the Cluster name
#   # ...............................................................................
#   cat( "\n \n")
#   cat( "\n##### Cluster", clusterName, "{.tabset}")
#   cat( "\n \n")
#   
#   
#   # Add a sub-tab with the enrichment of the positive markers in the GO BP terms
#   # ...............................................................................
#   cat( "\n \n")
#   cat( "\n###### Markers genes GO Biological Process")
#   cat( "\n \n")
#   
#   positive_ego <- enrichGO(gene          = markers[ which( markers$cluster == clusterName), "gene"],
#                            universe      = rownames( sc10x@assays$RNA),
#                            OrgDb         = org.Mm.eg.db,
#                            keyType       = "SYMBOL",
#                            ont           = "BP",
#                            pAdjustMethod = "BH",
#                            pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
#                            qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
#   
#   selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
#   
#   if( selected_term_number> 0){  
#     
#     print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of marker genes in BP GO terms"))
#     positive_ego_df = data.frame( positive_ego)
#     
#     positive_ego_simMatrix <- calculateSimMatrix( positive_ego_df$ID,
#                                                   orgdb="org.Mm.eg.db",
#                                                   ont="BP",
#                                                   method="Rel")
#     
#     positive_ego_scores <- setNames(-log10( positive_ego_df$qvalue), positive_ego_df$ID)
#     positive_ego_reducedTerms <- reduceSimMatrix(positive_ego_simMatrix,
#                                                  positive_ego_scores,
#                                                  threshold=0.7,
#                                                  orgdb="org.Mm.eg.db")
#     
#     print( scatterPlot( positive_ego_simMatrix, positive_ego_reducedTerms))
#     
#     treemapPlot( positive_ego_reducedTerms)
#   }
#   else{
#     cat("<BR>No result<BR>")
#   }
#   
#   
#   # Add a sub-tab with the enrichment of the positive markers in the GO MF terms
#   # ...............................................................................
#   cat( "\n \n")
#   cat( "\n###### Markers genes GO Molecular Function")
#   cat( "\n \n")
#   
#   positive_ego <- enrichGO(gene          = markers[ which( markers$cluster == clusterName), "gene"],
#                            universe      = rownames( sc10x@assays$RNA),
#                            OrgDb         = org.Mm.eg.db,
#                            keyType       = "SYMBOL",
#                            ont           = "MF",
#                            pAdjustMethod = "BH",
#                            pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
#                            qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
#   
#   selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
#   
#   if( selected_term_number> 0){  
#     
#     print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of positive markers in MF GO terms"))
#     positive_ego_df = data.frame( positive_ego)
#     
#     positive_ego_simMatrix <- calculateSimMatrix( positive_ego_df$ID,
#                                                   orgdb="org.Mm.eg.db",
#                                                   ont="MF",
#                                                   method="Rel")
#     
#     positive_ego_scores <- setNames(-log10( positive_ego_df$qvalue), positive_ego_df$ID)
#     positive_ego_reducedTerms <- reduceSimMatrix(positive_ego_simMatrix,
#                                                  positive_ego_scores,
#                                                  threshold=0.7,
#                                                  orgdb="org.Mm.eg.db")
#     
#     print( scatterPlot( positive_ego_simMatrix, positive_ego_reducedTerms))
#     
#     treemapPlot( positive_ego_reducedTerms)
#   }
#   else{
#     cat("<BR>No result<BR>")
#   }
# }
# 

for (clusterName in cell_types_keep) {
  
  cluster_markers <- markers[markers$cluster == clusterName & markers$p_val_adj < 0.05, ]
  genes_to_test   <- cluster_markers$gene
  
  cat("\n\n")
  cat("\n##### Cluster", clusterName, "{.tabset}")
  cat("\n\n")
  
  if (length(genes_to_test) >= 5) {
    
    # GO BP
    cat("\n\n")
    cat("\n###### Markers genes GO Biological Process")
    cat("\n\n")
    
    tryCatch({
      positive_ego_BP <- enrichGO(
        gene          = genes_to_test,
        universe      = rownames(sc10x@assays$RNA),
        OrgDb         = org.Mm.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
        qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF
      )
      
      if (!is.null(positive_ego_BP) && sum(positive_ego_BP@result$p.adjust <= positive_ego_BP@pvalueCutoff) > 0) {
        
        p_bp <- dotplot(positive_ego_BP, showCategory = 10) + 
          ggtitle(paste("BP Enrichment -", clusterName)) +
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
          theme(axis.text.y = element_text(size = 10, lineheight = 1.2))
        print(p_bp)
        
        bp_df        <- data.frame(positive_ego_BP)
        bp_simMatrix <- calculateSimMatrix(bp_df$ID, orgdb = "org.Mm.eg.db", ont = "BP", method = "Rel")
        bp_scores    <- setNames(-log10(bp_df$qvalue), bp_df$ID)
        bp_reduced   <- reduceSimMatrix(bp_simMatrix, bp_scores, threshold = 0.7, orgdb = "org.Mm.eg.db")
        
        tryCatch({
          print(scatterPlot(bp_simMatrix, bp_reduced))
          treemapPlot(bp_reduced)
        }, error = function(e) {
          cat("<BR>ScatterPlot/Treemap non disponible (trop peu de termes) :<BR>", conditionMessage(e), "<BR>")
        })
        
      } else {
        cat("<BR>No significant BP result for this cluster.<BR>")
      }
    }, error = function(e) {
      cat("<BR>Erreur GO BP pour le cluster", clusterName, ":", conditionMessage(e), "<BR>")
    })
    
    # GO MF
    cat("\n\n")
    cat("\n###### Markers genes GO Molecular Function")
    cat("\n\n")
    
    tryCatch({
      positive_ego_MF <- enrichGO(
        gene          = genes_to_test,
        universe      = rownames(sc10x@assays$RNA),
        OrgDb         = org.Mm.eg.db,
        keyType       = "SYMBOL",
        ont           = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
        qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF
      )
      
      if (!is.null(positive_ego_MF) && sum(positive_ego_MF@result$p.adjust <= positive_ego_MF@pvalueCutoff) > 0) {
        
        p_mf <- dotplot(positive_ego_MF, showCategory = 10) + 
          ggtitle(paste("MF Enrichment -", clusterName)) +
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
          theme(axis.text.y = element_text(size = 10, lineheight = 1.2))
        print(p_mf)
        
        mf_df        <- data.frame(positive_ego_MF)
        mf_simMatrix <- calculateSimMatrix(mf_df$ID, orgdb = "org.Mm.eg.db", ont = "MF", method = "Rel")
        mf_scores    <- setNames(-log10(mf_df$qvalue), mf_df$ID)
        mf_reduced   <- reduceSimMatrix(mf_simMatrix, mf_scores, threshold = 0.7, orgdb = "org.Mm.eg.db")
        
        tryCatch({
          print(scatterPlot(mf_simMatrix, mf_reduced))
          treemapPlot(mf_reduced)
        }, error = function(e) {
          cat("<BR>ScatterPlot/Treemap non disponible (trop peu de termes) :<BR>", conditionMessage(e), "<BR>")
        })
        
      } else {
        cat("<BR>No significant MF result for this cluster.<BR>")
      }
    }, error = function(e) {
      cat("<BR>Erreur GO MF pour le cluster", clusterName, ":", conditionMessage(e), "<BR>")
    })
    
  } else {
    cat("<BR>Pas assez de gènes différentiels significatifs (< 5) pour générer des analyses GO.<BR>")
  }
}