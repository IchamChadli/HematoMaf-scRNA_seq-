###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#

######## CONSTANTS ON THE ANALYSIS STEP IDENTITY
######## Modify the ANALYSIS_STEP_NAME and LITERAL_TITLE variable

# This is the literal title of the analysis step. It will be shown at the beginning
# of the HTML report
# Example : LITERAL_TITLE = "Quality Control, normalization and clustering"
LITERAL_TITLE = "Global heterogeneity HSC Wild Type vs Mafb KO"

# This is the name of the analysis step. Use the same name as the folder
# name of the analysis step
# Example : ANALYSIS_STEP_NAME = "02_GlobalHeterogeneity"
ANALYSIS_STEP_NAME = "03d_GlobalHeterogeneity_HemaScribe_EE_EnrichGo"

# Name of file
SAMPLE = "KO_NT_KO_MCSF"

# This is the path to the analysis step output folder. It will be automatically
# created at first analysis launch
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, SAMPLE)
PATH_ANALYSIS_EXTRA_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, SAMPLE, "extra")


######## CONSTANTS USED BY THE ANALYSIS STEP
######## Add here all the constants required for the analysis

SAMPLE_1 = c(file.path(PATH_EXPERIMENT_OUTPUT,SAMPLE_KO_NT))

SAMPLE_2 = c(file.path(PATH_EXPERIMENT_OUTPUT,SAMPLE_KO_MCSF))

SAMPLE_id1 = "KO_NT"

SAMPLE_id2 = "KO_MCSF"

# Path to the Cell cycle gene lists
CELL_CYCLE_SPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "01_CellCycle", "S_phase_genes.csv")))
CELL_CYCLE_G2MPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "01_CellCycle", "G2M_phase_genes.csv")))

# Path to Heat Shock stress genes list
#PATH_HS_STRESS_MARKER_GENES_TABLE_FILE = file.path( PATH_EXPERIMENT_REFERENCE, "02_HeatShock", "coregene_df-FALSE-v3.csv")

#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (for Seurat3 using 'future')
NBCORES = 4;

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 20000;



#### Filtering / Normalization

# Switch from QC exploration mode and QC filtering
# In exploration mode, the cells that would be filtered by QC threshold are not filtered but marked
# They are shown in the later analysis with a different color
QC_EXPLORATION_MODE = FALSE

# Filters for loading seurat object
LOAD_MIN_CELLS     = 3;    # Retain cells with at least this many features (annotations)
LOAD_MIN_FEATURES  = 200;  # Retain annotations appearing in at least this many cells

# Cells with number of UMIs outside the range will be excluded
FILTER_UMI_MIN     = 0;
FILTER_UMI_MAX     = 40000;

# Cells with number of genes outside the range will be excluded
FILTER_FEATURE_MIN = 0;
FILTER_FEATURE_MAX = 7000;

# Cells with percentage of mitochondrial genes above threshold will be excluded
FILTER_MITOPCT_MAX = 10;

# Cells with percentage of ribosomal genes below threshold will be excluded
FILTER_RIBOPCT_MIN = 5;

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL ;  # c("nCount_RNA") for UMIs (NULL to ignore)




#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report

# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION_RANGE = c( 0.1, 0.2, 0.4, 0.6, 0.8);
FINDCLUSTERS_RESOLUTION     = 0.8;
FINDCLUSTERS_USE_PCA_NBDIMS = 50;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden

# Annotation des clusters (progéniteurs)


#CLUSTER_ANNOTATION_KO_MCSF <- c(
#  "0"  = "LMPP",
#  "1"  = "HSC_like",
#  "2"  = "MPP",
#  "3"  = "G2M",
#  "4"  = "CLP",
#  "5"  = "pMono",
#  "6"  = "pNeu",
#  "7"  = "S_phase",
#  "8"  = "iHSC",
#  "9"  = "GMP",
#  "10" = "pEr",
#  "11" = "pNeuLate",
#  "12" = "cDC1_prog",
#  "13" = "PreB",
#  "14" = "pMast",
#  "15" = "ProB",
#  "16" = "pMk"
#)


# CLUSTER_ANNOTATION_WT_NT <- c(
#  "0"  = "MPP",
#  "1"  = "Myeloid_primed_prog",
#  "2"  = "HSC_like",
#  "3"  = "G2M",
#  "4"  = "GMP",
#  "5"  = "GMP",
#  "6"  = "iHSC",
#  "7"  = "pMono",
#  "8"  = "pNeu",
#  "9"  = "CLP",
#  "10" = "pNeu",
#  "11" = "pMono",
#  "12" = "pNeuLate",
#  "13" = "cDC1_prog",
#  "14" = "ProB",
#  "15" = "pMast")

#CLUSTER_ANNOTATION_WT_MCSF <- c(
#  "0"  = "MPP",
#  "1"  = "CMP",
#  "2"  = "CLP",
#  "3"  = "G2M",
#  "4"  = "S_phase",
#  "5"  = "GMP",
#  "6"  = "Myeloid_primed_prog",
#  "7"  = "pNeu",
#  "8"  = "S_phase",
#  "9"  = "pMono",
#  "10" = "pEr",
#  "11" = "pMast",
#  "12" = "cDC1_prog",
#  "13" = "Late_Ery",
#  "14" = "PreB",
#  "15" = "ProB"
#)

#CLUSTER_ANNOTATION_KO_NT <- c(
#  "0"  = "MPP",
#  "1"  = "CMP",
#  "2"  = "CLP",
#  "3"  = "G2M",
#  "4"  = "S_phase",
#  "5"  = "Myeloid_primed_prog",
#  "6"  = "GMP",
#  "7"  = "pNeu",
#  "8"  = "iHSC",
#  "9"  = "pEr",
#  "10" = "pMast",
#  "11" = "pNeu",
#  "12" = "pNeuLate",
#  "13" = "cDC1_prog",
#  "14" = "Late_Ery",
#  "15" = "pMast",
#  "16" = "PreB",
#  "17" = "ProB"
#)


# CLUSTER_ANNOTATION_WT_KO_NT <- c(
#   "0"  = "STHSC",
#   "1"  = "LTHSC",
#   "2"  = "G2M",
#   "3"  = "LMPP",
#   "4"  = "GMP",
#   "5"  = "pMk",
#   "6"  = "MPP3",
#   "7"  = "iHSC",
#   "8"  = "pNeu",
#   "9"  = "S",
#   "10" = "pEr2",
#   "11" = "MPP4", #Bank1 à verifier 
#   "12" = "pcDC1",
#   "13" = "pMono",
#   "14" = "pEr1",
#   "15" = "pMast", #pBaso
#   "16" = "pL1",
#   "17" = "pL2",
#   "18" = "pNk"
# )

# CLUSTER_ANNOTATION_WT_KO_MCSF <- c(
#   "0" = "LTHSC",
#   "1" = "STHSC",
#   "2" = "MPP4", #(CLP)
#   "3" = "G2M",
#   "4" = "pL2",
#   "5" = "pMono1",
#   "6" = "GMP",
#   "7" = "iHSC",
#   "8" = "S2",
#   "9" = "S1",
#   "10" = "pNeu1",
#   "11" = "pEr2",
#   "12" = "pMono2",
#   "13" = "pMk",
#   "14" = "pNeu2",
#   "15" = "pcDC1",
#   "16" = "pL1",
#   "17" = "pMast",
#   "18" = "pEr1",
#   "19" = "pNk"
# )

# CLUSTER_ANNOTATION_WT_NT_WT_MCSF <- c(
#   "0" = "STHSC",
#   "1" = "MPP4",
#   "2" = "pMk",
#   "3" = "G2M",
#   "4" = "LTHSC",
#   "5" = "GMP",
#   "6" = "S1",
#   "7" = "S2",
#   "8" = "MPP3",
#   "9" = "pNeu1",
#   "10" = "LMPP",
#   "11" = "pEr1",
#   "12" = "pNeu2",
#   "13" = "pcDC1",
#   "14" = "pL1",
#   "15" = "pBaso",
#   "16" = "pEr2",
#   "17" = "pL2",
#   "18" = "pNk"
#   )

# CLUSTER_ANNOTATION_KO_NT_KO_MCSF <- c(
#   "0" = "STHSC",
#   "1" = "LTHSC",
#   "2" = "LMPP",
#   "3" = "G2M",
#   "4" = "pMk1",
#   "5" = "MPP3",
#   "6" = "GMP",
#   "7" = "iHSC",  #(marqueurs d'histones = cycle ?)
#   "8" = "pNeu",
#   "9" = "S",
#   "10" = "pEr1",
#   "11" = "pMono",
#   "12" = "pEr2",
#   "13" = "pcDC1",
#   "14" = "pL",    #(ProB)
#   "15" = "pBaso",
#   "16" = "pMk2",
#   "17" = "pNk",
#   "18" = "pMonor"
#   )

# CLUSTERS_MARKERS <- list(
#   # HSC
#   STHSC = c("Cd34","Gata2","Mllt3"),
#   LTHSC = c("Hlf","Mecom","Fgd5"),
#   iHSC  = c("Gata2","Mecom","Procr"),
#   
#   # Prolifération
#   G2M = c("Mki67","Cenpf","Ccnb1"),
#   S   = c("Pcna","Mcm5","Tyms"),
#   
#   # Progéniteurs myéloïdes
#   CMP = c("Cebpa","Spi1","Cd34"),
#   GMP = c("Cebpa","Csf3r","Ms4a3"),
#   pMono = c("Cd14","Ccr2","S100a9"),
#   pNeu  = c("Mpo","Elane","S100a8"),
#   
#   # MPP
#   MPP3 = c("Cebpa","Ms4a3","Elane"),
#   MPP4 = c("Flt3","Il7r","Dntt"),
#   
#   # Mégacaryocyte
#   pMk = c("Itga2b","Pf4","Mpl"),
#   
#   # Érythroïde
#   pEr1 = c("Klf1","Gata1","Epor"),
#   pEr2 = c("Ahsp","Alas2","Hbb-bs"),
#   
#   # Dendritique
#   pcDC1 = c("Irf8","Batf3","Clec9a"),
#   
#   # Mastocytes
#   pMast = c("Kit","Cpa3","Ms4a2"),
#   
#   # Lymphoïde
#   pL1 = c("Il7r","Rag1","Dntt"),
#   pL2 = c("Il7r","Rag1","Cd79a")
# )
#cluster_annotation_hsc <- c(
#  "HSC_like"            = "LTHSC",
#  "iHSC"                = "STHSC",
#  "MPP"                 = "STHSC",
#  "Myeloid_primed_prog" = "MPP3",
#  "CLP"                 = "MPP4"
#)

#RENAME_CLUSTERS = list()
#RENAME_CLUSTERS[[ SAMPLE_NAME ]] = list()

#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "0"]]  = "MPP_LMPP"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "1"]]  = "Inflammatory_prog"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "2"]]  = "Myeloid_DC_prog"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "3"]]  = "G2M"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "4"]]  = "pMk"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "5"]]  = "pMy"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "6"]]  = "iHSC"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "7"]]  = "pMono"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "8"]]  = "S_phase"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "9"]]  = "pLY"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "10"]] = "pGMPn"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "11"]] = "pEr"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "12"]] = "pNeuLate"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "13"]] = "cDC1_prog"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "14"]] = "pB"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "15"]] = "Eos_Baso_prog"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "16"]] = "PreB_prog"
#RENAME_CLUSTERS[[ SAMPLE_NAME ]][[ "17"]] = "NK_like"
  

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.01;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP_TABLE     = 100;       # Number of marker genes to show in report tables (NULL for all)
FINDMARKERS_SHOWTOP_HEATMAP   = 5;       # Number of marker genes to show in repot heatmaps (NULL for all)

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05





#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 5;
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES=NULL


### Color list 

broad_cols = c("darkslategray3", "#BD0A36", "cadetblue", "palevioletred3", "#B6E4B3","#D1A391", "orange","#BFBB60", "#17BECF","coral3", "darkcyan", "#FB9A99", "#1170AA", "#BB7693","cornflowerblue", "palegreen3", 
           "sienna2","#AD8BC9", "#BCDBCC", "#9D7660",  "#9F8F12")


fine_cols = c("darkslategray3", "#449450", "cadetblue", "palevioletred3", "#B6E4B3","#D1A391", "orange","#BFBB60", "#D5D5D5","coral3", "sienna2", "#FB9A99", "#1170AA", "#BB7693","cornflowerblue", "palegreen3", 
              "sienna2","#AD8BC9", "#BCDBCC", "#9D7660",  "#9F8F12")


#combined_cols = c("#4BBAC3", "#358747", "#D15821", "#B6E4B3","#D1A391", "#B254A5","#BFBB60", "#17BECF","coral3", "#F4D166","#D2293D","#FB9A99", "#1170AA", "#BB7693","cornflowerblue", "palegreen3", 
#                           "sienna2","#AD8BC9", "#BCDBCC", "#9D7660",  "#9F8F12", "#BBE0DB", "#B3B0CE", "#E0B1D0", "#F38A77")


combined_cols <- c(
  # --- Cellules souches et progéniteurs immatures (bas/droite) ---
  "HSC"                 = "#D1A391",   # Beige rosé doux
  "STHSC"               = "#1170AA",   # Bleu roi intense (contraste fort avec le beige et le rouge)
  
  # --- Progéniteurs multipotents (Centre) ---
  "MPP2"                = "#F4D166",   # Jaune d'or lumineux
  "MPP3"                = "#D2293D",   # Rouge vif (très distinct du jaune MPP2 et du rose MPP4)
  "MPP4"                = "#FB9A99",   # Rose pastel
  
  # --- Lignée Mégalocaryocyte / Érythroïde (Haut et Droite) ---
  "EryP"                = "#D15821",   # Orange brique (pôle haut)
  "MkP"                 = "#9D7660",   # Brun terre (pôle droite, loin de HSC et MPP2)
  
  # --- Lignée Monocyte / Granulocyte (Pôle Gauche - Fort contraste requis) ---
  "mGMP"                = "#17BECF",   # Turquoise vibrant
  "cMoP"                = "#358747",   # Vert forêt foncé
  "GP"                  = "#BBE0DB",   # Menthe / Vert d'eau très clair
  "Immature_neutro"     = "#9F8F12",   # Olive / Jaune moutarde foncé
  
  # --- Lignée Lymphoïde (Bas Gauche) ---
  "CLP"                 = "#4BBAC3",   # Cyan
  "Immature_B_cell"     = "#B254A5",   # Violet pourpre (très distinct du vert/bleu au-dessus)
  
  # --- Autres populations rares ---
  "pDC"                 = "#AD8BC9",   # Mauve
  "cDC"                 = "#BB7693",   # Rose violacé foncé
  "Basophil"            = "#FC8D62"    # Corail
)


# Hematopoietic Tree

ordre_hematopoietique <- c(
  "HSC", "STHSC",                  # Cellules souches
  "MPP2", "MPP3", "MPP4",          # Progéniteurs multipotents
  "MkP", "EryP",                   # Lignée Mégalocaryocyte / Érythroïde
  "cMoP", "mGMP",                  # Lignée Monocyte / Granulocyte
  "Immature_neutro", "GP",     # Neutrophiles / Granulocytes progéniteurs
  "CLP", "Immature_B_cell", "pDC","cDC",
  "Basophil"         
)

## Genes monitored individually (tsv file, one column for each group of genes)
#MONITORED_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                   #               "03_MonitoredGenes",
                                  #                "Monitored.csv"),
                                 #      sep = ",",
                                    #   header = TRUE,
                                     #  stringsAsFactors = FALSE,
                                      # row.names = NULL, fill = TRUE));
#MONITORED_GENES = Map('[', MONITORED_GENES, lapply(MONITORED_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

## Genes monitored individually (tsv file, one column for each group of genes)
#CHOSEN_MARKER_GENES = read.table( file.path( PATH_EXPERIMENT_REFERENCE,
 #                                                 "03_MonitoredGenes",
  #                                                "ChosenMarkerGenes.csv"),
   #                                    sep = ",",
    #                                   header = TRUE,
     #                                  stringsAsFactors = FALSE,
      #                                 row.names = NULL, fill = TRUE);
#CHOSEN_MARKER_GENES = CHOSEN_MARKER_GENES$Gene

## Genes groups to analyze (csv file, one column for each group of genes)
#GROUP_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
 #                                             "04_GeneGroups",
  #                                            "GeneGroups.csv"),
   #                                sep = ",", quote = '"',
    #                               header = TRUE,
     #                              stringsAsFactors = FALSE,
      #                             row.names = NULL, fill = TRUE));
#GROUP_GENES = Map('[', GROUP_GENES, lapply(GROUP_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings


## Genes monitored as modules (tsv file, one column for each group of genes)
#MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
 #                                               "05_Modules",
  #                                              "Modules.csv"),
   #                                     sep = ",",
    #                                    header = TRUE,
     #                                   stringsAsFactors = FALSE,
      #                                  row.names = NULL, fill = TRUE));
#MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
#MODULES_GENES[[ "S_PHASE"]] = CELL_CYCLE_SPHASE_GENELIST
#MODULES_GENES[[ "G2M_PHASE"]] = CELL_CYCLE_G2MPHASE_GENELIST

