###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples and several analysis steps
###############################################################################


#### General
#### Variables that contains general information on the experiment

# The experiment name this file is in
# It must be the same name as the name of the experiment folder you defined
# Example : "10x_190712_m_moFluMemB"

EXPERIMENT_NAME = basename( dirname( dirname( WORKING_DIR)))


#### Automatic definition of path
#### You may not modify the following section where standard path are defined in
#### constants so to be used in your code. Have a look at the declared variable
#### and use it in your code.

#### Input / Output

SAMPLE_KO_MCSF = "01_Cellranger/mm10/outs/per_sample_outs/KO_MCSF/sample_filtered_feature_bc_matrix"
SAMPLE_KO_NT = "01_Cellranger/mm10/outs/per_sample_outs/KO_NT/sample_filtered_feature_bc_matrix"
SAMPLE_WT_MCSF = "01_Cellranger/mm10/outs/per_sample_outs/WT_MCSF/sample_filtered_feature_bc_matrix"
SAMPLE_WT_NT = "01_Cellranger/mm10/outs/per_sample_outs/WT_NT/sample_filtered_feature_bc_matrix"



# This is the path of the current experiment this file is in
PATH_EXPERIMENT = file.path( PATH_PROJECT, EXPERIMENT_NAME)

# Those are the path to the main folder : RAW data, REFERENCE data and Output of analysis
PATH_EXPERIMENT_RAWDATA   = file.path( PATH_EXPERIMENT, "00_Rawdata")
PATH_EXPERIMENT_REFERENCE = file.path( PATH_EXPERIMENT, "01_Reference")
PATH_EXPERIMENT_OUTPUT    = file.path( PATH_EXPERIMENT, "05_OUTPUT")


# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME
                            )



#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



