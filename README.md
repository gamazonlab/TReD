
# TReD
In this study, we proposed a novel drug reposition framework based on cosine similarity called “TReD” (Transcriptome-informed Reversal Distance) that embeds the disease signatures and drug response profiles into a high-dimensional normed space to quantify the reversal potential of candidate drugs in a disease-related cell-based screening.  
<img width="662" height="637" alt="image" src="https://github.com/user-attachments/assets/b8aee8df-2e52-4a3f-aec1-1e65007f6216" />


# Input Files
## Drug signature
### LINCS
Please download the compound signatures from LINCS portal:  
wget https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx  

## Population Disease signature
### TWAS
1. Covid19  
blood: data/disease_signature/Covid19/processed/twas_blood_new.Rdata
lung: data/disease_signature/Covid19/processed/twas_lung_new.Rdata  
lymphocytes: data/disease_signature/Covid19/processed/twas_lymphocytes_new.Rdata  
spleen: data/disease_signature/Covid19/processed/twas_spleen_new.Rdata

2. T2D  
Adipose_Visceral_Omentum: data/disease_signature/T2D/processed/twas_Adipose_Visceral_Omentum_new.Rdata  
Muscle_Skeletal: data/disease_signature/T2D/processed/twas_Adipose_Visceral_Omentum_new.Rdata  
Liver: data/disease_signature/T2D/processed/twas_Adipose_Visceral_Omentum_new.Rdata  
Pancreas: data/disease_signature/T2D/processed/twas_Adipose_Visceral_Omentum_new.Rdata  

3. AD  
Brain_Amygdala: data/disease_signature/T2D/processed/twas_Brain_Amygdala_new.Rdata  
Brain_Anterior_cingulate_cortex_BA24: data/disease_signature/T2D/processed/twas_Brain_Anterior_cingulate_cortex_BA24_new.Rdata  
Brain_Caudate_basal_ganglia: data/disease_signature/T2D/processed/twas_Brain_Caudate_basal_ganglia_new.Rdata  
Brain_Cerebellar_Hemisphere: data/disease_signature/T2D/processed/twas_Brain_Cerebellar_Hemisphere_new.Rdata  
Brain_Cerebellum: data/disease_signature/T2D/processed/twas_Brain_Cerebellum_new.Rdata  
Brain_Cortex: data/disease_signature/T2D/processed/twas_Brain_Cortex_new.Rdata  
Brain_Frontal_Cortex_BA9: data/disease_signature/T2D/processed/twas_Brain_Frontal_Cortex_BA9_new.Rdata  
Brain_Hippocampus: data/disease_signature/T2D/processed/twas_Brain_Hippocampus_new.Rdata  
Brain_Hypothalamus: data/disease_signature/T2D/processed/twas_Brain_Hypothalamus_new.Rdata  
Brain_Nucleus_accumbens_basal_ganglia: data/disease_signature/T2D/processed/twas_Brain_Nucleus_accumbens_basal_ganglia_new.Rdata  
Brain_Putamen_basal_ganglia: data/disease_signature/T2D/processed/twas_Brain_Putamen_basal_ganglia_new.Rdata  
Brain_Spinal_cord_cervical_c-1: data/disease_signature/T2D/processed/twas_Brain_Spinal_cord_cervical_c-1_new.Rdata  
Brain_Substantia_nigra: data/disease_signature/T2D/processed/twas_Brain_Substantia_nigra_new.Rdata  

### TWAS Tool
Predixcan (https://www.nature.com/articles/ng.3367, https://github.com/gamazonlab/MR-JTI/blob/master/model_training/predixcan/)

### DGE
1. Covid19  
ALV: data/disease_signature/Covid19/processed/dge_ALV_new.Rdata  
EXP: data/disease_signature/Covid19/processed/dge_EXP_new.Rdata  
BALF: data/disease_signature/Covid19/processed/dge_BALF_new.Rdata

2. T2D  
islets: data/disease_signature/T2D/processed/dge_islets_new.Rdata  
myoblasts: data/disease_signature/T2D/processed/dge_myoblasts_new.Rdata  
myotubes: data/disease_signature/T2D/processed/dge_myotubes_new.Rdata  

# Script
1. code/functions.R: providing shared functions
2. code/compute_d_equal_step_immune.R: calculating the reversal distance and performing the permutation test for each perturbation
3. code/rs_analysis_equal_step_immune.R: analyzing the screening results and searching for the potential drug for the disease
4. code/read_gctx/read_cmap_gctx.py: parsing the GCTX file and splitting the huge dataset into 100 small subsets
5. code/colocalization/covid19_result_for_coloc_leadsnp_200kbwindow.R: performing the colocalization analysis for Covid19
6. code/predixcan/predixcan_r.r: script for running the PrediXcan



