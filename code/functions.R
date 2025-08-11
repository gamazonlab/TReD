
xy_read = function(file_name){
    return(as.data.frame(data.table::fread(file_name)))
}


get_run_datasets = function(phenotype){
    if (phenotype == 'T2D'){
        twas_datasets = c('Adipose_Visceral_Omentum', 'Muscle_Skeletal', 'Liver', 'Pancreas')
        dge_datasets = c('islets', 'myoblasts', 'myotubes')
    } else if (phenotype == 'Covid19'){
        twas_datasets = c('blood', 'lung', 'lymphocytes', 'spleen')
        dge_datasets = c('ALV', 'EXP', 'BALF')
    } else if (phenotype == 'AD'){
        twas_datasets = c(
            "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", 
            "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", 
            "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", 
            "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra" 
        )
        dge_datasets = NULL
    } else{
        twas_datasets = NULL
        dge_datasets = NULL
    }
    return(list('TWAS'=twas_datasets, 'DGE'=dge_datasets))
}


getFocusedCellLines = function(dataset){
    focused_cells = c()
    if (dataset == 'Liver'){
        focused_cells = c('JHH7', 'SNU761', 'JHH5', 'HUH7', 'HEPG2', 'PHH', 'HUH751')
    } 
    
    if (dataset == 'Adipose_Visceral_Omentum'){
        focused_cells = c('ASC')
    }
    
    if (dataset == 'Muscle_Skeletal' | dataset == 'myoblasts' | dataset == 'myotubes'){
        focused_cells = c('A204', 'SKB')
    }
    
    if (dataset == 'Pancreas' | dataset == 'islets'){
        focused_cells = c('YAPC', 'DANG')
    }
    
    if (dataset %in% c('blood', 'lung', 'lymphocytes', 'spleen', 'ALV', 'EXP', 'BALF')){
        focused_cells = c(
            "KMS34", "U266B1", "OCILY3", "U937", "K562", "SUDHL4", "OCILY19", "MINO", "NALM6",
            "JURKAT", "WSUDLCL2", "OCILY10", "MICROGLIA-PSEN1", "BJAB", "TMD8", "HBL1"
        )
    }
    
    if (dataset %in% get_run_datasets('AD')$TWAS){
        focused_cells = c('SHSY5Y', 'SKNSH', 'LN229', 'U251MG', 'YH13', 'GI1', 'NEU', 'HIMG001', 'HIMG002')
    }
    
    return(focused_cells)
}
