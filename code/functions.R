
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
    } else{
        twas_datasets = NULL
        dge_datasets = NULL
    }
    return(list('TWAS'=twas_datasets, 'DGE'=dge_datasets))
}