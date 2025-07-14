##################################################____load packages____####################################################
library(dplyr)
library(ACAT)
library(ggplot2)
library(ggrepel)
library(readxl)

source('~/Rwork/covid19_drug_reposition/functions.R')

##################################################____define function____##################################################
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
    
    return(focused_cells)
}


##################################################____do analysis____######################################################

analyze_rs = function(phenotype_id = 1, method_id = 2){
    phenotype = c('T2D', 'Covid19')[phenotype_id]
    method = c('TReD', 'cmap_score')[method_id]
    work_dir = 'your_work_path/covid19_drug_reposition/compound_d_immune_equal_step/'
    annotation_df = xy_read("~/LINCS2_beta/annotation_df.txt")
    siginfo_beta = xy_read('your_LINCS_path/siginfo_beta.txt')
    datasets = as.character(unlist(get_run_datasets(phenotype)))
    
    # 'immune' means results of focused cell lines
    all_sub_immune_drug_instances_rs = list()
    all_sub_bind_drug_instances_rs = list()
    all_sub_immune_drug_rs_cell = list()
    for (i in 1:length(datasets)){
        dataset = datasets[i]
        print(dataset)
        
        focused_cells = getFocusedCellLines(dataset)
        rs_dir = paste0(work_dir, phenotype, '/', method)
        
        # get instances of focused cell lines ('immune')
        print('Binding all drug instances together...')
        immune_instances_rs = list()
        for (j in 0:99){
            rs_name = paste0(rs_dir, '/', dataset, '/d_', j, '.txt')
            if (file.exists(rs_name)){
                temp_rs = xy_read(rs_name)
                cell_name = as.character(sapply(temp_rs$pert_id, function(x) strsplit(x, "[_]")[[1]][2]))
                temp_rs = temp_rs[which(toupper(cell_name) %in% toupper(focused_cells)), ]
                if (nrow(temp_rs) == 0){
                    next
                }
                immune_instances_rs[[j+1]] = temp_rs
            }
        }
        immune_instances_rs = do.call(rbind, immune_instances_rs)
        
        immune_instances_rs = immune_instances_rs[, -which(colnames(immune_instances_rs) == 'V1')]
        colnames(immune_instances_rs)[which(colnames(immune_instances_rs) == 'pert_id')] = 'sig_id'
        immune_instances_rs = merge(immune_instances_rs, siginfo_beta[, c('sig_id', 'pert_id')], by = 'sig_id')
        # immune_instances_rs$permutation_p_immune_instances_adj = p.adjust(immune_instances_rs$permutation_p)
        immune_drug_instances_rs = merge(immune_instances_rs, unique(annotation_df[, c('pert_id', 'common_name')]))
        print(paste0('drug_instances: ', nrow(immune_drug_instances_rs)))
        if (!dir.exists(paste0(rs_dir, '/focused_cell_line_drug_instances/'))){
            dir.create(paste0(rs_dir, '/focused_cell_line_drug_instances/'))
        }
        write.table(immune_drug_instances_rs, file = paste0(rs_dir, '/focused_cell_line_drug_instances/', dataset, '.csv'), row.names = F, sep = ',')
        print('OK!')
        
        # immune_drug_instances_rs$permutation_p_immune_drug_instances_adj = p.adjust(immune_drug_instances_rs$permutation_p)
        print('Extracting the significant compound instances...')
        if (method == 'TReD'){
            ## subset d>mean+2sd & p<0.05
            sub_immune_drug_instances_rs = immune_drug_instances_rs[which(immune_drug_instances_rs$drug_d > (mean(immune_drug_instances_rs$drug_d) + 2*sd(immune_drug_instances_rs$drug_d))), ]
            sub_immune_drug_instances_rs = sub_immune_drug_instances_rs[which(sub_immune_drug_instances_rs$permutation_p < 0.05), ]
        } else{
            ## subset d<0 & p<0.05
            sub_immune_drug_instances_rs = immune_drug_instances_rs[which(immune_drug_instances_rs$drug_d < 0), ]
            sub_immune_drug_instances_rs = sub_immune_drug_instances_rs[which(sub_immune_drug_instances_rs$permutation_p < 0.05), ]
        }
        
        print(paste0('OK! No. extracted instances: ', dim(sub_immune_drug_instances_rs)[1]))
        if (nrow(sub_immune_drug_instances_rs) == 0){
            next
        }
        
        ## match pert_idose, pert_iname, cell_iname
        sub_immune_drug_instances_rs = merge(sub_immune_drug_instances_rs, siginfo_beta[, c('sig_id', 'pert_idose', 'pert_itime', 'cell_iname')], by = 'sig_id')
        
        ## whether effecticive under multi-cell or multi-idose
        print('Annotating multi-cell, multi-dose...')
        pert_ids = sub_immune_drug_instances_rs$pert_id
        cell_bool = c(); dose_bool = c(); time_bool = c();
        for (pert_id in pert_ids){
            temp_subrs = sub_immune_drug_instances_rs[which(pert_ids == pert_id), ]
            cell_bool = c(cell_bool, length(unique(temp_subrs$cell_iname)) >= 2)
            dose_bool = c(dose_bool, length(unique(temp_subrs$pert_idose)) >= 2)
            time_bool = c(time_bool, length(unique(temp_subrs$pert_itime)) >= 2)
        }
        sub_immune_drug_instances_rs$cell_bool = cell_bool
        sub_immune_drug_instances_rs$dose_bool = dose_bool
        sub_immune_drug_instances_rs$time_bool = time_bool
        
        # the instance with max reversal distance was select to represent the drug
        sub_bind_drug_instances_rs = sub_immune_drug_instances_rs %>% group_by(pert_id) %>% dplyr::slice(which.max(drug_d))
        sub_immune_drug_rs_cell = sub_immune_drug_instances_rs[sub_immune_drug_instances_rs$cell_bool, ] %>% group_by(pert_id) %>% dplyr::slice(which.max(drug_d)) %>% as.data.frame
        sub_immune_drug_rs_dose = sub_immune_drug_instances_rs[sub_immune_drug_instances_rs$dose_bool, ] %>% group_by(pert_id) %>% dplyr::slice(which.max(drug_d)) %>% as.data.frame
        sub_immune_drug_rs_time = sub_immune_drug_instances_rs[sub_immune_drug_instances_rs$time_bool, ] %>% group_by(pert_id) %>% dplyr::slice(which.max(drug_d)) %>% as.data.frame
        
        # save required files
        all_sub_immune_drug_instances_rs[[i]] = sub_immune_drug_instances_rs
        all_sub_bind_drug_instances_rs[[i]] = sub_bind_drug_instances_rs
        all_sub_immune_drug_rs_cell[[i]] = sub_immune_drug_rs_cell
        print('OK!')
    }
    print('Saving all results...')
    all_sub_immune_drug_instances_rs = do.call(rbind, all_sub_immune_drug_instances_rs)
    all_sub_bind_drug_instances_rs = do.call(rbind, all_sub_bind_drug_instances_rs)
    all_sub_immune_drug_rs_cell = do.call(rbind, all_sub_immune_drug_rs_cell)
    
    # sub_annotation_df = annotation_df[match(unique(all_sub_immune_drug_rs$pert_id), annotation_df$pert_id), ]
    # the final significant instances must be significant in more than 2 cell lines
    sub_annotation_df_cell = merge(annotation_df, data.frame(pert_id = unique(all_sub_immune_drug_rs_cell$pert_id)), by = 'pert_id')
    # sub_annotation_df_dose = merge(annotation_df, data.frame(pert_id = unique(all_sub_immune_drug_rs_dose$pert_id)), by = 'pert_id')
    # sub_annotation_df_time = merge(annotation_df, data.frame(pert_id = unique(all_sub_immune_drug_rs_time$pert_id)), by = 'pert_id')
    # write.table(sub_annotation_df_cell, file = 'sub_annotation_df_cell.csv', sep = ',', row.names = F)
    write.table(all_sub_immune_drug_instances_rs, file = paste0(rs_dir, '/all_sub_immune_drug_instances_rs.csv'), sep = ',', row.names = F)
    write.table(all_sub_bind_drug_instances_rs, file = paste0(rs_dir, '/all_sub_bind_drug_instances_rs.csv'), sep = ',', row.names = F)
    write.table(all_sub_immune_drug_rs_cell, file = paste0(rs_dir, '/all_sub_immune_drug_rs_cell.csv'), sep = ',', row.names = F)
    # write.table(all_sub_immune_drug_rs_dose, file = paste0(work_dir, 'all_sub_immune_drug_rs_dose.csv'), sep = ',', row.names = F)
    # write.table(all_sub_immune_drug_rs_time, file = paste0(work_dir, 'all_sub_immune_drug_rs_time.csv'), sep = ',', row.names = F)
    print('OK!')
}

