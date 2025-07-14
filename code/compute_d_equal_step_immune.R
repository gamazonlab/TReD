
##################################################____load packages____####################################################
library(readxl)
library(foreach)
library(doMC)
registerDoMC(cores=max(detectCores() - 1, 1))

source('~/code/functions.R')


##################################################____define function____##################################################

### compute_d: compute d between significant gene and pertubations
## twas_sig_rs: a dataframe consisting of twas sig genes' geneid and rank, this rank is a ratio which represents the rank
## drug_sig_rs: a dataframe of drug same as twas_sig_rs
## return: reversal distance, d
compute_d = function(twas_sig_rs, drug_sig_rs, method){
    if (method == 'TReD'){
        twas_rank = twas_sig_rs$twas_rank - 0.5 # vector of twas rank
        drug_rank = drug_sig_rs$drug_rank - 0.5 # vector of drug rank
        twas_rank_norm = sqrt(sum(twas_rank**2)) # norm of twas rank
        drug_rank_norm = sqrt(sum(drug_rank**2)) # norm of drug rank
        cosine_seta = sum(twas_rank*drug_rank) / (twas_rank_norm*drug_rank_norm) # cos(seta)
        d = drug_rank_norm*(-cosine_seta) # the d we need
        return(d)
    } else{
        ### calculates the reversal score between the input signature and the drug profile
        sig_up = subset(twas_sig_rs, twas_rank < 0.5, select="GeneID")
        sig_down = subset(twas_sig_rs, twas_rank >= 0.5, select="GeneID")
        drug_signature = drug_sig_rs
        
        num_genes <- nrow(drug_signature)
        ks_up <- 0
        ks_down <- 0
        connectivity_score <- 0
        
        drug_signature[, "rank"] <- rank(drug_signature[, "drug_rank"])
        up_tags_rank <- merge(drug_signature, sig_up, by.x = "GeneID", by.y = 1)
        down_tags_rank <- merge(drug_signature, sig_down, by.x = "GeneID", by.y = 1)
        up_tags_position <- sort(up_tags_rank$rank)
        down_tags_position <- sort(down_tags_rank$rank)
        num_tags_up <- length(up_tags_position)
        num_tags_down <- length(down_tags_position)
        
        if(num_tags_up > 1) {
            a_up <- 0; b_up <- 0;
            a_up <- max(sapply(1:num_tags_up, function(j) {
                j/num_tags_up - up_tags_position[j]/num_genes
            }))
            b_up <- max(sapply(1:num_tags_up, function(j) {
                up_tags_position[j]/num_genes - (j-1)/num_tags_up
            }))
            
            if(a_up > b_up) {
                ks_up <- a_up
            } else {
                ks_up <- -b_up
            }
        }else{
            ks_up <- 0
        }
        
        if (num_tags_down > 1){
            a_down <- 0; b_down <- 0;
            a_down <- max(sapply(1:num_tags_down,function(j) {
                j/num_tags_down - down_tags_position[j]/num_genes
            }))
            b_down <- max(sapply(1:num_tags_down,function(j) {
                down_tags_position[j]/num_genes - (j-1)/num_tags_down
            }))
            
            if(a_down > b_down) {
                ks_down <- a_down
            } else {
                ks_down <- -b_down
            }
        }else{
            ks_down <- 0
        }
        
        if (ks_up == 0 & ks_down != 0){ #only down gene inputed
            connectivity_score <- -ks_down
        }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
            connectivity_score <- ks_up
        }else if (sum(sign(c(ks_down,ks_up))) == 0) {
            connectivity_score <- ks_up - ks_down # different signs
        }
        return(connectivity_score)
    }
}


### pipeline: perform computing
## twas_sig_rs: a dataframe consisting of twas sig genes' geneid and rank, this rank is a ratio which represents the rank
## drug_sig_rs: a dataframe of drug same as twas_sig_rs
## return: results
pipeline = function(twas_sig_gene, compound_signatures, gene_list, method, dataset){
    ### start computing rand_drug_d  
    print('start computing rand_drug_d')
    N_PERMUTATIONS <- 1000
    rand_drug_d = foreach(drug_id = 1:ncol(compound_signatures)) %dopar% {
        ## get drug gene rank
        drug_rs = as.data.frame(subset(compound_signatures, select = drug_id))
        drug_rs[1] = rank(-drug_rs[1]) / (dim(drug_rs)[1])  #uniformly distributed
        drug_rs = cbind(gene_list, drug_rs)
        colnames(drug_rs) = c('GeneID', 'drug_rank')
        temp_sig_rs = merge(twas_sig_gene, drug_rs, by = 'GeneID')
        temp_sig_rs$drug_rank = rank(temp_sig_rs$drug_rank) / nrow(temp_sig_rs)
        
        ## get null distribution
        # print(paste0('Now we are computing null distribution of ', drug_id))
        null_drug_d_rand1000 = sapply(1:N_PERMUTATIONS, function(i){
            set.seed(i)
            temp_sig_rs$drug_rank = sample(temp_sig_rs$drug_rank, length(temp_sig_rs$drug_rank))
            twas_sig_rs = subset(temp_sig_rs, select = c('GeneID', 'twas_rank'))
            drug_sig_rs = subset(temp_sig_rs, select = c('GeneID', 'drug_rank'))
            
            return(compute_d(twas_sig_rs, drug_sig_rs, method))
        })
        return(null_drug_d_rand1000)
    }
    
    # ### start computing drug_d
    print('start computing drug_d')
    drug_d = sapply(1:ncol(compound_signatures), function(drug_id){
        # print(paste0('Now we are computing d of ', drug_id))
        drug_rs = as.data.frame(subset(compound_signatures, select = drug_id))
        drug_rs[1] = rank(-drug_rs[1]) / (dim(drug_rs)[1])
        drug_rs = cbind(gene_list, drug_rs)
        colnames(drug_rs) = c('GeneID', 'drug_rank')
        temp_sig_rs = merge(twas_sig_gene, drug_rs, by = 'GeneID')
        temp_sig_rs$drug_rank = rank(temp_sig_rs$drug_rank)/nrow(temp_sig_rs)

        # set.seed(2020)
        # temp_sig_rs$drug_rank = temp_sig_rs[sample(nrow(temp_sig_rs), nrow(temp_sig_rs), replace = F),'drug_rank']

        twas_sig_rs = subset(temp_sig_rs, select = c('GeneID', 'twas_rank'))
        drug_sig_rs = subset(temp_sig_rs, select = c('GeneID', 'drug_rank'))
        return(compute_d(twas_sig_rs, drug_sig_rs, method))
    })

    result = data.frame(
        pert_id = colnames(compound_signatures),
        drug_d = drug_d,
        data_set = dataset
    )
    if (method == 'TReD'){
        result$permutation_p = sapply(1:ncol(compound_signatures), function(i){
            d = result$drug_d[i]
            return(length(which(rand_drug_d[[i]] >= d)) / length(rand_drug_d[[i]]))
        })
        
        all_rand_drud_d = unlist(rand_drug_d)
        result$experiment_wise_permutation_p = sapply(1:ncol(compound_signatures), function(i){
            d = result$drug_d[i]
            return(length(which(all_rand_drud_d >= d)) / length(all_rand_drud_d))
        })
    } else{
        result$permutation_p = sapply(1:ncol(compound_signatures), function(i){
            d = result$drug_d[i]
            return(length(which(abs(rand_drug_d[[i]]) >= abs(d))) / length(rand_drug_d[[i]]))
        })
        
        all_rand_drud_d = unlist(rand_drug_d)
        result$experiment_wise_permutation_p = sapply(1:ncol(compound_signatures), function(i){
            d = result$drug_d[i]
            return(length(which(abs(all_rand_drud_d) >= abs(d))) / length(all_rand_drud_d))
        })
    }
    

    return(result)
}


main_fun = function(args, phenotype_id=1, method_id=1){
    phenotype = c('T2D', 'Covid19')[phenotype_id]
    method = c('TReD', 'cmap_score')[method_id]
    work_path = 'your_work_path/covid19_drug_reposition/'
    
    ## load compound_signatures
    gene_list = xy_read('your_LINCS_path/parse_gctx/row.csv')
    # please download LINCS data and processed it according to /read_gctx
    compound_signatures = xy_read(paste0('your_LINCS_path/parse_gctx/exp_mat_', args, '.csv'))
    row.names(compound_signatures) = gene_list$rid
    
    ## load cell name and select cell lines
    # cell_name = as.character(sapply(colnames(compound_signatures), function(x) strsplit(x, "[_]")[[1]][2]))
    # cell_focused = as.data.frame(read_excel('/data/g_gamazon_lab/zhoud2/covid19_drug_reposition/TWAS_rs/T2D/T2D_cmap_cellli_annot.xlsx'))
    # compound_signatures = compound_signatures[, which(toupper(cell_name) %in% toupper(cell_focused$`cell line name`))]
    
    ## datasets
    datasets = get_run_datasets(phenotype)
    all_data_sets = c(datasets$TWAS, datasets$DGE)
    
    for (dataset_index in 1:length(all_data_sets)){
        dataset = all_data_sets[dataset_index]
        print(dataset)
        
        output = paste0(work_path, 'compound_d_immune_equal_step/', phenotype, '/', method, '/', dataset, '/')
        if(!dir.exists(output)){
          dir.create(output)
        }
        rs_save_path = paste0(output, '/d_', args, '.txt')
        if (!file.exists(rs_save_path)){
            if (dataset %in% datasets$TWAS){
                load(paste0(work_path, 'disease_signature/', phenotype, '/processed/twas_', dataset, '_new.Rdata'))
                twas_rs = as.data.frame(twas_rs)
                
                # remove duplicated genes
                twas_rs = twas_rs[order(twas_rs$padj), ]
                twas_rs = twas_rs[!duplicated(twas_rs$gene), ]
                
                twas_rs$twas_rank = rank(-twas_rs$zscore) / (nrow(twas_rs)) 
                twas_sig_gene = subset(twas_rs, padj < 0.05, select = c('GeneID', 'twas_rank'))
            } else if (dataset %in% c('ALV', 'EXP', 'BALF')){
                ## load dge sig gene list, to reuse twas's code, so named twas_sig_gene
                load(paste0(work_path, 'disease_signature/', phenotype, '/processed/dge_', dataset, '_new.Rdata'))
                dge_rs = as.data.frame(dge_rs)
                
                if (dataset == "BALF"){
                    gene_name_col = 'symbol'
                    fold_change_col = 'log2FoldChange'
                } else{
                    gene_name_col = 'GeneName'
                    fold_change_col = 'log2FoldChange'
                }
                
                # remove duplicated genes
                dge_rs = dge_rs[order(abs(dge_rs[, fold_change_col]), decreasing = T), ]
                dge_rs = dge_rs[!duplicated(dge_rs[, gene_name_col]), ]
                
                # use up to 50 gene for positive and negative
                # positive = dge_rs[which(dge_rs[, fold_change_col]>0), ]
                # negative = dge_rs[which(dge_rs[, fold_change_col]<0), ]
                # min_size = min(nrow(positive), nrow(negative), 50)
                # dge_rs = rbind(positive[1:min_size, ], negative[1:min_size, ])
                # dge_rs = dge_rs[1:min(100, nrow(dge_rs)), ]
                
                twas_sig_gene = subset(dge_rs, select = c('GeneID', 'dge_rank'))
                colnames(twas_sig_gene)[2] = 'twas_rank'
            } else{
                load(paste0(work_path, 'disease_signature/', phenotype, '/processed/dge_', dataset, '_new.Rdata'))
                if (dataset == 'islets'){
                    dge_rs$twas_rank = rank(-as.numeric(dge_rs$`Fold D/N`)) / (nrow(dge_rs)) 
                    twas_sig_gene = subset(dge_rs, select = c('GeneID', 'twas_rank'))
                } else{
                    dge_rs$twas_rank = rank(-as.numeric(dge_rs$Ratio)) / (nrow(dge_rs)) 
                    twas_sig_gene = subset(dge_rs, select = c('GeneID', 'twas_rank'))
                }
            }
            
            # scale the rank to equal step series
            twas_sig_gene$twas_rank = rank(twas_sig_gene$twas_rank) / nrow(twas_sig_gene)
            
            result = pipeline(twas_sig_gene, compound_signatures, gene_list, method, dataset)
            print('Computing is over!!!')
            write.table(result, file = rs_save_path, sep = '\t')
            print('Computing d is OK!')
        }
    }
}


##################################################____run function____##################################################
args = as.numeric(commandArgs(TRUE))
main_fun(args, phenotype_id=2, method_id=2)


