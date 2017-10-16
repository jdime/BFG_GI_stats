devtools::use_package('dplyr')

add_fdrs <- function(gi_data,
                            fdr_positive = 0.05,
                            fdr_negative = 0.05,
                           nn_pair_type = 'broad'){
  gi_data_full <- gi_data
  gi_data <-
    dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')
  
  ddr_pairs <-
    gi_data$Type_of_gene_i == 'DNA_repair' &
    gi_data$Type_of_gene_j == 'DNA_repair'
  
  #Neutral-Neutral pairs
  if(nn_pair_type == 'broad'){
    #Can define all non DNA_repair - DNA_repair as neutral
    nn_pairs <-
    gi_data$Type_of_gene_i != 'DNA_repair' | gi_data$Type_of_gene_j != 'DNA_repair'
  }else if(nn_pair_type == 'narrow'){
    #Or just neutral-neutral pairs
    nn_pairs <-
      gi_data$Type_of_gene_i == 'Neutral' &
      gi_data$Type_of_gene_j == 'Neutral'
  }
  
  z_cols <-
    grep('Class',
         grep('^Z_GIS', colnames(gi_data), val = T),
         invert = T,
         val = T)
  
  z_scores <- gi_data[,z_cols]
  z_scores_full <- gi_data_full[,z_cols]
  
  nn_scores <- z_scores[nn_pairs,]
  non_nn_scores <- z_scores[ddr_pairs,]
  
  #ret_table <- c()
  fdr_matrix <- c()
  
  for(condition in colnames(z_scores)){
    nn_scores_cond <-  nn_scores[,condition]
    non_nn_scores_cond <- non_nn_scores[,condition]
    
    unequal_tests <- c(`>=`,`<=`)
    
    precision_list <- lapply(unequal_tests,function(unequal_test){
      sapply(1:length(non_nn_scores_cond),function(i){
        observed <- sum(unequal_test(non_nn_scores_cond, non_nn_scores_cond[i]))
        expected <- sum(unequal_test(nn_scores_cond, non_nn_scores_cond[i]))
        
        expected <- expected/length(nn_scores_cond)
        expected <- expected*length(non_nn_scores_cond)
        fdr <- (expected/observed)
        fdr <- min(fdr,1)
        return(fdr)
      })
    })
    
    pos_prec_vec <<- precision_list[[1]]
    neg_prec_vec <<- precision_list[[2]]
    
    pos_prec_func <- approxfun(non_nn_scores_cond,pos_prec_vec)
    neg_prec_func <- approxfun(non_nn_scores_cond,neg_prec_vec)
    
    fdrs_pos <- pos_prec_func(z_scores_full[,condition])
    fdrs_neg <- neg_prec_func(z_scores_full[,condition])
    
    #FDR is never 0
    fdrs_neg[fdrs_neg == 0] <- min(fdrs_neg[fdrs_neg > 0],na.rm=T)
    fdrs_pos[fdrs_pos == 0] <- min(fdrs_pos[fdrs_pos > 0],na.rm=T)
    
    #Extreme values from linked pairs get assigned lowest FDR
    fdrs_neg[z_scores_full[,condition] <= min(non_nn_scores_cond)] <- min(fdrs_neg,na.rm=T)
    fdrs_pos[z_scores_full[,condition] >= max(non_nn_scores_cond)] <- min(fdrs_pos,na.rm=T)
  
    #fdr_scores <- fdrs_neg
    comb_fdr <- cbind(fdrs_pos,fdrs_neg)
    
    fdr_scores <- sapply(1:nrow(comb_fdr),function(i){
      if(z_scores_full[i,condition] < 0){
        return(fdrs_neg[i])
      }else{
        return(fdrs_pos[i])
      }
    })
    

    fdr_matrix <- cbind(fdr_matrix,fdr_scores)
    
    
    
    #fdr_cutoff_neg <- max(non_nn_scores_cond[neg_prec_vec < fdr_negative])
    #fdr_cutoff_pos <- min(non_nn_scores_cond[pos_prec_vec < fdr_positive])
    
    #positive_inters <- sum(non_nn_scores_cond > fdr_cutoff_pos)
    #negative_inters <- sum(non_nn_scores_cond < fdr_cutoff_neg)
    
    #ret_table <- cbind(ret_table,c(fdr_cutoff_pos,fdr_cutoff_neg,positive_inters,negative_inters))
  }
  #rownames(ret_table) <- c('Z cutoff Positive',
  #                         'Z cutoff Negative',
  #                         'Positive Interactions at Cutoff',
  #                         'Negative Interactions at Cutoff')
  #ret_table <- as.data.frame(ret_table)
  drugs <- sapply(colnames(z_scores),function(name){
    strsplit(name,split='\\.')[[1]][2]
  })
  
  #colnames(ret_table) <- drugs
  
  fdr_names <- sapply(drugs,function(drug){
    paste(c('FDR.Internal_ij.',drug),collapse='')
  })
  
  colnames(fdr_matrix) <- fdr_names
  
  gi_data_full[,fdr_names] <- fdr_matrix
  
  return(gi_data_full)
  #return(ret_table)
  
}

update_calls <- function(gi_data,
                         z_table,
                         z_score_cols = NULL,
                         z_class_cols = NULL,
                         fdr_cols = NULL,
                         fdr_cutoff = 0.05){
  if(is.null(z_score_cols)){
    z_score_cols <-
      grep('Class',
           grep('^Z_GIS', colnames(gi_data), val = T),
           invert = T,
           val = T) 
  }
  if(is.null(z_class_cols)){
    z_class_cols <-
      grep('Class',
           grep('^Z_GIS', colnames(gi_data), val = T),
           invert = F,
           val = T)  
  }
  if(is.null(fdr_cols)){
    fdr_cols <- grep('^FDR', colnames(gi_data), val = T)
  }
  for(i in 1:length(z_score_cols)){
    gi_data[,z_class_cols[i]] <- "NEUTRAL"
    neg_crit <- gi_data[,fdr_cols[i]] <= fdr_cutoff & gi_data[,z_score_cols[i]] < 0
    pos_crit <- gi_data[,fdr_cols[i]] <= fdr_cutoff & gi_data[,z_score_cols[i]] > 0
    
    gi_data[neg_crit,z_class_cols[i]] <- "AGGRAVATING"
    gi_data[pos_crit,z_class_cols[i]] <- "ALLEVIATING"
  }
  
  return(gi_data)
}
