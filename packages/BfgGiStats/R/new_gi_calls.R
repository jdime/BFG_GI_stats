devtools::use_package('dplyr')


#' Adds FDR columns to the genetic interaction data
#'
#' @param gi_data input genetic interaction table
#' @param nn_pair_type "null distribution" of GI scores used to
#' compute FDR.Either 'broad' for all interactions with neutral
#' pairs or 'narrow' for only neutral-neutral pairs
#'
#' @return gi_data with 'FDR.internal' metrics added
add_fdrs <- function(gi_data,
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
  
  
  fdr_matrix <- c()
  
  for(condition in colnames(z_scores)){
    nn_scores_cond <-  nn_scores[,condition]
    non_nn_scores_cond <- non_nn_scores[,condition]
    non_nn_scores_cond_full <- z_scores_full[,condition]
    
    
    #Hacky way to make FDRs for both positive and negative
    #in a loop
    unequal_tests <- c(`>=`,`<=`)
    
    precision_list <- lapply(unequal_tests,function(unequal_test){
      mu <- mean(nn_scores_cond)
      sigma <- sd(nn_scores_cond)
      sapply(1:length(non_nn_scores_cond_full),function(i){
        observed <- sum(unequal_test(non_nn_scores_cond_full, non_nn_scores_cond_full[i]))
        if(all.equal(unequal_test,`<=`) == T){
          use_lower_tail <- T
        }else{
          use_lower_tail <- F
        }
        
        #Using normal distribuion to calculate expected
        expected <- pnorm(non_nn_scores_cond_full[i],mean=mu,sd=sigma,lower.tail=use_lower_tail)
        expected <- expected*length(non_nn_scores_cond_full)
        fdr <- (expected/observed)
        
        #No sense returning fdr estimates >100%
        return(min(fdr,1))
      })
    })
    
    fdrs_pos <- precision_list[[1]]
    fdrs_neg <- precision_list[[2]]
      comb_fdr <- cbind(fdrs_pos,fdrs_neg)
    
    #Return positive or negative FDRs based on the nominal sign
    #of the interaction
    fdr_scores <- sapply(1:nrow(comb_fdr),function(i){
      if(z_scores_full[i,condition] < 0){
        return(fdrs_neg[i])
      }else{
        return(fdrs_pos[i])
      }
    })
    

    fdr_matrix <- cbind(fdr_matrix,fdr_scores)
  }
  
  #Add the appropriate columns at the end and fill in the data
  drugs <- sapply(colnames(z_scores),function(name){
    strsplit(name,split='\\.')[[1]][2]
  })
  fdr_names <- sapply(drugs,function(drug){
    paste(c('FDR.Internal_ij.',drug),collapse='')
  })
  
  colnames(fdr_matrix) <- fdr_names
  
  gi_data_full[,fdr_names] <- fdr_matrix
  
  return(gi_data_full)
  
}

#' Update genetic interaction calls based on either Z metrics or internal FDR
#'
#' @param gi_data input genetic interaction table
#' @param z_score_cols which columns correspond to Z scores - default is to use grep to find automatically
#' @param z_class_cols which columns correspond to Z-based classifications - default is to use grep to find automatically
#' @param fdr_cols which columns correspond to internal FDRs - default is to use grep to find automatically
#' @param fdr_cutoff internal FDR cutoff to make GI calls - defaults to 0.05
#' @param use_z use Z scores instead of FDR to make GI calls? defaults to False
#' @param z_cutoff_neg if use_z is True, Z-cutoff for negative interactions
#' @param z_cutoff_pos  if use_z is True, Z-cutoff for positive interactions
#'
#' @return a version of gi_data with GI calls updated
update_calls <- function(gi_data,
                         z_score_cols = NULL,
                         fdr_cols = NULL,
                         fdr_cutoff = 0.05,
                         use_z = F,
                         z_cutoff_neg = NULL,
                         z_cutoff_pos = NULL){
  if(is.null(z_score_cols)){
    z_score_cols <-
      grep('Class',
           grep('^Z_GIS', colnames(gi_data), val = T),
           invert = T,
           val = T) 
  }
  
  z_scores <- gi_data[ , z_score_cols]
  colnames(z_scores) <-
    sapply(colnames(z_scores), function(name) {
      sprintf('%s_Class', name)
    })
  
  gi_data <- cbind(gi_data, z_scores)
  
  #if(is.null(z_class_cols)){
  z_class_cols <-
    grep('Class',
         grep('^Z_GIS', colnames(gi_data), val = T),
         invert = F,
         val = T)  
  #}
  if(is.null(fdr_cols)){
    fdr_cols <- grep('^FDR', colnames(gi_data), val = T)
  }
  for(i in 1:length(z_score_cols)){
    gi_data[,z_class_cols[i]] <- "NEUTRAL"
    
    if(use_z == F){
      neg_crit <- gi_data[,fdr_cols[i]] <= fdr_cutoff & gi_data[,z_score_cols[i]] < 0
      pos_crit <- gi_data[,fdr_cols[i]] <= fdr_cutoff & gi_data[,z_score_cols[i]] > 0
    }else{
      neg_crit <- gi_data[,z_score_cols[i]] < z_cutoff_neg
      pos_crit <- gi_data[,z_score_cols[i]] > z_cutoff_pos
    }
    gi_data[neg_crit,z_class_cols[i]] <- "AGGRAVATING"
    gi_data[pos_crit,z_class_cols[i]] <- "ALLEVIATING"
  }
  
  return(gi_data)
}
