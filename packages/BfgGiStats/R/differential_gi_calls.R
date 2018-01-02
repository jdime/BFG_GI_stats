differential_gi_analysis <- function(gi_data,
                                     fdr_cutoff = 0.05,
                                     require_sign_change = F,
                                     nn_pair_type = 'broad'){
  gi_data <-
    dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')
  
  
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
  
  ddr_pairs <-
    gi_data$Type_of_gene_i == 'DNA_repair' &
    gi_data$Type_of_gene_j == 'DNA_repair'
  
  
  z_cols <- grep('Class',grep('Z_GIS',colnames(gi_data),val=T),invert=T,val=T)
  
  
  conditions <- sapply(strsplit(grep('^GIS',colnames(gi_data),val=T),split='\\.'),function(x){x[2]})
  
  
  ret_list <- list()
  ret_df <- c()
  
  for(condition1 in conditions) {
    ret_list[[condition1]] <- list()
    for (condition2 in conditions) {
      if (condition1 > condition2) {
        
        z_name1 <- paste(c('Z_GIS_ij', condition1), collapse = '.')
        z_name2 <- paste(c('Z_GIS_ij', condition2), collapse = '.')
        
        z_class_name1 <- paste(c(z_name1,'Class'),collapse='_')
        z_class_name2 <- paste(c(z_name2,'Class'),collapse='_')
        
        nn_scores_cond <-
          gi_data[nn_pairs, z_name1] - gi_data[nn_pairs, z_name2]
        
        non_nn_scores_cond <-
          gi_data[ddr_pairs, z_name1] - gi_data[ddr_pairs, z_name2]
        
        unequal_tests <- c(`>=`, `<=`)
        precision_list <- lapply(unequal_tests, function(unequal_test) {
          mu <- mean(nn_scores_cond)
          sigma <- sd(nn_scores_cond)
          #sapply(1:length(non_nn_scores_cond), function(i) {
          sapply(1:length(non_nn_scores_cond), function(i) {
            observed <-
              sum(unequal_test(non_nn_scores_cond, non_nn_scores_cond[i]))
            if (all.equal(unequal_test, `<=`) == T) {
              use_lower_tail <- T
            } else{
              use_lower_tail <- F
            }
            
            #Using normal distribuion to calculate expected
            expected <-
              pnorm(
                non_nn_scores_cond[i],
                mean = mu,
                sd = sigma,
                lower.tail = use_lower_tail
              )
            
            
            #expected <- sum(unequal_test(nn_scores_cond, non_nn_scores_cond[i]))/length(non_nn_scores_cond)
            expected <- expected * length(non_nn_scores_cond)
            
            fdr <- (expected / observed)
            
            #No sense returning fdr estimates >100%
            return(min(fdr, 1))
          })
        })
        fdrs_pos <- precision_list[[1]]
        fdrs_neg <- precision_list[[2]]
        comb_fdr <- cbind(fdrs_pos, fdrs_neg)
        
        #Return positive or negative FDRs based on the nominal sign
        #of the interaction
        fdr_scores <- sapply(1:nrow(comb_fdr), function(i) {
          if (non_nn_scores_cond[i] < 0) {
            return(fdrs_neg[i])
          } else{
            return(fdrs_pos[i])
          }
        })
        
        if(require_sign_change){
          sig <- fdr_scores < fdr_cutoff & gi_data[ddr_pairs,z_class_name1] != gi_data[ddr_pairs,z_class_name2]
        }else{
          sig <- fdr_scores < fdr_cutoff
        }

        ddr_data <- gi_data[ddr_pairs,]
        
        sig_names <- apply(ddr_data[which(sig),1:2],1,function(x){paste(x,collapse='_')})
        
        ret_list[[condition1]][[condition2]] <- sig_names
        
        
        
        retvec <- cbind(ddr_data[which(sig),c("Barcode_i","Barcode_j")],
                        rep(condition1,sum(sig)),
                        rep(condition2,sum(sig)),
                        ddr_data[which(sig),c(z_name1,z_name2)],
                        ddr_data[which(sig),c(z_class_name1,z_class_name2)],
                        non_nn_scores_cond[which(sig)],
                        fdr_scores[which(sig)])
        
        colnames(retvec) <- c('Barcode_i','Barcode_j','Condition1','Condition2','Z_Condition1','Z_Condition2','Class_Condition1','Class_Condition2','DeltaZ','DeltaZ_FDR')
        
        ret_df <- rbind(ret_df,retvec)
      }
    }
  }
  return(ret_df)
}


