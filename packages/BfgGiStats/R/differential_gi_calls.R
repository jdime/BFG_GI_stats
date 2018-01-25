#' Differential Genetic Interaction analysis
#'
#' @param gi_data processed genetic interaction data
#' @param fdr_cutoff false discovery rate cutoff for calling significant GI changes
#' @param require_sign_change require different GI classifications to call a differential interaction; True or False
#' @param nn_pair_type null distribution" of GI scores used to
#' compute FDR.  Either 'broad' for all interactions with neutral
#' pairs or 'narrow' for only neutral-neutral pairs
#'
#' @return a data frame with differential genetic interaction comparisons
differential_gi_analysis <- function(gi_data,
                                     fdr_cutoff = 0.05,
                                     delta_gi_cutoff = 0,
                                     require_sign_change = T,
                                     nn_pair_type = 'broad',
                                     make_plots = F){
  gi_data <-
    dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')
  
  
  if(nn_pair_type == 'broad'){
    #Can define all non DNA_repair - DNA_repair as neutral
    nn_pairs <-
      gi_data$Type_of_gene_x != 'DNA_repair' | gi_data$Type_of_gene_y != 'DNA_repair'
  }else if(nn_pair_type == 'narrow'){
    #Or just neutral-neutral pairs
    nn_pairs <-
      gi_data$Type_of_gene_x == 'Neutral' &
      gi_data$Type_of_gene_y == 'Neutral'
  }
  
  ddr_pairs <-
    gi_data$Type_of_gene_x == 'DNA_repair' &
    gi_data$Type_of_gene_y == 'DNA_repair'
  
  
  z_cols <- grep('Class',grep('Z_GIS',colnames(gi_data),val=T),invert=T,val=T)
  
  
  conditions <- sapply(strsplit(grep('^GIS',colnames(gi_data),val=T),split='\\.'),function(x){x[2]})
  
  
  ret_list <- list()
  ret_df <- c()
  
  for(condition1 in conditions) {
    ret_list[[condition1]] <- list()
    for (condition2 in conditions) {
      if (condition1 > condition2) {
        
        
        
        z_name1 <- paste(c('Z_GIS_xy', condition1), collapse = '.')
        z_name2 <- paste(c('Z_GIS_xy', condition2), collapse = '.')
        
        gi_name1 <- paste(c('GIS_xy', condition1), collapse = '.')
        gi_name2 <- paste(c('GIS_xy', condition2), collapse = '.')
        
        gi_err_name1 <- paste(c('SE_GIS_xy', condition1), collapse = '.')
        gi_err_name2 <- paste(c('SE_GIS_xy', condition2), collapse = '.')
        
        
        z_class_name1 <- paste(c(z_name1,'Class'),collapse='_')
        z_class_name2 <- paste(c(z_name2,'Class'),collapse='_')
        
        #nn_pairs <- (gi_data[,z_class_name1] == gi_data[,z_class_name2])
        
        
        #nn_pairs <- (gi_data[,z_class_name1] == 'NEUTRAL' & gi_data[,z_class_name2] == 'NEUTRAL')# | (gi_data$Type_of_gene_i != 'DNA_repair' | gi_data$Type_of_gene_y != 'DNA_repair')
        nn_pairs <- (gi_data$Type_of_gene_x != 'DNA_repair' | gi_data$Type_of_gene_y != 'DNA_repair')
        
        nn_z_scores_cond <-
          (gi_data[nn_pairs, gi_name1] - gi_data[nn_pairs, gi_name2])/sqrt(gi_data[nn_pairs, gi_err_name1]^2 + gi_data[nn_pairs, gi_err_name2]^2)
        nn_gi_scores_cond <-
          gi_data[nn_pairs, gi_name1] - gi_data[nn_pairs, gi_name2]
        
         
        non_nn_z_scores_cond <-
          (gi_data[!nn_pairs, gi_name1] - gi_data[!nn_pairs, gi_name2])/sqrt(gi_data[!nn_pairs, gi_err_name1]^2 + gi_data[!nn_pairs, gi_err_name2]^2)
        non_nn_gi_scores_cond <-
          gi_data[!nn_pairs, gi_name1] - gi_data[!nn_pairs, gi_name2]
        
        
        all_z_scores_cond <-
          (gi_data[, gi_name1] - gi_data[, gi_name2])/sqrt(gi_data[, gi_err_name1]^2 + gi_data[, gi_err_name2]^2)
        all_gi_scores_cond <-
          gi_data[, gi_name1] - gi_data[, gi_name2]
        
        
        
        nn_scores_cond <- nn_z_scores_cond
        non_nn_scores_cond <- non_nn_z_scores_cond
        all_scores_cond <- all_z_scores_cond
        
        ##If plots of process desired
        if(make_plots){
          first_hist_test <-
            hist(
              all_z_scores_cond,
              breaks = 100,
              plot = F
            )
          second_hist_test <-
            hist(
              nn_z_scores_cond,
              breaks = first_hist_test$breaks,
              plot = F
            )
          
          hist(
            non_nn_z_scores_cond,
            breaks = first_hist_test$breaks,
            freq = F,
            border = NA,
            col = rgb(1, 0, 0, 0.5),
            xlab = '∆Z',
            ylab = '',
            main = paste(c(condition1,condition2),collapse=' - '),
            ylim = c(0, max(
              c(first_hist_test$density, second_hist_test$density)
            ))
          )
          hist(
            nn_z_scores_cond,
            breaks = first_hist_test$breaks,
            freq = F,
            border = NA,
            col = rgb(0, 0, 0, 0.7),
            add = T
          )
          mu <- mean(nn_z_scores_cond)
          sigma <- sd(nn_z_scores_cond)
          
          ranges <- c(non_nn_z_scores_cond,nn_z_scores_cond)
          line_range <- seq(min(ranges),max(ranges),by = 0.05)
          
          lines(line_range, dnorm(line_range, mean = mu, sd = sigma), lwd = 0.5)
        }
        
        
        unequal_tests <- c(`>=`, `<=`)
        precision_list <- lapply(unequal_tests, function(unequal_test) {
          mu <- mean(nn_scores_cond)
          sigma <- sd(nn_scores_cond)
          #sapply(1:length(non_nn_scores_cond), function(i) {
          sapply(1:length(all_scores_cond), function(i) {
            observed <-
              sum(unequal_test(all_scores_cond, all_scores_cond[i]))
            if (all.equal(unequal_test, `<=`) == T) {
              use_lower_tail <- T
            } else{
              use_lower_tail <- F
            }
            
            #Using normal distribuion to calculate expected
            expected <-
              pnorm(
                #non_nn_scores_cond[i],
                all_scores_cond[i],
                mean = mu,
                sd = sigma,
                lower.tail = use_lower_tail
              )#*length(all_scores_cond)
            
            #return(expected)
            
            #if(sign(all_scores_cond[i]) == 1){
            #  expected <- empPvals(all_scores_cond[i],nn_scores_cond)
            #}else{
            #  expected <- empPvals(-all_scores_cond[i],-nn_scores_cond)
            #}
            #expected <- empPvals(all_scores_cond[i],nn_scores_cond)
            
            
            #expected <- max(sum(unequal_test(nn_scores_cond, all_scores_cond[i])),1)
            #expected <- expected/length(nn_scores_cond)
            #expected <- expected*length(all_scores_cond)
            
            #/length(all_scores_cond)
            #expected <- max(expected, 1/length(all_scores_cond))
            
            #observed <- sum(unequal_test(all_scores_cond, all_scores_cond[i]))
            
            #fdr <- expected/observed
            
            #print(expected)
            #
            
            #
            #expected <- expected * length(all_scores_cond)
            
            return(min(c(expected,1 - expected))*2)
            
            #fdr <- (expected /(expected + observed))
            
            #fdr <- max(fdr, 1/length(nn_scores_cond))
            
            #return(fdr)
            #No sense returning fdr estimates >100%
            #return(min(fdr, 1))
          })
        })
        fdrs_pos <- precision_list[[1]]
        fdrs_neg <- precision_list[[2]]
        comb_fdr <- cbind(fdrs_pos, fdrs_neg)
        
        #Return positive or negative FDRs based on the nominal sign
        #of the interaction
        fdr_scores <- sapply(1:nrow(comb_fdr), function(i) {
          if (all_scores_cond[i] < 0) {
            return(fdrs_neg[i])
          } else{
            return(fdrs_pos[i])
          }
        })
        
       # print(fdrs_pos)
        
       hist(fdr_scores,main=paste(c(condition1,condition2),collapse=' - '))
       
       
       fdr_scores <- qvalue(fdr_scores)$q
        
        if(require_sign_change){
          sig <- (fdr_scores < fdr_cutoff) & (gi_data[,z_class_name1] != gi_data[,z_class_name2]) & (abs(all_gi_scores_cond) >= delta_gi_cutoff)
        }else{
          sig <- (fdr_scores < fdr_cutoff) & (abs(all_gi_scores_cond) >= delta_gi_cutoff)
        }

        ddr_data <- gi_data[,]#gi_data[ddr_pairs,]
        
        sig_names <- apply(ddr_data[which(sig),1:2],1,function(x){paste(x,collapse='_')})
        
        ret_list[[condition1]][[condition2]] <- sig_names
        
        
        
        retvec <- cbind(ddr_data[which(sig),c("Barcode_x","Barcode_y")],
                        rep(condition1,sum(sig)),
                        rep(condition2,sum(sig)),
                        ddr_data[which(sig),c("Type_of_gene_x","Type_of_gene_y")],
                        ddr_data[which(sig),c(z_name1,z_name2)],
                        ddr_data[which(sig),c(gi_name1,gi_name2)],
                        ddr_data[which(sig),c(z_class_name1,z_class_name2)],
                        all_z_scores_cond[which(sig)],
                        all_gi_scores_cond[which(sig)],
                        fdr_scores[which(sig)])
        
        colnames(retvec) <-
          c(
            'Barcode_x',
            'Barcode_y',
            'Condition1',
            'Condition2',
            'Type_of_gene_x',
            'Type_of_gene_y',
            'Z_Condition1',
            'Z_Condition2',
            'GI_Condition1',
            'GI_Condition2',
            'Class_Condition1',
            'Class_Condition2',
            'DeltaZ',
            'DeltaGI',
            'DeltaZ_FDR'
          )
        
        ret_df <- rbind(ret_df,retvec)
      }
    }
  }
  return(ret_df)
}

# 
# differential_calls_by_gene <- function(differential_calls, fdr_cutoff = 0.05) {
#   Gene1 <-
#     sapply(differential_calls$Barcode_x, function(name) {
#       strsplit(name, split = '_')[[1]][1]
#     })
#   
#   Gene2 <-
#     sapply(differential_calls$Barcode_y, function(name) {
#       strsplit(name, split = '_')[[1]][1]
#     })
#   
#   
#   differential_calls <- cbind(differential_calls, Gene1, Gene2)
#   
#   #differential_calls$Gene1 <- as.character(differential_calls$Gene1)
#   #differential_calls$Gene2 <- as.character(differential_calls$Gene2)
#   #rownames(differential_calls) <- NULL
#   
#   #grouped_calls <- differential_calls %>% group_by(Gene1,Gene2)#,Condition1,Condition2)
#   
#   #grouped_xndices <- attributes(grouped_calls)$indices
#   
#   split_calls <-
#     split(
#       differential_calls,
#       list(
#         Gene1,
#         Gene2,
#         differential_calls$Condition1,
#         differential_calls$Condition2
#       )
#     )
#   
#   ret_df <- c()
#   for (i in 1:length(split_calls)) {
#     #print(i)
#     gene_calls <- split_calls[[i]]
#     if (nrow(gene_calls) > 0) {
#       n_diff <-
#         sum(gene_calls$Class_Condition1 != gene_calls$Class_Condition2)
#       all_same_direction <-
#         length(unique(gene_calls$Class_Condition1) == 1) &
#         length(unique(gene_calls$Class_Condition2) == 1)
#       max_fdr <- max(gene_calls$DeltaZ_FDR)
#       mean_Z_Condition1 <- mean(gene_calls$Z_Condition1)
#       mean_Z_Condition2 <- mean(gene_calls$Z_Condition2)
#       mean_DeltaZ <- mean(gene_calls$DeltaZ)
#       
#       if (max_fdr < fdr_cutoff) {
#         call <- 'DIFFERENTIAL'
#       } else{
#         call <- 'NON_DIFFERENTIAL'
#       }
#       
#       
#       #Merge the directions by consensus
#       call_count_cond1 <- table(gene_calls$Class_Condition1)
#       call_count_cond2 <- table(gene_calls$Class_Condition2)
#       
#       m1 <- max(call_count_cond1)
#       m2 <- max(call_count_cond2)
#       
#       consensus_class1 <- names(which(call_count_cond1 == m1))
#       consensus_class2 <- names(which(call_count_cond2 == m2))
#       
#       if (length(consensus_class1) > 1) {
#         consensus_class1 <- 'AMBIGUOUS'
#       }
#       if (length(consensus_class2) > 1) {
#         consensus_class2 <- 'AMBIGUOUS'
#       }
#       
#       info <-
#         gene_calls[1, c(
#           'Gene1',
#           'Gene2',
#           'Condition1',
#           'Condition2',
#           'Class_Condition1',
#           'Class_Condition2'
#         )]
#       info <- apply(info, 2, as.character)
#       
#       info <- c(info, mean_Z_Condition1)
#       info <- c(info, mean_Z_Condition2)
#       info <- c(info, mean_DeltaZ)
#       info <- c(info, max_fdr)
#       info <- c(info, call)
#       
#       names(info)[(length(info) - 6):length(info)] <-
#         c(
#           'Consensus_Class_Condition1',
#           'Consensus_Class_Condition2',
#           'mean_Z_Condition1',
#           'mean_Z_Condition2',
#           'mean_DeltaZ',
#           'Max_FDR',
#           'Call'
#         )
#       info$Consensus_Class_Condition1 <- consensus_class1
#       info$Consensus_Class_Condition2 <- consensus_class2
#       
#       ret_df <- rbind(ret_df, unlist(info))
#     }
#   }
#   
#   ret_df <- as.data.frame(ret_df)
#   return(ret_df)
# }


differential_calls_histogram <- function(differential_calls){
  orig_hist <-
    hist(
      differential_calls$DeltaZ,
      breaks = 300,
      freq = F,
      border = NA,
      col = rgb(0, 0, 0, 0.7),
      xlim = c(-20, 20),
      main = '',
      xlab = '∆Z',
      ylim=c(0,0.3)
    )
  
  hist(
    as.matrix(
      filter(
        differential_calls,
        Class_Condition1 == "NEUTRAL" &
          Class_Condition2 == "NEUTRAL"# &
          #Type_of_gene_x != 'Neutral' &
          #Type_of_gene_y != 'Neutral'
      )$DeltaZ
    ),
    breaks = orig_hist$breaks,
    add = T,
    col = rgb(1, 0, 0, 0.5),
    freq = FALSE,
    border = NA
  )
  
  hist(
    as.matrix(
      filter(
        differential_calls,
        Type_of_gene_x == 'Neutral' |
          Type_of_gene_y == 'Neutral'
      )$DeltaZ
    ),
    add = T,
    col = rgb(0, 0, 1, 0.5),
    freq = F,
    border = NA,
    breaks = orig_hist$breaks
  )
  legend('topleft',legend=c('All pairs','non-DDR Pairs','neutral GI pairs'),col=c('grey30','blue','red'),pch=15)
}