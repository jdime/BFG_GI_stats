devtools::use_package('ROCR')
devtools::use_package('dplyr')


#Used for calculating AUC
rhombus_integration <- function(x, y) {
  sum <- 0
  for (i in 2:length(x)) {
    base <- x[i] - x[i - 1]
    mid_height <- mean(c(y[i], y[i - 1]))
    sum <- sum + base * mid_height
  }
  return(sum)
}

#Makes a precision 
precision_vs_stonge <- function(gi_data,
                                  control_name = "NoDrug",
                                  condition_name = "MMS",
                                  fdr_prefix = "FDR.Internal_ij",
                                  z_prefix = "Z_GIS_ij",
                                  fdr_cutoff = 0.05,
                                  xlims = c(-5,5),
                                  xlab = '-Log10 (Internal FDR)*sign(Z)') {
  
  for (condition in c(control_name, condition_name)) {
    if (condition == control_name) {
      st_onge_class <- 'SOJ_Class_NoMMS'
    }
    if (condition == condition_name) {
      st_onge_class <- 'SOJ_Class_MMS'
    }
    
    gi_data_filtered <-
      dplyr::filter(gi_data,
                    SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
    
    
    
    fdr_column <- paste(c(fdr_prefix,condition),collapse='.')
    z_column <- paste(c(z_prefix,condition),collapse='.')
    
    labels_pos <-
      gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
    labels_neg <-
      gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
    
    scores_cond <- sign(gi_data_filtered[,z_column])*-log10(gi_data_filtered[,fdr_column])
    
    pos_perf <- ROCR::performance(ROCR::prediction(scores_cond, labels_pos), 'prec')
    neg_perf <- ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), 'prec')
    
    pos_cutoff <- pos_perf@x.values[[1]]
    pos_precision <- pos_perf@y.values[[1]]
    
    neg_cutoff <- neg_perf@x.values[[1]]
    neg_precision <- neg_perf@y.values[[1]]
    
    if (condition == control_name) {
      if (is.null(xlims)) {
        xlims <-
          c(-1*neg_perf@x.values[[1]][2], pos_perf@x.values[[1]][2])
      }
      par(mar=c(4.5,5,3,1))
      plot(
        pos_cutoff[pos_cutoff > 0],
        (pos_precision*100)[pos_cutoff > 0],
        type = 'l',
        ylab = 'St.Onge Validation Rate (%)',
        xlab = xlab,
        main = 'GI Precision vs St. Onge',
        lwd = 1,
        xlim = xlims,
        ylim=c(0,100))
      lines(-neg_cutoff[neg_cutoff > 0],(neg_precision*100)[neg_cutoff > 0],lwd=1)
    } else{
      lines(pos_cutoff[pos_cutoff > 0],(pos_precision*100)[pos_cutoff > 0],lwd=1,col='red')
      lines(-neg_cutoff[neg_cutoff > 0],(neg_precision*100)[neg_cutoff > 0],lwd=1,col='red')
    }
  }
  
  abline(v=-log10(fdr_cutoff), lty = 3, lwd = 0.7)
  abline(v=log10(fdr_cutoff), lty = 3, lwd = 0.7)
  
  legend(xlims[1],50,legend=c(control_name, condition_name),fill=c('black','red'))
}






st_onge_comparison_plot <- function(gi_data_old,
                                    gi_data_new,
                                    control_name = "GIS_ij.DMSO",
                                    condition_name = "GIS_ij.MMS"){
  #par(mfrow = c(1, 2))
  
  gi_data_list <- list(gi_data_new, gi_data_old)
  names(gi_data_list) <- c('new', 'old')
  
  
  #Compare AUC
  #gi_scores <- gi_data_new[, grep('^GIS', colnames(gi_data_new))]
  for (condition in c(control_name, condition_name)) {
    for (i in 1:length(gi_data_list)) {
      gi_data <- gi_data_list[[i]]
      gi_data_filtered <-
        dplyr::filter(gi_data,
                      SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
      
      #gi_scores <- grep('^GIS', colnames(gi_data))]
      
      
      if (condition == control_name) {
        st_onge_class <- 'SOJ_Class_NoMMS'
      }
      if (condition == condition_name) {
        st_onge_class <- 'SOJ_Class_MMS'
      }
      
      labels_pos <- gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
      labels_neg <- gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
      z_scores_cond_pos <-
        gi_data_filtered[, condition]#, collapse = '')]
      z_scores_cond_neg <- -z_scores_cond_pos
      
      
      perf_pos <-
        ROCR::performance(ROCR::prediction(z_scores_cond_pos, labels_pos), 'sens', 'fpr')
      perf_neg <-
        ROCR::performance(ROCR::prediction(z_scores_cond_neg, labels_neg), 'sens', 'fpr')
      
      ROCR::plot(
        perf_pos,
        lwd = 2,
        col = 'blue',
        main = paste(names(gi_data_list)[i],condition)
      )
      lines(perf_neg@x.values[[1]],
            perf_neg@y.values[[1]],
            lwd = 2,
            col = 'red')
      
      auc_neg <-
        rhombus_integration(perf_neg@x.values[[1]], perf_neg@y.values[[1]])
      auc_pos <-
        rhombus_integration(perf_pos@x.values[[1]], perf_pos@y.values[[1]])
      
      text(0.8, 0.6, sprintf('AUC neg = %s', format(auc_neg, digits = 2)), col =
             'red')
      text(0.8, 0.4, sprintf('AUC pos = %s', format(auc_pos, digits = 2)), col =
             'blue')
      
    }
  }
  
  #Correlate with GI scores
  for (condition in c(control_name,condition_name)) {
    for (i in 1:length(gi_data_list)) {
      gi_data <- gi_data_list[[i]]
      if (condition == control_name) {
        st_onge_e <- 'SOJ_E_NoMMS'
      }
      if (condition == condition_name) {
        st_onge_e <- 'SOJ_E_MMS'
      }
      plot(gi_data[, st_onge_e], gi_data[, condition], pch = 16,
           main= paste(names(gi_data_list)[i],condition), col = rgb(0,0,0,0.3),
           xlab = 'St Onge E',ylab ='GIS')
      lines(smooth.spline(gi_data[is.finite(gi_data[, st_onge_e]), st_onge_e], 
                          gi_data[is.finite(gi_data[, st_onge_e]), condition],
                          penalty=10), col = 'blue')
      abline(c(0,1),lwd=2,col='red')
      correl <- cor(gi_data[, st_onge_e], gi_data[, condition],use='pair')
      
      text(0,0,sprintf('R = %s',format(correl,digits=2)),col='red',cex=1.5)
      #stop()
    }
    
    
  }

}






# 
# mcc_vs_stonge <- function(gi_data,
#                           score_cols = NULL,
#                           control_name = "GIS_ij.DMSO",
#                           condition_name = "GIS_ij.MMS",
#                           fdr_prefix = "FDR.Internal_ij",
#                           fdr_cutoff = 0.01,
#                           xlims = NULL) {
#   gene1 <- sapply(gi_data$Barcode_i, function(x) {
#     strsplit(x, split = '_')[[1]][1]
#   })
#   
#   gene2 <- sapply(gi_data$Barcode_j, function(x) {
#     strsplit(x, split = '_')[[1]][1]
#   })
#   
#   if (is.null(score_cols)) {
#     score_cols <- grep('^GIS', colnames(gi_data))
#   }
#   scores <- gi_data[, score_cols]
#   for (condition in c(control_name, condition_name)) {
#     if (condition == control_name) {
#       st_onge_class <- 'SOJ_Class_NoMMS'
#     }
#     if (condition == condition_name) {
#       st_onge_class <- 'SOJ_Class_MMS'
#     }
#     
#     gi_data_filtered <-
#       dplyr::filter(gi_data,
#                     SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
#     
#     scores_cond <- gi_data_filtered[, condition]
#     labels_pos <-
#       gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
#     labels_neg <-
#       gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
#     
#     pos_perf <- ROCR::performance(ROCR::prediction(scores_cond, labels_pos), 'mat')
#     neg_perf <- ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), 'mat')
#     
#     pos_cutoff <- pos_perf@x.values[[1]]
#     pos_precision <- pos_perf@y.values[[1]]
#     
#     neg_cutoff <- neg_perf@x.values[[1]]
#     neg_precision <- neg_perf@y.values[[1]]
#     
#     if (condition == control_name) {
#       if (is.null(xlims)) {
#         xlims <-
#           c(-1 * max(neg_perf@x.values[[1]]), max())
#       }
#       plot(
#         pos_cutoff[pos_cutoff > 0],
#         (pos_precision*100)[pos_cutoff > 0],
#         type = 'l',
#         ylab = 'MCC * 100',
#         xlab = 'Z Cutoff',
#         main = 'MCC vs St. Onge',
#         lwd = 2,
#         xlim = xlims,
#         ylim=c(0,100))
#       lines(-neg_cutoff[neg_cutoff > 0],(neg_precision*100)[neg_cutoff > 0],lwd=1)
#     } else{
#       lines(pos_cutoff[pos_cutoff > 0],(pos_precision*100)[pos_cutoff > 0],lwd=1,col='red')
#       lines(-neg_cutoff[neg_cutoff > 0],(neg_precision*100)[neg_cutoff > 0],lwd=1,col='red')
#     }
#   }
#   
# }
