devtools::use_package('ROCR')
devtools::use_package('dplyr')

rhombus_integration <- function(x, y) {
  sum <- 0
  for (i in 2:length(x)) {
    base <- x[i] - x[i - 1]
    mid_height <- mean(c(y[i], y[i - 1]))
    sum <- sum + base * mid_height
  }
  return(sum)
}


st_onge_comparison_plot <- function(gi_data_old,
                                    gi_data_new,
                                    control_name = "Z_GIS_ij.DMSO",
                                    condition_name = "Z_GIS_ij.MMS"){
  par(mfrow = c(1, 2))
  
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


prec_recall_vs_stonge <- function(gi_data,
                                  score_cols = NULL,
                                  control_name = "GIS_ij.DMSO",
                                  condition_name = "GIS_ij.MMS",
                                  fdr_prefix = "FDR.Internal_ij",
                                  fdr_cutoff = 0.05,
                                  xlims = NULL) {
  gene1 <- sapply(gi_data$Barcode_i, function(x) {
    strsplit(x, split = '_')[[1]][1]
  })
  
  gene2 <- sapply(gi_data$Barcode_j, function(x) {
    strsplit(x, split = '_')[[1]][1]
  })
  
  if (is.null(score_cols)) {
    score_cols <- grep('^GIS', colnames(gi_data))
  }
  scores <- gi_data[, score_cols]
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
    
    
    
    scores_cond <- gi_data_filtered[, condition]
    labels_pos <-
      gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
    labels_neg <-
      gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
    
    
    ##scores_cond_pos < - scores_cond[scores_cond > 0]
    score_ind_pos <- sort(scores_cond, index.return = T)$ix
    
    scores_cond_pos <- scores_cond[score_ind_pos]
    labels_pos <- labels_pos[score_ind_pos]
    
    
    
    prec_pos <- sapply(1:length(scores_cond_pos), function(i) {
      sum(labels_pos[i:length(labels_pos)]) / (length(labels_pos) - i + 1)
    })
    fdr_pos <- 1 - prec_pos
    
    
    
    score_ind_neg <- sort(scores_cond, index.return = T, decreasing = T)$ix
    
    scores_cond_neg <- scores_cond[score_ind_neg]
    labels_neg <- labels_neg[score_ind_neg]
     prec_neg <- sapply(1:length(scores_cond), function(i) {
      sum(labels_neg[i:length(labels_neg)]) / (length(labels_neg) - i + 1)
    })
    fdr_neg <- 1 - prec_neg
    
    drug <- strsplit(condition,split='\\.')[[1]][2]
    fdr_col <- paste(c(fdr_prefix,drug),collapse='.')
    
    if (condition == control_name) {
      if(is.null(xlims)){
        xlims <-c(min(scores_cond_pos), max(scores_cond_pos))
      }
      
      
      plot(
        scores_cond_pos[scores_cond_pos > 0],
        prec_pos[scores_cond_pos > 0]*100,
        type = 'l',
        ylab = 'St.Onge Validation Rate (%)',
        xlab = 'Z Cutoff',
        main = 'GI Precision vs St. Onge',
        #st_onge_class,
        lwd = 2,
        xlim = xlims)
        
        lines(scores_cond_neg[scores_cond_neg < 0],
              prec_neg[scores_cond_neg < 0]*100,
              lwd = 2)
        
        abline(v=min(scores_cond[scores_cond > 0 & gi_data_filtered[,fdr_col] <= fdr_cutoff]))
        abline(v=max(scores_cond[scores_cond < 0 & gi_data_filtered[,fdr_col] <= fdr_cutoff]))
        #abline(min(v=))
      
    } else{
      lines(scores_cond_pos[scores_cond_pos > 0],
            prec_pos[scores_cond_pos > 0]*100,
            col = 'red',
            lwd = 2)
      lines(scores_cond_neg[scores_cond_neg < 0],
            prec_neg[scores_cond_neg < 0]*100,
            col = 'red',
            lwd = 2)
      abline(v=min(scores_cond[scores_cond > 0 & gi_data_filtered[,fdr_col] <= fdr_cutoff]),col='red')
      abline(v=max(scores_cond[scores_cond < 0 & gi_data_filtered[,fdr_col] <= fdr_cutoff]),col='red')
      
    }
  }
}
  
