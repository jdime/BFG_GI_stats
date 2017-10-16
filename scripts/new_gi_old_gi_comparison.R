library(ROCR)

rhombus_integration <- function(x, y) {
  sum <- 0
  for (i in 2:length(x)) {
    base <- x[i] - x[i - 1]
    mid_height <- mean(c(y[i], y[i - 1]))
    sum <- sum + base * mid_height
  }
  return(sum)
}

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')

gi_data_new <- read.table('table_s1_new.tsv', head = T)
#gi_data_new <-
#  dplyr::filter(gi_data_new,
#                Remove_by_Chromosomal_distance_or_SameGene == 'no')
#gi_data_new <-
#  dplyr::filter(gi_data_new, C_ij.HetDipl > 100)

gi_data_old <- read.table('table_s1.tsv', head = T)
#gi_data_old <-
#  dplyr::filter(gi_data_old,
#                Remove_by_Chromosomal_distance_or_SameGene == 'no')
#gi_data_old <-
#  dplyr::filter(gi_data_old, C_ij.HetDipl > 100)
par(mfrow = c(1, 2))

gi_data_list <- list(gi_data_new, gi_data_old)
names(gi_data_list) <- c('new', 'old')


#Compare AUC
gi_scores <- gi_data_new[, grep('^GIS', colnames(gi_data_new))]
for (condition in colnames(gi_scores)[c(2, 3)]) {
  for (i in 1:length(gi_data_list)) {
    gi_data <- gi_data_list[[i]]
    gi_data_filtered <-
      dplyr::filter(gi_data,
                    SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
    
    gi_scores <- gi_data[, grep('^GIS', colnames(gi_data))]
    
    
    if (condition %in% c("GIS_ij.NoDrug", "GIS_ij.DMSO")) {
      st_onge_class <- 'SOJ_Class_NoMMS'
    }
    if (condition == "GIS_ij.MMS") {
      st_onge_class <- 'SOJ_Class_MMS'
    }
    
    labels_pos <- gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
    labels_neg <- gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
    z_scores_cond_pos <-
      gi_data_filtered[, condition]#, collapse = '')]
    z_scores_cond_neg <- -z_scores_cond_pos
    
    
    perf_pos <-
      performance(prediction(z_scores_cond_pos, labels_pos), 'sens', 'fpr')
    perf_neg <-
      performance(prediction(z_scores_cond_neg, labels_neg), 'sens', 'fpr')
    
    plot(
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
for (condition in colnames(gi_scores)[c(2,3)]) {
  for (i in 1:length(gi_data_list)) {
    gi_data <- gi_data_list[[i]]
    if (condition %in% c("GIS_ij.NoDrug", "GIS_ij.DMSO")) {
      st_onge_e <- 'SOJ_E_NoMMS'
    }
    if (condition == "GIS_ij.MMS") {
      st_onge_e <- 'SOJ_E_MMS'
    }
    plot(gi_data[, st_onge_e], gi_data[, condition], pch = 16,
         main= paste(names(gi_data_list)[i],condition), col = rgb(0,0,0,0.3),
         xlab = 'St Onge E',ylab ='GIS')
    lines(smooth.spline(gi_data[is.finite(gi_data[, st_onge_e]), st_onge_e], gi_data[is.finite(gi_data[, st_onge_e]), condition],penalty=10), col =
            'blue')
   correl <- cor(gi_data[, st_onge_e], gi_data[, condition],use='pair')
   
   text(0,0,sprintf('R = %s',format(correl,digits=2)),col='red',cex=1.5)
    #stop()
  }
  

}
