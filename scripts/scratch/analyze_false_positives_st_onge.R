library(dplyr)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')


gi_data <- read.table('table_s1_new.tsv', head = T, stringsAsFactors = F)
gene1 <- sapply(gi_data$Barcode_i, function(x) {
  strsplit(x, split = '_')[[1]][1]
})

gene2 <- sapply(gi_data$Barcode_j, function(x) {
  strsplit(x, split = '_')[[1]][1]
})
#gi_data <-
#  dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no', C_ij.HetDipl > 100)

precision_cutoff <- 0.8

#Neutral-Neutral pairs
nn_pairs <-
  gi_data$Type_of_gene_i == 'Neutral' &
  gi_data$Type_of_gene_j == 'Neutral'

#Alternatively, can define all non DNA_repair - DNA_repair as neutral
nn_pairs <-
  !(gi_data$Type_of_gene_i == 'DNA_repair' &
  gi_data$Type_of_gene_j == 'DNA_repair') &
  gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'



gi_scores <- gi_data[, grep('^GIS', colnames(gi_data))]
nn_scores <- gi_scores[nn_pairs, ]
non_nn_scores <- gi_scores[!(nn_pairs), ]
z_cols <-
  grep('Class',
       grep('^Z_GIS', colnames(gi_data), val = T),
       invert = T,
       val = T)
z_cols <- which(colnames(gi_data) %in% z_cols)

#Set false discovery ratea
fdr_cutoff <- 0.1
z_cutoff_positive <- 1
z_cutoff_negative <- -0.6

summary <- c()

prec_pos <- c()
prec_neg <- c()



for (analysis in c("gis_cutoff_plot")) {
  for (condition in colnames(gi_scores)[c(2, 3)]) {
    if (condition %in% c("GIS_ij.NoDrug", "GIS_ij.DMSO")) {
      st_onge_class <- 'SOJ_Class_NoMMS'
    }
    if (condition == "GIS_ij.MMS") {
      st_onge_class <- 'SOJ_Class_MMS'
    }
    
    gi_data_filtered <-
      dplyr::filter(gi_data,
                    SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
    
    
    
    z_scores_cond_pos <-
      gi_data_filtered[, paste(c('Z_', condition), collapse = '')]
    raw_scores_cond_pos <- gi_data_filtered[, condition]
    
    labels_pos <- gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
    z_score_ind <- sort(z_scores_cond_pos, index.return = T)$ix
    
    raw_scores_cond_pos <- raw_scores_cond_pos[z_score_ind]
    z_scores_cond_pos <- z_scores_cond_pos[z_score_ind]
    labels_pos <- labels_pos[z_score_ind]
    
    
    prec_pos <- sapply(1:length(z_scores_cond_pos), function(i) {
      sum(labels_pos[i:length(labels_pos)]) / (length(labels_pos) - i + 1)
    })
    fdr_pos <- 1 - prec_pos
    
    
    
    
    z_scores_cond_neg <-
      gi_data_filtered[, paste(c('Z_', condition), collapse = '')]
    raw_scores_cond_neg <- gi_data_filtered[, condition]
    
    labels_neg <- gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
    z_score_ind <-
      sort(z_scores_cond_neg,
           index.return = T,
           decreasing = T)$ix
    
    raw_scores_cond_neg <- raw_scores_cond_neg[z_score_ind]
    z_scores_cond_neg <- z_scores_cond_neg[z_score_ind]
    labels_neg <- labels_neg[z_score_ind]
    
    prec_neg <- sapply(1:length(z_scores_cond_neg), function(i) {
      sum(labels_neg[i:length(labels_neg)]) / (length(labels_neg) - i + 1)
    })
    fdr_neg <- 1 - prec_neg
    
    
    # if (analysis == "z_cutoff_positive_plot") {
    #   plot(
    #     z_scores_cond_pos,
    #     prec_pos,
    #     type = 'l',
    #     ylab = 'Accuracy (1 - FDR)',
    #     xlab = 'Z Score Cutoff',
    #     main = st_onge_class
    #   )
    #   lines(z_scores_cond_neg, prec_neg)
    #   abline(v = z_cutoff_positive)
    #   abline(v = z_cutoff_negative)
    #   
    #   prec_neg_z1 <-
    #     prec_neg[which.min(abs(z_scores_cond_neg - (z_cutoff_negative)))]
    #   prec_pos_z1 <-
    #     prec_pos[which.min(abs(z_scores_cond_pos - (z_cutoff_positive)))]
    #   
    #   text(
    #     z_cutoff_negative - 0.8,
    #     0.8,
    #     sprintf(
    #       "Acc at Z < %s\n= %s",
    #       z_cutoff_negative,
    #       format(prec_neg_z1, digits = 2)
    #     )
    #   )
    #   text(
    #     z_cutoff_positive - 0.8,
    #     0.8,
    #     sprintf(
    #       "Acc at Z > %s\n= %s",
    #       z_cutoff_positive,
    #       format(prec_pos_z1, digits = 2)
    #     )
    #   )
    # }
    
    
    if (analysis == "gis_cutoff_plot") {
      if (condition %in% c("GIS_ij.NoDrug", "GIS_ij.DMSO")) {
        plot(
          raw_scores_cond_pos[raw_scores_cond_pos > 0],
          prec_pos[raw_scores_cond_pos > 0],
          type = 'l',
          ylab = 'Precision (1 - FDR)',
          xlab = 'GIS Cutoff',
          main = 'GI Precision vs St. Onge',#st_onge_class,
          lwd = 2,
          xlim=c(min(raw_scores_cond_pos),max(raw_scores_cond_pos))
        )
        
      } else{
        lines(raw_scores_cond_pos[raw_scores_cond_pos > 0],
              prec_pos[raw_scores_cond_pos > 0],
              col = 'red',
              lwd = 2)
        pos_cut <- min(raw_scores_cond_pos[which(prec_pos >= precision_cutoff)])
        abline(v=pos_cut)
      }
      
      if (condition %in% c("GIS_ij.NoDrug", "GIS_ij.DMSO")) {
        lines(raw_scores_cond_neg[raw_scores_cond_neg < 0],
              prec_neg[raw_scores_cond_neg < 0],
              lwd = 2)
        
      } else{
        lines(raw_scores_cond_neg[raw_scores_cond_neg < 0],
              prec_neg[raw_scores_cond_neg < 0],
              col = 'red',
              lwd = 2)
        neg_cut <- max(raw_scores_cond_neg[which(prec_neg >= precision_cutoff)])
        abline(v=neg_cut)
      }
    }
    
  }
}