this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')

gi_data <- read.table('table_s1_new.tsv', head = T)

gi_data_filtered <-
  dplyr::filter(gi_data,
                SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
nn_pairs <-
  gi_data$Type_of_gene_i == 'Neutral' &
  gi_data$Type_of_gene_j == 'Neutral'

ddr_pairs <-
  gi_data$Type_of_gene_i == 'DNA_repair' &
  gi_data$Type_of_gene_j == 'DNA_repair'

#Alternatively, can define all non DNA_repair - DNA_repair as neutral
nn_pairs <-
  gi_data$Type_of_gene_i != 'DNA_repair' | gi_data$Type_of_gene_j != 'DNA_repair'

gi_scores <- gi_data[,grep('^GIS',colnames(gi_data))]

nn_scores <- gi_scores[nn_pairs,]
non_nn_scores <- gi_scores[ddr_pairs,]


for (condition in colnames(gi_scores)[c(3)]) {
  if (condition %in% c("GIS_ij.NoDrug", "GIS_ij.DMSO")) {
    st_onge_class <- 'SOJ_Class_NoMMS'
  }
  if (condition == "GIS_ij.MMS") {
    st_onge_class <- 'SOJ_Class_MMS'
  }
  
  gi_scores_filtered <- gi_data_filtered[,grep('^GIS',colnames(gi_data))]
  raw_scores_cond_filtered <- gi_scores_filtered[, condition]
  raw_scores_cond <- gi_scores[,condition]
  
  labels_pos <- gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
  labels_neg <- gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
  
  pos_prec_stonge <- sapply(raw_scores_cond_filtered,function(score){
    sum(labels_pos[raw_scores_cond_filtered >= score])/sum(raw_scores_cond_filtered >= score)
  })
  
  neg_prec_stonge <- sapply(raw_scores_cond_filtered,function(score){
    sum(labels_neg[raw_scores_cond_filtered <= score])/sum(raw_scores_cond_filtered <= score)
  })
  
  
  nn_scores_cond <-  unlist(as.vector(nn_scores))#nn_scores[,condition]#unlist(as.vector(nn_scores))#nn_scores[,condition]
  non_nn_scores_cond <- non_nn_scores[,condition]#unlist(as.vector(non_nn_scores))#
  
  pos_prec_internal <- sapply(1:length(raw_scores_cond),function(i){
    observed <- sum(non_nn_scores_cond >= raw_scores_cond[i])
    expected <- sum(nn_scores_cond >= raw_scores_cond[i])
    return(expected)
    #expected should never be 0
    expected <- max(expected,1)
    
    expected <- expected/length(nn_scores_cond)
    expected <- expected*length(non_nn_scores_cond)
    fdr <- (expected/observed)
    fdr <- min(fdr,1)
    return(1 - fdr)
  })
  
  neg_prec_internal <- sapply(1:length(raw_scores_cond),function(i){
    observed <- sum(non_nn_scores_cond <= raw_scores_cond[i])
    expected <- sum(nn_scores_cond <= raw_scores_cond[i])
    return(expected)
    #expected should never be 0
    #expected <- max(expected,1)
    
    #expected <- expected/length(nn_scores_cond)
    #expected <- expected*length(non_nn_scores_cond)
    
    #fdr <- (expected/observed)
    #fdr <- min(fdr,1)
    #return(1 - fdr)
  })
  
  relevant_internal_neg <- neg_prec_internal[which(gi_data$SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))]
  relevant_internal_pos <- pos_prec_internal[which(gi_data$SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))]
  #raw_scores_cond_neg <- raw_scores_cond_neg[z_score_ind]
  
  
  
  #predict
  #acc_pos <- sapply(1:length(z_scores_cond_pos), function(i) {
  #  sum(labels_pos[i:length(labels_pos)]) / (length(labels_pos) - i + 1)
  #})
  
  #acc_neg <- sapply(1:length(z_scores_cond_pos), function(i) {
  #  sum(labels_pos[i:length(labels_pos)]) / (length(labels_pos) - i + 1)
  #})
  
}