library(dplyr)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')

gi_data <- read.table('table_s1_new.tsv',head=T)
gi_data <-
  dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')

#Neutral-Neutral pairs
nn_pairs <-
  gi_data$Type_of_gene_i == 'Neutral' &
  gi_data$Type_of_gene_j == 'Neutral'

ddr_pairs <-
  gi_data$Type_of_gene_i == 'DNA_repair' &
  gi_data$Type_of_gene_j == 'DNA_repair'



#Alternatively, can define all non DNA_repair - DNA_repair as neutral
nn_pairs <-
  gi_data$Type_of_gene_i != 'DNA_repair' | gi_data$Type_of_gene_j != 'DNA_repair'
  
     #gi_data$Type_of_gene_j != 'DNA_repair')# |
  #(gi_data$Type_of_gene_i != 'DNA_repair' &
  #   gi_data$Type_of_gene_j == 'DNA_repair')


gi_scores <- gi_data[,grep('^GIS',colnames(gi_data))]
nn_scores <- gi_scores[nn_pairs,]
non_nn_scores <- gi_scores[ddr_pairs,]

#Set false discovery rate
fdr_cutoff <- 0.05
z_cutoff_positive <- 1.2
z_cutoff_negative <- -0.6

summary <- c()


for(condition in colnames(gi_scores)){
  #How many interactions are observed at a cutoff
  obs_vec_neg <- c()
  #How many false interactions are expected at a cutoff
  exp_vec_neg <- c()
  #Predicted true positive rate (1 - fdr_cutoff)
  fdr_vec_neg <- c()
  
  obs_vec_pos <- c()
  exp_vec_pos <- c()
  fdr_vec_pos <- c()
  
  
  nn_scores_cond <- sort(unlist(nn_scores[,condition])) 
  non_nn_scores_cond <- non_nn_scores[,condition]
  
  
  #nn_scores_cond <- nn_scores_cond# - median(nn_scores_cond)
  #non_nn_scores_cond <- non_nn_scores_cond# - median(nn_scores_cond)
  
  z_scores_cond <- gi_data[,paste(c('Z_',condition),collapse='')]#scale(gi_scores[,condition])
  z_based_cutoff_pos <-
    min(gi_scores[, condition][which(z_scores_cond > z_cutoff_positive)])
  z_based_cutoff_neg <-
    max(gi_scores[, condition][which(z_scores_cond < z_cutoff_negative)])
  
  #Test fdr_cutoff at all score cutoffs based on neutral-neutral scores
  for(i in 1:length(nn_scores_cond)){
    #Negative gis
    expected_false <- (i/length(nn_scores_cond))*length(non_nn_scores_cond)
    observed <- sum(non_nn_scores_cond < nn_scores_cond[i])
    
    
    obs_vec_neg <- c(obs_vec_neg,observed)
    exp_vec_neg <- c(exp_vec_neg,expected_false)
    #fdr <- 1 - (c(observed - expected_false)/observed)
    fdr <- expected_false/observed
    fdr <- min(fdr,1)
    fdr_vec_neg <- c(fdr_vec_neg,c(observed - expected_false)/observed)
    
    expected_false <- ((length(nn_scores_cond) - i)/length(nn_scores_cond))*length(non_nn_scores_cond)
    observed <- sum(non_nn_scores_cond > nn_scores_cond[i])
    
    obs_vec_pos <- c(obs_vec_pos,observed)
    exp_vec_pos <- c(exp_vec_pos,expected_false)
    
    fdr <- 1 - (c(observed - expected_false)/observed)
    fdr <- min(fdr,1)
    fdr_vec_pos <- c(fdr_vec_pos,fdr)
  }
  
  fdr_vec_neg <- 1 - (obs_vec_neg - exp_vec_neg)/obs_vec_neg
  fdr_vec_pos <- 1 - (obs_vec_pos - exp_vec_pos)/obs_vec_pos
  
  #Calculate fdr_cutoff at Z cutoff for positive and negative interactions
  expected_false_neg_gi_z <-
    (sum(nn_scores_cond < z_based_cutoff_neg) / length(nn_scores_cond)) * length(non_nn_scores_cond)
  
  expected_false_pos_gi_z <-
    (sum(nn_scores_cond > z_based_cutoff_pos) / length(nn_scores_cond)) * length(non_nn_scores_cond)
  
  observed_neg_gi_z <- 
    sum(non_nn_scores_cond < z_based_cutoff_neg)
  
  observed_pos_gi_z <- 
    sum(non_nn_scores_cond > z_based_cutoff_pos)
  
  
  fdr_z_pos <- 1 - (observed_pos_gi_z - expected_false_pos_gi_z)/observed_pos_gi_z
  fdr_z_neg <- 1 - (observed_neg_gi_z - expected_false_neg_gi_z)/observed_neg_gi_z
  
  
  #print(condition)
  #
  pos_fdr_cutoff <- min(which(fdr_vec_pos <= fdr_cutoff))
  obs_pos_at_fdr_cutoff <- obs_vec_pos[pos_fdr_cutoff]
  exp_false_pos_at_fdr_cutoff <- exp_vec_pos[pos_fdr_cutoff]
  pos_fdr <- fdr_vec_pos[pos_fdr_cutoff]
  
  
  #print(1 - fdr_vec_pos[pos_cutoff])
  #print(c(obs_vec_pos[pos_cutoff],exp_vec_pos[pos_cutoff]))
  
  neg_fdr_cutoff <- max(which(fdr_vec_neg <= fdr_cutoff))
  obs_neg_at_fdr_cutoff <- obs_vec_neg[neg_fdr_cutoff]
  exp_false_neg_at_fdr_cutoff <- exp_vec_neg[neg_fdr_cutoff]
  neg_fdr <- fdr_vec_neg[neg_fdr_cutoff]
  
  
  
  
  #summary <- rbind(summary,c(condition,fdr_z_pos,fdr_z_neg,observed_pos_gi_z,observed_neg_gi_z))
  #
  summary <-
    rbind(
      summary,
      c(
        condition,
        pos_fdr,
        neg_fdr,
        obs_pos_at_fdr_cutoff,
        obs_neg_at_fdr_cutoff,
        z_scores_cond[max(which(fdr_vec_neg <= fdr_cutoff))],
        z_scores_cond[min(which(fdr_vec_pos <= fdr_cutoff))]
      )
    )
  
  #print(1 - fdr_vec_neg[neg_cutoff])
  #print(c(obs_vec_neg[neg_cutoff],exp_vec_neg[neg_cutoff]))
  
  #print(c(observed_neg_gi_z,expected_false_neg_gi_z))
  #print(c(observed_pos_gi_z,expected_false_pos_gi_z))
  
  
  
  
  
}
