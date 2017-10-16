library(dplyr)
pos_cut <- 0.146
neg_cut <- -0.079

pos_cut <- 1#1e-10#0.089#0.146
neg_cut <- -1#e-10#0.01#-0.079

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')

gi_data <- read.table('table_s1_new.tsv',head=T,stringsAsFactors = F)

gi_data_unfiltered <- gi_data

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

#stop()

gi_scores <- gi_data[,grep('^GIS',colnames(gi_data))]
gi_scores_unfiltered <- gi_data_unfiltered[,grep('^GIS',colnames(gi_data))]
z_cols <-
  grep('Class',
       grep('^Z_GIS', colnames(gi_data), val = T),
       invert = T,
       val = T)

z_scores <- gi_data[,z_cols]

nn_scores <- gi_scores[nn_pairs,]
non_nn_scores <- gi_scores[ddr_pairs,]

plot(density(unlist(nn_scores)),lwd=2,xlab='GIS',main='')
lines(density(unlist(non_nn_scores)),lwd=2,col='red')

fdr_vec <- c()
result_table <- c()
for(condition in colnames(gi_scores)){
  all_scores <- gi_scores_unfiltered[,condition]
  nn_scores_cond <-  nn_scores[,condition]#unlist(as.vector(nn_scores))#nn_scores[,condition]
  non_nn_scores_cond <- non_nn_scores[,condition]#unlist(as.vector(non_nn_scores))#
  
  pos_prec_vec <- sapply(1:length(non_nn_scores_cond),function(i){
    observed <- sum(non_nn_scores_cond >= non_nn_scores_cond[i])
    expected <- sum(nn_scores_cond >= non_nn_scores_cond[i])
    #expected should never be 0
    expected <- max(expected,1)
    
    expected <- expected/length(nn_scores_cond)
    expected <- expected*length(non_nn_scores_cond)
    fdr <- (expected/observed)
    fdr <- min(fdr,1)
    return(1 - fdr)
  })
  
  neg_prec_vec <- sapply(1:length(non_nn_scores_cond),function(i){
    observed <- sum(non_nn_scores_cond <= non_nn_scores_cond[i])
    expected <- sum(nn_scores_cond <= non_nn_scores_cond[i])
    #expected should never be 0
    expected <- max(expected,1)
    
    expected <- expected/length(nn_scores_cond)
    expected <- expected*length(non_nn_scores_cond)
    
    fdr <- (expected/observed)
    fdr <- min(fdr,1)
    return(1 - fdr)
  })
  
  continuous_pos_prec_vec <- approx(non_nn_scores_cond,pos_prec_vec,n=10000)
  continuous_neg_prec_vec <- approx(non_nn_scores_cond,neg_prec_vec,n=10000)
  
  fdr_scores <- all_scores
  fdr_scores[all_scores < 0] <- 1 - sapply(all_scores[all_scores < 0],function(score){
    continuous_neg_prec_vec$y[which.min(abs(continuous_neg_prec_vec$x - score))]
  })
  fdr_scores[all_scores > 0] <- 1 - sapply(all_scores[all_scores > 0],function(score){
    continuous_pos_prec_vec$y[which.min(abs(continuous_pos_prec_vec$x - score))]
  })
  
  fdr_vec <- cbind(fdr_vec,fdr_scores)
  
  #print(condition)
  neg <- sum(non_nn_scores[,condition] < neg_cut)
  pos <- sum(non_nn_scores[,condition] > pos_cut)
  fp_neg_gis <- round((sum(nn_scores_cond < neg_cut)/length(nn_scores_cond))*length(non_nn_scores_cond))
  fp_pos_gis <- round((sum(nn_scores_cond > pos_cut)/length(nn_scores_cond))*length(non_nn_scores_cond))
  
  fdr_neg <- format(max(fp_neg_gis,1)/neg,digits=2)
  fdr_pos <- format(max(fp_pos_gis,1)/pos,digits=2)
  
  fdr_10_cutoff_neg <- max(non_nn_scores_cond[neg_prec_vec > 0.95])
  fdr_10_cutoff_pos <- min(non_nn_scores_cond[pos_prec_vec > 0.95])
  
  result_table <- cbind(result_table,c(pos,neg,fp_pos_gis,fp_neg_gis,fdr_pos,fdr_neg,fdr_10_cutoff_neg,fdr_10_cutoff_pos))
  
  #print(sum(non_nn_scores[,condition] > max(non_nn_scores_cond[neg_prec_vec > 0.98])))
  #print(sum(non_nn_scores[,condition] > min(non_nn_scores_cond[pos_prec_vec > 0.98])))
  
#  plot(non_nn_scores_cond,neg_prec_vec,type='p',xlim=c(-0.5,0.5),col='blue',pch=16,ylim=c(0.5,1),main=condition)
 # points(non_nn_scores_cond,pos_prec_vec,type='p',col='red',pch=16)
#  #abline(h=0.95)
#  abline(v=-0.17)
#  abline(v=0.10)
  
  
}
colnames(result_table) <- sapply(colnames(gi_scores),function(x){strsplit(x,split='\\.')[[1]][2]})
rownames(result_table) <-
  c(
    'Positive Interactions',
    'Negative Interactions',
    'Expected False Positives +GIs',
    'Expected False Positives -GIs',
    'Estimated FDR +GIs',
    'Estimated FDR -GIs',
    '10% FDR Cutoff +GIs',
    '10% FDR Cutoff -GIs'
  )
#for(i in 1:sort(non_nn_scores))

