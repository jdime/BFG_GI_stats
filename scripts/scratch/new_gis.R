this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')

gi_data <- read.table('table_s1.tsv', head = T, stringsAsFactors = F)
#gi_data <-
#  dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')

nn_pairs <-
  gi_data$Type_of_gene_i == 'Neutral' &
  gi_data$Type_of_gene_j == 'Neutral' &
  gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'

non_nn_unlinked_pairs <-
  gi_data$Type_of_gene_i != 'Neutral' &
  gi_data$Type_of_gene_j != 'Neutral' &
  gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'



count_data_real <- gi_data[, grep('^C_', colnames(gi_data))]

condition_sums_real <- apply(count_data_real, 2, sum)
freq_data_real <- count_data_real / condition_sums_real



dm_list <- list()
for(i in 1){
  print(i)
  
  if(i == 1){
    count_data <- sapply(1:length(condition_sums_real),function(i){
    sapply(1:nrow(freq_data_real),function(j){
      rbinom(1,condition_sums_real[i],freq_data_real[j,i])
    })
  }) + 1
  }
  else{
    count_data <- count_data_real + 10
  }
  condition_sums <- apply(count_data, 2, sum)
  freq_data <- count_data / condition_sums



#pseudo_counts <- 1 / condition_sums

#freq_data <- t(t(freq_data) + pseudo_counts)

norm_freq_data <- (freq_data/freq_data[,1])[,2:ncol(freq_data)]



#wt_fitness <- apply(norm_freq_data,2,function(x){
#  median(x[nn_pairs])
#})

#wt_norm_freq_data <- t(t(norm_freq_data)/wt_fitness)

genes <- unique(unlist(gi_data[, 1:2]))
neutral_genes <-
  unique(unlist(
    dplyr::filter(
      gi_data,
      Type_of_gene_i == 'Neutral' & Type_of_gene_j == 'Neutral'
    )[, 1:2]
  ))
ddr_genes <-
  unique(unlist(
    dplyr::filter(
      gi_data,
      Type_of_gene_i != 'Neutral' & Type_of_gene_j != 'Neutral'
    )[, 1:2]
  ))

single_fitness_scores <- t(sapply(genes,function(gene){
  criteria <- gi_data[,1] == gene & gi_data$Type_of_gene_j == 'Neutral' | gi_data[,2] == gene & gi_data$Type_of_gene_i == 'Neutral'
  criteria <- criteria & gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
  #wt_norm_freq_data[which(gi_data[,1]) == 1]
  return(apply(norm_freq_data[criteria,],2,mean))
}))

#stop()

#single_fitness_scores <- single_fitness_scores/single_fitness_scores[,1]
wt_fit <- apply(single_fitness_scores[neutral_genes,],2,median)



single_fitness_scores<- t(t(single_fitness_scores)/wt_fit)


norm_dm_fit <- t(t(norm_freq_data)/wt_fit)



bc1 <- gi_data$Barcode_i
bc2 <- gi_data$Barcode_j
dms <- t(sapply(1:nrow(norm_dm_fit),function(i){
  log2(norm_dm_fit[i,]/(single_fitness_scores[bc1[i],]*single_fitness_scores[bc2[i],]))
}))

dm_list[[length(dm_list) + 1]] <- dms
}



z_scores <- c()
for(i in 1:ncol(dm_list[[1]])){
  dm_vec <- c()
  for(j in 1:length(dm_list)){
    dm_vec <- cbind(dm_vec,dm_list[[j]][,i])
  }
  z_scores <- cbind(z_scores,apply(dm_vec,1,mean)/apply(dm_vec,1,sd))
}

sd_scores <- c()
for(i in 1:ncol(dm_list[[1]])){
  dm_vec <- c()
  for(j in 1:length(dm_list)){
    dm_vec <- cbind(dm_vec,dm_list[[j]][,i])
  }
  sd_scores <- cbind(sd_scores,apply(dm_vec,1,sd))#apply(dm_vec,1,mean)/apply(dm_vec,1,sd))
}

#max(nn)

#z_scores[z_scores < -1000] <- -1000
#z_scores[z_scores > 1000] <- -1000

#for(i in 1:length(dm_list)){
#  dm_vec_cond3 <- cbind(dm_vec_cond3,dm_list[[i]][,3])
#}



#z_scores <- apply(dm_vec_cond3,1,mean)/apply(dm_vec_cond3,1,sd)
gi_data[,grep('^GI',colnames(gi_data))] <- dm_list[[1]]#[,2:ncol(dms)]
z_cols <- grep('Class',grep('^Z_GIS',colnames(gi_data),val=T),invert=T,val=T)
gi_data[,z_cols] <- apply(dms,2,scale)

write.table(gi_data,quote=F,file='table_s1_new.tsv',sep='\t')

