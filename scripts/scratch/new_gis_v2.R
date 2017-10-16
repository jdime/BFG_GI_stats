this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')

gi_data <-
  read.table('table_s1.tsv', head = T, stringsAsFactors = F)
#gi_data <-
#  dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')

same_same <-
  sapply(gi_data$Barcode_i, function(x) {
    strsplit(x, split = '_')[[1]][1]
  }) == sapply(gi_data$Barcode_j, function(x) {
    strsplit(x, split = '_')[[1]][1]
  }) 


nn_pairs <-
  gi_data$Type_of_gene_i == 'Neutral' &
  gi_data$Type_of_gene_j == 'Neutral' &
  gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'

non_nn_unlinked_pairs <-
  gi_data$Type_of_gene_i != 'Neutral' &
  gi_data$Type_of_gene_j != 'Neutral' &
  gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'



count_data_real <- gi_data[, grep('^C_', colnames(gi_data))] + 1

condition_sums_real <- apply(count_data_real, 2, sum)
freq_data_real <- count_data_real / condition_sums_real




count_data <- count_data_real# + 1

condition_sums <- apply(count_data, 2, sum)
freq_data <- count_data / condition_sums



#pseudo_counts <- 1 / condition_sums

#freq_data <- t(t(freq_data) + pseudo_counts)

norm_freq_data <- (freq_data / freq_data[, 1])[, 2:ncol(freq_data)]



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

single_fitness_scores <- t(sapply(genes, function(gene) {
  criteria <-
    gi_data[, 1] == gene &
    gi_data$Type_of_gene_j == 'Neutral' |
    gi_data[, 2] == gene & gi_data$Type_of_gene_i == 'Neutral'
  criteria <-
    criteria &
    gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
  criteria <-
    criteria &
    gi_data$C_ij.HetDipl > 100
  #wt_norm_freq_data[which(gi_data[,1]) == 1]
  return(apply(norm_freq_data[criteria, ], 2, mean))
}))

#stop()

#single_fitness_scores <- single_fitness_scores/single_fitness_scores[,1]
wt_fit <- apply(single_fitness_scores[neutral_genes,], 2, mean)





single_fitness_scores <- t(t(single_fitness_scores) / wt_fit)


norm_dm_fit <- t(t(norm_freq_data) / wt_fit)

#stop()

#Use same-same gene pairs to infer number of generations
#wt_gens <-
#  abs(apply(log2(norm_dm_fit[same_same &
#                               gi_data$C_ij.HetDipl > 100, ]), 2, mean))#abs(log2(apply(norm_dm_fit,2,min)))

#stop()
#wt_gens <- 20#wt_gens*

dm_vec <- c()
gens_tested <- c(100:500)/50
for(wt_gens in 1){#gens_tested){#gens_tested){
#wt_fit_new <- t(t(log2(wt_fit))/wt_gens) + 2
  #wt_gens <- c(6.6,6.54,6.92,6.66,6.12,6.24,6.24,6.4,6.9,6.36)#rep(wt_gens,ncol(norm_dm_fit))
  wt_fit <-
    abs(apply(log2(norm_dm_fit[same_same &
                                 gi_data$C_ij.HetDipl > 100, ]), 2, mean))
  wt_gens <- wt_fit#/wt_gens
  wt_fit <- 1
  print(wt_fit)
  single_fitness_scores_new <- t(t(log2(single_fitness_scores))/wt_gens + wt_fit)
  norm_dm_fit_new <- t(t(log2(norm_dm_fit))/wt_gens + wt_fit)

#single_fitness_scores_new[single_fitness_scores_new < 0.05] <- 0.05
#norm_dm_fit_new[norm_dm_fit_new < 0.05] <- 0.05
single_fitness_scores_new[single_fitness_scores_new < 0] <- 0
norm_dm_fit_new[norm_dm_fit_new < 0] <- 0




bc1 <- gi_data$Barcode_i
bc2 <- gi_data$Barcode_j
dms <- t(sapply(1:nrow(norm_dm_fit), function(i) {
  #log2(norm_dm_fit_new[i, ] / (single_fitness_scores_new[bc1[i], ] * single_fitness_scores_new[bc2[i], ]))
  norm_dm_fit_new[i, ] - (single_fitness_scores_new[bc1[i], ] * single_fitness_scores_new[bc2[i], ])
}))

#Hacky fix
dms[!is.finite(dms)] <- min(dms[is.finite(dms)])


dm_vec <- c(dm_vec,max(abs(apply(dms[nn_pairs &
                                       gi_data$C_ij.HetDipl > 0,10,drop=F],2,mean))))
}

#stop()


gi_data[, grep('^GI', colnames(gi_data))] <-
  dms#[,2:ncol(dms)]
z_cols <-
  grep('Class',
       grep('^Z_GIS', colnames(gi_data), val = T),
       invert = T,
       val = T)
gi_data[, z_cols] <- apply(dms, 2, scale)

write.table(gi_data,
            quote = F,
            file = 'table_s1_new.tsv',
            sep = '\t')
