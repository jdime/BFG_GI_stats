this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')

gi_data <-
  read.table('table_s1.tsv', head = T, stringsAsFactors = F)

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


gis_list <- list()
nobs <- nrow(count_data_real)

for (i in 1) {
  print(i)
  if(i != 1){
  count_data <- sapply(1:ncol(count_data_real), function(i) {
    rbinom(nobs, condition_sums_real[i], freq_data_real[, i])
  }) + 1
  }
  
  else{
    count_data <- count_data_real# + 1
  }
  #count_data_real# + 1
  
  condition_sums <- apply(count_data, 2, sum)
  f_ij_data <- count_data / condition_sums
  
  r_ij_data <- (f_ij_data / f_ij_data[, 1])[, 2:ncol(f_ij_data)]
  
  
  #norm_freq_data <- (freq_data / freq_data[, 1])[, 2:ncol(freq_data)]
  
  #
  #
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
  
  
  r_ij_median_single_genes <- t(sapply(genes, function(gene) {
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
    return(apply(r_ij_data[criteria,], 2, median))
  }))
  
  
  
  r_ij_wt <- apply(r_ij_median_single_genes[neutral_genes, ], 2, median)
  
  rel_g_ij <- log2(t(t(r_ij_data) / r_ij_wt))
  
  #stop()
  
  g_wt <-
    apply(-(rel_g_ij[same_same & gi_data$C_ij.HetDipl > 100, ]), 2, median)
  
  g_ij <- t(t(rel_g_ij) + g_wt)
  
  w_ij_data <- t(t(g_ij) / g_wt)
  
  #Quick fix
  
  
  
  w_ij_single_genes <- t(sapply(genes, function(gene) {
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
    if(gene == 'RTT101_recipient_2'){
      print(w_ij_data[criteria,])
    }
    return(apply(w_ij_data[criteria,], 2, median))
  }))
  
  w_ij_error_single_genes <- t(sapply(genes, function(gene) {
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
     if(gene == 'RTT101_recipient_2'){
       print(w_ij_data[criteria,])
     }
     return(apply(w_ij_data[criteria,], 2, sd))
   }))
  
  #Quick fix
  w_ij_data[w_ij_data < 0] <- 0
  
  bc1 <- gi_data$Barcode_i
  bc2 <- gi_data$Barcode_j
  gis <- t(sapply(1:nrow(w_ij_data), function(i) {
    #log2(norm_dm_fit_new[i, ] / (single_fitness_scores_new[bc1[i], ] * single_fitness_scores_new[bc2[i], ]))
    #w_ij_data[i, ] - (w_ij_single_genes[bc1[i], ] * w_ij_single_genes[bc2[i], ])
    w_i <- w_ij_single_genes[bc1[i],]
    w_j <- w_ij_single_genes[bc2[i],]
    
    
    w_i_error <- w_ij_error_single_genes[bc1[i],]
    w_j_error <- w_ij_error_single_genes[bc2[i],]
    
    w_ij <- w_ij_data[i,]
    
    #Delta rule
    prod_uncertainty <-
      abs(w_i * w_j) * sqrt((w_i_error / w_i) ^ 2 + (w_j_error / w_j) ^ 2)
    #gi_uncertainty <- sqrt(fitness_uncertainty ^ 2 + prod_uncertainty ^ 2)
    gis <- (w_ij) - (w_i * w_j)
    
    
    return(gis/prod_uncertainty)
    #gis[abs(gis) < prod_uncertainty] <- 0
    #gis[gis < 0] <- 1#gis[gis < 0] + prod_uncertainty
    #gis[gis > 0] <- gis[gis > 0] - prod_uncertainty
    
    gis <- sapply(1:length(gis),function(i){
      print(prod_uncertainty[i])
      if(abs(gis[i]) < prod_uncertainty[i]){
        return(0)
      }else if(gis[i] > 0){
        return(gis[i] - prod_uncertainty[i])
      }else{
        return(gis[i] + prod_uncertainty[i])
      }
      #return(gis[i])
      #if(gis[i] < 0){
        
        #gis[i] <- gis[i] + prod_uncertainty[i]
      #}
      #else{
        #gis[i] <- gis[i] - prod_uncertainty[i]
      #}
    })
    
  
    
    
    
    return(gis)
    
    #if(abs(gis) < prod_uncertainty){
    #  return(0)
    #}
    
    #error_gis <- sapply(gis,function(gi){
    #  if(gi < 0){
    #   error_ratio <- gi + gi_uncertainty
    #   error_ratio <- min(error_ratio,0)
    #  }
    #  if(gi > 0){
    #   error_ratio <- gi - gi_uncertainty
    #   error_ratio <- max(error_ratio,0)
    #  }
    #  return(error_ratio)
    #})
    
    #return(error_gis)
    #prod <- wi
    
    
    # gis_ratio <- w_ij/(w_i*w_j)
    # ratio_uncertainty <- abs(gis_ratio)*sqrt((fitness_uncertainty/w_ij)^2+(prod_uncertainty/(w_i*w_j))^2)
    #
    
    #
    #
    # return(error_ratio)
    
    #
    
    
    
    
    
    #gis <- (w_ij + fitness_uncertainty) - (w_i + fitness_uncertainty * w_j + fitness_uncertainty)
    
    #gis <- ((w_ij + fitness_uncertainty)/(w_i*w_j + prod_uncertainty))
    
    #gis[abs(gis) < gi_uncertainty] <- 0
    #if(abs(gis) < gi_uncertainty){
    #  return(0)
    #}
    #return(gis)#_uncertainty)
    
    #return(gi_uncertainty)
    #log2(w_ij_data[i, ]/(w_ij_single_genes[bc1[i], ] * w_ij_single_genes[bc2[i], ]))
  }))
  gis_list[[length(gis_list) + 1]] <- gis
}

gis <- gis_list[[1]]
gi_uncertainty <- sapply(1:nrow(gis),function(i){
  sapply(1:ncol(gis),function(j){
    gi_vals <- sapply(gis_list,function(x){x[i,j]})
    return(sd(gi_vals))
  })
})


#stop()
#stop()


gi_data[, grep('^GI', colnames(gi_data))] <-
  gis#[,2:ncol(dms)]
z_cols <-
  grep('Class',
       grep('^Z_GIS', colnames(gi_data), val = T),
       invert = T,
       val = T)
gi_data[, z_cols] <- apply(gis, 2, scale)

write.table(gi_data,
            quote = F,
            file = 'table_s1_new.tsv',
            sep = '\t')