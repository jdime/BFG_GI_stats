devtools::use_package('dplyr')

update_gis <- function(gi_data,
                       well_measured_cutoff = 100,
                       g_wt_vec = c(12.62,8.34,8.44,7.04,7.7,7.84,7.5,7.76,6.94,6.28)) {
  #Define special pairs
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
  
  
  count_data <- gi_data[, grep('^C_', colnames(gi_data))] + 1
  condition_sums <- apply(count_data, 2, sum)
  freq_data <- count_data / condition_sums
  
  
  f_ij_data <- count_data / condition_sums
  r_ij_data <- (f_ij_data / f_ij_data[, 1])[, 2:ncol(f_ij_data)]
  
  
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
      gi_data$C_ij.HetDipl >= well_measured_cutoff
    #wt_norm_freq_data[which(gi_data[,1]) == 1]
    return(apply(r_ij_data[criteria, ], 2, median))
  }))
  
  
  
  r_ij_wt <-
    apply(r_ij_median_single_genes[neutral_genes,], 2, median)
  
  rel_g_ij <- log2(t(t(r_ij_data) / r_ij_wt))
  
  #g_wt <- -apply(rel_g_ij,2,function(x){
  #  
  #  x <- x[!(same_same) & gi_data$C_ij.HetDipl >= well_measured_cutoff]
  #  min(x[is.finite(x)])
  #  
  #})
  
  if(is.null(g_wt_vec)){
  g_wt <-
     apply(-(rel_g_ij[same_same &
                        gi_data$C_ij.HetDipl >= well_measured_cutoff, ]), 2, median)
  }else{
    g_wt <- g_wt_vec
  }
  g_ij <- t(t(rel_g_ij) + g_wt)
  
  g_ij[!is.finite(g_ij)] <- 0
  
  w_ij_data <- t(t(g_ij) / g_wt)
  
  
  
  #Using mean here instead of median because
  #otherwise standard deviation doesn't make sense
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
      gi_data$C_ij.HetDipl >= well_measured_cutoff
    
    return(apply(w_ij_data[criteria, ], 2, mean))
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
      gi_data$C_ij.HetDipl >= well_measured_cutoff
      return(apply(w_ij_data[criteria, ], 2, sd))
  }))
  
  #Quick fix
  w_ij_data[w_ij_data < 0] <- 0
  
  bc1 <- gi_data$Barcode_i
  bc2 <- gi_data$Barcode_j
  
  gis <- t(sapply(1:nrow(w_ij_data), function(i) {
    w_i <- w_ij_single_genes[bc1[i], ]
    w_j <- w_ij_single_genes[bc2[i], ]
    
    w_ij <- w_ij_data[i, ]
    
    gis <- (w_ij) - (w_i * w_j)
    return(gis)
  }))
  
  gi_uncertainty <- t(sapply(1:nrow(w_ij_data), function(i) {
    w_i <- w_ij_single_genes[bc1[i], ]
    w_j <- w_ij_single_genes[bc2[i], ]
    
    
    w_i_error <- w_ij_error_single_genes[bc1[i], ]
    w_j_error <- w_ij_error_single_genes[bc2[i], ]
    
    w_ij <- w_ij_data[i, ]
    
    
    #Delta rule
    prod_uncertainty <-
      abs(w_i * w_j) * sqrt((w_i_error / w_i) ^ 2 + (w_j_error / w_j) ^ 2)
    
    return(prod_uncertainty)
  }))
  
  
  
  
  gi_data[, grep('^GI', colnames(gi_data))] <-
    gis
  z_cols <-
    grep('Class',
         grep('^Z_GIS', colnames(gi_data), val = T),
         invert = T,
         val = T)
  gi_data[, z_cols] <- gis / gi_uncertainty
  
  return(gi_data)
}