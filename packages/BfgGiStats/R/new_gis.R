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
  
  
  count_data <- gi_data[, grep('^C_', colnames(gi_data))]
  condition_sums <- apply(count_data, 2, sum)
  freq_data <- count_data / condition_sums
  
  
  f_xy_data <- count_data / condition_sums
  r_xy_data <- (f_xy_data / f_xy_data[, 1])[, 2:ncol(f_xy_data)]
  g_xy_data <- t(t(log2(r_xy_data)) + g_wt_vec)
  g_xy_wt <- apply(g_xy_data[nn_pairs,],2,median)
  w_xy_data <- t(t(g_xy_data) / g_xy_wt)
  
  genes <- unique(unlist(gi_data[, 1:2]))

  #Quick fix
  w_xy_data[w_xy_data < 0] <- 0
  
  #Using mean here instead of median because
  #otherwise standard deviation doesn't make sense
  w_xy_single_genes <- t(sapply(genes, function(gene) {
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
    
    return(apply(w_xy_data[criteria, ], 2, mean))
  }))
  
  w_xy_error_single_genes <- t(sapply(genes, function(gene) {
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
      return(apply(w_xy_data[criteria, ], 2, function(x){sd(x)/sqrt(length(x))}))
  }))
  
  #Quick fix
  
  
  bc1 <- gi_data$Barcode_i
  bc2 <- gi_data$Barcode_j
  
  gis <- t(sapply(1:nrow(w_xy_data), function(i) {
    w_x <- w_xy_single_genes[bc1[i], ]
    w_y <- w_xy_single_genes[bc2[i], ]
    
    w_xy <- w_xy_data[i, ]
    
    gis <- (w_xy) - (w_x * w_y)
    
    #gis <- log2((w_xy + 0.01)/((w_x * w_y) + 0.01))
    
    if(sum(is.na(gis) > 0)){
      print(w_xy)
      print(w_x)
      print(w_y)
      #print(c(w_xy,w_x,w_y))
    }
    
    return(gis)
  }))
  
  gi_uncertainty <- t(sapply(1:nrow(w_xy_data), function(i) {
    w_x <- w_xy_single_genes[bc1[i], ]
    w_y <- w_xy_single_genes[bc2[i], ]
    
    
    w_x_error <- w_xy_error_single_genes[bc1[i], ]
    w_y_error <- w_xy_error_single_genes[bc2[i], ]
    
    w_xy <- w_xy_data[i, ]
    
    
    #Delta rule
    prod_uncertainty <-
      abs(w_x * w_y) * sqrt((w_x_error / w_x) ^ 2 + (w_y_error / w_y) ^ 2)
    
    #gis <- ((w_xy + 0.01)/((w_x * w_y) + 0.01))
    
    #prod_uncertainty <- (prod_uncertainty)/log(2)
    
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