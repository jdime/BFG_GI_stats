devtools::use_package('dplyr')

#' Update BFG-GI genetic interaction data
#'
#' @param gi_data input genetic interaction table
#' @param pseudocount integer, pseudocount to add in the calculations
#' @param count_data_grep_pattern grep pattern used to find count columns
#' @param g_wt_vec how many generations each pool grew, in the same order as the count columns
#'
#' @return an updated gi_data
update_gis <- function(gi_data,
                       #Old parameter, to use only well-measured pairs, no longer needed
                       #well_measured_cutoff = 100,
                       #This defines the number of generations each pool grew
                       pseudocount = 1,
                       count_data_grep_pattern = '^C_',
                       g_wt_vec = c(12.62,8.34,8.44,7.04,7.7,7.84,7.5,7.76,6.94,6.28)) {
  #Define special pairs
  same_same <-
    sapply(gi_data$Barcode_x, function(x) {
      strsplit(x, split = '_')[[1]][1]
    }) == sapply(gi_data$Barcode_y, function(x) {
      strsplit(x, split = '_')[[1]][1]
    })
  
  nn_pairs <-
    gi_data$Type_of_gene_x == 'Neutral' &
    gi_data$Type_of_gene_y == 'Neutral' &
    gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
  
  non_nn_unlinked_pairs <-
    gi_data$Type_of_gene_x != 'Neutral' &
    gi_data$Type_of_gene_y != 'Neutral' &
    gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
  
  #Add a pseudocount of 1
  count_data <- gi_data[, grep(count_data_grep_pattern, colnames(gi_data))] + pseudocount
  
  condition_sums <- apply(count_data, 2, sum)
  freq_data <- count_data / condition_sums
  
  
  #See methods for rationale of these metrics
  f_xy_data <- count_data / condition_sums
  r_xy_data <- (f_xy_data / f_xy_data[, 1])[, 2:ncol(f_xy_data)]
  g_xy_data <- t(t(log2(r_xy_data)) + g_wt_vec)
  
  #Count error estimated by the poisson distribution
  r_xy_error <-
    r_xy_data * condition_sums[2:length(condition_sums)]/condition_sums[1] *
    sqrt((sqrt(count_data[, 2:ncol(f_xy_data)]) /count_data[, 2:ncol(f_xy_data)])^2 +
           (sqrt(count_data[, 1]) / count_data[, 1])^2)
  
  #Error propagation to log
  g_xy_error <- abs(r_xy_error/(r_xy_data*log(2)))
  
  #wt fitness estimate and error
  g_xy_wt <- apply(g_xy_data[nn_pairs,],2,mean)
  g_wt_error <- apply(g_xy_data[nn_pairs,],2,function(x){sd(x)/sqrt(length(x))})
  
  #Normalize fitness by wildtype
  w_xy_data <- t(t(g_xy_data) / g_xy_wt)
  
  genes <- unique(unlist(gi_data[, 1:2]))

  
  #Using mean now instead of median because
  #otherwise standard deviation doesn't make sense
  w_xy_single_genes <- t(sapply(genes, function(gene) {
    criteria <-
      gi_data[, 1] == gene &
      gi_data$Type_of_gene_y == 'Neutral' |
      gi_data[, 2] == gene & gi_data$Type_of_gene_x == 'Neutral'
    criteria <-
      criteria &
      gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
    
    #Used to filter by counts, no need with error model
    #criteria <-
    #  criteria &
    #  gi_data$C_xy.HetDipl >= well_measured_cutoff
    
    return(apply(w_xy_data[criteria, ], 2, mean))
  }))
  
  g_xy_single_genes <- t(sapply(genes, function(gene) {
    criteria <-
      gi_data[, 1] == gene &
      gi_data$Type_of_gene_y == 'Neutral' |
      gi_data[, 2] == gene & gi_data$Type_of_gene_x == 'Neutral'
    criteria <-
      criteria &
      gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
    #Used to filter by counts, no need with error model
    #criteria <-
    #  criteria &
    #  gi_data$C_xy.HetDipl >= well_measured_cutoff
    
    return(apply(g_xy_data[criteria, ], 2, mean))
  }))
  
  
  #Using standard error of mean as error
  g_xy_error_single_genes <- t(sapply(genes, function(gene) {
    criteria <-
      gi_data[, 1] == gene &
      gi_data$Type_of_gene_y == 'Neutral' |
      gi_data[, 2] == gene & gi_data$Type_of_gene_x == 'Neutral'
    criteria <-
      criteria &
      gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
    #Used to filter by counts, no need with error model
    #criteria <-
    #  criteria &
    #  gi_data$C_xy.HetDipl >= well_measured_cutoff
    return(apply(g_xy_data[criteria, ], 2, function(x){sd(x)/sqrt(length(x))}))
  }))
  
  
  bc1 <- gi_data$Barcode_x
  bc2 <- gi_data$Barcode_y
  
  gis <- t(sapply(1:nrow(w_xy_data), function(i) {
    w_x <- w_xy_single_genes[bc1[i], ]
    w_y <- w_xy_single_genes[bc2[i], ]
    
    w_xy <- w_xy_data[i, ]
    
    #Quick fix
    #Set fitness to 0 if less than 0
    w_x[w_x < 0] <- 0
    w_y[w_y < 0] <- 0
    w_xy[w_xy < 0] <- 0
    
    gis <- (w_xy) - (w_x * w_y)
    
    #Alternative log-model of genetic interactions, not used
    #gis <- log2((w_xy + 0.01)/((w_x * w_y) + 0.01))
    
    #For debugging
    #if(sum(is.na(gis) > 0)){
    #  print(w_xy)
    #  print(w_x)
    #  print(w_y)
    #  #print(c(w_xy,w_x,w_y))
    #}
    
    return(gis)
  }))
  
  
  
  gi_uncertainty <- t(sapply(1:nrow(w_xy_data), function(i) {
    
    g_xy <- g_xy_data[i, ]
    g_x <- g_xy_single_genes[bc1[i], ]
    g_y <- g_xy_single_genes[bc2[i], ]
    g_wt <- g_xy_wt
    
    
    g_xy_error <- g_xy_error[i, ]
    g_x_error <- g_xy_error_single_genes[bc1[i], ]
    g_y_error <- g_xy_error_single_genes[bc2[i], ]
    g_wt_error <- g_wt_error
    
    gis <- (g_xy*g_wt)/(g_x*g_y)
    
    
    # numerator <- g_xy * g_wt
    # numerator_error <- numerator*sqrt((g_xy_error/g_xy)^2 + (g_wt_error/g_wt)^2)
    # 
    # denominator <- g_x * g_y
    # denominator_error <- denominator*sqrt((g_x_error/g_x)^2 + (g_y_error/g_y)^2)
    # 
    # ratio <- numerator / denominator
    # 
    # ratio_error <- ratio * sqrt((numerator_error/numerator)^2 + (denominator_error/denominator)^2)
    # 
    # log_ratio_error <- ratio_error/(ratio*log(2))
    # 
    # return(log_ratio_error)
    
    numerator <- g_xy - g_x - g_y

    gis <- numerator/g_wt

    numerator_error <- sqrt(g_xy_error^2 + g_x_error^2 + g_y_error^2)

    gi_error <- abs(gis)*sqrt((numerator_error/numerator)^2 + (g_wt_error/g_wt)^2)

    return(gi_error)
    
    
    
  }))
  
  
  
  
  
  #gis_global <<- gis
  #stop()
  
  #Modify data to replace old definition with new definitions
  colnames(gis) <- sapply(colnames(gis),function(name){
    gsub('^C_','GIS_',name)
  })
  
  
  gi_data <- cbind(gi_data,gis)
  
  
  z_scores <- gis / gi_uncertainty
  
  
  w_x_data <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_x <- g_xy_single_genes[bc1[i], ]
    g_wt <- g_xy_wt
    w_x <- g_x/g_wt
    w_x[w_x < 0] <- 0
    
    return(w_x)
  }))
  
  w_y_data <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_y <- g_xy_single_genes[bc2[i], ]
    g_wt <- g_xy_wt
    w_y <- g_y/g_wt
    w_y[w_y < 0] <- 0
    
    return(w_y)
  }))
  
  
  w_x_error <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_x <- g_xy_single_genes[bc1[i], ]
    g_wt <- g_xy_wt
    g_x_error <- g_xy_error_single_genes[bc1[i], ]
    g_wt_error <- g_wt_error
    
    w_x <- g_x/g_wt
    err <- abs(w_x)*sqrt((g_x_error/g_x)^2 + (g_wt_error/g_wt)^2)
    
    return(err)
  }))
  
  w_y_error <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_y <- g_xy_single_genes[bc2[i], ]
    g_wt <- g_xy_wt
    g_y_error <- g_xy_error_single_genes[bc2[i], ]
    g_wt_error <- g_wt_error
    
    
    
    
    w_y <- g_y/g_wt
    err <- abs(w_y)*sqrt((g_y_error/g_y)^2 + (g_wt_error/g_wt)^2)
    return(err)
  }))
  
  w_xy_error <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_xy <- g_xy_data[i, ]
    g_wt <- g_xy_wt
    
    g_xy_error <- g_xy_error[i, ]
    g_wt_error <- g_wt_error
    
    
    
    w_xy <- g_xy/g_wt
    
    err <- abs(w_xy)*sqrt((g_xy_error/g_xy)^2 + (g_wt_error/g_wt)^2)
    
    #err <- abs(w_y)*sqrt((g_y_error/g_y)^2 + (g_wt_error/g_wt)^2)
    return(err)
  }))
  
  
  
  print(colnames(w_xy_error))
  # Set Some column names
  colnames(w_x_data) <- colnames(w_x_data) %>% sapply(function(name){
    gsub('^C_xy.','W_x.', name)
  })
  colnames(w_y_data) <- colnames(w_y_data) %>% sapply(function(name){
    gsub('^C_xy.','W_y.', name)
  })
  colnames(w_xy_data) <- colnames(w_xy_data) %>% sapply(function(name){
    gsub('^C_xy.','W_xy.', name)
  })
  
  colnames(w_x_error) <- colnames(w_x_error) %>% sapply(function(name){
    gsub('^C_xy.','W_x_SE.', name)
  })
  colnames(w_y_error) <- colnames(w_y_error) %>% sapply(function(name){
    gsub('^C_xy.','W_y_SE.', name)
  })
  colnames(w_xy_error) <- colnames(w_x_error) %>% sapply(function(name){
    gsub('^W_x_SE.','W_xy_SE.', name)
  })
  
  
  
  
  
  
  #print(colnames(w_x_data))
  #
  #print(colnames(w_y_))
  print(colnames(w_xy_error))
  
  
  
  #colnames()
  
  
  colnames(z_scores) <- sapply(colnames(z_scores),function(name){
    gsub('^GIS_','Z_GIS_',name)
  })
  
  #Not returned in latest version
  colnames(gi_uncertainty) <- colnames(gis)
  colnames(gi_uncertainty) <- sapply(colnames(gi_uncertainty),function(name){
    gsub('^GIS_','SE_GIS_',name)
  })
  
  
  #Correct fitness
  w_xy_data[w_xy_data < 0] <- 0
  
  
  gi_data <- cbind(gi_data,
                   
                   w_x_data,
                   w_y_data,
                   w_xy_data,
                   
                   w_x_error,
                   w_y_error,
                   w_xy_error,
                   
                   gi_uncertainty,
                   z_scores)
  
  #Old way of updating the data
  #gi_data[, grep('^GI', colnames(gi_data))] <-
  #  gis
  #z_cols <-
  #  grep('Class',
  #       grep('^Z_GIS', colnames(gi_data), val = T),
  #       invert = T,
  #       val = T)
  #gi_data[, z_cols] <- gis / gi_uncertainty
  
  return(gi_data)
}