summarize_measurements <- function(gi_data_subset){
  ret_vec <- gi_data_subset[1,]
  
  print(nrow(gi_data_subset))
  
  gis_columns <- grep('^GIS',colnames(gi_data_subset))
  err_columns <- grep('^SE_GIS',colnames(gi_data_subset))
  z_columns <- grep('Class',grep('^Z_GIS',colnames(gi_data_subset),val=T),val=T,invert=T)
  z_class_columns <- grep('Class',grep('^Z_GIS',colnames(gi_data_subset),val=T),val=T,invert=F)
  
  
  #print('here')
  new_gis <- sapply(1:length(gis_columns),function(i){
    var_weights <- 1/(gi_data_subset[,err_columns[i]]^2)
    weighted.mean(x = gi_data_subset[,gis_columns[i]],w = var_weights)
    #mean(gi_data_subset[,gis_columns[i]])
  })
  #print('here2')
  
  new_errs <- sapply(1:length(gis_columns),function(i){
    var_weights <- 1/(gi_data_subset[,err_columns[i]]^2)
    sqrt(sum((gi_data_subset[,err_columns[i]]*var_weights)^2))/sum(var_weights)
    #sqrt(sum((gi_data_subset[,err_columns[i]])^2))/nrow(gi_data_subset)
  })
  #print('here3')
  
  new_zs <- sapply(1:length(z_columns),function(i){
    #old_zs <- gi_data_subset[,z_columns[i]]
    #old_zs[which.min(abs(old_zs))]
    #  mean(gi_data_subset[,z_columns[i]])
    #sum(gi_data_subset[,z_columns[i]])/sqrt(nrow(gi_data_subset))
    var_weights <- 1/(gi_data_subset[,err_columns[i]]^2)
    #sum(gi_data_subset[,z_columns[i]]*var_weights)/sqrt(sum(var_weights))
    sum(gi_data_subset[,z_columns[i]]*var_weights)/sqrt(sum(var_weights))
    #weighted.mean(x = gi_data_subset[,z_columns[i]],w = var_weights)
  })
  
  
  #new_zs <- new_gis/new_errs
  #print('here4')
  
  ret_vec[gis_columns] <- new_gis
  ret_vec[err_columns] <- new_errs
  ret_vec[z_columns] <- new_zs
  
  return(ret_vec)
}


Gene1 <-
  sapply(gi_data$Barcode_i, function(name) {
    strsplit(name, split = '_')[[1]][1]
  })

Gene2 <-
  sapply(gi_data$Barcode_j, function(name) {
    strsplit(name, split = '_')[[1]][1]
  })

gene_columns <- cbind(Gene1,Gene2)

gene_columns <- t(apply(gene_columns,1,sort))

Gene1 <- gene_columns[,1]
Gene2 <- gene_columns[,2]

pair <- apply(gene_columns,1,function(x){paste(x,collapse='_')})

preserved_colnames <- colnames(gi_data)
preserved_colnames <-
  preserved_colnames[grep('Barcode', preserved_colnames, invert = T)]
preserved_colnames <- c('Gene1', 'Gene2', preserved_colnames)

gi_data_grp <- cbind(gi_data, cbind(Gene1, Gene2))

split_gi_data_grp <- split(
  gi_data_grp,
  list(
    Gene1,
    Gene2
  )
)

n_measurements <- sapply(split_gi_data_grp,nrow)



split_gi_data_by_measurements <- list()
for(i in 1:max(n_measurements)){
  print(i)
  relevant_measurements <- split_gi_data_grp[n_measurements >= i]
  relevant_n <- n_measurements[n_measurements >= i]
  
  split_gi_data_by_measurements[[i]] <- list()
  
  for(j in 1:length(relevant_measurements)){
    n_fractions <- floor(relevant_n[j]/i)
    sub_scripts <- cbind((1:n_fractions)*i - (i-1),(1:n_fractions)*i)
    #print(sub_scripts)
    
    for(k in 1:nrow(sub_scripts)){
      sub_script <- sub_scripts[k,]
      
      relevant_measure <- relevant_measurements[[j]][sub_script[1]:sub_script[2],]
      
      
      split_gi_data_by_measurements[[i]] <- rbind(split_gi_data_by_measurements[[i]],summarize_measurements(relevant_measure))
    }
    #for(k in 1:n_fractions){
    #  sub_script <- ((i-1)*k +1):(i*k)
    #    
    #    #((i-1)*k + 1):(i*k[i])
    #  
    #
    #  split_gi_data_by_measurements[[i]] <- rbind(split_gi_data_by_measurements[[i]],relevant_measurements[[j]][sub_script,])
    #}
  }
  #n_fractions <- floor(relevant_n/i)
  
}


new_summary <- lapply(split_gi_data_by_measurements,function(gi_data){
  add_p_values(gi_data,z_grep_pattern = "^Z_GIS",nn_pair_type = 'broad')
  #non_nn_unlinked_pairs <-
  #  gi_data$Type_of_gene_i != 'Neutral' &
  #  gi_data$Type_of_gene_j != 'Neutral' &
  #  gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
  
  #gi_data <-
  #  dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no' & (Chromosomal_distance_bp > 1e05 | is.na(Chromosomal_distance_bp)))
  #nn_pairs <-
  #  gi_data$Type_of_gene_i != 'DNA_repair' | gi_data$Type_of_gene_j != 'DNA_repair'
  #hist(as.matrix(gi_data[nn_pairs,z_columns]),breaks=100)
})

pair_vals <- n_measurements[n_measurements > 1]
pair_names <- names(pair_vals)

ret_df <- c()
for(i in 1:length(pair_vals)){
  pair <- strsplit(pair_names[i],split='\\.')[[1]]
  val <- pair_vals[i]
  
  relevant_df <- new_summary[[val]]
  
  eyo <- dplyr::filter(relevant_df,Gene1 == pair[1] & Gene2 == pair[2])
  
  ret_df <- rbind(ret_df,eyo)
}