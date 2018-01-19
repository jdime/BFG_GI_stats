devtools::use_package('metap')

average_gi_data_by_gene <- function(gi_data,type_column_grep = 'Type_of_gene'){

  Gene1 <-
    sapply(gi_data$Barcode_i, function(name) {
      strsplit(name, split = '_')[[1]][1]
    })
  
  Gene2 <-
    sapply(gi_data$Barcode_j, function(name) {
      strsplit(name, split = '_')[[1]][1]
    })
  
  gene_columns <- cbind(Gene1,Gene2)
  
  #Alphabetically re-order genes for easier queries
  alpha_order <- t(apply(gene_columns,1,function(x){sort(x, index.return = T)$ix}))
  
  gene_columns <- t(sapply(1:nrow(gi_data),function(i){
    gene_columns[i,alpha_order[i,]]
  }))
  
  type_columns <- grep(type_column_grep,colnames(gi_data))
  
  gi_data[,type_columns] <- as.data.frame(t(sapply(1:nrow(gi_data),function(i){
    dat <- gi_data[i,type_columns]
    return(dat[alpha_order[i,]])
  })))
  
  Gene1 <- gene_columns[,1]
  Gene2 <- gene_columns[,2]
  
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
  
  
  ret_df <- c()
  for(i in 1:length(split_gi_data_grp)){
    testset <- split_gi_data_grp[[i]]
    
    if(nrow(testset) > 0){
      ret_vec <- testset[1,]
      
      gis_columns <- grep('^GIS',colnames(gi_data))
      err_columns <- grep('^SE_GIS',colnames(testset))
      z_columns <- grep('Class',grep('^Z_GIS',colnames(gi_data),val=T),val=T,invert=T)
      z_class_columns <- grep('Class',grep('^Z_GIS',colnames(gi_data),val=T),val=T,invert=F)
      fdr_columns <- grep('^FDR',colnames(gi_data))
      
      new_gis <- sapply(1:length(gis_columns),function(i){
        var_weights <- 1/(testset[,err_columns[i]]^2)
        weighted.mean(x = testset[,gis_columns[i]],w = var_weights)
        #mean(testset[,gis_columns[i]])
      })
      
      new_errs <- sapply(1:length(gis_columns),function(i){
        var_weights <- 1/(testset[,err_columns[i]]^2)
        sqrt(sum((testset[,err_columns[i]]*var_weights)^2))/sum(var_weights)
        #sqrt(sum((testset[,err_columns[i]])^2))/nrow(testset)
      })
      
      new_zs <- new_gis/new_errs
      
      new_fdrs <- sapply(1:length(fdr_columns),function(i){
        var_weights <- 1/(testset[,err_columns[i]]^2)

        p_vals <- testset[,fdr_columns[i]]
        if(length(p_vals) == 1){
          return(p_vals)
        }
        
        p_positive <- p_vals
        p_negative <- p_vals
        
        gi_scores <- testset[,gis_columns[i]]
        nominal_gi <- mean(testset[,gis_columns[i]])
        
        if(nominal_gi < 0){
          p_vals[gi_scores > 0] <- 1 - p_vals[gi_scores > 0]
        }else{
          p_vals[gi_scores < 0] <- 1 - p_vals[gi_scores < 0]
        }
        
        #Stops it from crashing
        p_vals[p_vals < 1e-200] <- 1e-200
        
        nomin_p <- metap::sumz(p_vals, weights = var_weights)$p[1]
        
        return(min(nomin_p, 1 - nomin_p)*2)

      })
    
      new_zs <- new_gis/new_errs
      
      ret_vec[gis_columns] <- new_gis
      ret_vec[z_columns] <- new_zs
      ret_vec[fdr_columns] <- new_fdrs
      ret_vec[err_columns] <- new_errs
      
      
      #Update w metrics
      w_xy_err_columns <- grep('W_ij_SE\\.',colnames(gi_data))       
      w_xy_columns <- grep('W_ij\\.',colnames(gi_data))       
      
      new_w_xy_errs <- sapply(1:length(w_xy_err_columns),function(i){
        var_weights <- 1/(testset[,w_xy_err_columns[i]]^2)
        sqrt(sum((testset[,w_xy_err_columns[i]]*var_weights)^2))/sum(var_weights)
      })
      
      new_w_xy <- sapply(1:length(w_xy_columns),function(i){
        var_weights <- 1/(testset[,w_xy_err_columns[i]]^2)
        weighted.mean(x = testset[,w_xy_columns[i]],w = var_weights)
        #mean(testset[,gis_columns[i]])
      })
      
      
      ret_vec[w_xy_err_columns] <- new_w_xy_errs
      ret_vec[w_xy_columns] <- new_w_xy
      
      
      ret_df <- rbind(ret_df,ret_vec)
    }
  }
  
  gi_data_grp <- ret_df
  
  gi_data_grp <- gi_data_grp[, preserved_colnames]
  
  gi_data_grp <- gi_data_grp[,grep('^C_ij',names(gi_data_grp),invert=T)]
  
  #stop()
  
  colnames(gi_data_grp)[1:2] <- colnames(gi_data[1:2])
  
  
  gi_data_grp[, 1] <- as.character(gi_data_grp[, 1])
  gi_data_grp[, 2] <- as.character(gi_data_grp[, 2])
  
  return(gi_data_grp)
}