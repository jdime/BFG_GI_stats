make_performance_matrix <- function(gi_data,
                                    control_name = "NoDrug",
                                    condition_name = "MMS",
                                    x_axis_grep = "FDR.Internal_ij",
                                    y_axis_grep = "^GIS_ij",
                                    metr = 'mat',
                                    x_cutoff_vec = c(0:50)/5,
                                    y_cutoff_vec = c(0:20)/200,
                                    st_onge_class_to_analyze = 'both')
                                    {
  negative_result_matrix <- matrix(nrow=length(x_cutoff_vec),ncol=length(y_cutoff_vec))
  positive_result_matrix <- negative_result_matrix
  
  if(st_onge_class_to_analyze == 'MMS'){
    conditions <- condition_name
  }
  if(st_onge_class_to_analyze == 'DMSO'){
    conditions <- control_name
  }
  if(st_onge_class_to_analyze == 'both'){
    conditions <- c(control_name,condition_name)
  }
  
  x_scores <- c()
  y_scores <- c()
  labels_pos <- c()
  labels_neg <- c()
  for (condition in conditions) {
    if (condition == control_name) {
      st_onge_class <- 'SOJ_Class_NoMMS'
    }
    if (condition == condition_name) {
      st_onge_class <- 'SOJ_Class_MMS'
    }
    
    gi_data_filtered <-
      dplyr::filter(gi_data,
                    SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
    
    relevant_x <- grep(x_axis_grep,colnames(gi_data_filtered),val=T)
    relevant_x <- grep(condition,relevant_x,val =T)
    
    relevant_y <- grep(y_axis_grep,colnames(gi_data_filtered),val=T)
    relevant_y <- grep(condition,relevant_y,val =T)
    
    
    
    x_vals <- c(gi_data_filtered[ , relevant_x])
    y_vals <- c(gi_data_filtered[ , relevant_y])
    
    x_scores <- c(x_scores, x_vals)
    y_scores <- c(y_scores, y_vals)
    
    labels_pos <-
      c(labels_pos, gi_data_filtered[, st_onge_class] == 'ALLEVIATING')
    labels_neg <-
      c(labels_neg, gi_data_filtered[, st_onge_class] == 'AGGRAVATING')
  }
  
  
  for (i in 1:length(x_cutoff_vec)) {
    for (j in 1:length(y_cutoff_vec)) {
      x_cutoff <- x_cutoff_vec[i]
      y_cutoff <- y_cutoff_vec[j]
      
      scores_cond_pos <-
        as.numeric(x_scores >= x_cutoff & y_scores >= y_cutoff)
      scores_cond_neg <-
        as.numeric(x_scores <= -x_cutoff & y_scores <= -y_cutoff)
      pos_perf <-
        ROCR::performance(ROCR::prediction(scores_cond_pos, labels_pos), metr)@y.values[[1]][2]
      neg_perf <-
        ROCR::performance(ROCR::prediction(scores_cond_neg, labels_neg), metr)@y.values[[1]][2]
      
      
      positive_result_matrix[i, j] <- pos_perf
      negative_result_matrix[i, j] <- neg_perf
    }
  }
  
  rownames(positive_result_matrix) <- x_cutoff_vec
  colnames(positive_result_matrix) <- y_cutoff_vec
  
  rownames(negative_result_matrix) <- x_cutoff_vec
  colnames(negative_result_matrix) <- y_cutoff_vec
  
  retval <- list(positive_result_matrix,negative_result_matrix)
  names(retval) <- c('positive_interactions','negative_interactions')
  
  return(retval)
}

performance_heatmap <- function(mat,
                                color_func = NULL,
                                min_val = NULL,
                                max_val = NULL,
                                xlab = NULL,
                                ylab = NULL,
                                main = NULL,
                                highlight_max_vals = T,
                                highlighted_co_ordinates = NULL,
                                add_legend = F,
                                legend_title = NULL
){
  
  map_color <- function(val,color_func,min_val,max_val,resolution=100){
    colors <- color_func(resolution)
    vals <- seq(min_val,max_val,length.out=resolution)
    return(colors[which.min(abs(val - vals))])
  }
  
  
  if(is.null(color_func)){
    # Master colour scale for heatmap-style plots
    my_color_list <- c(
      rgb(1, 0.45, 0.25),
      rgb(0.8, 0.25, 0.25),
      rgb(0, 0, 0),
      rgb(0.25, 0.45, 0.8),
      rgb(0.25, 0.75, 1)
    )
    color_func <- grDevices::colorRampPalette(my_color_list)
    
  }
  if(is.null(min_val)){
    min_val <- min(mat)  
  }
  if(is.null(max_val)){
    max_val <- max(mat)
  }
  x_vals <- as.numeric(rownames(mat))
  y_vals <- as.numeric(colnames(mat))
  
  start_x <- min(x_vals)
  end_x <- max(x_vals)
  
  start_y <- min(y_vals)
  end_y <- max(y_vals)
  
  print(y_vals)
  par(bty = 'n') 
  
  if(add_legend){
    layout(matrix(c(1,1,2,3), 2, 2, byrow = F),
           widths=c(3,1), heights=c(2,1))
  }
  par(mar=c(5,5,2,0))
  plot(NULL,xlim=c(start_x, end_x),ylim=c(start_y,end_y),xlab=xlab,ylab=ylab,main=main)
  
  
  row_num <- nrow(mat)
  col_num <- ncol(mat)
  
  max_list <- c()
  
  for(i in 1:row_num){
    for(j in 1:col_num){
      xleft <- ((i - 1)/row_num)*end_x + start_x
      xright <- (i/row_num)*end_x + start_x
      ytop <- (j/col_num)*end_y + start_y
      ybot <- ((j-1)/col_num)*end_y + start_y
      
      
      fill_col <- map_color(mat[i,j],color_func,min_val,max_val)
      
      rect(xleft,ybot,xright,ytop,border=NA,col=fill_col)
      
      if(mat[i,j] == max(mat)){
        max_list <- rbind(max_list,c(i,j))
      }
    }
  }
  
  if(highlight_max_vals){
    for (row in 1:nrow(max_list)) {
      i <- max_list[row, 1]
      j <- max_list[row, 2]
      
      xleft <- ((i - 1) / row_num) * end_x + start_x
      xright <- (i / row_num) * end_x + start_x
      ytop <- (j / col_num) * end_y + start_y
      ybot <- ((j - 1) / col_num) * end_y + start_y
      points(mean(c(xleft, xright)), mean(c(ytop, ybot)), col = 'red', pch =
               8)
    }
  }
  
  if(!is.null(nrow(highlighted_co_ordinates))){
    for(row in 1:nrow(highlighted_co_ordinates)){
      i <- highlighted_co_ordinates[row, 1]
      j <- highlighted_co_ordinates[row, 2]
      
      xleft <- ((i - 1) / row_num) * end_x + start_x
      xright <- (i / row_num) * end_x + start_x
      ytop <- (j / col_num) * end_y + start_y
      ybot <- ((j - 1) / col_num) * end_y + start_y
      points(mean(c(xleft, xright)), mean(c(ytop, ybot)), col = 'red', pch =
               8)
    }
  }
  
  if(add_legend){
    par(xpd = T)
    plot(NULL,xlim=c(0,1),ylim=c(0,1),main='',axes=F,xlab='',ylab='')
    
    seqs <- seq(min_val,max_val,length.out=100)
    
    for(i in 1:length(seqs)){
      fill_col <- map_color(seqs[i],color_func,min_val,max_val)
      lines(c(-0.5,0),c(i/length(seqs),i/length(seqs)),col=fill_col,lwd=5)
      
    }
    
    text(0,0,min_val,adj=-0.5)
    text(0,0.5,mean(c(min_val,max_val)),adj=-0.5)
    text(0,1,max_val,adj=-0.5)
    text(-0.25,1.1,legend_title,adj=c(0.5))
  }
  return(max_list)
}

performance_data <- gi_data
performance_data[, grep('^FDR', colnames(performance_data))] <-
  -log10(performance_data[, grep('^FDR', colnames(performance_data))]) *
  sign(performance_data[, grep('^GI', colnames(performance_data))])

performance_vs_st_onge <- make_performance_matrix(performance_data)
optim_pos <- performance_heatmap(
  performance_vs_st_onge$positive_interactions,
  xlab = '-Log10(FDR) Threshold',
  ylab = '|GIS| Threshold',
  min_val = 0.5,
  max_val = 0.7,
  highlight_max_vals = T,
  add_legend = T,
  legend_title = 'MCC',
  main = 'Performance vs St. Onge (positive GIs)'
)[1, , drop = F]

optim_neg <- performance_heatmap(
  performance_vs_st_onge$negative_interactions,
  xlab = '-Log10(FDR) Threshold',
  ylab = '|GIS| Threshold',
  min_val = 0.5,
  max_val = 0.7,
  highlight_max_vals = T,
  add_legend = T,
  legend_title = 'MCC',
  main = 'Performance vs St. Onge (negative GIs)'
)[1, , drop = F]

