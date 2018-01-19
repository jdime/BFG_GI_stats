pos_cut <- 0.05#0.089#0.146
neg_cut <- 0.05#0.01#-0.079
library(dplyr)
library(metap)

p_combine <- function(x){
  #if(max(x) == 1){
  #  return(1)
  #}
  #return(pnorm(sum(qnorm(x))/sqrt(length(x))) )
  return(mean(x))
  #if(length(x) < 2){
  #  return(x)
  #}
  #return(mean(x))
  #return(sumz(x)$p[1])
}

#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)
#setwd('../data')


#gi_data <-
#  read.table('table_s1_new.tsv', head = T, stringsAsFactors = F)
gi_data <-
  dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')#, C_ij.HetDipl > 100)

shu_complex <- c('CSM2','PSY3','SHU2','SHU1')



gi_data$gene1 <- sapply(gi_data$Barcode_i, function(x) {
    strsplit(x, split = '_')[[1]][1]
  })

gi_data$gene2 <- sapply(gi_data$Barcode_j, function(x) {
  strsplit(x, split = '_')[[1]][1]
})




gi_data$Pair <- sapply(1:length(gi_data$gene1),function(i){
  paste(sort(c(gi_data$gene1[i],gi_data$gene2[i])),collapse='_')
})
#gi_data$Barcode_j <- gene2


#gi_data <- gi_data %>% group_by(Pair) %>% select(c(1,2,grep('^GIS', colnames(gi_data)))) %>% summarize_each(funs(p_combine))

#stop()

shu_data <- filter(gi_data, gene1 %in% shu_complex & gene2 %in% shu_complex)

stop()

shu_data <- shu_data %>% group_by(Pair) %>% dplyr::select(grep('^FDR.Int', colnames(shu_data))) %>% summarize_all(funs(p_combine))
positive_conditions_shu <- apply(shu_data[,grep('^FDR.Int', colnames(shu_data))],1,function(x){names(which(x > pos_cut))})
names(positive_conditions_shu) <- shu_data$Pair

mag1_data <-
  dplyr::filter(
    gi_data,
    gene1 %in% c(shu_complex, 'SLX4') &
      gene2 == 'MAG1' |
      gene2 %in% c(shu_complex, 'SLX4') & gene1 == 'MAG1'
  )

stop()
mag1_data <-
  mag1_data %>% dplyr::group_by(Pair) %>% dplyr::select(grep('^FDR.Int', colnames(gi_data))) %>% summarize_all(funs(p_combine))
negative_conditions_mag1 <-
  apply(mag1_data[, grep('^FDR.Int', colnames(mag1_data))], 1, function(x) {
    names(which(x < neg_cut))
  })
names(negative_conditions_mag1) <- mag1_data$Pair

slx4_data <-
  filter(
    gi_data,
    gene1 %in% c(shu_complex, 'MAG1') &
      gene2 == 'SLX4' |
      gene2 %in% c(shu_complex, 'MAG1') & gene1 == 'SLX4'
  )
slx4_data <-
  slx4_data %>% group_by(Pair) %>% dplyr::select(grep('^FDR.Int', colnames(gi_data))) %>% summarize_all(funs(p_combine))
negative_conditions_slx4 <-
  apply(slx4_data[, grep('^FDR.Int', colnames(slx4_data))], 1, function(x) {
    names(which(x < neg_cut))
  })
names(negative_conditions_slx4) <- slx4_data$Pair

#mag1_data <- filter(gi_data, gene1 %in% shu_complex & gene2 == 'MAG1' | )


