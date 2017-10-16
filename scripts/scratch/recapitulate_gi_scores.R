this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../data')

gi_data <- read.table('table_s1.tsv', head = T, stringsAsFactors = F)
gi_data <-
  dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')

count_data <- gi_data[, grep('^C_', colnames(gi_data))]

condition_sums <- apply(count_data, 2, sum)



freq_data <- count_data / condition_sums

pseudo_counts <- 1 / condition_sums


genes <- unique(unlist(gi_data[, 1:2]))
neutral_genes <-
  unique(unlist(
    dplyr::filter(
      gi_data,
      Type_of_gene_i == 'Neutral' & Type_of_gene_j == 'Neutral'
    )[, 1:2]
  ))


marginal_counts_original <- sapply(genes, function(gene) {
  apply(freq_data[which(gi_data[, 1] == gene |
                          gi_data[, 2] == gene), ], 2, sum)
})

marginal_counts <- t(marginal_counts_original + pseudo_counts)

relative_marginals <-
  marginal_counts / marginal_counts[, 'C_ij.HetDipl']

freq_data_dm <- t(t(freq_data) + pseudo_counts)
relative_dms <- freq_data_dm / freq_data_dm[, 'C_ij.HetDipl']


neutral_freq <- freq_data[gi_data$Type_of_gene_i == 'Neutral' & gi_data$Type_of_gene_j == 'Neutral',]
neutral_freq <- apply(neutral_freq,2,sum) + pseudo_counts
wt_norm <- neutral_freq/neutral_freq[1]




wt_normalized_relative_marginals <-
  t(apply(relative_marginals, 1, function(x) {
    x / wt_norm
  }))


wt_normalized_relative_dms <- t(apply(relative_dms, 1, function(x) {
  x / wt_norm
}))


#sm_wt <- mean(wt_normalized_relative_marginals[neutral_genes, 2:ncol(wt_normalized_relative_marginals)])

wt_normalized_relative_marginals <- wt_normalized_relative_marginals#/sm_wt
wt_normalized_relative_dms <- wt_normalized_relative_dms#/sm_wt
#wt_relative_marginals <-


#new_gis <- cxzc


new_gis <- t(sapply(1:nrow(gi_data), function(i) {
  log2(wt_normalized_relative_dms[i, ]/(wt_normalized_relative_marginals[gi_data[i, 1], ] *
    wt_normalized_relative_marginals[gi_data[i, 2], ]))
}))

#new_gis <- t(sapply(1:nrow(gi_data), function(i) {
#  wt_normalized_relative_dms[i, ] - (wt_normalized_relative_marginals[gi_data[i, 1], ] *
#    wt_normalized_relative_marginals[gi_data[i, 2], ])
#}))

g1 <- as.vector(new_gis[,2:ncol(new_gis)])
g2 <- as.vector(unlist(gi_data[,grep('^GI',colnames(gi_data))]))


smoothScatter(g1,g2,col=rgb(0,0,0,0.1),pch=16,xlab='Recreated GIS',ylab='Original GIS')
abline(c(0,1),col='red',lwd=3)


#hist(wt_normalized_relative_marginals[neutral_genes,2:ncol(wt_normalized_relative_marginals)])

#hist(wt_normalized_relative_marginals[neutral_genes,2:ncol(wt_normalized_relative_marginals)],xlab='Neutral-Neutral Single Mutant Fitness')

sm_wt <- wt_normalized_relative_marginals[neutral_genes, 2:ncol(wt_normalized_relative_marginals)][,'C_ij.MMS']
hist(
  sm_wt,
  xlab = 'Neutral-Neutral Single Mutant Fitness',
  main = '',
  col = 'grey'
)
abline(
  v = mean(sm_wt),
  col = 'red',
  lwd = 2
)

dm_wt <- wt_normalized_relative_dms[gi_data$Type_of_gene_i == 'Neutral' &
                                      gi_data$Type_of_gene_j == 'Neutral', 2:ncol(wt_normalized_relative_dms)][,'C_ij.MMS']

hist(dm_wt,
     xlab = 'Neutral-Neutral Double Mutant Fitness',
     main = '',
     col='grey'
     )
abline(
  v = mean(dm_wt),
  col = 'red',
  lwd = 2
)

gi_data[,grep('^GI',colnames(gi_data))] <- new_gis[,2:ncol(new_gis)]
z_cols <- grep('Class',grep('^Z_GIS',colnames(gi_data),val=T),invert=T,val=T)
gi_data[,z_cols] <- apply(new_gis[,2:ncol(new_gis)],2,scale)

write.table(gi_data,quote=F,file='table_s1_new.tsv',sep='\t')

