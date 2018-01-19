################################################################
################################################################
################################################################
########                                                ########
########                                                ########
########     Pipeline for update of BFG-GI data         ########
########                                                ########
########                                                ########
########                                                ########
################################################################
################################################################
################################################################
require(Cairo)
require(qvalue)
require(dplyr)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

devtools::load_all('../packages/BfgGiStats')
devtools::document('../packages/BfgGiStats')


setwd('../data')

#Read in genetic interaction data from two tables, using sum column from second table
gi_data <-
  read.table('table_s1_input.tsv', head = T, stringsAsFactors = F)
#Remove count data from first table
gi_data <- gi_data[,-grep('^C_',colnames(gi_data))]

#Update removal criteria to 100kb
gi_data[gi_data$Chromosomal_distance_bp <= 75000 &
          !is.na(gi_data$Chromosomal_distance_bp), 'Remove_by_Chromosomal_distance_or_SameGene'] <-
  'YES'

#Can change this to look at only one of the two technical replicates, analyzes the sum by default
#Change  to 'R1' or 'R2' to look at one of the two only
to_analyze <- 'sum'

grep_pattern1 <- sprintf('^C_ij.*%s',to_analyze)
grep_pattern2 <- sprintf('.%s',to_analyze)

gi_data_2 <-
  read.table('Cij_Table_WithReplicates.txt', head = T, stringsAsFactors = F)
gi_data <- cbind(gi_data,gi_data_2[,grep(grep_pattern1,colnames(gi_data_2))])
#Fix the names
colnames(gi_data) <- gsub(grep_pattern2,'',colnames(gi_data))

#Remove non well-measured strains
well_measured <- gi_data[,grep('HetDipl',colnames(gi_data))] >= 100
gi_data <- gi_data[well_measured, ]

#For troubleshooting/checking
#gi_data_original <- gi_data

#Change genetic interactions and significance calls
gi_data <- update_gis(gi_data,
                      #These values are empirically determined and have to be
                      #in the same order as in gi_data count columns
                      g_wt_vec = c(12.62,
                                   8.34,
                                   8.44,
                                   7.04,
                                   7.7,
                                   7.84,
                                   7.5,
                                   7.76,
                                   6.94,
                                   6.28))

#Add p values columns
setwd(this.dir)
setwd('../results')
Cairo::CairoPDF(file='Z_distribution.pdf',width=12,height=5)
par(mfrow=c(2,5))
par(mar=c(5,4,2,0))
gi_data <- add_p_values(gi_data,z_grep_pattern = "^Z_GIS",nn_pair_type = 'broad', make_plots = T)
dev.off()



#Preserve at this step for differential gi_calls
gi_data_old <- gi_data


#Make AUC plot
setwd(this.dir)
setwd('../results')

Cairo::CairoPDF(file='auc_vs_st_onge_new.pdf',width=7.5,height=4)
par(mfrow=c(1,2))
st_onge_auc_plot(gi_data)
dev.off()

#Score comparison scatterplot - by barcode
Cairo::CairoPDF(file='st_onge_scatterplot_barcodewise.pdf',width=7.5,height=3.5)
par(mfrow=c(1,2))
st_onge_scatterplot(gi_data)
dev.off()

########
#Calculate a per-gene GI score, Z score, and p value
########
gi_data <- average_gi_data_by_gene(gi_data)

stop()

Cairo::CairoPDF(file='gene_averaged_p_value_histogram.pdf',width=5.5,height=4.5)
hist(as.matrix(gi_data[,grep('^FDR',colnames(gi_data))]),main='',xlab='p-value',col='grey30')
dev.off()


#Calculate FDR from per-gene p values
gi_data[,grep('^FDR',colnames(gi_data))] <- 
  apply(gi_data[,grep('^FDR',colnames(gi_data))],2,function(x){qvalue(x)$q})


#Make St Onge MCC and precision plot
setwd(this.dir)
setwd('../results')


performance_data <- gi_data
performance_data[, grep('^FDR', colnames(performance_data))] <-
  -log10(performance_data[, grep('^FDR', colnames(performance_data))]) *
  sign(performance_data[, grep('^GI', colnames(performance_data))])
performance_vs_st_onge <- make_performance_matrix(performance_data)

Cairo::CairoPDF(file='pos_perf_vs_st_onge.pdf',width=5,height=4)
#Plots and stores which values are best MCC
optim_pos <- performance_heatmap(
  performance_vs_st_onge$positive_interactions,
  xlab = '-Log10(FDR) Threshold',
  ylab = 'GIS Threshold',
  min_val = 0.5,
  max_val = 0.7,
  highlight_max_vals = T,
  add_legend = T,
  legend_title = 'MCC',
  main = 'MCC vs St. Onge (positive GIs)'
)[1, , drop = F]
dev.off()

Cairo::CairoPDF(file='neg_perf_vs_st_onge.pdf',width=5,height=4)
optim_neg <- performance_heatmap(
  performance_vs_st_onge$negative_interactions,
  xlab = '-Log10(FDR) Threshold',
  ylab = 'GIS Threshold',
  min_val = 0.5,
  max_val = 0.7,
  highlight_max_vals = T,
  add_legend = T,
  legend_title = 'MCC',
  main = 'MCC vs St. Onge (negative GIs)'
)[1, , drop = F]
dev.off()


##Get cutoffs from performance matrix
log_fdr_cutoff_neg <- as.numeric(rownames(performance_vs_st_onge$negative_interactions)[optim_neg[1]])
gi_cutoff_neg <- as.numeric(colnames(performance_vs_st_onge$negative_interactions)[optim_neg[2]])

log_fdr_cutoff_pos <- as.numeric(rownames(performance_vs_st_onge$positive_interactions)[optim_pos[1]])
gi_cutoff_pos <- as.numeric(colnames(performance_vs_st_onge$positive_interactions)[optim_pos[2]])




setwd(this.dir)
setwd('../results')
Cairo::CairoPDF(file='gi_mcc_vs_st_onge.pdf',width=4.5,height=4)
optim_cutoffs <- precision_vs_stonge(gi_data,fdr_cutoff = 0.05,metr='mat',xlims=c(-10,10),ylab = "Matthew's Correlation Coefficient")
dev.off()

Cairo::CairoPDF(file='gi_prec_vs_st_onge.pdf',width=4.5,height=4)
precision_vs_stonge(gi_data,fdr_cutoff = 0.05,xlims=c(-10,10),cutoffs_drawn = optim_cutoffs)
dev.off()



#Update calls based on optimal cutoffs
gi_data <- update_calls(gi_data,
                        fdr_cutoff_pos = 10^-optim_cutoffs[2],#log_fdr_cutoff_pos,
                        gi_cutoff_pos = 0,#-gi_cutoff_pos,
                        fdr_cutoff_neg = 10^-optim_cutoffs[1],#log_fdr_cutoff_neg,
                        gi_cutoff_neg = 0)#-gi_cutoff_neg)

#Make gene-wise scatterplot
Cairo::CairoPDF(file='st_onge_scatterplot_genewise.pdf',width=7.5,height=3.5)
par(mfrow=c(1,2))
st_onge_scatterplot(gi_data)
dev.off()


###Write data
setwd(this.dir)
setwd('../data')
write.table(gi_data_old,row.names=F,quote=F,sep='\t',file='table_s1.tsv')
write.table(gi_data,row.names=F,quote=F,sep='\t',file='table_s1_genewise.tsv')

##Make/write differential calls
##Make differential calls, reporting everything
par(mfrow=c(7,7))
par(mar=c(3,2,4,1))
differential_calls <- differential_gi_analysis(gi_data,
                                               #Report everything
                                               fdr_cutoff = 1.1,
                                               require_sign_change = F,
                                               make_plots = T)

write.table(differential_calls,row.names=F,quote=F,sep='\t',file='table_s2_all.tsv')

#Plot some null distributions
setwd(this.dir)
setwd('../results')
Cairo::CairoPDF(file='delta_gi_distribution.pdf',width=5,height=3)
par(mar=c(4,4,1,1))
differential_calls_histogram(differential_calls)
dev.off()


#Filter for significance
differential_calls_sig <- differential_gi_analysis(gi_data,
                                               fdr_cutoff = 0.05,
                                               delta_gi_cutoff = 0,
                                               require_sign_change = T)
write.table(differential_calls_sig,row.names=F,quote=F,sep='\t',file='table_s2_significant.tsv')

