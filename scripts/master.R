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

#Can change this to look at only one of the two technical replicates, analyzes the sum by default
#Change grep pattern to '^C_ij.*R1' or '^C_ij.*R2' to look at one of the two only

gi_data_2 <-
  read.table('Cij_Table_WithReplicates.txt', head = T, stringsAsFactors = F)
gi_data <- cbind(gi_data,gi_data_2[,grep('^C_ij.*sum',colnames(gi_data_2))])
#Fix the names
colnames(gi_data) <- gsub('.sum','',colnames(gi_data))

#Remove non well-measured strains
well_measured <- gi_data[,grep('HetDipl',colnames(gi_data))] >= 100
gi_data <- gi_data[well_measured, ]


gi_data_old <- gi_data


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


#Add FDR columns
gi_data <- add_fdrs(gi_data)


#Update calls based on FDR columns
gi_data <- update_calls(gi_data,use_z = F,fdr_cutoff = 0.05)


##Making several plots
setwd(this.dir)
setwd('../results')
Cairo::CairoPDF(file='gi_fdr_prec_vs_st_onge.pdf',width=4.5,height=4)
precision_vs_stonge(gi_data)
dev.off()

Cairo::CairoPDF(file='auc_vs_st_onge_new.pdf',width=7.5,height=4)
par(mfrow=c(1,2))
st_onge_auc_plot(gi_data)
dev.off()

#Cairo::CairoPDF(file='auc_vs_st_onge_old.pdf',width=7.5,height=4)
#par(mfrow=c(1,2))
#st_onge_auc_plot(gi_data_old,old_data=T)
#dev.off()

Cairo::CairoPDF(file='st_onge_scatterplot.pdf',width=7.5,height=3.5)
par(mfrow=c(1,2))
st_onge_scatterplot(gi_data)
dev.off()

#Update the table
setwd(this.dir)
setwd('../data')
write.table(gi_data,row.names=F,quote=F,sep='\t',file='table_s1.tsv')

#Make differential calls
differential_calls <- differential_gi_analysis(gi_data,
                                               #Report everything
                                               fdr_cutoff = 2,
                                               require_sign_change = F)
write.table(differential_calls,row.names=F,quote=F,sep='\t',file='table_s2.tsv')