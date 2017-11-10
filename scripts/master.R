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

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

devtools::load_all('../packages/BfgGiStats')
devtools::document('../packages/BfgGiStats')

setwd('../data')

#Read in genetic interaction data from two tables, using sum column from second table
gi_data <-
  read.table('table_s1.tsv', head = T, stringsAsFactors = F)

gi_data <- gi_data[,-grep('^C_',colnames(gi_data))]

#Can change this to look at only one of the two technical replicates, analyzes the sum by default
gi_data_2 <-
  read.table('TableForAlbi_Cij_WithReplicates.txt', head = T, stringsAsFactors = F)
gi_data <- cbind(gi_data,gi_data_2[,grep('^C_ij.*sum',colnames(gi_data_2))])


#Used for various functions
z_cols <-
  grep('Class',
       grep('^Z_GIS', colnames(gi_data), val = T),
       invert = T,
       val = T)
z_cols <- which(colnames(gi_data) %in% z_cols)

gi_data_old <- gi_data


#Change genetic interactions and significance calls
gi_data <- update_gis(gi_data,g_wt_vec = c(12.62,8.34,8.44,7.04,7.7,7.84,7.5,7.76,6.94,6.28))

#Add FDR columns
gi_data <- add_fdrs(gi_data)

#Update calls based on FDR columns
gi_data <- update_calls(gi_data,use_z = F,fdr_cutoff = 0.05)

##Making a few plots
setwd(this.dir)
setwd('../results')
Cairo::CairoPDF(file='gi_fdr_prec_vs_st_onge.pdf',width=4.5,height=4)
precision_vs_stonge(gi_data)
dev.off()

#Update the table
setwd(this.dir)
setwd('../data')
write.table(gi_data,row.names=F,quote=F,sep='\t',file='table_s1_new.tsv')
