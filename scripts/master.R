this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

devtools::load_all('../packages/BfgGiStats')
devtools::document('../packages/BfgGiStats')

setwd('../data')
gi_data <-
  read.table('table_s1.tsv', head = T, stringsAsFactors = F)

z_cols <-
  grep('Class',
       grep('^Z_GIS', colnames(gi_data), val = T),
       invert = T,
       val = T)
z_cols <- which(colnames(gi_data) %in% z_cols)

gi_data_old <- gi_data

#stop()

gi_data <- update_gis(gi_data,g_wt_vec = c(12.62,8.34,8.44,7.04,7.7,7.84,7.5,7.76,6.94,6.28))
gi_data <- add_fdrs(gi_data)
gi_data <- update_calls(gi_data,use_z = T,
                        z_cutoff_neg = -5,
                        z_cutoff_pos = 6)

#gi_data[,grep('^GIS',colnames(gi_data))] <- 2^gi_data[,grep('^GIS',colnames(gi_data))]

st_onge_comparison_plot(gi_data_old,gi_data)
par(mfrow=c(1,1))
prec_recall_vs_stonge(gi_data,xlims=c(-0.2,0.2))
prec_recall_vs_stonge(
  gi_data,
  #score_cols = z_cols,
  control_name = "Z_GIS_ij.NoDrug",
  condition_name = "Z_GIS_ij.MMS",
  xlims=c(-20,10),
  fdr_cutoff = 0
)