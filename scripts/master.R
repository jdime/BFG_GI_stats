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

gi_data <- update_gis(gi_data)
gi_data <- add_fdrs(gi_data)
gi_data <- update_calls(gi_data)

st_onge_comparison_plot(gi_data_old,gi_data)
par(mfrow=c(1,1))
prec_recall_vs_stonge(gi_data)
prec_recall_vs_stonge(
  gi_data,
  score_cols = z_cols,
  control_name = "Z_GIS_ij.DMSO",
  condition_name = "Z_GIS_ij.MMS"
)