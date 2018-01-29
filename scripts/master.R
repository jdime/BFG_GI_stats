##################################################################
##################################################################
##################################################################
########                                                  ########
########                                                  ########
######## Pipeline to calculate Genetic Interaction Scores ########
######## from Barcode Fusion Genetics experiments         ########
######## and stastical evaluation of condition-dependent  ########
######## genetic interactions                             ########         
########                                                  ########
########                                                  ########
########                                                  ########
##################################################################
##################################################################
##################################################################
require(qvalue)
require(dplyr)

### #Needs to be run from source, not Rscript or this line won't work
### this_dir <- function() {
###   cmdArgs <- commandArgs(trailingOnly = FALSE)
###   needle <- "--file="
###   match <- grep(needle, cmdArgs)
###   if (length(match) > 0) {
###     # Rscript
###     return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
###   } else {
###     # 'source'd via R console
###     return(dirname(sys.frame(1)$ofile))
###   }
### }
### this.dir <- this_dir()

print(this.dir)

### Get directory where 'master.R' is contained
this.dir <- getwd() ### this must be the path to ~/.../BFG_GI_stats-master/scripts directory (where 'master.R' is contained)
setwd(this.dir)

### Create output and results directories
resultsDir<-(sub("scripts", "results",this.dir))
outputDir<-(sub("scripts", "output",this.dir))

if (dir.exists(resultsDir) == T){
## drectory already exists
} else {
    dir.create(resultsDir,recursive = T)

}

if (dir.exists(outputDir) == T){
## drectory already exists
} else {
    dir.create(outputDir,recursive = T)

}




devtools::load_all('../packages/BfgGiStats')
devtools::document('../packages/BfgGiStats')


setwd('../data')

#Read in genetic interaction data from two tables, using sum column from second table
gi_data <-
  read.table('table_s1_input.tsv',
             head = T,
             stringsAsFactors = F)
#Remove count data from first table
gi_data <- gi_data[, -grep('^C_', colnames(gi_data))]


#Update removal criteria to 100kb
gi_data[gi_data$Chromosomal_distance_bp <= 75000 &
          !is.na(gi_data$Chromosomal_distance_bp), 'Remove_by_Chromosomal_distance_or_SameGene'] <-
  'YES'

#Can change this to look at only one of the two technical replicates, analyzes the sum by default
#Change  to 'R1' or 'R2' to look at one of the two only
to_analyze <- 'R1'

grep_pattern1 <- sprintf('^C_xy.*%s', to_analyze)
grep_pattern2 <- sprintf('.%s', to_analyze)

gi_data_2 <-
  read.table('Cxy_Table_WithReplicates.txt',
             head = T,
             stringsAsFactors = F)
gi_data <-
  cbind(gi_data, gi_data_2[, grep(grep_pattern1, colnames(gi_data_2))])
#Fix the names
colnames(gi_data) <- gsub(grep_pattern2, '', colnames(gi_data))




#Remove non well-measured strains
well_measured <- gi_data[, grep('HetDipl', colnames(gi_data))] >= 100
gi_data <- gi_data[well_measured,]


#For troubleshooting/checking
#gi_data_original <- gi_data

#Change genetic interactions and significance calls
gi_data <- update_gis(gi_data,
                         #These values are empirically determined and have to be
                         #in the same order as in gi_data count columns
                         g_wt_vec = c('NoDrug' = 12.62,
                                      'DMSO' = 8.34,
                                      'MMS' = 7.84,
                                      '4NQO' = 7.5,
                                      'BLMC' = 6.94,
                                      'ZEOC' = 6.28,
                                      'HYDX' = 7.76,
                                      'DXRB' = 7.04,
                                      'CMPT' = 7.7,
                                      'CSPL' = 8.44))


#Analyze linkage patterns to justify above removal criteria
setwd(this.dir)
setwd('../results')
pdf(file = 'GIS_NoDrug_vs_distance.pdf',
                width = 4,
                height = 4)
par(mar=c(4.5,4.5,1,1))
plot(log10(gi_data$Chromosomal_distance_bp + 1),
     gi_data$GIS_xy.NoDrug,
     xlab = expression(Log[10](Chromosomal~distance~+1)),
     ylab = expression(GIS[xy]),
     pch = 16,
     col=rgb(0,0,0,0.3))
abline(v= log10(75001),col='red',lty=3,lwd=2)
dev.off()

#Add p values columns
setwd(this.dir)
setwd('../results')
pdf(file = 'Z_distribution.pdf',
                width = 12,
                height = 5)
par(mfrow = c(2, 5))
par(mar = c(5, 4, 2, 0))
gi_data <-
  add_p_values(
    gi_data,
    z_grep_pattern = "^Z_GIS",
    nn_pair_type = 'broad',
    make_plots = T
  )
dev.off()

#Remove camptothecin
gi_data <- gi_data[, -grep('CMPT', colnames(gi_data))]



#Preserve at this step for differential gi_calls
gi_data_old <- gi_data

## Reviewer plot, have to temporarily allow non well-measured
## strains for it to work
# setwd(this.dir)
# setwd('../results')
# pdf(file = 'gis_well_measured_vs_non.pdf',
#                 width = 4,
#                 height = 3.5)
# par(mar=c(4,4,3,1))
# plot(density(as.matrix(
#   gi_data_old %>% dplyr::filter(C_xy.HetDipl >= 100, Chromosomal_distance_bp == 0) %>% dplyr::select(grep('^GIS_xy.', colnames(gi_data_old)))
# )), xlim = c(-1.1, 0.5),lwd=2,col='red',xlab='GIS',main='GIS of Same-Same Pairs')
# lines(density(as.matrix(
#   gi_data_old %>% dplyr::filter(C_xy.HetDipl < 100, Chromosomal_distance_bp == 0) %>% dplyr::select(grep('^GIS_xy.', colnames(gi_data_old)))
# )),
# col='blue',
# lwd='2'
# )
# legend(-0.1,3,
#        legend=c(expression(C[xy]>=100),
#                 expression(C[xy]<100)),
#        col=c('red','blue'),
#        pch=15)
# dev.off()


#Make AUC plot
setwd(this.dir)
setwd('../results')

pdf(file = 'auc_vs_st_onge_new.pdf',
                width = 7.5,
                height = 4)
par(mfrow = c(1, 2))
st_onge_auc_plot(gi_data)
dev.off()

#Score comparison scatterplot - by barcode
pdf(file = 'st_onge_scatterplot_barcodewise.pdf',
                width = 7.5,
                height = 3.5)
par(mfrow = c(1, 2))
st_onge_scatterplot(gi_data)
dev.off()


########
#Calculate a per-gene GI score, Z score, and p value
########
gi_data <- average_gi_data_by_gene(gi_data)

#Make sure p-value distribution looks sane
pdf(file = 'gene_averaged_p_value_histogram.pdf',
                width = 5.5,
                height = 4.5)
hist(as.matrix(gi_data[, grep('^P.neutral', colnames(gi_data))]),
     main = '',
     xlab = 'p-value',
     col = 'grey30')
dev.off()



#Calculate FDR from per-gene p values
gi_data[, grep('^P.neutral', colnames(gi_data))] <-
  apply(gi_data[, grep('^P.neutral', colnames(gi_data))], 2, function(x) {
    qvalue(x)$q
  })

fdr_colnames <-
  sapply(grep('P.neut', colnames(gi_data), val = T), function(x) {
    gsub('^P.', 'FDR.', x)
  })

colnames(gi_data)[grep('^P.neutral', colnames(gi_data))] <- fdr_colnames


setwd(this.dir)
setwd('../results')
pdf(file = 'gi_prec_vs_st_onge.pdf',
                width = 4.5,
                height = 4)
precision_vs_stonge(
  gi_data,
  fdr_cutoff = 0.05,
  xlims = c(-4, 4),
  cutoffs_drawn = c(2,2)
)
dev.off()


#Update calls based on 1% FDR cutoffs
#No effect sizes used, but option exists
gi_data <- update_calls(
  gi_data,
  fdr_cutoff_pos = 0.01,
  gi_cutoff_pos = 0,
  fdr_cutoff_neg = 0.01,
  gi_cutoff_neg = 0
)



#Make gene-wise scatterplot
#Manuscript uses barcode-wise plot
pdf(file = 'st_onge_scatterplot_genewise.pdf',
                width = 7.5,
                height = 3.5)
par(mfrow = c(1, 2))
st_onge_scatterplot(gi_data)
dev.off()


###Write some data
setwd(this.dir)
setwd('../data')
dir.create('output', showWarnings = FALSE)
setwd('output')

write.table(
  gi_data_old,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s1.tsv'
)
write.table(
  gi_data,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s1_genewise.tsv'
)

##Make/write differential calls
##Make differential calls, reporting everything
par(mfrow = c(7, 7))
par(mar = c(3, 2, 4, 1))
differential_calls <- differential_gi_analysis(
  gi_data,
  #Report everything
  fdr_cutoff = 1.1,
  require_sign_change = F
  #make_plots = T
)

setwd(this.dir)
setwd('../data/output')
write.table(
  differential_calls,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s2_all.tsv'
)


#Filter for significance and sign change
differential_calls_sig <- differential_gi_analysis(
  gi_data,
  fdr_cutoff = 0.01,
  delta_gi_cutoff = 0,
  require_sign_change = T
)

#Can compare what happens if no sign change enforced
#differential_calls_sig_no_sign <- differential_gi_analysis(
#  gi_data,
#  fdr_cutoff = 0.01,
#  delta_gi_cutoff = 0,
#  require_sign_change = F
#)



#Write some data
setwd(this.dir)
setwd('../data/output')
write.table(
  differential_calls_sig,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s2_significant.tsv'
)


#Looks at frequency of differetial calls by gene
setwd(this.dir)
setwd('../results')
ddr_data <- dplyr::filter(gi_data, Type_of_gene_x != "Neutral", Type_of_gene_y != "Neutral")
genes <- unique(c(ddr_data$Gene_x,ddr_data$Gene_y))
differential_count <-
  sapply(genes, function(gene) {
    return(nrow(
      dplyr::filter(differential_calls_sig, Gene_x == gene |
               Gene_y == gene)
    ))
  })
pdf(file = 'differential_inters_by_gene.pdf',
                width = 4,
                height = 4)
hist(
  differential_count,
  breaks = 20,
  col = rgb(0, 0, 0, 0.7),
  xlab = 'Number of Differential Pairs',
  ylab = 'Frequency (number of genes)',
  main = ''
)
dev.off() 


#Look at only sign reversals
reversal_count <-
  sapply(genes, function(gene) {
    return(nrow(
      dplyr::filter(differential_calls_sig, Gene_x == gene |
               Gene_y == gene, Class_Condition1 != "Expected", Class_Condition2 != "Expected")
    ))
  })

pdf(file = 'sign_reversals_by_gene.pdf',
                width = 4,
                height = 4)
hist(
  reversal_count,
  breaks = 20,
  col = rgb(0, 0, 0, 0.7),
  xlab = 'Number of Sign-Reversed Pairs',
  ylab = 'Frequency (number of genes)',
  main = ''
)
dev.off()
