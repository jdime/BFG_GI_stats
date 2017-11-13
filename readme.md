# Genetic interaction definitions for BFG-GI

Updated statistical analysis of genetic interactions mapped using Barcode Fusion Genetics

## Getting Started

To update the genetic interaction data from the old definition, just run master.R in the scripts folder.  The only data used is table_s1.tsv and TableForAlbi_Cij_WithReplicates.txt in the data folder.  All functions this script depends on are in the packages folder.  

```
Rscript master.R
```

To run the analysis on one of two counts instead of both replicates change this line in master.R

```
Rgi_data <- cbind(gi_data,gi_data_2[,grep('^C_ij.*sum',colnames(gi_data_2))])
```

to
```
gi_data <- cbind(gi_data,gi_data_2[,grep('^C_ij.*R1',colnames(gi_data_2))])

```
or
```
gi_data <- cbind(gi_data,gi_data_2[,grep('^C_ij.*R2',colnames(gi_data_2))])
```

## Folder Contents
**scripts/** - contains the main script (master.R), as well as prototypes in the scratch folder

**data/** - contains the input data, as well as the modified input data (as table_s1_new.tsv)

**doc/** - contains an explanation of the new genetic interaction definitions and the significance calling algorithm

**packages/** - contains all scripts that master.R depends on in the BfgGiStats subfolder, in an R package format

**results/** - all resulting plots from master.R

## Authors

* **Albi Celaj**

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.