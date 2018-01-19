library(parallel)
library(reshape2)

st_onge_data <- read.table('st_onge_function_link_table.tsv',sep='\t',head=T)
st_onge_data[,1:2] <- t(apply(st_onge_data[,1:2],1,function(x){
  sort(x)
}))


num_cores <- detectCores()
cl <- makeCluster(num_cores)

gi_data[,1:2] <- t(apply(gi_data[,1:2],1,function(x){
  sort(x)
}))

gi_data <- gi_data[,grep('CMPT',colnames(gi_data),invert=T)]


merged_data <- merge(gi_data,st_onge_data,by=c('Barcode_i','Barcode_j'))

conditions <- gsub('GIS_ij.','',grep('^GIS',colnames(merged_data),val=T))

w_i <- merged_data[,grep('W_i\\.',colnames(merged_data),val=T)]
w_j <- merged_data[,grep('W_j\\.',colnames(merged_data),val=T)]
w_ij <- merged_data[,grep('W_ij\\.',colnames(merged_data),val=T)]
w_i_err <- merged_data[,grep('W_i_SE\\.',colnames(merged_data),val=T)]
w_j_err <- merged_data[,grep('W_j_SE\\.',colnames(merged_data),val=T)]
w_ij_err <- merged_data[,grep('W_ij_SE\\.',colnames(merged_data),val=T)]

z_i <- (w_ij - w_i)/sqrt(w_i_err^2 + w_ij_err^2)
z_j <- (w_ij - w_j)/sqrt(w_j_err^2 + w_ij_err^2)

colnames(z_i) <- sapply(conditions,function(condition){
  paste(c('A_i.',condition),collapse='')
})

colnames(z_j) <- sapply(conditions,function(condition){
  paste(c('A_j.',condition),collapse='')
})



calculate_congruence <- function(cl,merged_data,gi_name){
  parApply(cl, merged_data, 1, function(x) {
    gene1 <- x['Barcode_i']
    gene2 <- x['Barcode_j']
    
    gene1_data <-
      dplyr::filter(merged_data, Barcode_i == gene1 | Barcode_j == gene1)
    gene2_data <-
      dplyr::filter(merged_data, Barcode_i == gene2 | Barcode_j == gene2)
    
    gene1_partners <- apply(gene1_data, 1, function(x) {
      x <- x[c('Barcode_i', 'Barcode_j')]
      return(x[(x != gene1)])
    })
    gene2_partners <- apply(gene2_data, 1, function(x) {
      x <- x[c('Barcode_i', 'Barcode_j')]
      return(x[(x != gene2)])
    })
    
    common_partners <- intersect(gene1_partners, gene2_partners)
    
    gene1_data <- unlist(sapply(common_partners, function(partner) {
      filter(gene1_data, (Barcode_i == partner |
                            Barcode_j == partner))[,gi_name]
    }))
    gene2_data <- unlist(sapply(common_partners, function(partner) {
      filter(gene2_data, (Barcode_i == partner |
                            Barcode_j == partner))[,gi_name]
    }))
    
    return(cor(gene1_data, gene2_data))
  })
}


#new_congruence_data <- sapply(conditions[1:2], function(condition) {
#  gi_name <- paste(c('GIS_ij.',condition),collapse='')
#  calculate_congruence(cl,merged_data,gi_name)
#})

gi_names <- sapply(conditions,function(condition){paste(c('GIS_ij.',condition),collapse='')})
clusterExport(cl,'gi_names')
clusterExport(cl,'merged_data')
clusterEvalQ(cl, library("dplyr"))

#lul <- calculate_congruence(cl,merged_data,gi_names)
#lul2 <- sapply(gi_names,function(gi_name){calculate_congruence(cl,merged_data,gi_name)})

z_cols <- grep('^Z_',colnames(merged_data),val=T)
gi_cols <- grep('^GIS_',colnames(merged_data),val=T)

all_data <- as.data.frame(cbind(merged_data$Sig_GO_link,lul,lul2,merged_data[,c(gi_cols)],z_i,z_j))
colnames(all_data)[1] <- 'func'

mms_data <- all_data[,c('func',grep('MMS',colnames(all_data),val=T))]
mms_data_st_onge <- merged_data[,c('Sig_GO_link','z.Sxy.Sx.','z.Sxy.Sy.','Epsilon','Congruence')]
colnames(mms_data_st_onge)[1] <- 'func'

  

pred <- glm(func ~., data = mms_data, family=binomial(link='logit'))
plot(performance(prediction(pred$fitted.values,all_data$func),'tpr','fpr'))

pred_st_onge <- glm(func ~., data = mms_data_st_onge, family=binomial(link='logit'))
plot(performance(prediction(pred_st_onge$fitted.values,all_data$func),'tpr','fpr'))

pred_full <- glm(func ~., data = all_data, family=binomial(link='logit'))
plot(performance(prediction(pred_full$fitted.values,all_data$func),'tpr','fpr'))


cross_validated_glm <- function(merged_data,all_data){
  eks <- as.matrix(all_data[,-1])
  why <- as.matrix(all_data[,1])
  eks <- apply(eks,2,scale)
  
  clusterExport(cl,'all_data')
  clusterEvalQ(cl, library("glmnet"))
  
  parSapply(cl,1:nrow(merged_data),function(i){
    print(i)
    #gene1 <- merged_data$Barcode_i[i]
    #gene2 <- merged_data$Barcode_i[i]
    
    #query <- !(merged_data$Barcode_i %in% c(gene1,gene2) | merged_data$Barcode_j %in% c(gene1,gene2))
    
    #cut_data <- all_data[(1:nrow(all_data))[-i],]
    
    
    cut_eks <- eks[(1:nrow(all_data))[-i],]
    cut_why <- why[(1:nrow(all_data))[-i],]
    
    glm_pred <- glmnet(cut_eks,cut_why,lambda=cv.glmnet(cut_eks,cut_why,family="binomial")$lambda.min,family="binomial")
    
    #left_out <- all_data[i,-1]
    ey <- predict(glm_pred,eks[i,,drop=F])
    return(ey)
    
    #pred_full <- glm(func ~., data = cut_data, family=binomial(link='logit'))
    #return(predict(pred_full,all_data[i,]))
  })
}


#eks <- as.matrix(all_data[,-1])
#eks <- eks[,names(which(apply(all_data[,2:ncol(all_data)],2,function(i){cor.test(i,all_data$func)$p.v}) < 0.05))]
#why <- as.matrix(all_data[,1])
#glm_pred <- glmnet(eks,why,lambda=cv.glmnet(eks,why,family="binomial")$lambda.min,family="binomial")
#ey <- predict(glm_pred,eks)

#stop()

#for(i in 2:20){
#top10 <- names(which(apply(all_data[,2:ncol(all_data)],2,function(i){cor.test(i,all_data$func)$p.v}) < 0.05/ncol(all_data)))
top10 <- names(sort(apply(all_data[,2:ncol(all_data)],2,function(i){cor.test(i,all_data$func)$p.v}))[1:10])
  #names(sort(apply(all_data[,2:ncol(all_data)],2,function(i){cor.test(i,all_data$func)$p.v}))[1:10])


#  preds_cv_all <- cross_validated_glm(merged_data,all_data[,c('func',top5)])
#  print(max(performance(prediction(preds_cv_all,all_data$func),'mat')@y.values[[1]],na.rm=T))
#}

preds_cv_all <- cross_validated_glm(merged_data,all_data)#cbind(all_data[,c('func',top10)]))#all_data[,c('func',names(cols))])
preds_cv_mms <- cross_validated_glm(merged_data,mms_data)
preds_cv_mms_st_onge <- cross_validated_glm(merged_data,mms_data_st_onge)
preds_cv_mms_both <- cross_validated_glm(merged_data,cbind(mms_data,mms_data_st_onge[,2:ncol(mms_data_st_onge)]))

p1 <- performance(prediction(preds_cv_all,all_data$func),'prec','rec')
p2 <- performance(prediction(preds_cv_mms,all_data$func),'prec','rec')
p3 <- performance(prediction(preds_cv_mms_st_onge,all_data$func),'prec','rec')
p4 <- performance(prediction(preds_cv_mms_both,all_data$func),'prec','rec')

#distro <- sapply(1:nrow(all_data),function(i){
#  sampvec <- sample(1:nrow(all_data),replace=T)
#  perf <- performance(prediction(preds_cv_mms[sampvec],all_data$func[sampvec]),'mat')
#  max(perf@y.values[[1]],na.rm=T)
#})


#stopCluster(cl)