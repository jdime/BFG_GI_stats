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


merged_data <- merge(gi_data,st_onge_data,by=c('Barcode_x','Barcode_y'))

conditions <- gsub('GIS_xy.','',grep('^GIS',colnames(merged_data),val=T))

w_x <- merged_data[,grep('W_x\\.',colnames(merged_data),val=T)]
w_y <- merged_data[,grep('W_y\\.',colnames(merged_data),val=T)]
w_xy <- merged_data[,grep('W_xy\\.',colnames(merged_data),val=T)]
w_x_err <- merged_data[,grep('W_x_SE\\.',colnames(merged_data),val=T)]
w_y_err <- merged_data[,grep('W_y_SE\\.',colnames(merged_data),val=T)]
w_xy_err <- merged_data[,grep('W_xy_SE\\.',colnames(merged_data),val=T)]

z_x <- (w_xy - w_x)/sqrt(w_x_err^2 + w_xy_err^2)
z_y <- (w_xy - w_y)/sqrt(w_y_err^2 + w_xy_err^2)

colnames(z_x) <- sapply(conditions,function(condition){
  paste(c('A_x.',condition),collapse='')
})

colnames(z_y) <- sapply(conditions,function(condition){
  paste(c('A_y.',condition),collapse='')
})



calculate_congruence <- function(cl,merged_data,gi_name){
  parApply(cl, merged_data, 1, function(x) {
    gene1 <- x['Barcode_x']
    gene2 <- x['Barcode_y']
    
    gene1_data <-
      dplyr::filter(merged_data, Barcode_x == gene1 | Barcode_y == gene1)
    gene2_data <-
      dplyr::filter(merged_data, Barcode_x == gene2 | Barcode_y == gene2)
    
    gene1_partners <- apply(gene1_data, 1, function(x) {
      x <- x[c('Barcode_x', 'Barcode_y')]
      return(x[(x != gene1)])
    })
    gene2_partners <- apply(gene2_data, 1, function(x) {
      x <- x[c('Barcode_x', 'Barcode_y')]
      return(x[(x != gene2)])
    })
    
    common_partners <- intersect(gene1_partners, gene2_partners)
    
    gene1_data <- unlist(sapply(common_partners, function(partner) {
      filter(gene1_data, (Barcode_x == partner |
                            Barcode_y == partner))[,gi_name]
    }))
    gene2_data <- unlist(sapply(common_partners, function(partner) {
      filter(gene2_data, (Barcode_x == partner |
                            Barcode_y == partner))[,gi_name]
    }))
    
    return(cor(gene1_data, gene2_data))
  })
}


#new_congruence_data <- sapply(conditions[1:2], function(condition) {
#  gi_name <- paste(c('GIS_xy.',condition),collapse='')
#  calculate_congruence(cl,merged_data,gi_name)
#})

gi_names <- sapply(conditions,function(condition){paste(c('GIS_xy.',condition),collapse='')})
clusterExport(cl,'gi_names')
clusterExport(cl,'merged_data')
#clusterExport(cl,'calculate_congruence')
clusterEvalQ(cl, library("dplyr"))




#lul <- calculate_congruence(cl,merged_data,gi_names)
lul2 <- sapply(gi_names,function(gi_name){calculate_congruence(cl,merged_data,gi_name)})

#merged_data[,grep('FDR',colnames(merged_data))] <- -log10(merged_data[,grep('FDR',colnames(merged_data))])

z_cols <- grep('^Z_',colnames(merged_data),val=T)
z_cols <- grep('Class',z_cols,invert=T,val=T)
gi_cols <- grep('^GIS_',colnames(merged_data),val=T)
fdr_cols <- grep('^FDR',colnames(merged_data),val=T)
w_cols <- grep('^W_[xy]{1,}',colnames(merged_data),val=T)

all_data <-
  as.data.frame(cbind(merged_data$Sig_GO_link, lul2, merged_data[, c(z_cols, w_cols, gi_cols, fdr_cols)], z_x, z_y))
colnames(all_data)[1] <- 'func'

stop()  

mms_data <- all_data[,c('func',grep('MMS',colnames(all_data),val=T))]
mms_data_st_onge <- merged_data[,c('Sig_GO_link','z.Sxy.Sx.','z.Sxy.Sy.','Epsilon','Congruence')]
mms_data_st_onge[,2:ncol(mms_data_st_onge)] <- apply(mms_data_st_onge[,2:ncol(mms_data_st_onge)],2,scale)
colnames(mms_data_st_onge)[1] <- 'func'


pred <- glm(func ~., data = mms_data, family=binomial(link='logit'))
plot(performance(prediction(pred$fitted.values,all_data$func),'tpr','fpr'))

pred_st_onge <- glm(func ~., data = mms_data_st_onge, family=binomial(link='logit'))
plot(performance(prediction(pred_st_onge$fitted.values,all_data$func),'prec','rec'))

pred_full <- glm(func ~., data = all_data, family=binomial(link='logit'))
plot(performance(prediction(pred_full$fitted.values,all_data$func),'tpr','fpr'))


cross_validated_glm <- function(merged_data,all_data,conditions){
  eks <- as.matrix(all_data[,-1])
  why <- as.matrix(all_data[,1])
  bof <- all_data
  
  clusterExport(cl,'all_data')
  clusterExport(cl,'conditions')
  clusterEvalQ(cl, library("glmnet"))
  clusterEvalQ(cl, library("randomForest"))
  clusterEvalQ(cl, library("fastAdaboost"))
  
  eks <- parSapply(cl,1:nrow(merged_data),function(i){
    
    print(i)
    cut_eks <- eks[(1:nrow(all_data))[-i],]
    cut_why <- why[(1:nrow(all_data))[-i],]
    cut_bof <- bof[(1:nrow(all_data))[-i],]
    
    left_out <- eks[i,,drop=F]
    
    
    #ey <- fastAdaboost::adaboost(func~.,data=cut_bof,nIter=100)
    
    #return(predict(ey,left_out)$prob[,2])
    
    #Balance
    #crit <- which(cut_why == 1)
    #crit <- c(crit,which(cut_why == 0)[1:length(crit)])
    
    #cut_eks <- cut_eks[crit,]
    #cut_why <- cut_why[crit]
    
    initial_predictors <- lapply(conditions,function(condition){
      #cut_eks_cond <- cut_eks[,grep(condition,colnames(cut_eks))]
      
      crit <- grep(condition,colnames(cut_eks), val= T)#colnames(cut_eks) %in% c(condition)
      cut_eks_cond <- cut_eks[,crit]
      cut_bof_cond <- cut_bof[,c('func',crit)]
      
      return(fastAdaboost::adaboost(func~.,data=cut_bof_cond,nIter=10))
      
      #optim_lambda <- cv.glmnet(cut_eks_cond,cut_why,family="binomial")$lambda.min
      #return(glmnet(cut_eks_cond,cut_why,lambda=optim_lambda,family="binomial"))
      
      #return(randomForest(cut_eks_cond,as.factor(cut_why)))
    })
    
    
    initial_predictors <<- initial_predictors
    #print(initial_predictors)
    
    initial_prediction_results <- sapply(1:length(conditions),function(i){
      #print(i)
      #cut_eks_cond <- cut_eks[,grep(conditions[i],colnames(cut_eks))]
      crit <- grep(conditions[i],colnames(cut_eks))#colnames(cut_eks) %in% c(conditions[i])
      cut_eks_cond <- cut_eks[,crit]
      
      #as.numeric(as.vector(predict(initial_predictors[[i]],cut_eks_cond)))
      #predict(initial_predictors[[i]],cut_eks_cond,type='prob')[,2]
      predict(initial_predictors[[i]],cut_eks_cond)$prob[,2]
    })
    
    for_ada <- cbind(cut_why,initial_prediction_results)
    colnames(for_ada)[1] <- 'func'
    for_ada <- as.data.frame(for_ada)
    
    #meta_predictor_optim_lambda <- cv.glmnet(initial_prediction_results,cut_why,family="binomial")$lambda.min
    #meta_predictor <<- glmnet(initial_prediction_results,cut_why,lambda=meta_predictor_optim_lambda,family="binomial")
    
    meta_predictor <- fastAdaboost::adaboost(func~.,data=for_ada,nIter=10)
    
    #return(5)
    
    left_out <- eks[i,]
    #left_out <<- left_out
    
    cv_initial_prediction_results <- sapply(1:length(conditions),function(i){
      crit <- grep(conditions[i],names(left_out))#names(left_out) %in% c(conditions[i])
      cut_eks_cond <- left_out[crit]#grep(conditions[i],names(left_out))]
      
      #predict(initial_predictors[[i]],t(cut_eks_cond))
      predict(initial_predictors[[i]],t(cut_eks_cond))$prob[,2]
    })
    
    #return(5)
    return(predict(meta_predictor,t(cv_initial_prediction_results))$prob[,2])
    
    #ey <- predict(glm_pred,eks[i,,drop=F])
    
    
    #return(ey)
    
    #pred_full <- glm(func ~., data = cut_data, family=binomial(link='logit'))
    #return(predict(pred_full,all_data[i,]))
  })
  
  return(eks)
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

crit <- which(all_data$func == 1)
crit <- c(crit,which(all_data$func == 0)[1:length(crit)])

preds_cv_all <- cross_validated_glm(merged_data,all_data,conditions)#cbind(all_data[,c('func',top10)]))#all_data[,c('func',names(cols))])
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