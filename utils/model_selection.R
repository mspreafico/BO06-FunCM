library(data.table)
library(survival)
library(tdROC)
library(ipred)


best_AIC_BIC_models = function(data, fpca.list1, fpca.list2, fpca.list3){
  K1=max(max(which(fpca.list1$cumPVE<0.99)+1),5)
  K2=max(max(which(fpca.list2$cumPVE<0.99)+1),5)
  K3=max(max(which(fpca.list3$cumPVE<0.99)+1),5)
  
  listAIC = c('i'=NA, 'j'=NA, 'k'=NA)
  listBIC = c('i'=NA, 'j'=NA, 'k'=NA)
  aic = bic = Inf
  pb = txtProgressBar(min = 0, max = K1*K2*K3, style = 3)
  barindex = 0
  for(i in 1:K1){
    current_data_i = merge(data, fpca.list1$PCscores[,1:(i+1)], by='patid')
    for(j in 1:K2){
      current_data_ij = merge(current_data_i, fpca.list2$PCscores[,1:(j+1)], by='patid')
      for(k in 1:K3){
        barindex = barindex+1
        setTxtProgressBar(pb, barindex)
        current_data_ijk = merge(current_data_ij, fpca.list3$PCscores[,1:(k+1)], by='patid')
        model = coxph(Surv(timeOUT_new,death) ~ . + strata(age_group) - age_group, data=current_data_ijk[,-1])
        if(AIC(model)<aic){
          listAIC = c('i'=i, 'j'=j, 'k'=k)
          aic = AIC(model)
        }
        if(BIC(model)<bic){
          listBIC = c('i'=i, 'j'=j, 'k'=k)
          bic = BIC(model)
        }
      }
    }
  }
  return(list('bestAIC' = listAIC, 'bestBIC' = listBIC))
}


cond.prob = function(model, newdata, Tstart, Tpred){
  risk.Tstart = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tstart)$surv)
  risk.Tpred = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tpred)$surv)
  return(risk.Tpred/risk.Tstart)
}


td_comparison = function(data, models, Tstart=0, u, n.folds=5, status='status', time='time'){
  
  set.seed(12345)
  data = data[sample(1:dim(data)[1]),]
  folds = cut(seq(1,dim(data)[1]),breaks=n.folds,labels=FALSE)
  
  scores = list('AUC' = matrix(NA, nrow=length(models), ncol=length(u)),
                 'Brier' = matrix(NA, nrow=length(models), ncol=length(u)),
                 'AUC_sd' = matrix(NA, nrow=length(models), ncol=length(u)),
                 'Brier_sd' = matrix(NA, nrow=length(models), ncol=length(u)),
                 'iAUC' = matrix(NA, nrow=1, ncol=length(models))
                 
  )
  pb = txtProgressBar(min = 0, max = length(models)*length(u), style = 3)
  
  for(j in 1:length(models)){
    for(t in 1:length(u)){
      setTxtProgressBar(pb, t+(j-1)*length(u))
      auc = c()
      ibs = c()
      for(i in 1:n.folds){
        validIndexes = which(folds==i,arr.ind=TRUE)
        valid = data[validIndexes,]
        train = data[-validIndexes,]
        
        fit = coxph(as.formula(models[[j]]), data=train)
        
        pred = cond.prob(fit, newdata = valid, Tstart, u[t])
        roc = tdROC( X = 1-pred, Y = valid[,get(time)], 
                     delta = valid[,get(status)], tau=u[t])
        auc = c(auc, roc$AUC2$value)
        ibs = c(ibs, sbrier(obj = Surv(valid[,get(time)], valid[,get(status)]), pred=pred, btime=u[t]))
        
      }
      scores$AUC[[j,t]] = mean(auc, na.rm=T)
      scores$Brier[[j,t]] = mean(ibs)
      scores$AUC_sd[[j,t]] = sd(auc, na.rm = T)
      scores$Brier_sd[[j,t]] = sd(ibs)
    }
  }
  
  if(!is.null(names(models))){
    rownames(scores$AUC)=names(models)
    rownames(scores$Brier)=names(models)
    rownames(scores$AUC_sd)=names(models)
    rownames(scores$Brier_sd)=names(models)
  }
  
  surv.prob = summary( survfit(Surv(timeOUT_new,death)~1, data=data), times=u)$surv
  for(i in 1:length(models)){
    scores$iAUC[1,i] = IntegrateAUC(scores$AUC[i,], u, surv.prob, tmax=u[7])
  }
  
  colnames(scores$AUC)=paste0('t=',u)
  colnames(scores$Brier)=paste0('t=',u)
  colnames(scores$AUC_sd)=paste0('t=',u)
  colnames(scores$Brier_sd)=paste0('t=',u)
  colnames(scores$iAUC)=names(models)
  
  return(scores)
  
}


best_models_list=function(models, u=NULL){
  best_models = list('formula' = NULL,
                      'AUC' = matrix(NA, nrow=8, ncol=length(u)),
                      'Brier' = matrix(NA, nrow=8, ncol=length(u)),
                      'iAUC' = matrix(NA, nrow=1, ncol=length(models))
                      
  )
  nM=1:length(models)
  for(i in nM){
    best_models$formula[[i]] = models[[i]]$formula
    best_models$AUC[i,] = models[[i]]$AUC
    best_models$Brier[i,] = models[[i]]$Brier
    best_models$iAUC[,i] = models[[i]]$iAUC
  }
  
  if(!is.null(u)){
    colnames(best_models$AUC)=paste0('t=',u)
    colnames(best_models$Brier)=paste0('t=',u)
  }
  
  rownames(best_models$AUC)=paste0('model',nM)
  rownames(best_models$Brier)=paste0('model',nM)
  colnames(best_models$iAUC)=paste0('model',nM)

  return(best_models)
}

