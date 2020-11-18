library(data.table)
library(dplyr)

fpca_prcomp <- function(fd_list, name=NULL, deriv=F){
  
  # Reformat fdata list in a matrix format
  if(deriv==T){
    evals_matrix = fd_list$deriv1
  }else{
    evals_matrix = fd_list$evals
  }
  patient_ids = unique(fd_list$IDs)
  evals_matrix = bind_cols(evals_matrix)
  colnames(evals_matrix) = as.character(patient_ids)
  
  # Get the mean curves
  mean_evals = rowMeans(evals_matrix)
  
  # Perform FPCA using standard prcomp function as explained in 
  # Ramsay and Silverman (2005) Functional Data Analysis, Chaper 8, Section 8.4.1
  nc <- dim(evals_matrix)[2] # number of columns: number of patients
  m <- dim(evals_matrix)[1] # number of rows: time
  h <- (fd_list$grid[m] - fd_list$grid[1])/(m - 1)
  
  # FPCA: Eigenvalues and eigenfunctions
  pca <- prcomp(t(evals_matrix))
  efuncs <- pca$rotation*sqrt(1/h)
  evalues <- pca$sdev^2*h
  
  # Eigenfunctions
  # Sign chosen to have positive integral (Optional)
  for(i in 1:m) {
    tempfun <- efuncs[,i]
    tempsign <- sum(tempfun)
    efuncs[,i] <- ifelse(tempsign<0, -1,1) * tempfun
  }
  
  # Principal Components Scores
  scores <- t(evals_matrix-mean_evals) %*% efuncs*h
  if(!is.null(name)){
    colnames(scores) <- paste0(colnames(scores), '_', name)
  }
  scores <- cbind("patid" = patient_ids, scores)
  
  # Create list
  fpca_list <- list("grid" = fd_list$grid,
                    "mean_evals" = mean_evals,
                    "efuncs" = efuncs,
                    "evalues" = evalues,
                    "PVE" = evalues/sum(evalues),
                    "cumPVE" = cumsum(evalues/sum(evalues)),
                    "PCscores" = scores)
  
  return(fpca_list)
  
}


fpca_plot <- function(fpca_list, K, xlab='Time', ylim2, byrow=T, cost=NULL){
  
  if(is.null(cost)){
    cost = rep(1,K)
  }
  if(length(cost)!=K){
    stop("Vector cost must be of lenght K")
  }
  time = fpca_list$grid
  x11()
  if(byrow==T){
    par(mfrow=c(K,2), mar=c(2.5,2,2.5,1))
  }else{
    par(mfcol=c(2,K), mar=c(2.5,2,2.5,1))
  }
  ylim <- range(fpca_list$efuncs[,1:K])
  for(i in 1:K) {
    v1 <- fpca_list$efuncs[,i]
    plot(time, v1, type="l", ylab="", ylim=ylim, lwd=2, 
         main=paste("Component ", i,", ",100*round(fpca_list$PVE[i],3), "%", sep=""))
    v2 <- v1 * cost[i] * sqrt(fpca_list$evalues[i])
    plot(time, fpca_list$mean_evals, type="l", ylim=ylim2, lwd=3)
    points(time, fpca_list$mean_evals + v2, lwd=3, pch="+", col='red')
    points(time, fpca_list$mean_evals - v2, lwd=3, pch="-", col='blue')
    
  }
}

merge_list <- function(data, fpca_list, varID='patid'){
  K <- max(which(fpca_list$cumPVE<0.99)+1)
  data <- merge(data, fpca_list$PCscores[,1:(K+1)], by=varID)
  return(data)
}
