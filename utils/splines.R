library(data.table)
library(fda)
library(splines2)

get_times_curves <- function(data, varid='patid', timename='time', varname){
  
  data<-data.table(data)
  
  obs = list()
  times = list()
  IDs = unlist(unique(data[, ..varid]),use.names=FALSE)
  for(i in 1:length(IDs)){
    obs[[i]] = unlist(data[get(varid)==IDs[i], ..varname],use.names=FALSE)
    times[[i]] = unlist(data[get(varid)==IDs[i], ..timename],use.names=FALSE)
  }
  return(list("times" = times, "obs" = obs, "IDs" = IDs))
}

fit_splines <- function(lista,m=4,nbasis=5,constrained=F,lb=0,up=0,n.points=50,gridlim=NULL){
  
  times = lista$times
  obs = lista$obs
  
  if(is.null(gridlim)){
    grid = seq(min(unlist(times)),max(unlist(times)),length.out=n.points)
  }else{
    grid = seq(gridlim[1],gridlim[2],length.out=n.points)
  }
  
  basis = create.bspline.basis(c(min(grid),max(grid)), nbasis, m)
  
  evals = list()
  deriv1 = list()
  
  
  for (i in 1:length(times)){
    print(i)
    basismat  = eval.basis(times[[i]], basis)
    
    if(constrained==T){
      coeff <-lsfit(basismat, log((obs[[i]]-lb)/(up-obs[[i]])), intercept=F)$coef
    }else{
      coeff <-lsfit(basismat, obs[[i]], intercept=F)$coef
    }
    
    if(constrained==T){
      esp <- eval.basis(grid, basis)
      esp1 <- eval.basis(grid, basis, Lfdobj=1)
      f <- lb+ up*exp(esp %*% coeff)
      g <- 1+exp(esp %*% coeff)
      evals[[i]] <- f/g
      deriv1[[i]] <- (up*exp(esp %*% coeff)*(esp1 %*% coeff)*g - f*exp(esp %*% coeff)*(esp1 %*% coeff)) / g^2
    }else{
      evals[[i]] <- t(eval.basis(grid, basis) %*% coeff)
      deriv1[[i]] <- t(eval.basis(grid, basis, Lfdobj=1) %*% coeff)
    }
  }
  
  return((list("grid" = grid, "evals" = evals, "deriv1" = deriv1, "IDs" = lista$IDs)))
}


fit_monotone_splines <- function(lista, m=5, nbasis=5, constrained=F, lb=0, up=0, lambda=10^(-1), conv=0.001, n.points=50){
  
  times = lista$times
  obs = lista$obs
  
  basis = create.bspline.basis(c(min(unlist(times)),max(unlist(times))), nbasis, m)
  
  #  starting values for coefficient
  cvec0 <- matrix(0,nbasis,1)
  Wfd0  <- fd(cvec0, basis)
  growfdPar <- fdPar(fdobj=Wfd0,lambda=lambda)
  
  grid=seq(min(unlist(times)),max(unlist(times)),length.out=n.points)
  evals=list() 
  deriv1=list()
  for (i in 1:length(times)){
    
    if(constrained==T){
      result <- smooth.monotone(times[[i]], log((obs[[i]]-lb)/(up-obs[[i]])), growfdPar,conv=conv)
    }else{
      result <- smooth.monotone(times[[i]], obs[[i]], growfdPar,conv=conv)
    }
    Wfd  <- result$Wfdobj
    beta <- result$beta
    print(i)
    if(constrained==T){
      
      esp <- eval.monfd(grid, Wfd)
      esp1 <- eval.monfd(grid, Wfd, Lfdobj=1)
      f <- (lb+ up*exp(beta[1] + beta[2]*esp))
      g <- (1+exp(beta[1] + beta[2]*esp))
      evals[[i]] <- f/g
      
      deriv1[[i]] <- ( (up*exp(beta[1] + beta[2]*esp)*beta[2]*esp1)*g - (exp(beta[1] + beta[2]*esp)*beta[2]*esp1)*f ) / g^2
      
    }else{
      evals[[i]] <- beta[1] + beta[2]*eval.monfd(grid, Wfd)
      deriv1[[i]] <- beta[2]*eval.monfd(grid, Wfd, Lfdobj=1)
    }
  }
  
  return((list("grid" = grid, "evals" = evals, "deriv1" = deriv1, "IDs" = lista$IDs)))
}