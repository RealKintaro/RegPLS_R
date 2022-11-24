library(Hmisc)
library(MASS)
library(data.table)

PLS<-function(data, alpha, beta){
  T1_vars_indices <- function(data, alpha){
    verified <- c()
    i <- 1
    rc <- rcorr(data, type=c("pearson","spearman"))
    for(x in rc$r[,"Y"]){
      if(abs(x) > alpha && x != 1){
        verified = append(verified, i)
      }
      i <- i+1
    }
    return(verified)
  }
  # Calcule de risdu
  residu<- function(DT_scale, verified){
    rn <- as.data.frame(matrix(nrow = length(DT_scale[,"Y_cn"]), ncol = length(verified)))
    colnames(rn) <- paste(verified, "res", sep="_")
    for(x in verified){
      f <- paste(paste(x), " ~ 0 +", paste(colnames(T), collapse = " + "))
      lm_<-lm(formula = f, data=DT_scale)
      r <- resid(lm_)
      rn[,paste(x, "res", sep="_")] <- r
    }
    return(rn)
  }
  # residu normé
  residu_norm <- function(r){
    temp <- r
    for(x in colnames(temp)){
      temp[,x] <- temp[,x] / var(temp[,x]) 
    }
    colnames(temp) <- colnames(temp) <- paste(verified, "resn", sep="_")
    return(temp)
  }
  
  
  #==========================================================
  
  
  #dt <- fread("Data_Cornell.csv")
  
  names <- c()
  for (i in 1:(length(data)-1) ){
    names <- append(names, paste("X", i, sep=""))
  }
  
  f <- paste("Y ~ 0 + ", paste(names, collapse=" + "),collapse='+')
  
  names <- append(names,"Y")
  colnames(data) <- names
  
  
  lm <- lm(formula = f, data = data)
  dtm<-as.matrix(data)
  rcorr(dtm, type=c("pearson","spearman"))
  dtms<-scale(dtm)
  
  colnames(dtms)<- paste(names, "cn", sep='_')
  DT_scale <- cbind(dtm,dtms)
  cov_m <- cov(dtms)
  
  
  verified <- T1_vars_indices(data = dtm, alpha)
  #========================================================
  #calculate T1
  # first term in T1
  A <- 0
  for (x in colnames(dtms)[verified]){
    A <- A + cov_m["Y_cn",x] ^ 2
  }
  # second term in T1
  B <- 0
  for (x in colnames(dtms)[verified]){
    B <- B + dtms[,x] * cov_m["Y_cn",x]
  }
  # T1
  T1 = B / sqrt(A)
  T = as.data.frame(T1)
  DT_scale <- cbind(DT_scale, T)
  # probably loop start here
  
  k <- 2
  repeat{
    
    names<-c()
    for (i in 1:length(T)){
      names<-append(names, paste("T", i, sep=""))
    }
    names
    colnames(T) <- names
    
    DT_scale<- as.data.frame(dtm)
    DT_scale <- cbind(DT_scale, dtms)
    DT_scale <- cbind(DT_scale, T)
    
    verified <- c()
    for(i in 1:(length(colnames(dtms)) - 1) ){
      f <- paste("Y_cn ~ 0 + ", paste(paste(colnames(T), collapse=" + "), "+" ),paste(colnames(dtms)[i], collapse=" + ") )
      lm_<-lm(formula = f, data=DT_scale)
      pval <- summary(lm_)$coefficients[,4][colnames(dtms)[i]]
      if (pval < beta){
        verified <- append(verified, colnames(dtms)[i])
      }
    }
    if(length(verified) == 0){
      break
    }
    
    r <- residu(DT_scale, verified)
    rn <- residu_norm(r)
    DT_scale <- cbind(DT_scale, r)
    DT_scale <- cbind(DT_scale, rn)
    
    coef_ <- c()
    verified_ <- c()
    for(x in colnames(rn)){
      f<- paste("Y_cn ~ 0 + ", paste(paste(colnames(T), collapse=" + "), "+" ),paste(x, collapse=" + ") )
      lm_<-lm(formula = f, data=DT_scale)
      c <- substr(x,1,nchar(x)-1)
      verified_ <- append(verified_, c)
      coef_ <- append(coef_, summary(lm_)$coefficients[,1][x])
      
    }
    verified_
    coef_
    # first term in Ti
    A <- 0
    for (x in coef_){
      A <- A + (x * x)
    }
    A <- sqrt(A)
    # second term in Ti
    B <- 0
    for (i in 1:length(verified_)){
      B <- B + (coef_[i] * r[,verified_[i]]) 
    }
    assign(paste("T",k, sep=""), B / A)
    T <- cbind(T,get(paste("T",k, sep="")))
    
    k<-k+1
  }
  return(DT_scale)
}
