library(Hmisc)
library(MASS)
library(data.table)

PLS<-function(data,Y,pvalue){

  # pour extraire les indices des variable pour la construction du T1
  T1_vars_indices <- function(data, pvalue){
    verified <- c()
    i <- 1
    # extraction des correlations et les P values
    rc <- rcorr(data, type=c("pearson","spearman"))
    # pour chaque p value on verifie si elle est inferieur a la valeur de P value entrée
    for(x in rc$P[1:nrow(rc$P) - 1,"Y"]){
      if(abs(x) < pvalue && x != 1){
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

  # Lecture des données
  data = as.data.frame(data)
  # extraction de la variable Y
  y = data[,paste(Y)]
  data[Y] = NULL

  # Renommer les variables
  names <- c()
  for (i in 1:(length(data)) ){
    names <- append(names, paste("X", i, sep=""))
  }
  colnames(data) <- names
  data[,'Y'] <- as.numeric(y)

  f <- paste("Y ~ 0 + ", paste(names, collapse=" + "),collapse='+')

  # Premiere regression
  lm <- lm(formula = f, data = data)

  # Extraction des corrélations et les pvalues
  dtm<-as.matrix(data)
  rcorr(dtm, type=c("pearson","spearman"))

  # Normalisation des données
  dtms<-scale(dtm)
  
  colnames(dtms)<- paste(colnames(data), "cn", sep='_')
  DT_scale <- cbind(dtm,dtms)
  cov_m <- cov(dtms)
  
  
  verified <- T1_vars_indices(data = dtm, pvalue)
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
  
  # Construction des composantes T
  #========================================================
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
    
    # Extraction des variables pour la construction de Tk
    verified <- c()
    for(i in 1:(length(colnames(dtms)) - 1) ){
      f <- paste("Y_cn ~ 0 + ", paste(paste(colnames(T), collapse=" + "), "+" ),paste(colnames(dtms)[i], collapse=" + ") )
      lm_<-lm(formula = f, data=DT_scale)
      pval <- summary(lm_)$coefficients[,4][colnames(dtms)[i]]
      if(! is.na(pval)){
        if (pval < pvalue){
          verified <- append(verified, colnames(dtms)[i])
        }
      }
    }

    #========================================================
    # Condition d'arret
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

# get of this file path

path2data<-file.path("c:","Users","makch","STUDIES","s3","ADD","regPLS")

# read data
DM.dt <- fread(file.path(path2data, "Data_Cornell.csv"))
mtcars <- fread("mtcars.csv")

Tcor = PLS(DM.dt,'Y',0.2)
Tmtcar = PLS(mtcars[,-c(1)],'mpg',0.2)
