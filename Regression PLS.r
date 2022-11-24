# Regression PLS

# read data function
read_data <- function(path2data, file_name) {
  data <- read.csv(file.path(path2data, file_name), header = TRUE, sep = ";")
  return(data)
}

# linear regression function
linear_regression <- function(data, formula) {
  model <- lm(formula = formula, data = data)
  return(model)
}

# correlation matrix function
correlation_matrix <- function(data, type) {
    data_matrix <- as.matrix(data)
  corr <- cor(data, method = type)
  return(corr)
}

# scale data function
scale_data <- function(data) {
  data_scaled <- scale(data)
  return(data_scaled)
}

# PLS algorithm function
pls_algorithm <- function(data, ncomp) {
  pls <- plsr(data, ncomp = ncomp)
  return(pls)
}

