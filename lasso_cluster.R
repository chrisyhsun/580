# 1 day lag
lag <- 1
lagged_dataset <- make_lagged_dataset(lag)
x <- lagged_dataset[[1]]
y <- lagged_dataset[[2]]

window_length <- 90
n_windows <- nrow(x) - window_length + 1
coef_arr <- c()
# For each realized covariance, train a LASSO model
for (i in 1:31) {
  coefs <- c()
  for (j in 1:n_windows) {
    x_data <- x[j:(j+window_length-1),]
    y_data <- unlist(y[j:(j+window_length-1), ..i])
    
    lasso <- ic.glmnet(x_data, y_data, crit="aic", intercept = FALSE)
    # size of coefs: no. of covariates * no. of windows
    # Stores the lasso coefficients for every windows
    coefs <- cbind(coefs, coef(lasso))
  }
  # coef_arr stores coefs for each of the 31 dependent variables
  coef_arr[[i]] <- coefs
}
# Calculate the proportion of beta assigned to each predictor
coef_arr_norm <- c()
for (i in 1:31) {
  coefs <- coef_arr[[i]]
  # normalized_coefs is also of size no. of covariates * no. of windows
  # except that every column has been normalized to add up to 1
  normalized_coefs <- abs(coefs) %*% diag(1/colSums(abs(coefs)))
  # Take average across all windows to get vector of length no. of covariates
  # representing average proportion of beta on each covariate.
  coef_arr_norm[[i]] <- rowMeans(normalized_coefs)
}

## SCRATCH
index_names <- substr(symbols, 2, nchar(symbols))
index_countries <- c("Amsterdam", "Australia", "Belgium", "India1", "Portugal", "Brazil", "USA1", "France", "Italy", "UK1", "Germany", "Canada", "Hong Kong", "Spain", "USA2", "Korea", "Pakistan", "Mexico", "Japan", "India2", "Denmark", "Finland", "Sweden", "Norway", "UK2", "Spain2", "USA3", "China", "Switzerland", "Singapore", "EU")

dd  <-  as.data.frame(matrix(unlist(coef_arr_norm), nrow=length(unlist(coef_arr_norm[1]))))
# first row all zeros bc no intercept
dd <- dd[-1,]

#intuitive starts (India1, USA2, China, UK1)
starts <- c(4, 15, 28, 10)
start_mat <- c()
for(s in starts) {
  start_mat <- cbind(start_mat, dd[,s])
}

colnames(dd) <- index_countries
rownames(dd) <- paste(index_countries, "_rv", sep="")

tmp <- Kmeans(t(dd), centers = t(start_mat), method = "maximum")
sort(tmp$cluster)
