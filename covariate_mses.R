closeAllConnections()
rm(list=ls())

library('fasttime')
library('data.table')
library('dplyr')
library('matrixStats')
library('bit64')
library('ggplot2')
library('gdata')
library('naniar')
library('mice')
library('forecast')
library('imputeTS')
library('zoo')
library('png')
library('grid')
library('HDeconometrics')
library('tictoc')
library('amap')
library('gtools')
library(gtrendsR)
library(reshape2)


setwd("D:/Spring 2020/FIN 580")

#############################################################################
add_cross_terms <- function(x, indices, itself = FALSE) {
  toRet <- x
  combos <- combinations(length(indices), 2, indices)
  for(r_idx in seq(1, dim(combos)[1])) {
    c_idx1 <- combos[r_idx,][1]
    c_idx2 <- combos[r_idx,][2]
    toAdd <- x[,..c_idx1] * x[,..c_idx2]
    toRet <- cbind(toRet, toAdd)
  }
  if(itself) {
    for(c_idx in indices) {
      toAdd <- x[,..c_idx] * x[,..c_idx]
      toRet <- cbind(toRet, toAdd)
    }
  }
  return(toRet)
}

make_lagged_dataset <- function(lag) {
  n_rows <- nrow(realized_covariances) - lag
  
  x <- c()
  y <- realized_covariances[(lag+1):nrow(realized_covariances), ]
  
  for (i in 1:lag) {
    x <- cbind(realized_covariances[(lag+1-i):(nrow(realized_covariances)-i), ], x)   
  }
  
  return (list(x, y)) 
}

get_word_data <- function(words, countries, start, end) {
  all_word_data <- c()
  start_date <- start
  while(start_date < end) {
    end_date <- min(start_date + 90, end)
    time = paste(as.character(start_date), as.character(end_date))
    #print(time) #debugging
    
    toReplace <- 0.5 # numeric chosen to represent "<1" in data
    
    datalist = list()
    i <- 1
    for(w_idx in seq(length(words))) {
      for(c_idx in seq(length(countries))) {
        google.trends <- gtrends(words[w_idx], geo = countries[c_idx], 
                                 gprop = "web", time = paste(as.character(start_date), as.character(end_date)))$interest_over_time
        if(!is.null(google.trends)) {
          #print(colnames(google.trends)) #debugging
          google.trends <- dcast(google.trends, date ~ keyword + geo, value.var = "hits")
          rownames(google.trends) <- google.trends$date
          google.trends$date <- NULL
          datalist[[i]] <- google.trends
          i <- i + 1
        }
      }
    }
    
    big_data <- do.call(cbind, datalist)
    big_data[big_data == "<1"] <- toReplace
    big_data$date <- seq(as.Date(start_date), by = "day", length.out = nrow(big_data))
    
    all_word_data <- rbind.fill(all_word_data, big_data)
    start_date <- start_date + 91
  }
  rownames(all_word_data) <- all_word_data$date
  all_word_data$date <- NULL
  return(all_word_data)
}

#############################

realized_covariances_orig <- fread('data/final_data.csv')
all_dates <- realized_covariances_orig$Date
realized_covariances_extra <- realized_covariances_orig[, -c("V1", "Date", "covars")]
realized_covariances <- realized_covariances_extra[, c(seq(3, 93, 3))]

tau=1

#############################################################################

get_mses <- function(group_idx, lags_to_try = c(1), sq = FALSE, cross = FALSE, cross_itself = FALSE, 
                     gtrends = NULL, debug = FALSE) {
  mse_mat <- matrix(0, nrow = length(lags_to_try), ncol = length(group_idx))
  mse_rw_mat <- matrix(0, nrow = length(lags_to_try), ncol = length(group_idx))
  
  for(l_idx in 1:length(lags_to_try)) {
    lag <- lags_to_try[l_idx]
    lagged_dataset <- make_lagged_dataset(lag)
    x <- lagged_dataset[[1]]
    if(sq) { x <- cbind(x, x*x)}
    if(cross) { x <- add_cross_terms(x, group_idx, itself = cross_itself) }
    if(!is.null(gtrends)) {
      gtrends <- gtrends[1:(dim(gtrends)[1] - lag),]
      x <- cbind(x, gtrends)
    }
    x <- data.matrix(x)
    y <- lagged_dataset[[2]]
    
    window_length <- 90
    n_windows <- nrow(x) - window_length
    coef_arr <- c()
    mse_grp <- 0
    mse_rw <- 0
    
    tmp_cnt <- 1
    
    for (i in group_idx) {
      coefs <- c()
      err <- c()
      err_rw <- c()
  
      for (j in 1:n_windows) {
        x_data <- x[j:(j+window_length-1),]
        y_data <- unlist(y[j:(j+window_length-1), ..i])
        
        lasso <- ic.glmnet(x_data, y_data, crit="aic", intercept = FALSE, maxit = 1e+07)
        first.step.coef <- coef(lasso)[-1]
        pf <- (abs(first.step.coef)+1/sqrt(abs(dim(x_data)[1])))^(-tau)
        adalasso <- ic.glmnet(x_data,y_data,crit="aic",penalty.factor=pf, intercept = FALSE)
        rv_pred <- predict(adalasso, x[(j+window_length),])[1][1]
        rv_actual_val <-unlist(y[(j+window_length),..i])
        # check error only when the actual value is present
        if (is.na(rv_actual_val) == FALSE){
          err <- c(err, (rv_pred - rv_actual_val)^2)
          err_rw <- c(err_rw, (unlist(y[(j+window_length - 1),..i]) - rv_actual_val)^2)
        }
        
        # size of coefs: no. of covariates * no. of windows
        # Stores the lasso coefficients for every windows
        coefs <- cbind(coefs, coef(adalasso))
      }
      # coef_arr stores coefs for each of the 31 dependent variables
      coef_arr[[tmp_cnt]] <- coefs
      tmp_cnt <- tmp_cnt + 1
      mse_grp <- c(mse_grp, mean(err))
      mse_rw <- c(mse_rw, mean(err_rw))
    }
    
    mse_mat[l_idx,] <- mse_grp[2:length(mse_grp)]
    mse_rw_mat[l_idx,] <- mse_rw[2:length(mse_rw)]
  }
  if(debug) {
    toRet <- list()
    toRet[[1]] <- mse_mat
    toRet[[2]] <- mse_rw_mat
    return(toRet)
  }
  return(mse_mat)
}

countries <- c("US", "IT", "CA", "IN", "CN", 
               "MX", "SG", "CH", "ES", "GB", "SE", "DK", "JP", 
               "PK", "KR", "HK", "DE", "AU", "BE", 
               "NL", "PT", "BR", "FR", "FI")
# US, Italy, Canada, India, China, 
# Mexico, Singapore, Switzerland, Spain, United Kingdom, Stockholm, Copenhagen, 
# Japan, Pakistan, South Korea, Hong Kong, Germany, Australia, Belgium, 
# Netherlands, Portugal, Brazil, France, Finland

group3_idx <- c(18, 28)
group2_idx <- c(1, 2, 6, 7, 11, 12, 15, 25, 27, 29, 31)
lags_to_try <- c(1, 2, 3, 4, 5)

words <- c("coronavirus", "COVID-19")
start <- as.Date(all_dates[1])
end <- as.Date(Sys.Date())

countries2 <- c("AU", "BR", "US", "DE", "CA", "GB", "CH")
gtrends2 <- get_word_data(words, countries2, start, end)
gtrends2 <- gtrends2[which(rownames(gtrends2) %in% all_dates),]
gtrends2[is.na(gtrends2)] <- 0

countries3 <- c("MX", "CN")
gtrends3 <- get_word_data(words, countries3, start, end)
gtrends3 <- gtrends3[which(rownames(gtrends3) %in% all_dates),]
gtrends3[is.na(gtrends3)] <- 0


grp2_res_lag <- get_mses(group2_idx, lags_to_try)
grp2_cluster_mses_lag <- rowMeans(grp2_res_lag)
grp2_res_sq <- get_mses(group2_idx, sq = TRUE)
grp2_cluster_mses_sq <- mean(grp2_res_sq)
grp2_res_cross <- get_mses(group2_idx, cross = TRUE, cross_itself = TRUE)
grp2_cluster_mses_cross <- mean(grp2_res_cross)
grp2_res_gtrends <- get_mses(group2_idx, gtrends = gtrends2)
grp2_cluster_mses_gtrends <- mean(grp2_res_gtrends) 


grp3_res_lag <- get_mses(group3_idx, lags_to_try = lags_to_try)
grp3_cluster_mses_lag <- rowMeans(grp3_res_lag)
grp3_res_sq <- get_mses(group3_idx, sq = TRUE)
grp3_cluster_mses_sq <- mean(grp3_res_sq)
grp3_res_cross <- get_mses(group3_idx, cross = TRUE, cross_itself = TRUE)
grp3_cluster_mses_cross <- mean(grp3_res_cross) 
grp3_res_gtrends <- get_mses(group3_idx, gtrends = gtrends3)
grp3_cluster_mses_gtrends <- mean(grp3_res_gtrends) 