closeAllConnections()
rm(list=ls())

library('fasttime')
library('data.table')
library(plyr)
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
library(pracma)

setwd("D:/Spring 2020/FIN 580")

#############################################################################
add_sq_terms <- function(x, indices) {
  toRet <- x
  for(c_idx in indices) {
    toAdd <- x[,..c_idx] * x[,c_idx]
    toRet <- cbind(toRet, toAdd)
  }
  return(toRet)
}

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

## Clustering Analysis ##

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
        
        lasso <- ic.glmnet(x_data, y_data, crit="bic", intercept = FALSE, maxit = 1e+07)
        first.step.coef <- coef(lasso)[-1]
        pf <- (abs(first.step.coef)+1/sqrt(abs(dim(x_data)[1])))^(-tau)
        adalasso <- ic.glmnet(x_data,y_data,crit="bic",penalty.factor=pf, intercept = FALSE)
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
               "MX", "SG", "CH", "ES", "GB", "SE", "DK", 
               "JP", "PK", "KR", "HK", "DE", "AU", "BE", 
               "NL", "PT", "BR", "FR", "FI")
# US, Italy, Canada, India, China, 
# Mexico, Singapore, Switzerland, Spain, United Kingdom, Stockholm, Copenhagen, 
# Japan, Pakistan, South Korea, Hong Kong, Germany, Australia, Belgium, 
# Netherlands, Portugal, Brazil, France, Finland

group1_idx <- c(3, 4, 5, 8, 9, 20, 21, 22, 23, 24, 26)
group3_idx <- c(18, 28)
group2_idx <- c(1, 2, 6, 7, 11, 12, 15, 25, 27, 29, 31)
group4_idx <- c(10, 13, 14, 16, 17, 19, 30)
lags_to_try <- c(1, 2, 3, 4, 5)

words <- c("coronavirus", "COVID-19")
start <- as.Date(all_dates[1])
end <- as.Date(Sys.Date())

countries1 <- c("BE", "IN", "PT", "FR", "IT", "FI", "DK", "ES")
gtrends1 <- get_word_data(words, countries1, start, end)
gtrends1 <- gtrends1[which(rownames(gtrends1) %in% all_dates),]
gtrends1[is.na(gtrends1)] <- 0

countries2 <- c("AU", "BR", "US", "DE", "CA", "GB", "CH")
gtrends2 <- get_word_data(words, countries2, start, end)
gtrends2 <- gtrends2[which(rownames(gtrends2) %in% all_dates),]
gtrends2[is.na(gtrends2)] <- 0

countries3 <- c("MX", "CN")
gtrends3 <- get_word_data(words, countries3, start, end)
gtrends3 <- gtrends3[which(rownames(gtrends3) %in% all_dates),]
gtrends3[is.na(gtrends3)] <- 0

countries4 <- c("GB", "HK", "ES", "KR", "PK", "JP", "SG")
gtrends4 <- get_word_data(words, countries4, start, end)
gtrends4 <- gtrends4[which(rownames(gtrends4) %in% all_dates),]
gtrends4[is.na(gtrends4)] <- 0


labels <- c('Lag 1', 'Lag 2', 'Lag 3', 'Lag 4', 'Lag 5',
            'Squared', 'Cross', 'Google trends')

grp1_res_lag <- get_mses(group1_idx, lags_to_try)
grp1_cluster_mses_lag <- rowMeans(grp1_res_lag)
grp1_res_sq <- get_mses(group1_idx, sq = TRUE)
grp1_cluster_mses_sq <- mean(grp1_res_sq)
grp1_res_cross <- get_mses(group1_idx, cross = TRUE, cross_itself = TRUE)
grp1_cluster_mses_cross <- mean(grp1_res_cross)
grp1_res_gtrends <- get_mses(group1_idx, gtrends = gtrends1)
grp1_cluster_mses_gtrends <- mean(grp1_res_gtrends) 

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

grp4_res_lag <- get_mses(group4_idx, lags_to_try)
grp4_cluster_mses_lag <- rowMeans(grp4_res_lag)
grp4_res_sq <- get_mses(group4_idx, sq = TRUE)
grp4_cluster_mses_sq <- mean(grp4_res_sq)
grp4_res_cross <- get_mses(group4_idx, cross = TRUE, cross_itself = TRUE)
grp4_cluster_mses_cross <- mean(grp4_res_cross)
grp4_res_gtrends <- get_mses(group4_idx, gtrends = gtrends4)
grp4_cluster_mses_gtrends <- mean(grp4_res_gtrends) 

grp1_all_mses <- c(grp1_cluster_mses_lag, grp1_cluster_mses_sq, grp1_cluster_mses_cross, grp1_cluster_mses_gtrends)
grp2_all_mses <- c(grp2_cluster_mses_lag, grp2_cluster_mses_sq, grp2_cluster_mses_cross, grp2_cluster_mses_gtrends)
grp3_all_mses <- c(grp3_cluster_mses_lag, grp3_cluster_mses_sq, grp3_cluster_mses_cross, grp3_cluster_mses_gtrends)
grp4_all_mses <- c(grp4_cluster_mses_lag, grp4_cluster_mses_sq, grp4_cluster_mses_cross, grp4_cluster_mses_gtrends)
res <- cbind(grp1_all_mses, grp2_all_mses, grp3_all_mses, grp4_all_mses)

rownames(res) <- labels
res

#############################################################################

## Modified HAR (Time and Log-scale) ##

get_mses_har <- function(group_idx, har_days = c(5, 22)) {
  lag <- 1
  lagged_dataset <- make_lagged_dataset(lag)
  x <- lagged_dataset[[1]]
  idx_names <- colnames(x)
  for(har_day in har_days) {
    toAdd <- as.data.frame(rollmean(zoo(x[, 1:31]), har_day))
    colnames(toAdd) <- paste0(idx_names, "_rm", har_day)
    pad <- as.data.frame(matrix(NA, nrow = (har_day - 1), ncol = ncol(toAdd)))
    colnames(pad) <- paste0(idx_names, "_rm", har_day)
    toAdd <- rbind(pad, toAdd)
    x <- cbind(x, toAdd) 
  }
  x <- x[(max(har_days):nrow(x)),]
  x <- data.matrix(x)
  y <- lagged_dataset[[2]]
  y <- tail(y, nrow(x))
  
  window_length <- 90
  n_windows <- nrow(x) - window_length
  coef_arr <- c()
  mse_grp <- 0

  tmp_cnt <- 1
  
  for (i in group_idx) {
    coefs <- c()
    err <- c()

    for (j in 1:n_windows) {
      x_data <- x[j:(j+window_length),]
      rel_cols <- seq(0, length(har_days)) * 31 + i
      x_data <- data.frame(x_data) %>% select(rel_cols)
      y_data <- unlist(y[j:(j+window_length-1), ..i])
      
      train_x <- x_data[1:(nrow(x_data) - 1),]
      test_x <- x_data[nrow(x_data),]
      
      linear_mod <- lm(y_data ~ ., data = data.frame(train_x))
      rv_pred <- sum(linear_mod$coefficients * c(1, as.numeric(test_x)), na.rm = TRUE)

      rv_actual_val <-unlist(y[(j+window_length),..i])
      
      # check error only when the actual value is present
      if (is.na(rv_actual_val) == FALSE){
        err <- c(err, (rv_pred - rv_actual_val)^2)
      }
      
      # size of coefs: no. of covariates * no. of windows
      # Stores the lasso coefficients for every windows
      coefs <- cbind(coefs, linear_mod$coefficients)
    }
    # coef_arr stores coefs for each of the 31 dependent variables
    coef_arr[[tmp_cnt]] <- coefs
    tmp_cnt <- tmp_cnt + 1
    mse_grp <- cbind(mse_grp, err)
  }
  return(mse_grp)[, 2:(ncol(mse_grp))]
}

get_mses_har_ln <- function(group_idx, har_days = c(5, 22)) {
  lag <- 1
  lagged_dataset <- make_lagged_dataset(lag)
  x <- lagged_dataset[[1]]
  x <- data.frame(x)
  x[x == 0] <- 1e-15
  x <- log(x)
  idx_names <- colnames(x)
  for(har_day in har_days) {
    toAdd <- as.data.frame(rollmean(zoo(x[, 1:31]), har_day))
    colnames(toAdd) <- paste0(idx_names, "_rm", har_day)
    pad <- as.data.frame(matrix(NA, nrow = (har_day - 1), ncol = ncol(toAdd)))
    colnames(pad) <- paste0(idx_names, "_rm", har_day)
    toAdd <- rbind(pad, toAdd)
    x <- cbind(x, toAdd) 
  }
  x <- x[(max(har_days):nrow(x)),]
  y <- lagged_dataset[[2]]
  y <- tail(y, nrow(x))
  
  window_length <- 90
  n_windows <- nrow(x) - window_length
  coef_arr <- c()
  mse_grp <- 0

  tmp_cnt <- 1
  
  for (i in group_idx) {
    coefs <- c()
    err <- c()

    for (j in 1:n_windows) {
      x_data <- x[j:(j+window_length),]
      rel_cols <- seq(0, length(har_days)) * 31 + i
      x_data <- data.frame(x_data) %>% select(rel_cols)
      y_data <- unlist(y[j:(j+window_length-1), ..i])
      y_data[y_data == 0] <- 1e-15
      y_data <- log(y_data)
      
      train_x <- x_data[1:(nrow(x_data) - 1),]
      test_x <- x_data[nrow(x_data),]
      
      linear_mod <- lm(y_data ~ ., data = data.frame(train_x))
      rv_pred <- exp(sum(linear_mod$coefficients * c(1, as.numeric(test_x)), na.rm = TRUE))

      rv_actual_val <-unlist(y[(j+window_length),..i])
      
      # check error only when the actual value is present
      if (is.na(rv_actual_val) == FALSE){
        err <- c(err, (rv_pred - rv_actual_val)^2)
      }
      
      # size of coefs: no. of covariates * no. of windows
      # Stores the lasso coefficients for every windows
      coefs <- cbind(coefs, linear_mod$coefficients)
    }
    # coef_arr stores coefs for each of the 31 dependent variables
    coef_arr[[tmp_cnt]] <- coefs
    tmp_cnt <- tmp_cnt + 1
    mse_grp <- cbind(mse_grp, err)
  }
  return(mse_grp)[, 2:(ncol(mse_grp))]
}

#############################################################################

split_mses <- function(mses, max_har_days, start_idx = 2) {
  split_idx <- which(all_dates == '2020-02-01') - 90 - max_har_days
  post_feb <- mses[split_idx:(nrow(mses)),]
  toRet <- cbind(colMeans(post_feb), colMeans(mses))
  return(toRet[start_idx:nrow(toRet),])
}

#############################################################################

h_days_1 <- c(5, 22)
grp1_res_har_5_22 <- get_mses_har(group_idx = group1_idx, har_days = h_days_1)
grp1_res_har_ln_5_22 <- get_mses_har_ln(group_idx = group1_idx, har_days = h_days_1)
grp1_5_22 <- cbind(split_mses(grp1_res_har_5_22, 22), split_mses(grp1_res_har_ln_5_22, 22))
grp2_res_har_5_22 <- get_mses_har(group_idx = group2_idx, har_days = h_days_1)
grp2_res_har_ln_5_22 <- get_mses_har_ln(group_idx = group2_idx, har_days = h_days_1)
grp2_5_22 <- cbind(split_mses(grp2_res_har_5_22, 22), split_mses(grp2_res_har_ln_5_22, 22))
grp3_res_har_5_22 <- get_mses_har(group_idx = group3_idx, har_days = h_days_1)
grp3_res_har_ln_5_22 <- get_mses_har_ln(group_idx = group3_idx, har_days = h_days_1)
grp3_5_22 <- cbind(split_mses(grp3_res_har_5_22, 22), split_mses(grp3_res_har_ln_5_22, 22))
grp4_res_har_5_22 <- get_mses_har(group_idx = group4_idx, har_days = h_days_1)
grp4_res_har_ln_5_22 <- get_mses_har_ln(group_idx = group4_idx, har_days = h_days_1)
grp4_5_22 <- cbind(split_mses(grp4_res_har_5_22, 22), split_mses(grp4_res_har_ln_5_22, 22))

h_days_2 <- c(3, 5)
grp1_res_har_3_5 <- get_mses_har(group_idx = group1_idx, har_days = h_days_2)
grp1_res_har_ln_3_5 <- get_mses_har_ln(group_idx = group1_idx, har_days = h_days_2)
grp1_3_5 <- cbind(split_mses(grp1_res_har_3_5, 5), split_mses(grp1_res_har_ln_3_5, 5))
grp2_res_har_3_5 <- get_mses_har(group_idx = group2_idx, har_days = h_days_2)
grp2_res_har_ln_3_5 <- get_mses_har_ln(group_idx = group2_idx, har_days = h_days_2)
grp2_3_5 <- cbind(split_mses(grp2_res_har_3_5, 5), split_mses(grp2_res_har_ln_3_5, 5))
grp3_res_har_3_5 <- get_mses_har(group_idx = group3_idx, har_days = h_days_2)
grp3_res_har_ln_3_5 <- get_mses_har_ln(group_idx = group3_idx, har_days = h_days_2)
grp3_3_5 <- cbind(split_mses(grp3_res_har_3_5, 5), split_mses(grp3_res_har_ln_3_5, 5))
grp4_res_har_3_5 <- get_mses_har(group_idx = group4_idx, har_days = h_days_2)
grp4_res_har_ln_3_5 <- get_mses_har_ln(group_idx = group4_idx, har_days = h_days_2)
grp4_3_5 <- cbind(split_mses(grp4_res_har_3_5, 5), split_mses(grp4_res_har_ln_3_5, 5))

h_days_3 <- c(5, 10, 22)
grp1_res_har_5_10_22 <- get_mses_har(group_idx = group1_idx, har_days = h_days_3)
grp1_res_har_ln_5_10_22 <- get_mses_har_ln(group_idx = group1_idx, har_days = h_days_3)
grp1_5_10_22 <- cbind(split_mses(grp1_res_har_5_10_22, 22), split_mses(grp1_res_har_ln_5_10_22, 22))
grp2_res_har_5_10_22 <- get_mses_har(group_idx = group2_idx, har_days = h_days_3)
grp2_res_har_ln_5_10_22 <- get_mses_har_ln(group_idx = group2_idx, har_days = h_days_3)
grp2_5_10_22 <- cbind(split_mses(grp2_res_har_5_10_22, 22), split_mses(grp2_res_har_ln_5_10_22, 22))
grp3_res_har_5_10_22 <- get_mses_har(group_idx = group3_idx, har_days = h_days_3)
grp3_res_har_ln_5_10_22 <- get_mses_har_ln(group_idx = group3_idx, har_days = h_days_3)
grp3_5_10_22 <- cbind(split_mses(grp3_res_har_5_10_22, 22), split_mses(grp3_res_har_ln_5_10_22, 22))
grp4_res_har_ln_5_10_22 <- get_mses_har_ln(group_idx = group4_idx, har_days = h_days_3)
grp4_res_har_5_10_22 <- get_mses_har(group_idx = group4_idx, har_days = h_days_3)
grp4_5_10_22 <- cbind(split_mses(grp4_res_har_5_10_22, 22), split_mses(grp4_res_har_ln_5_10_22, 22))

h_days_4 <- c(3, 5, 10, 22)
grp1_res_har_3_5_10_22 <- get_mses_har(group_idx = group1_idx, har_days = h_days_4)
grp1_res_har_ln_3_5_10_22 <- get_mses_har_ln(group_idx = group1_idx, har_days = h_days_4)
grp1_3_5_10_22 <- cbind(split_mses(grp1_res_har_3_5_10_22, 22), split_mses(grp1_res_har_ln_3_5_10_22, 22))
grp2_res_har_3_5_10_22 <- get_mses_har(group_idx = group2_idx, har_days = h_days_4)
grp2_res_har_ln_3_5_10_22 <- get_mses_har_ln(group_idx = group2_idx, har_days = h_days_4)
grp2_3_5_10_22 <- cbind(split_mses(grp2_res_har_3_5_10_22, 22), split_mses(grp2_res_har_ln_3_5_10_22, 22))
grp3_res_har_3_5_10_22 <- get_mses_har(group_idx = group3_idx, har_days = h_days_4)
grp3_res_har_ln_3_5_10_22 <- get_mses_har_ln(group_idx = group3_idx, har_days = h_days_4)
grp3_3_5_10_22 <- cbind(split_mses(grp3_res_har_3_5_10_22, 22), split_mses(grp3_res_har_ln_3_5_10_22, 22))
grp4_res_har_3_5_10_22 <- get_mses_har(group_idx = group4_idx, har_days = h_days_4)
grp4_res_har_ln_3_5_10_22 <- get_mses_har_ln(group_idx = group4_idx, har_days = h_days_4)
grp4_3_5_10_22 <- cbind(split_mses(grp4_res_har_3_5_10_22, 22), split_mses(grp4_res_har_ln_3_5_10_22, 22))

allgrp1_har <- cbind(grp1_5_22, grp1_3_5, grp1_5_10_22, grp1_3_5_10_22)
write.csv(data.frame(allgrp1_har), "grp1.csv", row.names = FALSE)
allgrp2_har <- cbind(grp2_5_22, grp2_3_5, grp2_5_10_22, grp2_3_5_10_22)
write.csv(data.frame(allgrp2_har), "grp2.csv", row.names = FALSE)
allgrp3_har <- cbind(grp3_5_22, grp3_3_5, grp3_5_10_22, grp3_3_5_10_22)
write.csv(data.frame(allgrp3_har), "grp3.csv", row.names = FALSE)
allgrp4_har <- cbind(grp4_5_22, grp4_3_5, grp4_5_10_22, grp4_3_5_10_22)
write.csv(data.frame(allgrp4_har), "grp4.csv", row.names = FALSE)

#############################################################################

## Modified HAR (Residuals) ##

get_mses_har_ln_resid <- function(group_idx, har_days = c(3, 5)) {
  lag <- 1
  lagged_dataset <- make_lagged_dataset(lag)
  x <- lagged_dataset[[1]]
  x <- data.frame(x)
  x[x == 0] <- 1e-15
  x <- log(x)
  idx_names <- colnames(x)
  for(har_day in har_days) {
    toAdd <- as.data.frame(rollmean(zoo(x[, 1:31]), har_day))
    colnames(toAdd) <- paste0(idx_names, "_rm", har_day)
    pad <- as.data.frame(matrix(NA, nrow = (har_day - 1), ncol = ncol(toAdd)))
    colnames(pad) <- paste0(idx_names, "_rm", har_day)
    toAdd <- rbind(pad, toAdd)
    x <- cbind(x, toAdd) 
  }
  x <- x[(max(har_days):nrow(x)),]
  y <- lagged_dataset[[2]]
  y <- tail(y, nrow(x))
  
  window_length <- 90
  n_windows <- nrow(x) - window_length
  coef_arr <- c()
  mse_grp <- 0

  tmp_cnt <- 1
  
  for (i in group_idx) {
    coefs <- c()
    err <- c()
    err_rw <- c()
    
    for (j in 1:n_windows) {
      x_data <- x[j:(j+window_length),]
      rel_cols <- seq(0, length(har_days)) * 31 + i
      x_data1 <- data.frame(x_data) %>% select(rel_cols)
      x_data2 <- select(data.frame(x_data), -rel_cols)
      y_data <- unlist(y[j:(j+window_length-1), ..i])
      y_data1 <- y_data
      y_data1[y_data1 == 0] <- 1e-15
      y_data1 <- log(y_data1)
      
      train_x1 <- x_data1[1:(nrow(x_data1) - 1),]
      test_x1 <- x_data1[nrow(x_data1),]
      linear_mod <- lm(y_data1 ~ ., data = data.frame(train_x1))
      resid <- y_data - exp(linear_mod$fitted.values)
      train_x2 <- x_data2[1:(nrow(x_data2) - 1),]
      test_x2 <- x_data2[nrow(x_data2),]
      linear_mod_resid <- lm(resid ~ ., data = data.frame(train_x2))
      
      rv_pred <- exp(sum(linear_mod$coefficients * c(1, as.numeric(test_x1)), na.rm = TRUE))
                     + sum(linear_mod_resid$coefficients * c(1, as.numeric(test_x2)), na.rm = TRUE)
      
      rv_actual_val <-unlist(y[(j+window_length),..i])
      
      # check error only when the actual value is present
      if (is.na(rv_actual_val) == FALSE){
        err <- c(err, (rv_pred - rv_actual_val)^2)
      }
      
      # size of coefs: no. of covariates * no. of windows
      # Stores the lasso coefficients for every windows
      coefs <- cbind(coefs, linear_mod$coefficients)
    }
    # coef_arr stores coefs for each of the 31 dependent variables
    coef_arr[[tmp_cnt]] <- coefs
    tmp_cnt <- tmp_cnt + 1
    mse_grp <- cbind(mse_grp, err)
  }
  return(mse_grp)[, 2:(ncol(mse_grp))]
}

get_mses_har_ln_resid2 <- function(group_idx, har_days = c(3, 5)) {
  lag <- 1
  lagged_dataset <- make_lagged_dataset(lag)
  x <- lagged_dataset[[1]]
  x <- data.frame(x)
  
  x[x == 0] <- 1e-15
  x <- log(x)
  idx_names <- colnames(x)
  for(har_day in har_days) {
    toAdd <- as.data.frame(rollmean(zoo(x[, 1:31]), har_day))
    colnames(toAdd) <- paste0(idx_names, "_rm", har_day)
    pad <- as.data.frame(matrix(NA, nrow = (har_day - 1), ncol = ncol(toAdd)))
    colnames(pad) <- paste0(idx_names, "_rm", har_day)
    toAdd <- rbind(pad, toAdd)
    x <- cbind(x, toAdd) 
  }
  x <- x[(max(har_days):nrow(x)),]

  y <- lagged_dataset[[2]]
  y <- tail(y, nrow(x))
  
  window_length <- 90
  n_windows <- nrow(x) - window_length
  coef_arr <- c()
  mse_grp <- 0

  tmp_cnt <- 1
  
  for (i in group_idx) {
    coefs <- c()
    err <- c()

    for (j in 1:n_windows) {
      x_data <- x[j:(j+window_length),]
      rel_cols <- seq(0, length(har_days)) * 31 + i
      x_data1 <- data.frame(x_data) %>% select(rel_cols)
      y_data <- unlist(y[j:(j+window_length-1), ..i])
      y_data1 <- y_data
      y_data1[y_data1 == 0] <- 1e-15
      y_data1 <- log(y_data1)
      
      
      #print(head(x_data))
      #toRet <- list()
      #toRet[[1]] <- x_data
      #toRet[[2]] <- y_data
      #return(toRet)
      
      train_x1 <- x_data1[1:(nrow(x_data1) - 1),]
      test_x1 <- x_data1[nrow(x_data1),]
      linear_mod <- lm(y_data1 ~ ., data = data.frame(train_x1))
      
      resid <- y_data - exp(linear_mod$fitted.values)
      train_x2 <- train_x1^2
      test_x2 <- test_x1^2
      linear_mod_resid <- lm(resid ~ ., data = data.frame(train_x2))
      
      rv_pred <- exp(sum(linear_mod$coefficients * c(1, as.numeric(test_x1)), na.rm = TRUE))
      + sum(linear_mod_resid$coefficients * c(1, as.numeric(test_x2)), na.rm = TRUE)

      rv_actual_val <-unlist(y[(j+window_length),..i])
      
      # check error only when the actual value is present
      if (is.na(rv_actual_val) == FALSE){
        err <- c(err, (rv_pred - rv_actual_val)^2)
      }
      
      # size of coefs: no. of covariates * no. of windows
      # Stores the lasso coefficients for every windows
      coefs <- cbind(coefs, linear_mod$coefficients)
    }
    # coef_arr stores coefs for each of the 31 dependent variables
    coef_arr[[tmp_cnt]] <- coefs
    tmp_cnt <- tmp_cnt + 1
    mse_grp <- cbind(mse_grp, err)
  }
  return(mse_grp)[, 2:(ncol(mse_grp))]
}

#############################################################################

grp1_res_har_ln_resid <- get_mses_har_ln_resid(group1_idx)
grp1_res_har_ln_resid2 <- get_mses_har_ln_resid2(group1_idx)
grp1_resid <- cbind(split_mses(grp1_res_har_ln_resid, 5), split_mses(grp1_res_har_ln_resid2, 5))
grp2_res_har_ln_resid <- get_mses_har_ln_resid(group2_idx)
grp2_res_har_ln_resid2 <- get_mses_har_ln_resid2(group2_idx)
grp2_resid <- cbind(split_mses(grp2_res_har_ln_resid, 5), split_mses(grp2_res_har_ln_resid2, 5))
grp3_res_har_ln_resid <- get_mses_har_ln_resid(group3_idx)
grp3_res_har_ln_resid2 <- get_mses_har_ln_resid2(group3_idx)
grp3_resid <- cbind(split_mses(grp3_res_har_ln_resid, 5), split_mses(grp3_res_har_ln_resid2, 5))
grp4_res_har_ln_resid <- get_mses_har_ln_resid(group4_idx)
grp4_res_har_ln_resid2 <- get_mses_har_ln_resid2(group4_idx)
grp4_resid <- cbind(split_mses(grp4_res_har_ln_resid, 5), split_mses(grp4_res_har_ln_resid2, 5))

write.csv(data.frame(grp1_resid), "grp1_resid.csv", row.names = FALSE)
write.csv(data.frame(grp2_resid), "grp2_resid.csv", row.names = FALSE)
write.csv(data.frame(grp3_resid), "grp3_resid.csv", row.names = FALSE)
write.csv(data.frame(grp4_resid), "grp4_resid.csv", row.names = FALSE)

#############################################################################

get_mses_rw <- function(group_idx) {
  lag <- 1
  lagged_dataset <- make_lagged_dataset(lag)
  x <- lagged_dataset[[1]]
  y <- lagged_dataset[[2]]
  
  window_length <- 90
  n_windows <- nrow(x) - window_length
  toRet <- c()
  
  for (i in group_idx) {
    rv_pred <- y[window_length:(nrow(y)-1), ..i]
    rv_actual <- y[(window_length + 1):(nrow(y)), ..i]
    toRet <- cbind(toRet, (rv_actual - rv_pred)^2)
  }
  return(toRet)
}

## Bar Plot vs. RW ##

mses_rw <- get_mses_rw(seq(1, 31))
mses_mod_har <- get_mses_har_ln(group_idx = seq(1, 31), har_days = h_days_2)
toPlot <- cbind(split_mses(mses_mod_har, 5), split_mses(mses_rw, 0, start_idx = 1))
toPlot <- toPlot[, c(1, 3)]
colnames(toPlot) <- c('Har(1, 3, 5) with log scale', 'Random Walk')
rownames(toPlot) <- colnames(realized_covariances)
prev_mar <- par()$mar
prev_mgp <- par()$mgp
par(mar = c(6.3, 6.3, 4.1, 2.1))
par(mgp = c(5, 1, 0))
barplot(t(toPlot), beside=TRUE, ylab="MSE", 
        cex.names=0.8, cex.lab = .8, las=2, col=c("black","lightgrey"))
box(bty="l")
par(mar = prev_mar)
par(mgp = prev_mgp)

## Table Values ##
grp1_rw <- mean(split_mses(get_mses_rw(group1_idx), 1, start_idx=1)[,1])
grp2_rw <- mean(split_mses(get_mses_rw(group2_idx), 1, start_idx=1)[,1])
grp3_rw <- mean(split_mses(get_mses_rw(group3_idx), 1, start_idx=1)[,1])
grp4_rw <- mean(split_mses(get_mses_rw(group4_idx), 1, start_idx=1)[,1])
grp1_postfeb <- colMeans(allgrp1_har[,seq(1, ncol(allgrp1_har), 2)])
grp2_postfeb <- colMeans(allgrp2_har[,seq(1, ncol(allgrp2_har), 2)])
grp3_postfeb <- colMeans(allgrp3_har[,seq(1, ncol(allgrp3_har), 2)])
grp4_postfeb <- colMeans(allgrp4_har[,seq(1, ncol(allgrp4_har), 2)])
grp1_postfeb_resid <- colMeans(grp1_resid[,seq(1, ncol(grp1_resid), 2)])
grp2_postfeb_resid <- colMeans(grp2_resid[,seq(1, ncol(grp2_resid), 2)])
grp3_postfeb_resid <- colMeans(grp3_resid[,seq(1, ncol(grp3_resid), 2)])
grp4_postfeb_resid <- colMeans(grp4_resid[,seq(1, ncol(grp4_resid), 2)])
tab_vals <- cbind(c(grp1_rw, grp1_postfeb, grp1_postfeb_resid),
                  c(grp2_rw, grp2_postfeb, grp2_postfeb_resid),
                  c(grp3_rw, grp3_postfeb, grp3_postfeb_resid),
                  c(grp4_rw, grp4_postfeb, grp4_postfeb_resid))
colnames(tab_vals) <- c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4')
rownames(tab_vals) <- c('Random Walk', 
                        'HAR(1, 5, 22)', 'HAR(1, 5, 22), log', 
                        'HAR(1, 3, 5)', 'HAR(1, 3, 5), log', 
                        'HAR(1, 5, 10, 22)', 'HAR(1, 5, 10, 22), log',
                        'HAR(1, 3, 5, 10, 22)', 'HAR(1, 3, 5, 10, 22), log',
                        'HAR(1, 3, 5), log, lm on others',
                        'HAR(1, 3, 5), log, lm on squares')

#############################################################################

get_preds_har <- function(group_idx, har_days = c(5, 22)) {
  lag <- 1
  lagged_dataset <- make_lagged_dataset(lag)
  x <- lagged_dataset[[1]]
  idx_names <- colnames(x)
  for(har_day in har_days) {
    toAdd <- as.data.frame(rollmean(zoo(x[, 1:31]), har_day))
    colnames(toAdd) <- paste0(idx_names, "_rm", har_day)
    pad <- as.data.frame(matrix(NA, nrow = (har_day - 1), ncol = ncol(toAdd)))
    colnames(pad) <- paste0(idx_names, "_rm", har_day)
    toAdd <- rbind(pad, toAdd)
    x <- cbind(x, toAdd) 
  }
  x <- x[(max(har_days):nrow(x)),]
  x <- data.matrix(x)
  y <- lagged_dataset[[2]]
  y <- tail(y, nrow(x))
  
  window_length <- 90
  n_windows <- nrow(x) - window_length
  pred_grp <- c()
  
  tmp_cnt <- 1
  
  for (i in group_idx) {
    preds <- c()
    
    for (j in 1:n_windows) {
      x_data <- x[j:(j+window_length),]
      rel_cols <- seq(0, length(har_days)) * 31 + i
      x_data <- data.frame(x_data) %>% select(rel_cols)
      y_data <- unlist(y[j:(j+window_length-1), ..i])
      
      train_x <- x_data[1:(nrow(x_data) - 1),]
      test_x <- x_data[nrow(x_data),]
      
      linear_mod <- lm(y_data ~ ., data = data.frame(train_x))
      rv_pred <- sum(linear_mod$coefficients * c(1, as.numeric(test_x)), na.rm = TRUE)
      preds <- c(preds, rv_pred)
    }
    pred_grp <- cbind(pred_grp, preds)
  }
  return(pred_grp)[, 2:(ncol(pred_grp))]
}

get_preds_har_ln <- function(group_idx, har_days = c(5, 22)) {
  lag <- 1
  lagged_dataset <- make_lagged_dataset(lag)
  x <- lagged_dataset[[1]]
  x <- data.frame(x)
  x[x == 0] <- 1e-15
  x <- log(x)
  idx_names <- colnames(x)
  for(har_day in har_days) {
    toAdd <- as.data.frame(rollmean(zoo(x[, 1:31]), har_day))
    colnames(toAdd) <- paste0(idx_names, "_rm", har_day)
    pad <- as.data.frame(matrix(NA, nrow = (har_day - 1), ncol = ncol(toAdd)))
    colnames(pad) <- paste0(idx_names, "_rm", har_day)
    toAdd <- rbind(pad, toAdd)
    x <- cbind(x, toAdd) 
  }
  x <- x[(max(har_days):nrow(x)),]
  y <- lagged_dataset[[2]]
  y <- tail(y, nrow(x))
  
  window_length <- 90
  n_windows <- nrow(x) - window_length
  pred_grp <- c()
  
  for (i in group_idx) {
    preds <- c()
    
    for (j in 1:n_windows) {
      x_data <- x[j:(j+window_length),]
      rel_cols <- seq(0, length(har_days)) * 31 + i
      x_data <- data.frame(x_data) %>% select(rel_cols)
      y_data <- unlist(y[j:(j+window_length-1), ..i])
      y_data[y_data == 0] <- 1e-15
      y_data <- log(y_data)
      
      train_x <- x_data[1:(nrow(x_data) - 1),]
      test_x <- x_data[nrow(x_data),]
      
      linear_mod <- lm(y_data ~ ., data = data.frame(train_x))
      rv_pred <- exp(sum(linear_mod$coefficients * c(1, as.numeric(test_x)), na.rm = TRUE))
      preds <- c(preds, rv_pred)
    }
    pred_grp <- cbind(pred_grp, preds)
  }
  return(pred_grp)[, 2:(ncol(pred_grp))]
}

grp1_preds <- get_preds_har(group1_idx, har_days = c(5, 22))
grp2_preds <- get_preds_har(group2_idx, har_days = c(5, 22))
grp3_preds <- get_preds_har(group3_idx, har_days = c(5, 22))
grp4_preds <- get_preds_har(group4_idx, har_days = c(5, 22))
all_preds <- get_preds_har(seq(1, 31), har_days = c(5, 22))
toWrite <- data.frame(all_preds)
colnames(toWrite) <- colnames(realized_covariances)
write.csv(toWrite, "har_preds.csv", row.names = FALSE)

grp1_preds_ln <- get_preds_har_ln(group1_idx, har_days = c(3, 5))
grp2_preds_ln <- get_preds_har_ln(group2_idx, har_days = c(3, 5))
grp3_preds_ln <- get_preds_har_ln(group3_idx, har_days = c(3, 5))
grp4_preds_ln <- get_preds_har_ln(group4_idx, har_days = c(3, 5))
all_preds_ln <- get_preds_har_ln(seq(1, 31), har_days = c(3, 5))
toWrite_ln <- data.frame(all_preds_ln)
colnames(toWrite_ln) <- colnames(realized_covariances)
write.csv(toWrite_ln, "modhar_preds.csv", row.names = FALSE)

#############################################################################
