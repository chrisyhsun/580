## install libraries first ##
library(gtrendsR)
library(reshape2)

## Args ##
words <- c("coronavirus", "COVID-19")
countries <- c("US", "IT", "CN")
start_date <- "2020-01-01"
end_date <- "2020-04-11"
toReplace <- 0.5 # numeric chosen to represent "<1" in data

datalist = list()
i <- 1
for(w_idx in seq(length(words))) {
  for(c_idx in seq(length(countries))) {
    google.trends <- gtrends(words[w_idx], geo = countries[c_idx], 
                            gprop = "web", time = paste(start_date, end_date))[[1]]
    google.trends <- dcast(google.trends, date ~ keyword + geo, value.var = "hits")
    rownames(google.trends) <- google.trends$date
    google.trends$date <- NULL
    datalist[[i]] <- google.trends
    i <- i + 1
  }
}

big_data <- do.call(cbind, datalist)
big_data[big_data == "<1"] <- toReplace