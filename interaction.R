library('gtools')

# x: x from make_lagged_dataset
# indices: indices for which you want all combinations of cross terms
# itself: whether to multiply a vector by itself; default is to not include these
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