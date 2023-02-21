#'@description calculate pooled standard deviation for equal sample sizes
#'@references https://www.statisticshowto.com/pooled-standard-deviation/
#'@author Rodrigo Souto-Veiga

fun_pooled_sd_equal_sample_sizes <- function(x) {
  sqrt(sum(x^2) / length(x))
}

fun_pooled_sd_two_samples <- function(x) {
  sqrt(sum(x^2) / 2)
}

fun_pooled_sd_unequal_samples <- function(x) {
  sqrt(sum((length(x) - 1) * x^2) / length(x))
}