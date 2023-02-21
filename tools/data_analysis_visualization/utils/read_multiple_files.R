#'@description read multiple csv files and create dataframe
#'@author Rodrigo Souto-Veiga

# load libraries
library(plyr)
library(dplyr)

fun_read_files <- function(file_path, sim_type) {
  
  # names of the .csv files
  g_files <- base::list.files(
    path = file_path,
    pattern = sim_type,
    full.names = TRUE
  )
  
  # create df
  plyr::ldply(
    g_files,
    utils::read.csv,
    sep = ",",
    fill = TRUE,
    header = TRUE
  )
}
