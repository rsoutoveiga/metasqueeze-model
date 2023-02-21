###############################################################################
# Visualisation of the results of Experiment 4 (Supplementary 1) using the 
# MetaSqueeze model in the following article:
# 
# **Souto-Veiga, R., Groeneveld, J., Enright, N. J., Fontaine, J. B., &
# Jeltsch, F. (20XX). Climate change may shift metapopulations towards
# unstable source-sink dynamics in a fire-killed, serotinous shrub**
# 
# Corresponding author:
# 
# Rodrigo Souto-Veiga^1, 2^
# 
# [rsoutoveiga\@uni-potsdam.de](rsoutoveiga@uni-potsdam.de)
# 
# ORCID: [0000-0001-8639-620X](https://orcid.org/0000-0001-8639-620X)
# 
# ^1^ Plant Ecology and Nature Conservation, University of Potsdam, Am
# Mühlenberg 3, 14476, Potsdam, Germany
# 
# ^2^ Environmental and Conservation Sciences, Murdoch University, Murdoch
# 6150, WA, Australia
###############################################################################

# load libraries
library(dplyr)
library(viridis)
library(plotly)

library(reticulate)
reticulate::py_run_string("import sys")
reticulate::py_run_string("sys.path.insert(0, '/usr/lib/python3/dist-packages')")

# load utils
source("utils/pooled_sd.R")
source("utils/read_multiple_files.R")

# input files
g_path <- "data/sim_outputs/experiments/04_suppl_1_experiment"
g_sim_type <- "*metapop.csv"

# output figure
g_figure_name <- "04_suppl_1_experiment"
plot_text <- ""

# read simulation parameters
g_params <- base::list.files(
  path = base::paste0(g_path, "/in/sim"),
  pattern = "*_simfile.csv",
  full.names = TRUE
)

g_params <- plyr::ldply(
  g_params,
  utils::read.table,
  sep = " ",
  fill = TRUE,
  header = TRUE
)

# take the parameter values for plotting
g_params <- g_params[, c("sim_id", "val4", "val5")]

# rename columns
base::colnames(g_params) <- c("sim_id", "fire_interval", "num_dunes")

# sort by fire interval and num of dunes
g_params <- g_params[base::order(g_params$fire_interval, g_params$num_dunes), ]

# change values of number of dunes
g_params$num_dunes <- 0:37

# read simulation experiment output files and create df
g_df <- fun_read_files(paste0(g_path, "/raw"), g_sim_type)

# calculate metapopulation persistence
g_df_summ <- g_df %>%
  dplyr::group_by(sim_id) %>%
  dplyr::summarise(
    mean_persistence = base::mean(year_max)
  )

g_df_summ <- base::merge(x = g_params, y = g_df_summ, by = "sim_id")

# nrow refers to fire_mean from 5 years to 30 years in 1 year steps
# ncol refers to number of high dunes from 0 to 37 in 1 steps
matrix_mean_year_max <- base::matrix(
  g_df_summ$mean_persistence,
  nrow = 26,
  ncol = 38
)

#-----------------------------------------------------------------------------
# 3d plot parameters

t <- list(
  #  family = "mono",
  size = 17,
  family = "arial"
)

axx <- list(
  title = "<b>'High' dunes (#)</b>",
  font = t,
  tickvals = seq(0, 37, 5),
  # autotick = TRUE,
  # nticks = 7,
  # domain = c(25, 175),
  showgrid = T,
  tickmode = "array",
  titlefont = list(size = 26),
  tickcolor = "grey",
  ticklen = 0.1,
  tickwidth = 0.4
)

axy <- list(
  title = "<b>Fire interval (yr)</b>",
  font = t,
  zeroline = T,
  autotick = T,
  nticks = 10,
  start = 5,
  domain = c(5, 30),
  range = c(5, 33),
  titlefont = list(size = 26)
)

axz <- list(
  title = "<b>Persistence (yr)</b>",
  nticks = 6,
  domain = c(0, 501),
  range = c(0, 501),
  titlefont = list(size = 26)
)

p1 <- plotly::plot_ly(
  x = seq(0, 37, by = 1),
  y = seq(5, 30, by = 1),
  z = ~matrix_mean_year_max,
  scene = "scene2",
  showscale = FALSE
) %>%
  add_surface() %>%
  layout(
    # title = list(
    #   text = base::paste0("<br><b>  ", plot_text, "</b>"),
    #   font = list(
    #     size = 38,
    #     family = "arial"
    #   ),
    #   x = 0,
    #   xref = "paper",
    #   yref = "paper"
    # ),
    margin = list(
      l = 0,
      r = 0,
      t = 30,
      b = 0
    ),
    font = t,
    scene2 = list(
      camera = list(
        eye = list(
          x = -1.25, y = -0.75, z = 1.85
        )
      ),
      xaxis = axx,
      yaxis = axy,
      zaxis = axz,
      aspectmode = "cube", #this string can be 'data', 'cube', 'auto', 'manual'
      aspectratio = list(
        x = 1, y = 1, z = 1
      )
    )
  )

plotly::save_image(
  p = p1,
  file = base::paste0(g_path, "/plots/", g_figure_name, "_persistence",".png"),
  height = 700,
  width = 700
)

plotly::save_image(
  p = p1,
  file = base::paste0(g_path, "/plots/", g_figure_name, "_persistence",".svg"),
  height = 700,
  width = 700
)
