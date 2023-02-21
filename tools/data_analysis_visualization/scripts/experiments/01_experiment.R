###############################################################################
# Visualisation of the results of Experiment 1 using the MetaSqueeze model in 
# the following article:
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

# the tree lines below are needed to use plotly::save_image
library(reticulate)
reticulate::py_run_string("import sys")
reticulate::py_run_string("sys.path.insert(0, '/usr/lib/python3/dist-packages')")

# load utils
source("utils/pooled_sd.R")
source("utils/read_multiple_files.R")

# input files
g_path <- "data/sim_outputs/experiments/01_01_experiment"
g_sim_type <- "*metapop.csv"

# output figure
g_figure_name <- "01_01_experiment"
plot_text <- "(a)"

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
g_params <- g_params[, c("sim_id", "val3", "val4")]

# rename columns
base::colnames(g_params) <- c("sim_id", "fire_interval", "fire_size")

# change fire size values to percentages
g_params <- g_params %>%
  dplyr::mutate(
    fire_size = dplyr::case_when(
      fire_size == 530.33 ~ 25,
      fire_size == 1060.66 ~ 50,
      fire_size == 1590.99 ~ 75,
      fire_size == 2121.32 ~ 100,
      fire_size == 2651.65 ~ 125,
      fire_size == 3181.98 ~ 150,
      fire_size == 3712.31 ~ 175
    )
  )

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
# ncol refers to fire_size from 25% to 175% in 25% steps
matrix_mean_year_max <- base::matrix(
  g_df_summ$mean_persistence,
  #nrow = 26,
  nrow = 31,
  ncol = 7
)

#-----------------------------------------------------------------------------
# 3d plot parameters

t <- list(
  #  family = "mono",
  size = 17,
  family = "arial"
)

axx <- list(
  title = "<b>Fire size (%)</b>",
  font = t,
  tickvals = c(25, 50, 75, 100, 125, 150, 175),
  # autotick = TRUE,
  nticks = 7,
  domain = c(25, 175),
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
  domain = c(5, 39),
  range = c(5, 39),
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
  x = seq(25, 175, by = 25),
  y = seq(5, 35, by = 1),
  z = ~matrix_mean_year_max,
  scene = "scene2",
  showscale = FALSE
) %>%
  add_surface() %>%
  layout(
    title = list(
      text = base::paste0("<br><b>  ", plot_text, "</b>"),
      font = list(
        size = 38,
        # family = "sans serif"
        family = "arial"
      ),
      x = 0,
      xref = "paper",
      yref = "paper"
    ),
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

# plotly::save_image(
#   p = p1,
#   file = base::paste0(g_path, "/plots/", g_figure_name, "_persistence_TEST_RETICULATE",".png"),
#   height = 700,
#   width = 700
# )
# 
# plotly::save_image(
#   p = p1,
#   file = base::paste0(g_path, "/plots/", g_figure_name, "_persistence",".svg"),
#   height = 700,
#   width = 700
# )



#############################################

# calculate metapopulation persistence in on fire interval scenario
g_df_one_fi <- g_df %>%
  dplyr::group_by(sim_id) %>%
  dplyr::summarise(
    mean_persistence = base::mean(year_max)
  )

g_df_one_fi <- g_df %>%
  dplyr::group_by(sim_id, study_num) %>%
  dplyr::summarise(
    num_runs = dplyr::n(),
    mean_year_max = base::mean(year_max, na.rm = TRUE),
    sd_year_max = stats::sd(year_max, na.rm = TRUE)
  )

g_df_one_fi <- g_df_one_fi %>%
  dplyr::group_by(sim_id) %>%
  dplyr::summarise(
    num_study = dplyr::n(),
    grand_mean_year_max = base::mean(mean_year_max, na.rm = TRUE),
    grand_sd_year_max = fun_pooled_sd_equal_sample_sizes(sd_year_max),
    plus_sd_year_max = dplyr::if_else(
      (grand_mean_year_max + grand_sd_year_max) > 500,
      500,
      grand_mean_year_max + grand_sd_year_max
    ),
    minus_sd_year_max = dplyr::if_else(
      (grand_mean_year_max - grand_sd_year_max) < 0,
      0,
      grand_mean_year_max - grand_sd_year_max
    )
  )

g_df_one_fi <- base::merge(x = g_params, y = g_df_one_fi, by = "sim_id")

g_df_one_fi <- g_df_one_fi %>% 
  dplyr::filter(fire_interval == 5)

p2 <- ggplot2::ggplot(data = g_df_one_fi) +
  ggplot2::geom_line(
    mapping = aes(
      x = fire_size,
      y = grand_mean_year_max
    ),
    size = 0.8) + 
  ggplot2::geom_ribbon(
    mapping = aes(
      x = fire_size,
      ymin = minus_sd_year_max,
      ymax = plus_sd_year_max
    ),
    alpha = 0.3) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  ggplot2::xlab("Fire size (%)") + 
  ggplot2::ylab("Persistence (years)") + 
  ggplot2::scale_x_continuous(
    breaks = base::seq(from = 25, to = 175, by = 25)
  )

# save plot
ggplot2::ggsave(
  plot = p2,
  filename = base::paste0(g_figure_name, "one_fire_interval_dpi300", ".tiff"),
  path = base::paste0(g_path, "/plots"),
  width = 3.31 * 2,
  height = 4,
  bg = "white",
  dpi = 300
)
