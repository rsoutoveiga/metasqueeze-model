###############################################################################
# Visualisation of the results of Experiment 4 using the MetaSqueeze model in 
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
library(ggplot2)

# load utils
source("utils/pooled_sd.R")
source("utils/read_multiple_files.R")

# input files
g_root_path <- "data/sim_outputs/experiments"

g_folder_1 <- "/04_01_experiment"
g_folder_2 <- "/04_02_experiment"

# simulation output type
g_sim_type <- "*metapop.csv"

# name output figure
g_figure_name <- "04_experiment"

g_df_all <- base::data.frame()
g_folders <- c(g_folder_1, g_folder_2)


for (i in g_folders) {
  g_path <- base::paste0(g_root_path, i)
  
  # read simulation experiment output files and create df
  g_df <- fun_read_files(base::paste0(g_path, "/raw"), g_sim_type)
  
  # calculate metapopulation persistence
  g_df_summ <- g_df %>%
    dplyr::group_by(sim_id, study_num) %>%
    dplyr::summarise(
      num_runs = dplyr::n(),
      mean_year_max = base::mean(year_max, na.rm = TRUE),
      sd_year_max = stats::sd(year_max, na.rm = TRUE),
    )
  
  g_df_summ <- g_df_summ %>%
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
      ),
    )
  
  
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
  g_params <- g_params[, c("sim_id", "val5")]
  
  # rename columns
  base::colnames(g_params) <- c("sim_id", "num_pops")
  
  g_params$num_pops <- 0:37
  g_params$sort_size <- i
  
  g_df_summ <- base::merge(x = g_params, y = g_df_summ, by = "sim_id")
  

  g_df_all <- base::rbind(g_df_all, g_df_summ)
}

# plot
p1 <- ggplot2::ggplot(data = g_df_all) +
  ggplot2::geom_line(
    mapping = aes(
      x = num_pops,
      y = grand_mean_year_max,
      group = base::as.factor(sort_size), 
      colour = base::as.factor(sort_size)),
    size = 0.8) + 
  ggplot2::geom_ribbon(
    mapping = aes(
      x = num_pops,
      ymin = minus_sd_year_max,
      ymax = plus_sd_year_max,
      group = base::as.factor(sort_size), 
      fill = base::as.factor(sort_size)),
    alpha = 0.30) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  ggplot2::xlab("Dunes as 'high' habitat quality (#)") + 
  ggplot2::ylab("Persistence (years)") + 
  ggplot2::scale_y_continuous(
    breaks = base::seq(from = 0, to = 500, by = 100)
  ) +
  ggplot2::scale_x_continuous(
    breaks = base::seq(from = 0, to = 37, by = 5)
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 500)) +
  ggplot2::scale_color_manual(
    name = "Inter-dune sort",
    labels = c("small to large", "large to small"),
    values = c("#0073C2FF", "#EFC000FF")
    # values = c("#56B4E9", "#E69F00")
  ) +
  ggplot2::scale_fill_manual(
    name = "Inter-dune sort",
    labels = c("small to large", "large to small"),
    values = c("#0073C2FF", "#EFC000FF")
    # values = c("#56B4E9", "#E69F00")
  ) +
  ggplot2::geom_vline(
    mapping = aes(
      xintercept = 18, 
      linetype = "dunes occupied\n at the start of \n the simulations"), 
    colour = "black") +
  ggplot2::scale_linetype_manual(name = "", values = c("dashed")) 


ggplot2::ggsave(
  plot = p1,
  filename = base::paste0(g_figure_name, "_dpi300.tiff"),
  path = base::paste0(g_root_path, g_folder_1, "/plots"),
  width = 3.31 * 2,
  height = 3.31,
  bg = "white",
  dpi = 300
)

# save the same image in all folders of the experiment
# ggplot2::ggsave(
#   plot = p1,
#   filename = base::paste0(g_figure_name, ".png"),
#   path = base::paste0(g_root_path, g_folder_1, "/plots"),
#   width = 7,
#   height = 4,
#   bg = "white"
# )
# 
# 
# ggplot2::ggsave(
#   plot = p1,
#   filename = base::paste0(g_figure_name, ".png"),
#   path = base::paste0(g_root_path, g_folder_2, "/plots"),
#   width = 7,
#   height = 4,
#   bg = "white"
# )
# 
# 
# ggplot2::ggsave(
#   plot = p1,
#   filename = base::paste0(g_figure_name, "dpi100.png"),
#   path = base::paste0(g_root_path, g_folder_1, "/plots"),
#   width = 3.31 * 2,
#   height = 3.31,
#   bg = "white",
#   dpi = 100
# )
# 
# ggplot2::ggsave(
#   plot = p1,
#   filename = base::paste0(g_figure_name, "dpi1200.png"),
#   path = base::paste0(g_root_path, g_folder_1, "/plots"),
#   width = 3.31 * 2,
#   height = 3.31,
#   bg = "white",
#   dpi = 1200
# )
# 
# ggplot2::ggsave(
#   plot = p1,
#   filename = base::paste0(g_figure_name, "dpi1200.tiff"),
#   path = base::paste0(g_root_path, g_folder_1, "/plots"),
#   width = 3.31 * 2,
#   height = 3.31,
#   bg = "white",
#   dpi = 1200
# )
