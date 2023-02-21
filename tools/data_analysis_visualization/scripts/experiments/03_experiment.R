###############################################################################
# Visualisation of the results of Experiment 3 using the MetaSqueeze model in 
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
library(ggplot2)

# load utils
source("utils/pooled_sd.R")
source("utils/read_multiple_files.R")

# input files
g_path <- "./data/sim_outputs/experiments/03_experiment"
g_sim_type <- "*metapop.csv"

# output figure
g_figure_name <- "03_experiment"
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
base::colnames(g_params) <- c("sim_id", "habitat", "survival")

# rename parameter values
g_params <- g_params %>%
  dplyr::mutate(
    survival = dplyr::case_when(
      survival == "plant_types.csv" ~ 0,
      survival == "plant_types_3.csv" ~ 3,
      survival == "plant_types_6.csv" ~ 6,
      survival == "plant_types_9.csv" ~ 9,
      survival == "plant_types_12.csv" ~ 12,
      survival == "plant_types_15.csv" ~ 15,
      survival == "plant_types_18.csv" ~ 18,
      survival == "plant_types_21.csv" ~ 21,
      survival == "plant_types_24.csv" ~ 24,
      survival == "plant_types_27.csv" ~ 27
    )
  )

g_params <- g_params %>%
  dplyr::mutate(
    habitat = dplyr::case_when(
      habitat == "metapop_baseline.csv" ~ 1,
      habitat == "metapop_baseline_medium.csv" ~ 2,
      habitat == "metapop_baseline_high.csv" ~ 3
    )
  )

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

g_df_summ <- base::merge(x = g_params, y = g_df_summ, by = "sim_id")

p1 <- ggplot2::ggplot(data = g_df_summ) +
  ggplot2::geom_line(
    mapping = aes(
      x = survival,
      y = grand_mean_year_max,
      group = base::as.factor(habitat), 
      colour = base::as.factor(habitat)),
    size = 0.8) + 
  ggplot2::geom_ribbon(
    mapping = aes(
      x = survival,
      ymin = minus_sd_year_max,
      ymax = plus_sd_year_max,
      group = base::as.factor(habitat), 
      fill = base::as.factor(habitat)),
    alpha = 0.3) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  ggplot2::xlab("Survival sum (%)") + 
  ggplot2::ylab("Persistence (years)") + 
  ggplot2::scale_y_continuous(
    breaks = base::seq(from = 0, to = 500, by = 100)
  ) +
  ggplot2::scale_x_continuous(
    breaks = base::seq(from = 0, to = 27, by = 3)
  ) +
  ggplot2::scale_color_manual(
    name = "Dune habitat quality",
    breaks = c(1, 2, 3),
    labels = c("low", "moderate", "high"),
    values = c("#868686FF", "#0073C2FF", "#EFC000FF")
    # values = c("#56B4E9", "#E69F00")
  ) +
  ggplot2::scale_fill_manual(
    name = "Dune habitat quality",
    breaks = c(1, 2, 3),
    labels = c("low", "moderate", "high"),
    values = c("#868686FF", "#0073C2FF", "#EFC000FF")
    # values = c("#56B4E9", "#E69F00")
  )

# ggplot2::ggsave(
#   plot = p1,
#   filename = base::paste0(g_figure_name, "_dpi100.png"),
#   path = base::paste0(g_path, "/plots"),
#   width = 3.31 * 2,
#   height = 3.31,
#   bg = "white",
#   dpi = 100
# )

ggplot2::ggsave(
  plot = p1,
  filename = base::paste0(g_figure_name, "_dpi300.tiff"),
  path = base::paste0(g_path, "/plots"),
  width = 3.31 * 2,
  height = 3.31,
  bg = "white",
  dpi = 300
)