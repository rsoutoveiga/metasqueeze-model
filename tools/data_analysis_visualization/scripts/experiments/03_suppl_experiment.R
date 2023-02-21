###############################################################################
# Visualisation of the results of Experiment 3 (Suppl.) using the MetaSqueeze  
# model in the following article:
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
g_path <- "data/sim_outputs/experiments/03_suppl_experiment"
g_sim_type <- "*metapop.csv"

# output figure
g_figure_name <- "03_suppl_experiment"
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
g_params <- g_params[, c("sim_id", "val2", "val3")]

# rename columns
base::colnames(g_params) <- c("sim_id", "fire_interval", "model")

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

g_df_summ$model <- base::as.factor(g_df_summ$model)


# calculate persistence difference between mortality scenarios
g_df_model_0 <- g_df_summ %>% 
  dplyr::filter(model == 0)

g_df_model_1 <- g_df_summ %>% 
  dplyr::filter(model == 1)

g_df_model_1$difference <- g_df_model_1$grand_mean_year_max - g_df_model_0$grand_mean_year_max

base::min(g_df_model_1$difference)
base::max(g_df_model_1$difference)


# plot
p1 <- ggplot2::ggplot(data = g_df_summ) +
  ggplot2::geom_line(
    mapping = aes(
      x = fire_interval,
      y = grand_mean_year_max,
      group = model, 
      colour = model),
    size = 0.8) + 
  ggplot2::geom_ribbon(
    mapping = aes(
      x = fire_interval,
      ymin = minus_sd_year_max,
      ymax = plus_sd_year_max,
      group = model, 
      fill = model),
    alpha = 0.3) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  ggplot2::xlab("Fire interval (years)") + 
  ggplot2::ylab("Persistence (years)") + 
  ggplot2::scale_y_continuous(
    breaks = base::seq(from = 50, to = 200, by = 25)
  ) +
  ggplot2::scale_x_continuous(
    breaks = base::seq(from = 5, to = 35, by = 5)
  ) +
  ggplot2::coord_cartesian(ylim = c(50, 200)) +
  ggplot2::scale_color_manual(
    name = "Intra-dune variability",
    breaks = c(0, 1),
    labels = c("none", "flower"),
    values = c("#0073C2FF", "#EFC000FF")
    # values = c("#56B4E9", "#E69F00")
  ) + 
  ggplot2::scale_fill_manual(
    name = "Intra-dune variability",
    breaks = c(0, 1),
    labels = c("none", "flower"),
    values = c("#0073C2FF", "#EFC000FF")
    # values = c("#56B4E9", "#E69F00")
  )


ggplot2::ggsave(
  plot = p1,
  filename = base::paste0(g_figure_name, "_dpi300.tiff"),
  path = base::paste0(g_path, "/plots"),
  width = 3.31 * 2,
  height = 3.31,
  bg = "white",
  dpi = 300
)

