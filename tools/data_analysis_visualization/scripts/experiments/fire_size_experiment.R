###############################################################################
# Visualisation of the results of fire experiment using the MetaSqueeze model 
# in the following article:
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
g_path <- "data/sim_outputs/experiments/fire_size_experiment"
g_sim_type <- "*_fire.csv"

# output figure
g_figure_name <- "fire_size_experiment"
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
g_params <- g_params[, c("sim_id", "val1")]

# rename columns
base::colnames(g_params) <- c("sim_id", "fire_size")

# rename values of fire size from meter to percentage
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
g_df <- fun_read_files(base::paste0(g_path, "/raw"), g_sim_type)


g_df <- base::merge(x = g_params, y = g_df, by = "sim_id")


# calculate the minimum burned area to achieve optimal survival
# (data used in Discussion)
g_df_75 <- g_df %>% dplyr::filter(fire_size == 75)
g_df_100 <- g_df %>% dplyr::filter(fire_size == 100)

median_75 <- stats::median(g_df_75$burned_area_perct)
median_100 <- stats::median(g_df_100$burned_area_perct)


# plot 
p1 <- ggplot(data = g_df) +
  geom_boxplot(
    mapping = aes(
      x = fire_size, 
      y = burned_pops_perct, 
      group = fire_size
    ),
    outlier.size = 0.4
  ) +
  xlab("Fire size (%)") +
  ylab("Burnt dunes (%)") +
  scale_x_continuous(breaks = seq(25, 175, 25)) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  ggtitle("(a)")

p2 <- ggplot(data = g_df) +
  geom_boxplot(
    mapping = aes(
      x = fire_size, 
      y = burned_area, 
      group = fire_size
    )
  ) +
  xlab("Fire size (%)") +
  ylab("Burnt suitable area (ha)") +
  scale_x_continuous(breaks = seq(25, 175, 25)) +
  scale_y_continuous(breaks = seq(0, 550, 100)) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  ggtitle("(b)")


p3 <- ggplot(data = g_df) +
  geom_boxplot(
    mapping = aes(
      x = fire_size, 
      y = burned_area_perct, 
      group = fire_size
    ),
    outlier.size = 0.4
  ) +
  xlab("Fire size (%)") +
  ylab("Burnt suitable area (%)") +
  scale_x_continuous(breaks = seq(25, 175, 25)) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  ggtitle("(b)")

p4 <- ggplot(data = g_df) +
  geom_boxplot(
    mapping = aes(
      x = fire_size, 
      y = burned_area_perct, 
      group = fire_size
    ),
    outlier.size = 0.01
  ) +
  xlab("Fire size (%)") +
  ylab("Burnt suitable area (%)") +
  scale_x_continuous(breaks = seq(25, 175, 25)) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  )

# classify outliers as true in column outlier
g_df <-
  g_df %>%
  dplyr::group_by(fire_size) %>%
  dplyr::mutate(
    outlier = (burned_area_perct > stats::median(burned_area_perct) + 
      stats::IQR(burned_area_perct) * 1.5) | 
      (burned_area_perct < stats::median(burned_area_perct) -
      stats::IQR(burned_area_perct) * 1.5) 
  ) %>%
  dplyr::ungroup(fire_size)

# g_df_outliers <- g_df %>% 
#   dplyr::filter(base::isTRUE(outlier))

p4 <- ggplot() +
  geom_jitter(
    data = g_df %>% 
      dplyr::filter(outlier == TRUE),
    mapping = aes(
      x = fire_size, 
      y = burned_area_perct, 
      group = fire_size),
    alpha = 0.3,
    size = 0.1
  ) +
  geom_boxplot(
    data = g_df, 
    mapping = aes(
      x = fire_size, 
      y = burned_area_perct, 
      group = fire_size
    ),
    color = "#0073C2FF",
    outlier.shape = NA
  ) +
  ggtitle("(d)") +
  xlab("Fire size (%)") +
  ylab("% suitable habitat burned") +
  scale_x_continuous(breaks = seq(25, 175, 25)) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  )

# save plot
ggplot2::ggsave(
  plot = p4,
  filename = base::paste0(g_figure_name, ".png"),
  path = base::paste0(g_path, "/plots"),
  width = 3.31,
  height = 3.31,
  bg = "white"
)
