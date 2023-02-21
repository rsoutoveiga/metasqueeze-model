###############################################################################
# Visualisation of landscape generated in the MetaSqueeze model. Used in the 
# following article:
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

# load necessary libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

source("./scripts/experiments/fire_size_experiment.R")

g_path <- "./data/sim_outputs/study_site_eneabba/"
g_file_name <- "seed_13_study.csv"
g_figure_name <- "patchy_fires"

landscape = utils::read.table(base::paste0(g_path, g_file_name), sep = ",")

landscape <- base::as.matrix(landscape)
base::colnames(landscape) <- NULL

landscape_raster <- landscape

landscape = reshape2::melt(landscape)

base::colnames(landscape)[1] <- "row"
base::colnames(landscape)[2] <- "col"

# landscape$value[landscape$value == 0] <- NA

landscape <- landscape %>% 
  dplyr::mutate(value = base::ifelse(value > 1, 1, value))

# p1 <- ggplot(landscape, aes(x=col, y=row, fill=value)) +
#   geom_raster(colour="grey80") +
#   scale_fill_gradientn(colours = terrain.colors(18)) +
#   coord_equal()
# 
# p2 <- ggplot(landscape, aes(x=col, y=row, fill=value)) +
#   geom_tile(colour="grey20") +
#   scale_fill_viridis_c(option = "cividis", direction = 1) +
#   coord_equal()
# 
# 
pland1 <- ggplot(data = landscape) +
  geom_tile(aes(x = col, y = row, fill = as.factor(value)), color = "grey80") +
  scale_fill_manual(
    name = "Habitat:",
    breaks = c(0, 1),
    labels = c("unsuitable", "suitable"),
    values = c("white", "#009E73")) +
  coord_equal() +
  theme_void() +
  ggtitle("(a)\n") +
  ggplot2::theme(
    text = element_text(size = 12, face = "bold"),
    # axis.text = element_text(size = 14),
    # axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(20, 20, 20, 20))

pland2 <- ggplot(data = landscape) +
  geom_tile(aes(x = col, y = row, fill = as.factor(value)), colour = "grey80") +
  scale_fill_manual(
    name = "Habitat:",
    breaks = c(0, 1),
    labels = c("unsuitable", "suitable"),
    values = c("white", "#009E73")) +
  coord_equal() +
  theme_void() +
  ggtitle("(b)\n") +
  ggplot2::theme(
    text = element_text(size = 12, face = "bold"),
    # axis.text = element_text(size = 14),
    # axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(20, 20, 20, 20))

####

g_path <- "./data/sim_outputs/study_site_eneabba/"
g_file_name <- "seed_13_study_small.csv"

landscape = utils::read.table(base::paste0(g_path, g_file_name), sep = ",")

landscape <- base::as.matrix(landscape)
base::colnames(landscape) <- NULL

landscape_raster <- landscape

landscape = reshape2::melt(landscape)

base::colnames(landscape)[1] <- "row"
base::colnames(landscape)[2] <- "col"


pland3 <- ggplot(data = landscape) +
  geom_tile(aes(x = col, y = row, fill = as.factor(value)), colour = "grey50") +
  scale_fill_manual(
    name = "Fire feedback:",
    breaks = c(0, 1, 2, 3),
    labels = c("unsuitable", "unburnt", "insufficient\nfuel load", "burnt"),
    values = c("white", "#009E73", "#0073C2FF", "#EFC000FF")) +
  coord_equal() +
  theme_void() +
  ggtitle("(c)") +
  # theme(legend.position = "bottom")
  ggplot2::theme(
    text = element_text(size = 12, face = "bold"),
    # axis.text = element_text(size = 14),
    # axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom"
  )

parr1 <- ggpubr::ggarrange(pland1, pland2, ncol = 2, common.legend = TRUE, legend = "bottom")

parr2 <- ggpubr::ggarrange(pland3, p4, ncol = 2, common.legend = TRUE, legend = "bottom")

parr3 <- ggpubr::ggarrange(parr1, parr2, ncol = 1, nrow = 2)



ggplot2::ggsave(
  plot = parr3,
  filename = base::paste0(g_figure_name, "closer_new_new.tiff"),
  path = base::paste0(g_path, "/plots"),
  width = 3.31 * 2,
  height = 3.31 * 2,
  bg = "white",
  dpi = 1200
)
