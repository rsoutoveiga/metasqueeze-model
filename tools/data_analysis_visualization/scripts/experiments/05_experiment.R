###############################################################################
# Visualisation of the results of Experiment 5 using the MetaSqueeze model in 
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
library(ggpubr)
library(ggrepel)

# load utils
source("utils/pooled_sd.R")
source("utils/read_multiple_files.R")

# input files
g_path <- "data/sim_outputs/experiments/06_experiment"

# simulation output type
g_sim_type <- "*_metapop_pops.csv"

# name output figure
g_figure_name <- "06_experiment"

# sources patches (dunes)
g_source_patches <- c(3, 4, 6, 9, 1, 18, 15)

# read simulation experiment output files and create df
g_df <- fun_read_files(base::paste0(g_path, "/raw"), g_sim_type)

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
base::colnames(g_params) <- c("sim_id", "metapop_file")


g_df_metapop <- plyr::ldply(
  base::paste0(g_path, "/in/metapopulation/", g_params$metapop_file),
  utils::read.table,
  sep = " ",
  fill = TRUE,
  header = TRUE
)

g_df_metapop <- g_df_metapop %>% 
  dplyr::mutate(is_occupied = base::ifelse(init_individuals == 0, 0, 1))

g_df_metapop$is_occupied <- base::as.factor(g_df_metapop$is_occupied)

g_df_metapop <- g_df_metapop[, c("code", "is_occupied")]

base::colnames(g_df_metapop)[1] <- "pop_id"

g_df <- g_df %>% 
  dplyr::left_join(g_df_metapop, by = "pop_id")

# classification of patches (dunes) into sources and sinks
g_df <- g_df %>% 
  dplyr::mutate(
    is_source = dplyr::if_else(
      pop_id %in% g_source_patches,
      "source",
      "sink"
    )
  )

# calculations
g_df_summ <- g_df %>%
  dplyr::group_by(sim_id, study_num, is_source, is_occupied) %>%
  dplyr::summarise(
    mean_year_max = base::mean(year_max, na.rm = TRUE),
    sd_year_max = stats::sd(year_max, na.rm = TRUE),
    mean_persistence = base::mean(persistence, na.rm = TRUE),
    sd_persistence = stats::sd(persistence, na.rm = TRUE),
    mean_generations = base::mean(generations, na.rm = TRUE),
    sd_generations = stats::sd(generations, na.rm = TRUE),
    mean_burned = base::mean(burned, na.rm = TRUE),
    sd_burned = stats::sd(burned, na.rm = TRUE)
  )

g_df_summ <- g_df_summ %>%
  dplyr::group_by(sim_id, is_source, is_occupied) %>%
  dplyr::summarise(
    num_study = dplyr::n(),
    grand_mean_year_max = base::mean(mean_year_max, na.rm = TRUE),
    grand_sd_year_max = fun_pooled_sd_equal_sample_sizes(sd_year_max),
    grand_mean_persistence = base::mean(mean_persistence, na.rm = TRUE),
    grand_sd_persistence = fun_pooled_sd_equal_sample_sizes(sd_persistence),
    grand_mean_generations = base::mean(mean_generations, na.rm = TRUE),
    grand_sd_generations = fun_pooled_sd_equal_sample_sizes(sd_generations),
    grand_mean_burned = base::mean(mean_burned, na.rm = TRUE),
    grand_sd_burned = fun_pooled_sd_equal_sample_sizes(sd_burned)
  )

#display factor levels for region
g_df_summ$sim_id <- base::as.factor(g_df_summ$sim_id)
base::levels(g_df_summ$sim_id) <- base::c("baseline", "current")

# plot
# subplot (d): persistence time
pd <- ggplot2::ggplot(data = g_df_summ, mapping = ggplot2::aes(x = sim_id)) +
  ggplot2::geom_line(
    mapping = aes(
      y = grand_mean_persistence,
      group = base::interaction(is_source, is_occupied),
      colour = base::interaction(is_source, is_occupied),
    ),
    position = ggplot2::position_dodge(width = 0.3),
    size = 1,
    alpha = 0.4
  ) +
  ggplot2::geom_point(
    mapping = aes(
      y = grand_mean_persistence,
      group = base::interaction(is_source, is_occupied),
      colour = base::interaction(is_source, is_occupied),
      shape = base::interaction(is_source, is_occupied)
    ),
    position = ggplot2::position_dodge(width = 0.3),
    size = 3
  ) +
  ggplot2::geom_errorbar(
    mapping = aes(
      ymin = dplyr::if_else(
        (grand_mean_persistence - grand_sd_persistence) < 0,
        0, 
        grand_mean_persistence - grand_sd_persistence
      ),
      ymax = dplyr::if_else(
        (grand_mean_persistence + grand_sd_persistence) > 500,
        500, 
        grand_mean_persistence + grand_sd_persistence
      ),
      group = base::interaction(is_source, is_occupied),
      colour = base::interaction(is_source, is_occupied),
    ),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 500)) +
  ggplot2::scale_y_continuous(limits = c(0, 500)) +
  ggplot2::theme_classic() +
  ggplot2::ggtitle("(d)") +
  ggplot2::xlab("Scenario") +
  ggplot2::ylab("Persistence (years)") +
  ggplot2::scale_color_manual(
    name = "Dune type:",
    values = c("#868686FF", "#0073C2FF", "#EFC000FF"),
    labels = c("sink unoccupied", "sink occupied", "source")
  ) +
  ggplot2::scale_shape_manual(
    name = "Dune type:",
    values = c(18, 17, 15),
    labels = c("sink unoccupied", "sink occupied", "source")
  ) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 12, face = "bold"),
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

###########################################################################
### subplots (a) and (b): rescue effects and recolonizations ##############
###########################################################################

# simulation output type
g_sim_type <- "*_rescue.csv"

# read simulation experiment output files and create df
g_df <- fun_read_files(base::paste0(g_path, "/raw"), g_sim_type)

# add the initial empty patches/dunes as extinct dunes
# should be this counted in the model (c++ code)??
g_df$extinctions <- g_df$extinctions + 19

# calculations
g_df_summ <- g_df %>%
  dplyr::group_by(sim_id, study_num) %>%
  dplyr::summarise(
    mean_year_max = base::mean(year_max, na.rm = TRUE),
    sd_year_max = stats::sd(year_max, na.rm = TRUE),
    
    mean_generations = base::mean(generations, na.rm = TRUE),
    sd_generations = stats::sd(generations, na.rm = TRUE),
    
    mean_immi_source = base::mean(immi_source, na.rm = TRUE),
    sd_immi_source = stats::sd(immi_source, na.rm = TRUE),
    
    mean_immi_adult = base::mean(immi_adult, na.rm = TRUE),
    sd_immi_adult = stats::sd(immi_adult, na.rm = TRUE),
    
    mean_immi_seedling = base::mean(immi_seedling, na.rm = TRUE),
    sd_immi_seedling = stats::sd(immi_seedling, na.rm = TRUE),
    
    mean_immi_germi = base::mean(immi_germi, na.rm = TRUE),
    sd_immi_germi = stats::sd(immi_germi, na.rm = TRUE),
    
    
    mean_extinctions = base::mean(extinctions, na.rm = TRUE),
    sd_extinctions = stats::sd(extinctions, na.rm = TRUE),
    
    mean_recol_source = base::mean(recol_source, na.rm = TRUE),
    sd_recol_source = stats::sd(recol_source, na.rm = TRUE),
    
    mean_recol_adult = base::mean(recol_adult, na.rm = TRUE),
    sd_recol_adult = stats::sd(recol_adult, na.rm = TRUE),
    
    mean_recol_seedling = base::mean(recol_seedling, na.rm = TRUE),
    sd_recol_seedling = stats::sd(recol_seedling, na.rm = TRUE),
    
    mean_recol_germi = base::mean(recol_germi, na.rm = TRUE),
    sd_recol_germi = stats::sd(recol_germi, na.rm = TRUE)
  )

g_df_summ <- g_df_summ %>%
  dplyr::group_by(sim_id) %>%
  dplyr::summarise(
    num_study = dplyr::n(),
    
    grand_mean_year_max = base::mean(mean_year_max, na.rm = TRUE),
    grand_sd_year_max = base::sqrt(
      base::sum(sd_year_max ^ 2, na.rm = TRUE) / num_study),
    
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
    
    
    ##########
    grand_mean_generations = base::mean(mean_generations, na.rm = TRUE),
    grand_sd_generations = base::sqrt(
      base::sum(sd_generations ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_generations = dplyr::if_else(
      (grand_mean_generations - grand_sd_generations) < 0,
      0,
      grand_mean_generations - grand_sd_generations
    ),
    
    plus_sd_generations = grand_mean_generations + grand_sd_generations,
    
    grand_mean_immi_source = base::mean(mean_immi_source, na.rm = TRUE),
    grand_sd_immi_source = base::sqrt(
      base::sum(sd_immi_source ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_immi_source = dplyr::if_else(
      (grand_mean_immi_source - grand_sd_immi_source) < 0,
      0,
      grand_mean_immi_source - grand_sd_immi_source
    ),
    
    plus_sd_immi_source = grand_mean_immi_source + grand_sd_immi_source,
    
    grand_mean_immi_adult = base::mean(mean_immi_adult, na.rm = TRUE),
    grand_sd_immi_adult = base::sqrt(
      base::sum(sd_immi_adult ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_immi_adult = dplyr::if_else(
      (grand_mean_immi_adult - grand_sd_immi_adult) < 0,
      0,
      grand_mean_immi_adult - grand_sd_immi_adult
    ),
    
    plus_sd_immi_adult = grand_mean_immi_adult + grand_sd_immi_adult,
    
    grand_mean_immi_seedling = base::mean(mean_immi_seedling, na.rm = TRUE),
    grand_sd_immi_seedling = base::sqrt(
      base::sum(sd_immi_seedling ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_immi_seedling = dplyr::if_else(
      (grand_mean_immi_seedling - grand_sd_immi_seedling) < 0,
      0,
      grand_mean_immi_seedling - grand_sd_immi_seedling
    ),
    
    plus_sd_immi_seedling = grand_mean_immi_seedling + grand_sd_immi_seedling,
    
    grand_mean_immi_germi = base::mean(mean_immi_germi, na.rm = TRUE),
    grand_sd_immi_germi = base::sqrt(
      base::sum(sd_immi_germi ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_immi_germi = dplyr::if_else(
      (grand_mean_immi_germi - grand_sd_immi_germi) < 0,
      0,
      grand_mean_immi_germi - grand_sd_immi_germi
    ),
    
    plus_sd_immi_germi = grand_mean_immi_germi + grand_sd_immi_germi,
    #########
    ##########
    grand_mean_extinctions = base::mean(mean_extinctions, na.rm = TRUE),
    grand_sd_extinctions = base::sqrt(
      base::sum(sd_extinctions ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_extinctions = dplyr::if_else(
      (grand_mean_extinctions - grand_sd_extinctions) < 0,
      0,
      grand_mean_extinctions - grand_sd_extinctions
    ),
    
    plus_sd_extinctions = grand_mean_extinctions + grand_sd_extinctions,
    
    grand_mean_recol_source = base::mean(mean_recol_source, na.rm = TRUE),
    grand_sd_recol_source = base::sqrt(
      base::sum(sd_recol_source ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_recol_source = dplyr::if_else(
      (grand_mean_recol_source - grand_sd_recol_source) < 0,
      0,
      grand_mean_recol_source - grand_sd_recol_source
    ),
    
    plus_sd_recol_source = grand_mean_recol_source + grand_sd_recol_source,
    
    grand_mean_recol_adult = base::mean(mean_recol_adult, na.rm = TRUE),
    grand_sd_recol_adult = base::sqrt(
      base::sum(sd_recol_adult ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_recol_adult = dplyr::if_else(
      (grand_mean_recol_adult - grand_sd_recol_adult) < 0,
      0,
      grand_mean_recol_adult - grand_sd_recol_adult
    ),
    
    plus_sd_recol_adult = grand_mean_recol_adult + grand_sd_recol_adult,
    
    grand_mean_recol_seedling = base::mean(mean_recol_seedling, na.rm = TRUE),
    grand_sd_recol_seedling = base::sqrt(
      base::sum(sd_recol_seedling ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_recol_seedling = dplyr::if_else(
      (grand_mean_recol_seedling - grand_sd_recol_seedling) < 0,
      0,
      grand_mean_recol_seedling - grand_sd_recol_seedling
    ),
    
    plus_sd_recol_seedling = grand_mean_recol_seedling + grand_sd_recol_seedling,
    
    grand_mean_recol_germi = base::mean(mean_recol_germi, na.rm = TRUE),
    grand_sd_recol_germi = base::sqrt(
      base::sum(sd_recol_germi ^ 2, na.rm = TRUE) / num_study),
    
    minus_sd_recol_germi = dplyr::if_else(
      (grand_mean_recol_germi - grand_sd_recol_germi) < 0,
      0,
      grand_mean_recol_germi - grand_sd_recol_germi
    ),
    
    plus_sd_recol_germi = grand_mean_recol_germi + grand_sd_recol_germi
  )

g_df_summ <- g_df_summ %>%
  dplyr::group_by(sim_id) %>% 
  dplyr::mutate(
    rescue_source = base::abs(
      ((grand_mean_generations - grand_mean_immi_source) / 
         grand_mean_generations) - 1) * 100.0,
    rescue_source_ymax = base::abs(
      ((plus_sd_generations - plus_sd_immi_source) / 
         plus_sd_generations) - 1) * 100.0,
    rescue_source_ymin = base::abs(
      ((minus_sd_generations - minus_sd_immi_source) / 
         minus_sd_generations) - 1) * 100.0,
    rescue_adult = base::abs( 
      ((grand_mean_generations - grand_mean_immi_adult) / 
         grand_mean_generations) - 1) * 100.0,
    rescue_adult_ymax = base::abs(
      ((plus_sd_generations - plus_sd_immi_adult) / 
         plus_sd_generations) - 1) * 100.0,
    rescue_adult_ymin = base::abs(
      ((minus_sd_generations - minus_sd_immi_adult) / 
         minus_sd_generations) - 1) * 100.0,
    rescue_seedling = base::abs(
      ((grand_mean_generations - grand_mean_immi_seedling) / 
         grand_mean_generations) - 1) * 100.0,
    rescue_seedling_ymax = base::abs(
      ((plus_sd_generations - plus_sd_immi_seedling) / 
         plus_sd_generations) - 1) * 100.0,
    rescue_seedling_ymin = base::abs(
      ((minus_sd_generations - minus_sd_immi_seedling) / 
         minus_sd_generations) - 1) * 100.0,
    rescue_germi = base::abs(
      ((grand_mean_generations - grand_mean_immi_germi) / 
         grand_mean_generations) - 1) * 100.0,
    rescue_germi_ymax = base::abs(
      ((plus_sd_generations - plus_sd_immi_germi) / 
         plus_sd_generations) - 1) * 100.0,
    rescue_germi_ymin = base::abs(
      ((minus_sd_generations - minus_sd_immi_germi) / 
         minus_sd_generations) - 1) * 100.0
  )

g_df_summ <- g_df_summ %>%
  dplyr::group_by(sim_id) %>% 
  dplyr::mutate(
    recol_source = base::abs(
      ((grand_mean_extinctions - grand_mean_recol_source) / 
         grand_mean_extinctions) - 1) * 100.0,
    recol_adult = base::abs(
      ((grand_mean_extinctions - grand_mean_recol_adult) / 
         grand_mean_extinctions) - 1) * 100.0,
    recol_seedling = base::abs(
      ((grand_mean_extinctions - grand_mean_recol_seedling) / 
         grand_mean_extinctions) - 1) * 100.0,
    recol_germi = base::abs(
      ((grand_mean_extinctions - grand_mean_recol_germi) / 
         grand_mean_extinctions) - 1) * 100.0
  )


g_df_summ$sim_id <- base::as.factor(g_df_summ$sim_id)
base::levels(g_df_summ$sim_id) <- c("baseline", "current")
# g_df_summ$sim_id <- base::factor(g_df_summ$sim_id, levels = c("baseline", "current"))

# plot
pa <- ggplot(data = g_df_summ, aes(x = sim_id)) +
  geom_line(
    mapping = aes(
      y = recol_germi,
      group = 1,
      colour = "seed"
    ),
    size = 1,
    alpha = 0.5
  ) +
  geom_line(
    mapping = aes(
      y = recol_seedling,
      group = 1,
      colour = "seedling"
    ),
    size = 1,
    alpha = 0.5
  ) +
  geom_line(
    mapping = aes(
      y = recol_adult,
      group = 1,
      colour = "adult"
    ),
    size = 1,
    alpha = 0.5
  ) +
  geom_line(
    mapping = aes(
      y = recol_source,
      group = 1,
      colour = "source"
    ),
    size = 1,
    alpha = 0.5
  ) +
  geom_point(
    mapping = aes(
      y = recol_germi,
      colour = "seed",
      shape = "seed"
    ),
    size = 3
  ) +
  geom_point(
    mapping = aes(
      y = recol_seedling,
      group = sim_id,
      colour = "seedling",
      shape = "seedling"
    ),
    size = 3
  ) +
  geom_point(
    mapping = aes(
      y = recol_adult,
      colour = "adult",
      shape = "adult",
      group = sim_id
    ),
    size = 3
  ) +
  geom_point(
    mapping = aes(
      y = recol_source,
      colour = "source",
      shape = "source",
      group = sim_id
    ),
    size = 3
  ) +
  theme_classic() +
  ggtitle("(a)") +
  xlab("Scenario") +
  ylab("Recolonizations (%)") +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_color_manual(
    name = "Plant stage:",
    breaks = c("seed", "seedling", "adult", "source"),
    values = c("#999999", "#009E73", "#56B4E9", "#E69F00")
  ) +
  scale_shape_manual(
    name = "Plant stage:",
    breaks = c("seed", "seedling", "adult", "source"),
    values = c(2, 3, 4, 5)
  ) +
  theme(
    text         = element_text(size = 12, face = "bold"),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )


pb <- ggplot(data = g_df_summ, aes(x = sim_id)) +
  geom_line(
    mapping = aes(
      y = rescue_germi,
      group = 1,
      colour = "seed"
    ),
    size = 1,
    alpha = 0.5
  ) +
  geom_line(
    mapping = aes(
      y = rescue_seedling,
      group = 1,
      colour = "seedling"
    ),
    size = 1,
    alpha = 0.5
  ) +
  geom_line(
    mapping = aes(
      y = rescue_adult,
      group = 1,
      colour = "adult"
    ),
    size = 1,
    alpha = 0.5
  ) +
  geom_line(
    mapping = aes(
      y = rescue_source,
      group = 1,
      colour = "source"
    ),
    size = 1,
    alpha = 0.5
  ) +
  geom_point(
    mapping = aes(
      y = rescue_germi,
      # group = sim_id,
      colour = "seed",
      shape = "seed"
    ),
    size = 3
    # width = 0.1,
    # position=position_dodge(width=0.3)
  ) +
  geom_point(
    mapping = aes(
      y = rescue_seedling,
      group = sim_id,
      colour = "seedling",
      shape = "seedling"
    ),
    size = 3
    # width = 0.1,
    # position=position_dodge(width=2)
  ) +
  geom_point(
    mapping = aes(
      y = rescue_adult,
      colour = "adult",
      shape = "adult",
      group = sim_id
    ),
    size = 3
    # width = 0.1,
    # position=position_dodge(width=2)
  ) +
  geom_point(
    mapping = aes(
      y = rescue_source,
      # ymin = rescue_source_ymin,
      # ymax = rescue_source_ymax,
      # group = sim_id,
      colour = "source",
      shape = "source",
      group = sim_id
    ),
    size = 3
  ) +
  theme_classic() +
  ggtitle("(b)") +
  xlab("Scenario") +
  ylab("Rescue effects (%)") +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_color_manual(
    name = "Plant stage:",
    breaks = c("seed", "seedling", "adult", "source"),
    values = c("#999999", "#009E73", "#56B4E9", "#E69F00")
  ) +
  scale_shape_manual(
    name = "Plant stage:",
    breaks = c("seed", "seedling", "adult", "source"),
    values = c(2, 3, 4, 5)
  ) +
  theme(
    text         = element_text(size = 12, face = "bold"),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

###########################################################################
### subplot (c): rescue sink-to-source ####################################
###########################################################################

# simulation output type
g_sim_type <- "*_plants_per_pop.csv"

# read simulation experiment output files and create df
g_df <- fun_read_files(base::paste0(g_path, "/raw"), g_sim_type)

g_df_metapop <- plyr::ldply(
  base::paste0(g_path, "/in/metapopulation/", g_params$metapop_file),
  utils::read.table,
  sep = " ",
  fill = TRUE,
  header = TRUE
)

# classification of patches (dunes) into sources and sinks
g_df <- g_df %>% 
  dplyr::mutate(
    is_source = dplyr::if_else(
      pop_id %in% g_source_patches,
      "source",
      "sink"
    )
  )

g_df <- g_df %>% 
  dplyr::filter(is_source == "source")

g_df <- g_df %>% 
  dplyr::mutate(
    is_sink_to_source = dplyr::if_else(
      pop_source_id %in% g_source_patches,
      "no",
      "yes"
    )
  )

g_df <- g_df %>% 
  dplyr::filter(is_sink_to_source == "yes")

g_df_summ <- g_df %>% 
  dplyr::group_by(sim_id, study_num) %>% 
  dplyr::summarise(
    sum_mean_germination = base::sum(germination_mean),
    sum_sd_germination = base::sqrt(base::sum(germination_sd^2)),
    sum_mean_seedling = base::sum(seedling_mean),
    sum_sd_seedling = base::sqrt(base::sum(seedling_sd^2)),
    sum_mean_adult = base::sum(adult_mean),
    sum_sd_adult = base::sqrt(base::sum(adult_sd^2)),
    sum_mean_source = base::sum(source_mean),
    sum_sd_source = base::sqrt(base::sum(source_sd^2))
  )


g_df_summ <- g_df_summ %>%
  dplyr::group_by(sim_id) %>%
  dplyr::summarise(
    num_study = dplyr::n(),
    grand_mean_germination = base::mean(sum_mean_germination, na.rm = TRUE),
    grand_sd_germination = fun_pooled_sd_equal_sample_sizes(sum_sd_germination),
    grand_mean_seedling = base::mean(sum_mean_seedling, na.rm = TRUE),
    grand_sd_seedling = fun_pooled_sd_equal_sample_sizes(sum_sd_seedling),
    grand_mean_adult = base::mean(sum_mean_adult, na.rm = TRUE),
    grand_sd_adult = fun_pooled_sd_equal_sample_sizes(sum_sd_adult),
    grand_mean_source = base::mean(sum_mean_source, na.rm = TRUE),
    grand_sd_source = fun_pooled_sd_equal_sample_sizes(sum_sd_source)
  )

# display factor levels for region
g_df_summ$sim_id <- base::as.factor(g_df_summ$sim_id)
base::levels(g_df$sim_id) <- base::c("baseline", "current")

# plot
pc <- ggplot(data = g_df_summ) +
    geom_errorbar(
    mapping = aes(
      x = sim_id,
      ymin = dplyr::if_else(
        (grand_mean_source - grand_sd_source) < 0,
        0, 
        grand_mean_source - grand_sd_source
      ),
      ymax = grand_mean_source + grand_sd_source
    ),
    width = 0.2
  ) +

  geom_point(
    mapping = aes(
      x = sim_id,
      y = grand_mean_source
    ),
    size = 3
  ) +
  ggrepel::geom_label_repel(
    mapping = aes(
      label = base::round(grand_mean_source, 0),
      x = sim_id,
      y = grand_mean_source,
    ),
    hjust = 0,
    size = 3.5, 
    segment.size = 0.05,
    nudge_x = 0.4,
    nudge_y = 300,
    direction = "both",
    label.size = NA
  ) +
  # scale_y_continuous(limits = c(0, 3000)) +
  theme_classic() +
  ggtitle("(c)") +
  xlab("Scenario") +
  ylab("Rescue sink-to-source (# plants)") +
  scale_x_discrete(labels = c("baseline", "current")) +
  theme(
    text         = element_text(size = 12, face = "bold"),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )


# combine subplots
parr <- ggpubr::ggarrange(pa, pb, ncol = 2, common.legend = TRUE, legend = "bottom")
parr2 <- ggpubr::ggarrange(pc, pd, ncol = 2, common.legend = TRUE, legend = "bottom")
parr3 <- ggpubr::ggarrange(parr, parr2, ncol = 1)

# save figure
ggplot2::ggsave(
  plot = parr3,
  filename = base::paste0(g_figure_name, "_dpi300.tiff"),
  path = base::paste0(g_path, "/plots"),
  width = 3.31 * 2,
  height = 3.31 * 2,
  bg = "white",
  dpi = 300
)
# 
# ggplot2::ggsave(
#   plot = parr3,
#   filename = base::paste0(g_figure_name, "_dpi100.png"),
#   path = base::paste0(g_path, "/plots"),
#   width = 3.31 * 2,
#   height = 3.31 * 2,
#   bg = "white",
#   dpi = 100
# )
# 
# ggplot2::ggsave(
#   plot = parr3,
#   filename = base::paste0(g_figure_name, "_dpi1200.png"),
#   path = base::paste0(g_path, "/plots"),
#   width = 3.31 * 2,
#   height = 3.31 * 2,
#   bg = "white",
#   dpi = 1200
# )
# 
# ggplot2::ggsave(
#   plot = parr3,
#   filename = base::paste0(g_figure_name, "_dpi100_by_4.png"),
#   path = base::paste0(g_path, "/plots"),
#   width = 3.31 * 2,
#   height = 4 * 2,
#   bg = "white",
#   dpi = 100
# )

