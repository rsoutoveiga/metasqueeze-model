###############################################################################
# Sensitivity Analysis Generator for the MetaSqueeze model. Used in the 
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

# load necessary packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, plyr, readr)

# g_args <- base::commandArgs(trailingOnly = TRUE)

# if (base::is.null(g_args) | g_args == "") {
#   base::stop(
#     base::paste0(
#       "no input file name specified.\n",
#       "the current Working Directory is:\n",
#       base::getwd()
#     )
#   )
# }

# g_var_all <- readr::read_file(g_args[1])
g_var_all <- readr::read_file("./sims/sensitivity_02_current_ref.txt")

g_calculate_scenarios <- function(a_param, a_var) {
  
  g_args_experiment <- base::list()
  
  g_var_split <- base::unlist(
    base::strsplit(a_var, split = ",")
  )
  
  for (i in 1:base::length(a_param)) {
    
    g_args_experiment_tmp <- g_var_split[base::length(a_param) + i]
    
    # argument values are range values if they are between square brackets
    g_is_range <- base::grepl("\\[.*?\\]", g_args_experiment_tmp)
    
    # remove parenthesis and square brackets
    g_args_experiment_tmp <- base::gsub('[(]', '', g_args_experiment_tmp)
    g_args_experiment_tmp <- base::gsub('[)]', '', g_args_experiment_tmp)
    g_args_experiment_tmp <- base::gsub('[[]', '', g_args_experiment_tmp)
    g_args_experiment_tmp <- base::gsub('[]]', '', g_args_experiment_tmp)
    
    g_args_experiment_tmp <- base::unlist(
      base::strsplit(g_args_experiment_tmp, split = " ")
    )
    
    g_args_experiment_tmp <- base::gsub(' ', '', g_args_experiment_tmp, fixed = TRUE)  # remove white spaces
    g_args_experiment_tmp <- stringi::stri_remove_empty(g_args_experiment_tmp, na_empty = FALSE)  # remove empty chars
    
    # if arguments are within brackets, then, values are range arguments
    # i.e. [from by to]
    if (g_is_range) {
      if (base::length(g_args_experiment_tmp) == 3) {
        g_args_experiment_tmp <- base::as.numeric(g_args_experiment_tmp)
        g_args_experiment_tmp <- base::seq(from = g_args_experiment_tmp[1], to = g_args_experiment_tmp[2], by = g_args_experiment_tmp[3])
      } else {
        base::stop(
          base::paste0(
            "there must be three range arguments:\n",
            "[from by to]\n e.g. [100 50 1000]"
          )
        )
      }
    }
    g_args_experiment <- base::append(g_args_experiment, base::list(g_args_experiment_tmp))
  }
  g_args_experiment
}

######################
# parameter sweeping
# using same/similar methods to generate simulation inputs as netlogo (Behavespace)

g_sim_id <- 1
g_is_study_rep <- 1
g_study_rep <- "study_rep.csv"
g_sim_rep <- 10
g_run_time <- 500
g_is_study <- 0
g_is_fire_output <- 0
g_is_dynamics <- 0
g_is_metapop <- 0
g_is_metapop_full <- 0
g_is_metapop_dynamics <- 0
g_is_metapop_pops <- 0
g_is_recol_rescue <- 0
g_is_plants_pop <- 0
g_is_dispersal <- 0
g_scenario_file <- "scenario_base.csv"


g_sim_name <- base::unlist(
  stringr::str_extract_all(
    string = g_var_all, 
    pattern = "(?<=\\#).*(?=\\#)"
  )
)

g_sim_name <- base::gsub('[\"]', '', g_sim_name)
g_sim_name <- base::gsub('[)]', '', g_sim_name)
g_sim_name <- base::gsub('[(]', '', g_sim_name)
g_sim_name <- base::gsub('[)]', '', g_sim_name)
g_sim_name <- base::gsub('[[]', '', g_sim_name)
g_sim_name <- base::gsub('[]]', '', g_sim_name)
g_sim_name <- base::gsub(' ', '', g_sim_name, fixed = TRUE)  # remove white spaces
g_sim_name <- stringi::stri_remove_empty(g_sim_name, na_empty = FALSE)  # remove empty chars


g_perc <- base::unlist(
  stringr::str_extract_all(
    string = g_var_all, 
    pattern = "(?<=\\$).*(?=\\$)"
  )
)

g_var_sim <- base::unlist(
  stringr::str_extract_all(
    string = g_var_all, 
    pattern = "(?<=\\@).*(?=\\@)"
  )
)


g_var_fix <- base::unlist(
  stringr::str_extract_all(
    string = g_var_all, 
    pattern = "(?<=\\{).*(?=\\})"
  )
)

g_var_sets <- base::unlist(
  stringr::str_extract_all(
    string = g_var_all, 
    pattern = "(?<=\\<).*(?=\\>)"
  )
)

g_var_perc <- NULL

for (i in g_perc) {
  
  g_var_split <- base::unlist(
    base::strsplit(i, split = " ")
  )
  
  g_var_perc <- c(g_var_perc, g_var_split)
}

g_var_sensitivity <- NULL

for (i in g_var_sets) {
  
  g_var_split <- base::unlist(
    base::strsplit(i, split = " ")
  )
  
  g_var_sensitivity <- c(g_var_sensitivity, g_var_split)
}

# get simulation settings variables
g_var_fix_split <- base::unlist(
  base::strsplit(g_var_fix, split = ",")
)

g_var_fix_split <- base::gsub('[\"]', '', g_var_fix_split)
g_var_fix_split <- base::gsub('[)]', '', g_var_fix_split)
g_var_fix_split <- base::gsub('[(]', '', g_var_fix_split)
g_var_fix_split <- base::gsub('[)]', '', g_var_fix_split)
g_var_fix_split <- base::gsub('[[]', '', g_var_fix_split)
g_var_fix_split <- base::gsub('[]]', '', g_var_fix_split)
g_var_fix_split <- base::gsub(' ', '', g_var_fix_split, fixed = TRUE)  # remove white spaces
g_var_fix_split <- stringi::stri_remove_empty(g_var_fix_split, na_empty = FALSE)  # remove empty chars

for (i in g_var_sim) {
  g_sim_split <- base::unlist(
    base::strsplit(i, split = ",")
  )

    if (g_sim_split[1] == "is_study_rep" | g_sim_split[1] == "\"is_study_rep\"") {
    g_is_study_rep <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "study_rep" | g_sim_split[1] == "\"study_rep\"") {
    g_study_rep <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"sim_rep\"" | g_sim_split[1] == "sim_rep") {
    g_sim_rep <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"run_time\"" | g_sim_split[1] == "run_time") {
    g_run_time <- g_sim_split[2]
    next
  }
  
    if (g_sim_split[1] == "\"is_study\"" | g_sim_split[1] == "is_study") {
    g_is_study <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"is_fire_output\"" | g_sim_split[1] == "is_fire_output") {
    g_is_fire_output <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"is_dynamics\"" | g_sim_split[1] == "is_dynamics") {
    g_is_dynamics <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"is_metapop\"" | g_sim_split[1] == "is_metapop") {
    g_is_metapop <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"is_metapop_full\"" | g_sim_split[1] == "is_metapop_full") {
    g_is_metapop_full <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"is_metapop_dynamics\"" | g_sim_split[1] == "is_metapop_dynamics") {
    g_is_metapop_dynamics <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"is_metapop_pops\"" | g_sim_split[1] == "is_metapop_pops") {
    g_is_metapop_pops <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"is_recol_rescue\"" | g_sim_split[1] == "is_recol_rescue") {
    g_is_recol_rescue <- g_sim_split[2]
    next
  }

  if (g_sim_split[1] == "\"is_plants_pop\"" | g_sim_split[1] == "is_plants_pop") {
    g_is_plants_pop <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"is_dispersal\"" | g_sim_split[1] == "is_dispersal") {
    g_is_plants_pop <- g_sim_split[2]
    next
  }
  
  if (g_sim_split[1] == "\"scenario_file\"" | g_sim_split[1] == "scenario_file") {
    g_scenario_file <- g_sim_split[2]
    next
  }
}

# this vector will be the same for all simulations in the experiment
g_sim_settings <- base::c(
  g_sim_id,
  g_is_study_rep,
  g_study_rep,
  g_sim_rep,
  g_run_time,
  g_is_study,
  g_is_fire_output,
  g_is_dynamics,
  g_is_metapop,
  g_is_metapop_full,
  g_is_metapop_dynamics,
  g_is_metapop_pops,
  g_is_recol_rescue,
  g_is_plants_pop,
  g_is_dispersal,
  g_scenario_file
)

# remove parenthesis and square brackets
g_sim_settings <- base::gsub('[(]', '', g_sim_settings)
g_sim_settings <- base::gsub('[)]', '', g_sim_settings)
g_sim_settings <- base::gsub('[[]', '', g_sim_settings)
g_sim_settings <- base::gsub('[]]', '', g_sim_settings)
g_sim_settings <- base::gsub(' ', '', g_sim_settings, fixed = TRUE)  # remove white spaces
g_sim_settings <- stringi::stri_remove_empty(g_sim_settings, na_empty = FALSE)  # remove empty chars

g_sim_settings_with_fix <- base::append(g_sim_settings, g_var_fix_split)

# created vector with 5 characters
g_columns <- base::c(
  "sim_id",
  "is_study_rep",
  "study_rep",
  "sim_rep",
  "run_time",
  "is_study",
  "is_fire_output",
  "is_dynamics", 
  "is_metapop",
  "is_metapop_full",
  "is_metapop_dynamics",
  "is_metapop_pops",
  "is_recol_rescue",
  "is_plants_pop",
  "is_dispersal",
  "scenario_file"
)

# input files
g_path <- "../../data/in/"

# read simulation parameters
g_df_scenario <- plyr::ldply(
  base::paste0(g_path, "scenarios/",g_sim_settings[16]),
  utils::read.table,
  sep = " ",
  fill = TRUE,
  header = TRUE
)

g_df_scenario$fire_file[1]

g_df_fire <- plyr::ldply(
  base::paste0(g_path, "fire/", g_df_scenario$fire_file[1]),
  utils::read.table,
  sep = " ",
  fill = TRUE,
  header = TRUE
)

g_df_dispersal <- plyr::ldply(
  base::paste0(g_path, "dispersal/", g_df_scenario$dispersal_file[1]),
  utils::read.table,
  sep = " ",
  fill = TRUE,
  header = TRUE
)

g_df_species <- plyr::ldply(
  base::paste0(g_path, "species/", g_df_scenario$species_file[1]),
  utils::read.table,
  sep = " ",
  fill = TRUE,
  header = TRUE
)

# combine the four files: scenario, fire, dispersal, species
g_df_all <- base::cbind(
  g_df_scenario,
  g_df_fire,
  g_df_dispersal,
  g_df_species
)

# name columns of changed variables
g_max_param <- 1 + base::length(g_var_fix)

g_changed_param <- base::c()

for (i in 1:g_max_param) {
  g_changed_param <- base::append(
    g_changed_param,
    base::c(base::paste0("param", i), base::paste0("val", i))
  )
}

# add extra columns names for the changed parameters (fix and var_sets)
# this is important for parallel!
g_columns <- base::append(
  g_columns,
  g_changed_param
)

g_sim_rows <- base::list()
g_sim_rows <- base::append(g_sim_rows, base::list(g_columns))


# save simulations from negative to positive (e.g. -10%, -5%, +5%, +10%)
g_all_changed <- NULL

g_var_perc <- base::as.numeric(g_var_perc)

for (i in g_var_sensitivity) {
  g_base_value <- g_df_all[1, i]

  # negative percentage
  g_var_perc <- base::sort(g_var_perc, decreasing = TRUE)
  
  for (j in g_var_perc) {
    g_all_changed <- base::append(
      g_all_changed, 
      base::paste(i, g_base_value * (100 - j) / 100)
    )
  }
  
  # positive percentage
  g_var_perc <- base::sort(g_var_perc, decreasing = FALSE)
  
  for (j in g_var_perc) {
    g_all_changed <- base::append(
      g_all_changed, 
      base::paste(i, g_base_value * (100 + j) / 100)
    )
  }
}

for (i in g_all_changed) {
  
  g_split <- base::unlist(
    base::strsplit(i, split = " ")
  )
  g_sim_rows <- base::append(
    g_sim_rows, base::list(base::append(g_sim_settings_with_fix, g_split))
  )
}

# update simulation number at the end of the lists
for (i in 2:base::length(g_sim_rows)) {
  g_sim_rows[[i]][1] <- (i - 1)
}

# save simulations in file
for (i in 1:base::length(g_sim_rows)) {
  df <- base::t(base::as.data.frame(g_sim_rows[[i]]))
  utils::write.table(
    df,
    file = base::paste0("./sims/generated_sims/", g_sim_name, ".csv"),
    quote = FALSE,
    append = TRUE,
    col.names = FALSE,
    row.names = FALSE
  )
}