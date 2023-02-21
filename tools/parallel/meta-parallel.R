require(foreach)
require(doParallel)

g_args <- base::commandArgs(trailingOnly = TRUE)

#my_working_directory <- "../../"
my_working_directory <- "~/metasqueeze/"

base::message(
  base::paste0(
    "the current Working Directory is:/n",
    base::getwd()
    )
)

# if (base::is.null(g_args) | g_args == "") {
# if (g_args == "") {
if (base::is.null(g_args)) {
  base::stop(
    "no arguments specified:/n
    argument 1: 'build type' (choose between: 'debug' or 'release')/n
    argument 2: 'number of cores'/n
    argument 3: 'name of the input folder' (located in ~/data/in/sim/
    argument 4: 'name of the output folder' (optional)"
  )
}
# read arguments
build_type    <- g_args[1]    # build type: debug or release
num_cores     <- g_args[2]    # number of cores
input_folder  <- g_args[3]    # input folder name must be in "./data/in"
output_folder <- g_args[4]    # this argument is optional

#build_type <- "release"
#num_cores <- 20
#input_folder <- "test" 

if (base::is.na(build_type)) {
  base::stop(
    "no arguments specified:/n
    argument 1: 'build type' (choose between: 'debug' or 'release')/n
    argument 2: 'number of cores'/n
    argument 3: 'name of the input folder' (located in ~/data/in/sim/)
    argument 4: 'name of the output folder' (optional)"
  )
}


# build_type    <- "debug"    # build type: debug or release
# num_cores     <- 1    # number of cores
# input_folder  <- "default"    # input folder name must be in "./data/in"
# output_folder <- NA    # this argument is optional


if (build_type != "debug" & build_type != "release") {
  base::stop(
    base::paste0(
      "build type is not valid:\n",
      "  1. debug\n",
      "  2. release\n"
    )
  )
}

num_cores <- base::as.numeric(num_cores)

if (base::is.numeric(num_cores)) {
  num_cores <- base::round(num_cores)
} else {
  base::stop("num cores is not numeric")
}

max_cores <- parallel::detectCores() - 1

if (num_cores > max_cores) {
  num_cores <- max_cores
  base::warning(
    base::paste0(
      "num_cores argument is too high.\n",
      "The num_cores is set to ",
      num_cores, " cores"
    )
  )
}

if (base::is.na(output_folder)) {
  base::message(
    base::paste0(
      "output_folder argument is NA.\n",
      "output_folder name is equal to input_folder name"
    )
  )
  output_folder <- input_folder
}

# set name of output folder in ~/data/out
output_time <- format(Sys.time(), "%Y_%m_%d_%H%M%S")
output_folder <- base::paste0(output_time, "_", output_folder)
base::remove(output_time)

# create folder output
tmp_folder <- base::paste0(
  my_working_directory,
  "data/out/",
  output_folder
)

if (base::dir.exists(base::paste0(tmp_folder, "/in/sim/"))) {
  base::stop("The output folder already exists")
} else {
  base::dir.create(base::paste0(tmp_folder, "/in/sim/"), recursive = TRUE)
}

my_split_sim_files <- base::paste0(
  my_working_directory, 
  "data/in/sim/", 
  input_folder, 
  "/split/"
) 

# create the 'split' folder in input_folder. Below each row of the simulation
# file will be save as independent .csv simulation file
if (base::dir.exists(my_split_sim_files)) {
  base::warning(
    "The sim folder in output folder already exists. All .csv files were deleted"
  )
  base::unlink(base::paste0(my_split_sim_files, "*.csv"))
} else {
  base::dir.create(my_split_sim_files, recursive = TRUE)
}

# create sim folder in output folder
# tmp_folder_sim_split <- base::paste0(
#   tmp_folder,
#   "/in/sim/split/"
# )


# read all simulation files (normally, there should be only one file, but
# there can be more)
my_files <- base::list.files(
  path = base::paste0(my_working_directory, "data/in/sim/", input_folder),
  pattern = "*.csv",
  full.names = TRUE
)

# merge all simulation input files in one dataframe
sim_params <- plyr::ldply(
  my_files,
  utils::read.table,
  sep = " ",
  fill = TRUE,
  header = TRUE
)

# rewrite the column 'sim_id" to avoid duplicates
sim_params$sim_id <- 1:nrow(sim_params)

# set to empty all cells as NA. 
# These are for the extra parameters included at the end of the simulation row
sim_params[base::is.na(sim_params)] <- "" 

# save simulation input file in out folder
# this is to have all the parameters used with the simulation results
utils::write.table(
  x = sim_params, 
  file = base::paste0(tmp_folder, "/in/sim/", output_folder, "_simfile", ".csv"), 
  sep = " ",
  row.names = FALSE,
  quote = FALSE
)

# save .csv file for each row of dataframe (i.e. split input files)
# these files will be submitted in foreach using parallel, that is, 
# each split file (i.e. row) is assigned to one core
for(i in 1:nrow(sim_params)) {
  utils::write.table(
    sim_params[i, ], 
    base::paste0(my_split_sim_files, input_folder, "_split_", i, ".csv"),
    sep = " ",
    row.names = FALSE,
    quote = FALSE
  )
}

# # copy simulation files in sim folder
# sim_files_full <- base::list.files(
#   path = base::paste0(my_working_directory, "data/in/sim/", input_folder),
#   pattern = "*",
#   full.names = TRUE
# )
# 
# for (my_file in sim_files_full) {
#   base::file.copy(my_file, tmp_folder_sim)
# }

# base::remove(sim_files_full)

# copy the rest of input files
# this save all input simulation files with the simulation results.
# this might be lead too heave folders in the case there are too many simulation
# files.
# NOTE the "sim" file is not included
# TODO find the way to copy only those files that are used in the simulations
my_folders_in <- c(
  "climate",
  "dispersal",
  "fire",
  "flower_distributions",
  "fuzzy_sets",
  "habitat_quality",
  "metapopulation",
  "plant_types",
  "scenarios",
  "species",
  "study_rep",
  "study_size"
)

for (my_folder in my_folders_in) {
  base::dir.create(
    base::paste0(
      tmp_folder,
      "/in/",
      my_folder
    )
  )
  
  base::file.copy(
    from = base::paste0(my_working_directory, "data/in/", my_folder),
    to = base::paste0(tmp_folder, "/in/"),
    recursive = TRUE
  )
}

# base::remove(tmp_folder, tmp_folder_in, max_cores)

# read split sim files to use in foreach parallel
sim_files <- base::list.files(
  path = my_split_sim_files,
  pattern = "*.csv"
)

base::print("LIST OF SPLIT FILES: ")
base::print(sim_files)
##############################

 
# this is the number of simulation files
MC_max <- base::length(sim_files)

# create the cluster
my_cluster <- parallel::makeCluster(num_cores)

# check cluster definition (optional)
base::message(my_cluster)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my_cluster)

# check if it is registered (optional)
foreach::getDoParRegistered()

input_folder <- base::paste0(input_folder, "/split")

base::print("MY_WORKING_DIRECTORY: ")
base::print(my_working_directory)
base::print("INPUT FOLDER")
base::print(input_folder)
base::print("OUTPUT FOLDER")
base::print(output_folder)
base::print("MAX NUMBER OF CORES")
base::print(MC_max)
base::print("NUMBER OF CORES TAKEN")
base::print(num_cores)

# start the simulations in parallel
foreach::foreach(MC = 1:MC_max) %dopar% {
  # this should be equal to what you type in the command line to start your program.
  # MC is my counter for the repetition
  base::system(
    base::paste(
      base::paste0(
	my_working_directory,
       	"build/",
       	build_type,
       	"/metapop"),
      input_folder,
      sim_files[MC],
      output_folder,
      sep = " "
    ),
    intern=T,
    show.output.on.console = F
  )
}

# at the end stop the cluster
parallel::stopCluster(cl = my_cluster)
  
