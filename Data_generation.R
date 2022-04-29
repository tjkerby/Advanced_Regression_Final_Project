library(tidyverse)
library(dyngen)

set.seed(42)

backbone <- backbone_bifurcating()

config <-
  initialise_model(
    backbone = backbone,
    num_cells = 1000,
    num_tfs = nrow(backbone$module_info),
    num_targets = 1000,
    num_hks = 4000,
    simulation_params = simulation_default(
      census_interval = 10, 
      ssa_algorithm = ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(num_simulations = 100),
      compute_cellwise_grn = TRUE
    )
  )

model_common <-
  config %>%
  generate_tf_network() %>%
  generate_feature_network() %>% 
  generate_kinetics() %>%
  generate_gold_standard()


plot_backbone_modulenet(model_common)

#####################
##### KO for B2 #####
#####################

b2_genes <- model_common$feature_info %>% filter(module_id == "B2") %>% pull(feature_id)
b2_ko_model <- model_common
b2_ko_model$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 100L,
  timepoint = 0, 
  genes = b2_genes,
  num_genes = length(b2_genes),
  multiplier = 0
)

b2_ko_model <- b2_ko_model %>%
  generate_cells() %>% 
  generate_experiment()

b2_ko_dataset <- as_dyno(b2_ko_model)
b2_ko_data <- list(cell_info = as.data.frame(b2_ko_dataset$cell_info), expression = as.data.frame(as.matrix(b2_ko_dataset$expression)))
b2_ko_data[['data_raw/bulk_grn']] <- as.data.frame(b2_ko_dataset$regulatory_network)
b2_ko_data[['data_raw/cellwise_grn']] <- as.data.frame(b2_ko_dataset$regulatory_network_sc)

plot_gold_mappings(b2_ko_model, do_facet = FALSE)

#####################
##### KO for B3 #####
#####################

b3_genes <- model_common$feature_info %>% filter(module_id == "B3") %>% pull(feature_id)
b3_ko_model <- model_common
b3_ko_model$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 100L,
  timepoint = 0, 
  genes = b3_genes,
  num_genes = length(b3_genes),
  multiplier = 0
)

b3_ko_model <- b3_ko_model %>%
  generate_cells() %>% 
  generate_experiment()

b3_ko_dataset <- as_dyno(b3_ko_model)
b3_ko_data <- list(cell_info = as.data.frame(b3_ko_dataset$cell_info), expression = as.data.frame(as.matrix(b3_ko_dataset$expression)))
b3_ko_data[['data_raw/bulk_grn']] <- as.data.frame(b3_ko_dataset$regulatory_network)
b3_ko_data[['data_raw/cellwise_grn']] <- as.data.frame(b3_ko_dataset$regulatory_network_sc)

plot_gold_mappings(b3_ko_model, do_facet = FALSE)

############################
##### COMBINE DATASETS #####
############################

full_data <- list()

b2_temp <- cbind(b2_ko_data$expression, rep("D", nrow(b2_ko_data$expression)))
colnames(b2_temp)[ncol(b2_temp)] <- "label"
b3_temp <- cbind(b3_ko_data$expression, rep("C", nrow(b3_ko_data$expression)))
colnames(b3_temp)[ncol(b3_temp)] <- "label"

full_data[["expression"]] <- rbind(b2_temp, b3_temp)
full_data[["b2_cellwise_grn"]] <- b2_ko_data$bulk_grn
full_data[["b2_bulk_grn"]] <- b2_ko_data$bulk_grn
full_data[["b2_cell_info"]] <- b2_ko_data$cell_info
full_data[["b3_cellwise_grn"]] <- b3_ko_data$bulk_grn
full_data[["b3_bulk_grn"]] <- b3_ko_data$bulk_grn
full_data[["b3_cell_info"]] <- b3_ko_data$cell_info

lapply(names(full_data), function(df) write.csv(full_data[[df]], file=paste0(df, ".csv")))