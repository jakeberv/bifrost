# Read in example simulation data
sim_share <- readRDS(file='sim_share.RDS')


# Start the R profiler
Rprof(memory.profiling = TRUE, event = "elapsed")
# Run the model with simulated data
searchOptimalConfiguration(
  baseline_tree = paintSubTree(
    reorder(as.phylo(sim_share[[1]]$paintedTree), order = 'postorder'),
    node = length(sim_share[[1]]$paintedTree$tip.label) + 1,
    state = 0
  ),
  trait_data = sim_share[[1]]$simulatedData,
  formula = 'trait_data ~ 1',
  min_descendant_tips = 10,
  num_cores = 8,
  shift_acceptance_threshold = 20,
  plot = T,
  IC = 'GIC',
  store_model_fit_history = FALSE,
  method = 'LL'
)
# Stop the profiler
Rprof(NULL)
# Look at the results
summaryRprof(memory = "both")
mem_res <- summaryRprof(memory = "both")
write.csv(x = as_tibble(x = mem_res$by.total, .name_repair = "universal"), file = "mem_profile_total.txt", row.names = FALSE, quote = FALSE)
write.csv(x = mem_res$by.self, file = "mem_profile_self.txt", quote = FALSE)

# Read in a more nicely formatted memory profile and view it
mem_res_self <- read_csv("mem_profile_self.txt")
mem_res_self %>% View()
