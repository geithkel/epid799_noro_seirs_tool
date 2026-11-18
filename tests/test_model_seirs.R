# Basic smoke tests for the norovirus SEIRS model.
# Run from the repository root with: Rscript tests/test_model_seirs.R

source("tool_code.R")

results <- run_seirs(
  population = 1000,
  initial_exposed = 5,
  initial_infectious = 5,
  initial_recovered = 0,
  r0 = 2,
  incubation_period = 1.2,
  infectious_period = 2,
  immunity_duration = 180,
  simulation_days = 30
)

total <- results$S + results$E + results$I + results$R
stopifnot(max(abs(total - 1000)) < 1e-4)
stopifnot(all(results$S >= -1e-6))
stopifnot(all(results$E >= -1e-6))
stopifnot(all(results$I >= -1e-6))
stopifnot(all(results$R >= -1e-6))

message("All SEIRS model smoke tests passed.")
