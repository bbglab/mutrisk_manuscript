# these functions are to model HSPC expansion directly as a population instead of different lineages
library(progress)
source("code/0_functions/analysis_variables.R")

# function to model expansions (grow a cell over a specific time interval)
grow_pop = function(pop, division_rate, dt, s = 1 ) {
  pop_division_rate = pop * division_rate * dt
  n_div <- rpois(1, pop_division_rate)
  d  = 1 + 1 + s # denominator
  probs = c(s/d, 1/d, 1/d)
  if (n_div > 0) {
    # For each division event, choose an outcome:
    # "expansion" gives +1 cell, "no_expansion" gives 0, "death" gives -1.
    outcomes <- sample(c(1, 0, -1),
                       size = n_div,
                       replace = TRUE,
                       prob = probs)

    # Update population: no_expansion events do not change the count.
    pop <- pop + sum(outcomes)
  }
  pop[pop < 0] = 0 # reset population to 0 every time there is a 'negative event'
  pop
}


sim_mut_expansion = function(years = 100, ncells, division_rate = 1.3, mut_rate, mut_start_rate, s = 1 , dt = 1 ) {

  n_steps <- ceiling(years / dt)

  pop = ncells
  pop_history = numeric(n_steps)
  pop_history[1] = pop

  mut_rate_dt = mut_rate * dt
  start_rate_pop = pop * mut_start_rate
  mut_cells = rpois(1, start_rate_pop)  # model the starting chance of having a mutation from birth
  pop = pop - mut_cells

  mut_cells_time = vector("numeric", n_steps)
  mut_cells_time[1] = mut_cells
  for (step in 2:n_steps) {
    # model the risk for R882H mutation rates
    mrate_pop = mut_rate_dt * pop
    mut = rpois(1, mrate_pop)
    pop = pop - mut # remove the mutated cell from the normal population
    mut_cells = mut_cells +  mut

    if (mut_cells > 0) {
        mut_cells = grow_pop(mut_cells, division_rate, dt, s)
    }
    mut_cells_time[step] = mut_cells
  }

  if (sum(mut_cells_time > 0)) {
    return(mut_cells_time)
  } else { return(NULL)}
}

# model mutation rates in large cohorts with low mutation rates:
# in this case the blood
mut_start_rate = 9.705839e-08
mut_rate = 1.585985e-08
years = 100
dt = 0.1
n_year_bins = 100/dt
n_simulations = 50

grid = expand.grid(c(1e5, 1e6, 1e7, 1e8, 1e9), c(1, 1.5, 2, 5,10, 20))

grid_list = list()
for (i in 1:nrow(grid)) {

  print(grid[i,])
  ncells = grid[i, 1]
  division_rate = grid[i, 2]
  plot_list = list()

  plot_list = as.list(lapply(1:n_simulations, \(x) {
    muts = sim_mut_expansion(ncells = ncells, mut_rate = mut_rate, mut_start_rate =  mut_start_rate, s = 1, dt = dt, division_rate = division_rate)
    if (!is.null(muts)) {
      m = matrix(unlist(muts), ncol = length(unlist(muts))/n_year_bins)
      rowSums(m)}
    }))

  plot_list = plot_list[!sapply(plot_list, is.null)]
  if (length(plot_list) > 0) {
    names(plot_list) = paste0("rep", 1:length(plot_list))
    plot_list  = as.data.frame(plot_list) |>
      mutate(age_seq = seq(dt, 100, dt)) |>
      pivot_longer(-age_seq) |>
      mutate(ncells = ncells,
             division_rate = division_rate,
             fixed_muts = age_seq * (ncells * mut_rate) + mut_start_rate * ncells)
  }

 grid_list[[i]] = plot_list
}


df = data.table::rbindlist(grid_list) |>
  mutate(division_rate_strip = paste("yearly division rate:", division_rate),
         ncells_strip = paste("# cells modeled:", ncells)) |>
  mutate(division_rate_strip = fct_reorder(division_rate_strip, as.numeric(str_extract(division_rate_strip, "\\d+"))))


df_line = df |>
  select(age_seq, ncells_strip, division_rate_strip, fixed_muts) |>
  distinct()


labs = labs(y = "number of mutated cells", x = "Age (years)",
              subtitle = paste0("n simulations: ", n_simulations,
                                "\nmut rate birth:",format(mut_start_rate, digits = 2),
                                "\nyearly mut_rate: ", format(mut_rate, digits =2)))


plot = df |>
  ggplot(aes(x = age_seq)) +
  geom_line(aes(y = value, group = name), color = "black", alpha = 0.1) +
  geom_line(data = df_line, aes(y = fixed_muts), color = "black") +
  cowplot::theme_cowplot() +
  ggh4x::facet_grid2(ncells_strip ~ division_rate_strip, scales = "free", independent = "all") +
  scale_y_continuous(labels = scales::comma) + labs

plot_log = plot +
  scale_y_log10()

output_dir = "manuscript/Supplementary_notes/Supplementary_Note_5/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
ggsave("manuscript/Supplementary_notes/Supplementary_Note_5/SN5_Figure_2.png", plot, width = 18, height = 15, bg = "white")
ggsave("manuscript/Supplementary_notes/Supplementary_Note_5/SN5_Figure_3.png", plot_log, width = 18, height = 15, bg = "white")

# further exploration: Check if this is worth the effort for the manuscript
# go for a bigger run: model 200.000 individuals
n_simulations = 2e5
grid = expand.grid(c(1e5), c(1.3))
clone_expansion_list = list()
for (i in 1:nrow(grid)) {

  print(grid[i,])
  ncells = grid[i, 1]
  division_rate = grid[i, 2]

  sim_muts = lapply(1:n_year_bins, \(x) list())
  pb <- progress_bar$new(total = n_simulations)
  for (j in 1:n_simulations) {
    muts = sim_mut_expansion(ncells = ncells, mut_rate = mut_rate, mut_start_rate =  mut_start_rate, s = 1, dt = dt, division_rate = division_rate)

    if (!is.null(muts)) {

      for (index in which(muts > 0)) {
        sim_muts[[index]] = c(sim_muts[[index]], muts[index])
      }
    }
    pb$tick()
  }

}

list_mut_clones = lapply(sim_muts, unlist)
mean_fraction_mutated = sapply(list_mut_clones, \(x) mean(c(x, rep(0, n_simulations - length(x)))))
df_200k_sims = data.frame(dynamic = mean_fraction_mutated) |>
  mutate(age_seq = seq(dt, 100, dt),
         fixed = age_seq * (ncells * mut_rate) + mut_start_rate * ncells)

# fraction of individuals with the mutation
plot_2e5_pats  = df_200k_sims |>
  pivot_longer(cols = c(dynamic, fixed), names_to = "model") |>
  ggplot(aes(x = age_seq)) +
  geom_line(aes(y = value, color = model)) +
  theme_cowplot() +
  theme(legend.position = "inside", legend.position.inside = c(0.07, 0.8)) +
  labs(x = 'age (Years)', y = 'fraction of cells with mutation', subtitle = "Simulation of DNMT3A mtuations",
       title = "simulation of DNMT3A R882H\nin 200,000 individuals")
plot_2e5_pats
ggsave("manuscript/Supplementary_notes/Supplementary_Note_5/SN5_Figure_4.png", plot_2e5_pats, width = 7, height = 4, bg = "white")
# todo: model the mutation rates with different initial mutation rates
