
# This script runs all simulations for different types of interventions

# The following code remains the same for all vaccination scenarios - please run 
# these first

# Libraries required for graphs
library(tidyverse)
library(cowplot)
library(incidence2) 
source("CI_prop_binomial.R")

# Load in the posterior distribution
sampled_post <- read.csv("Posterior_Distrib_MVD.csv")

# Check posterior distribution of R values
sampled_post %>%
  tibble() %>%
  mutate(R_before = R_basic,
         R_after = R_before * (1 - intervention_efficacy)) %>% 
  ggplot(aes(x = R_before, y = R_after)) + 
  geom_density_2d_filled()

## Define parameters that stay constant for all simulations
n_sims <-  1000

# Specify the serial interval mean and s.d.
SI_mean <- 9.2
SI_sd <- 4.4

# Number of days simulated
max_duration <- 365

# Time required for vaccine efficacy to peak
vaccination_time_to_max <- 7

# Specify the maximum vaccination delay (equal to the length of one cluster)
max_delay <- 30

# Break point (this is the time at which interventions start)
break_point <- 90

#############################################################################
## Simulate the no vaccination (control) scenario 

source("branching_basics_no_vaccination.R")

sim_params <- expand.grid(
  r_intro = c(0.002, 0.06)
) %>% mutate(label = paste(
  "E:0",
  "/",
  "C:0",
  sep = ""
))

list_results <- list()

for (j in seq_len(nrow(sim_params))) {
  n_case <- integer(n_sims)
  #Also record number of intros
  n_intros <- integer(n_sims)
  controlled <- logical(n_sims)
  R_before <- numeric(n_sims)
  R_after <- numeric(n_sims)
  last_force_infection <- numeric(n_sims)
  R_basic <- current_R_basic <- numeric(n_sims)
  intervention_efficacy <- numeric(n_sims)
  last_intro <- integer(n_sims)
  break_point_used <- integer(n_sims)
  
  for (i in 1:n_sims) {
    # select posterior sample to use
    current_idx <- sample.int(nrow(sampled_post), 1)
    current_R_basic <- sampled_post$R_basic[current_idx]
    current_intervention_efficacy <- sampled_post$intervention_efficacy[current_idx]
    
    out <- branching_process_model_no_vacc(
      R_basic = current_R_basic,
      intervention_efficacy = current_intervention_efficacy,
      serial_interval_mean = SI_mean,
      serial_interval_sd = SI_sd,
      r_daily_intro = sim_params$r_intro[j],
      max_duration = max_duration,
      break_point = break_point
    )
    
    # extract summaries from the simulation:
    # number of cases, whether the outbreak was controlled, and mean values of R
    # before and after intervention
    n_case[i] <- nrow(out) 
    controlled[i] <- attr(out, "controlled")
    n_intros[i] <- attr(out, "n_intros")
    break_point_used[i] <- attr(out, "break_point")
    before <- out$date_onset <= break_point
    R_before[i] <- mean(out$R[before])
    R_after[i] <- mean(out$R[!before])
    last_force_infection[i] <- attr(out, "last_force_infection")
    R_basic[i] <- current_R_basic
    intervention_efficacy[i] <- current_intervention_efficacy
    last_intro[i] <- out %>% 
      filter(is_intro) %>% 
      pull(date_onset) %>% 
      max()
    
    gc()
  }
  
  list_results[[j]] <- tibble(
    # keep track of parameter values
    label = sim_params$label[j],
    r_intro = sim_params$r_intro[j],
    n_case,
    n_intros,
    break_point_used,
    R_before,
    R_after, 
    R_basic,
    intervention_efficacy, 
    controlled,
    last_intro,
    last_force_infection
  )
}

# make a tibble for this scenario
res_no_vaccine <- bind_rows(list_results) %>%
  mutate(scenario = "ring") %>%
  select(scenario, everything())

######################################################################

## Simulate Ring Vaccination
source("branching_basics_ring.R")

## These parameters are relevant to ring/reactive

# Fraction of cases that are reported after intervention
reporting <- 0.9

# The amount of time between vaccination and infection
# For ring vaccination, a reactive strategy, we sample from the slope of
# the logistic curve
mean_delay_vaccination_infection <- 9
sd_delay_vaccination_infection <- 4

## Define the parameters values for different simulations
# sim_params <- expand.grid(
#   max_vaccine_efficacy = c(1, .9, .7),
#   vaccine_coverage = c(.9, .7, .5),
#   r_intro = c(0.002, 0.06)
# ) %>% mutate(label = paste(
#   "E:",
#   round(max_vaccine_efficacy * 100),
#   "/",
#   "C:",
#   round(vaccine_coverage * 100),
#   sep = ""
# ))
sim_params <- expand.grid(
  max_vaccine_efficacy = c(.9),
  vaccine_coverage = c(.9),
  r_intro = c(0.002, 0.06)
) %>% mutate(label = paste(
  "E:",
  round(max_vaccine_efficacy * 100),
  "/",
  "C:",
  round(vaccine_coverage * 100),
  sep = ""
))

# # Include the no vaccine scenario
# NV_low_intros <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.002,
#   label = "E:0/C:0"
# )
# NV_high_intros <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.06,
#   label = "E:0/C:0"
# )
# sim_params <- bind_rows(sim_params, NV_low_intros, NV_high_intros)

list_results <- list()

for (j in seq_len(nrow(sim_params))) {
  n_case <- integer(n_sims)
  #Also record number of intros
  n_intros <- integer(n_sims)
  controlled <- logical(n_sims)
  R_before <- numeric(n_sims)
  R_after <- numeric(n_sims)
  last_force_infection <- numeric(n_sims)
  R_basic <- current_R_basic <- numeric(n_sims)
  intervention_efficacy <- numeric(n_sims)
  last_intro <- integer(n_sims)
  break_point_used <- integer(n_sims)
  
  for (i in 1:n_sims) {
    # select posterior sample to use
    current_idx <- sample.int(nrow(sampled_post), 1)
    current_R_basic <- sampled_post$R_basic[current_idx]
    current_intervention_efficacy <- sampled_post$intervention_efficacy[current_idx]
    
    out <- branching_process_model_ring(
      R_basic = current_R_basic,
      intervention_efficacy = current_intervention_efficacy,
      serial_interval_mean = SI_mean,
      serial_interval_sd = SI_sd,
      r_daily_intro = sim_params$r_intro[j],
      max_duration = max_duration,
      reporting = reporting,
      vaccination_coverage = sim_params$vaccine_coverage[j],
      max_vaccination_efficacy = sim_params$max_vaccine_efficacy[j],
      vaccination_time_to_max = vaccination_time_to_max,
      mean_delay_vaccination_to_infection = mean_delay_vaccination_infection,
      sd_delay_vaccination_to_infection = sd_delay_vaccination_infection,
      max_delay = max_delay,
      break_point = break_point
    )
    
    # extract summaries from the simulation:
    # number of cases, whether the outbreak was controlled, and mean values of R
    # before and after intervention
    n_case[i] <- nrow(out) 
    controlled[i] <- attr(out, "controlled")
    n_intros[i] <- attr(out, "n_intros")
    break_point_used[i] <- attr(out, "break_point")
    before <- out$date_onset <= break_point
    R_before[i] <- mean(out$R[before])
    R_after[i] <- mean(out$R[!before])
    last_force_infection[i] <- attr(out, "last_force_infection")
    R_basic[i] <- current_R_basic
    intervention_efficacy[i] <- current_intervention_efficacy
    last_intro[i] <- out %>% 
      filter(is_intro) %>% 
      pull(date_onset) %>% 
      max()
    
    gc()
  }
  
  list_results[[j]] <- tibble(
    # keep track of parameter values
    label = sim_params$label[j],
    r_intro = sim_params$r_intro[j],
    n_case,
    n_intros,
    break_point_used,
    R_before,
    R_after, 
    R_basic,
    intervention_efficacy, 
    controlled,
    last_intro,
    last_force_infection
  )
}

# make a tibble for this scenario
res_ring <- bind_rows(list_results) %>%
  mutate(scenario = "ring") %>%
  select(scenario, everything())


#########################################################################

# Simulate prophylactic targeted vaccination

source("branching_basics_prophylactic_targetted.R")

## These parameters are only relevant to targeted vaccination

# Proportion of infected cases (historically) who were HCWs
prop_hcw = 0.06

# The amount of time between vaccination and infection
# For ring vaccination, a reactive strategy, we sample from the top of the
# logistic curve
mean_delay_vaccination_infection <- 20
sd_delay_vaccination_infection <- 5

## Define the parameters values for different simulations
# sim_params <- expand.grid(
#   max_vaccine_efficacy = c(1, .9, .7),
#   vaccine_coverage = c(.9, .7, .5),
#   r_intro = c(0.002, 0.06)
# ) %>% mutate(label = paste(
#   "E:",
#   round(max_vaccine_efficacy * 100),
#   "/",
#   "C:",
#   round(vaccine_coverage * 100),
#   sep = ""
# ))

sim_params <- expand.grid(
  max_vaccine_efficacy = c(.9),
  vaccine_coverage = c(.9),
  r_intro = c(0.002, 0.06)
) %>% mutate(label = paste(
  "E:",
  round(max_vaccine_efficacy * 100),
  "/",
  "C:",
  round(vaccine_coverage * 100),
  sep = ""
))

# #Include the no vaccine scenario
# NV_low_intros  <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.002,
#   label = "E:0/C:0"
# )
# NV_high_intros <-  data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.06,
#   label = "E:0/C:0"
# )
# sim_params <- bind_rows(sim_params, NV_low_intros, NV_high_intros)


list_results <- list()

for (j in seq_len(nrow(sim_params))) {
  n_case <- integer(n_sims)
  #Also record number of intros
  n_intros <- integer(n_sims)
  controlled <- logical(n_sims)
  R_before <- numeric(n_sims)
  R_after <- numeric(n_sims)
  last_force_infection <- numeric(n_sims)
  R_basic <- current_R_basic <- numeric(n_sims)
  intervention_efficacy <- numeric(n_sims)
  last_intro <- integer(n_sims)
  break_point_used <- integer(n_sims)
  
  for (i in 1:n_sims) {
    # select posterior sample to use
    current_idx <- sample.int(nrow(sampled_post), 1)
    current_R_basic <- sampled_post$R_basic[current_idx]
    current_intervention_efficacy <- sampled_post$intervention_efficacy[current_idx]
    
    out <- branching_process_model_target(
      R_basic = current_R_basic,
      intervention_efficacy = current_intervention_efficacy,
      serial_interval_mean = SI_mean,
      serial_interval_sd = SI_sd,
      r_daily_intro = sim_params$r_intro[j],
      max_duration = max_duration,
      vaccination_coverage = sim_params$vaccine_coverage[j],
      max_vaccination_efficacy = sim_params$max_vaccine_efficacy[j],
      vaccination_time_to_max = vaccination_time_to_max,
      mean_delay_vaccination_to_infection = mean_delay_vaccination_infection,
      sd_delay_vaccination_to_infection = sd_delay_vaccination_infection,
      max_delay = max_delay,
      prop_hcw = prop_hcw,
      break_point = break_point
    )
    
    # extract summaries from the simulation:
    # number of cases, whether the outbreak was controlled, and mean values of R
    # before and after intervention
    n_case[i] <- nrow(out) 
    controlled[i] <- attr(out, "controlled")
    n_intros[i] <- attr(out, "n_intros")
    break_point_used[i] <- attr(out, "break_point")
    before <- out$date_onset <= break_point
    R_before[i] <- mean(out$R[before])
    R_after[i] <- mean(out$R[!before])
    last_force_infection[i] <- attr(out, "last_force_infection")
    R_basic[i] <- current_R_basic
    intervention_efficacy[i] <- current_intervention_efficacy
    last_intro[i] <- out %>% 
      filter(is_intro) %>% 
      pull(date_onset) %>% 
      max()
    
    gc()
  }
  
  list_results[[j]] <- tibble(
    # keep track of parameter values
    label = sim_params$label[j],
    r_intro = sim_params$r_intro[j],
    n_case,
    n_intros,
    break_point_used,
    R_before,
    R_after, 
    R_basic,
    intervention_efficacy, 
    controlled,
    last_intro,
    last_force_infection
  )
}
  
# make a tibble for this scenario
res_prophylactic_target <- bind_rows(list_results) %>%
  mutate(scenario = "prophylactic targeted") %>%
  select(scenario, everything())

#########################################################################

# Simulate prophylactic mass vaccination

source("branching_basics_prophylactic_mass.R")

## These parameters are only relevant to targeted/prophylactic vaccination

# The amount of time between vaccination and infection
# For ring vaccination, a reactive strategy, we sample from the top of the
# logistic curve
mean_delay_vaccination_infection <- 20
sd_delay_vaccination_infection <- 5

## Define the parameters values for different simulations
# sim_params <- expand.grid(
#   max_vaccine_efficacy = c(1, .9, .7),
#   vaccine_coverage = c(.9, .7, .5),
#   r_intro = c(0.002, 0.06)
# ) %>% mutate(label = paste(
#   "E:",
#   round(max_vaccine_efficacy * 100),
#   "/",
#   "C:",
#   round(vaccine_coverage * 100),
#   sep = ""
# ))

sim_params <- expand.grid(
  max_vaccine_efficacy = c(.9),
  vaccine_coverage = c(.9),
  r_intro = c(0.002, 0.06)
) %>% mutate(label = paste(
  "E:",
  round(max_vaccine_efficacy * 100),
  "/",
  "C:",
  round(vaccine_coverage * 100),
  sep = ""
))

# #Include the no vaccine scenario
# NV_low_intros  <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.002,
#   label = "E:0/C:0"
# )
# NV_high_intros <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.06,
#   label = "E:0/C:0"
# )
# sim_params <- bind_rows(sim_params, NV_low_intros, NV_high_intros)

list_results <- list()

for (j in seq_len(nrow(sim_params))) {
  n_case <- integer(n_sims)
  #Also record number of intros
  n_intros <- integer(n_sims)
  controlled <- logical(n_sims)
  R_before <- numeric(n_sims)
  R_after <- numeric(n_sims)
  last_force_infection <- numeric(n_sims)
  R_basic <- current_R_basic <- numeric(n_sims)
  intervention_efficacy <- numeric(n_sims)
  last_intro <- integer(n_sims)
  break_point_used <- integer(n_sims)
  
  for (i in 1:n_sims) {
    # select posterior sample to use
    current_idx <- sample.int(nrow(sampled_post), 1)
    current_R_basic <- sampled_post$R_basic[current_idx]
    current_intervention_efficacy <- sampled_post$intervention_efficacy[current_idx]

    out <- branching_process_model_mass(
      R_basic = current_R_basic,
      intervention_efficacy = current_intervention_efficacy,
      serial_interval_mean = SI_mean,
      serial_interval_sd = SI_sd,
      r_daily_intro = sim_params$r_intro[j],
      max_duration = max_duration,
      vaccination_coverage = sim_params$vaccine_coverage[j],
      max_vaccination_efficacy = sim_params$max_vaccine_efficacy[j],
      vaccination_time_to_max = vaccination_time_to_max,
      mean_delay_vaccination_to_infection = mean_delay_vaccination_infection,
      sd_delay_vaccination_to_infection = sd_delay_vaccination_infection,
      max_delay = max_delay,
      break_point = break_point
    )
    
    # extract summaries from the simulation:
    # number of cases, whether the outbreak was controlled, and mean values of R
    # before and after intervention
    n_case[i] <- nrow(out) 
    controlled[i] <- attr(out, "controlled")
    n_intros[i] <- attr(out, "n_intros")
    break_point_used[i] <- attr(out, "break_point")
    before <- out$date_onset <= break_point
    R_before[i] <- mean(out$R[before])
    R_after[i] <- mean(out$R[!before])
    last_force_infection[i] <- attr(out, "last_force_infection")
    R_basic[i] <- current_R_basic
    intervention_efficacy[i] <- current_intervention_efficacy
    last_intro[i] <- out %>% 
      filter(is_intro) %>% 
      pull(date_onset) %>% 
      max()
    
    gc()
  }
  
  list_results[[j]] <- tibble(
    # keep track of parameter values
    label = sim_params$label[j],
    r_intro = sim_params$r_intro[j],
    n_case,
    n_intros,
    break_point_used,
    R_before,
    R_after, 
    R_basic,
    intervention_efficacy, 
    controlled,
    last_intro,
    last_force_infection
  )
  
}

# make a tibble for this scenario
res_prophylactic_mass <- bind_rows(list_results) %>%
  mutate(scenario = "prophylactic mass") %>%
  select(scenario, everything())

#################################################################

# Reactive Mass Vaccination
source("branching_basics_reactive_mass.R")

## These parameters are relevant to ring/reactive vaccination schemes

# The amount of time between vaccination and infection
# For ring vaccination, a reactive strategy, we sample from the slope of
# the logistic curve
mean_delay_vaccination_infection <- 9
sd_delay_vaccination_infection <- 4

## Define the parameters values for different simulations
# sim_params <- expand.grid(
#   max_vaccine_efficacy = c(1, .9, .7),
#   vaccine_coverage = c(.9, .7, .5),
#   r_intro = c(0.002, 0.06)
# ) %>% mutate(label = paste(
#   "E:",
#   round(max_vaccine_efficacy * 100),
#   "/",
#   "C:",
#   round(vaccine_coverage * 100),
#   sep = ""
# ))

sim_params <- expand.grid(
  max_vaccine_efficacy = c(.9),
  vaccine_coverage = c(.9),
  r_intro = c(0.002, 0.06)
) %>% mutate(label = paste(
  "E:",
  round(max_vaccine_efficacy * 100),
  "/",
  "C:",
  round(vaccine_coverage * 100),
  sep = ""
))

# #Include the no vaccine scenario
# NV_low_intros  <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.002,
#   label = "E:0/C:0"
# )
# NV_high_intros <-  data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.06,
#   label = "E:0/C:0"
# )
# sim_params <- bind_rows(sim_params, NV_low_intros, NV_high_intros)

list_results <- list()

for (j in seq_len(nrow(sim_params))) {
  n_case <- integer(n_sims)
  #Also record number of intros
  n_intros <- integer(n_sims)
  controlled <- logical(n_sims)
  R_before <- numeric(n_sims)
  R_after <- numeric(n_sims)
  last_force_infection <- numeric(n_sims)
  R_basic <- current_R_basic <- numeric(n_sims)
  intervention_efficacy <- numeric(n_sims)
  last_intro <- integer(n_sims)
  break_point_used <- integer(n_sims)
  
  for (i in 1:n_sims) {
    # select posterior sample to use
    current_idx <- sample.int(nrow(sampled_post), 1)
    current_R_basic <- sampled_post$R_basic[current_idx]
    current_intervention_efficacy <- sampled_post$intervention_efficacy[current_idx]
    
    out <- branching_process_model_reactive_mass(
      R_basic = current_R_basic,
      intervention_efficacy = current_intervention_efficacy,
      serial_interval_mean = SI_mean,
      serial_interval_sd = SI_sd,
      r_daily_intro = sim_params$r_intro[j],
      max_duration = max_duration,
      vaccination_coverage = sim_params$vaccine_coverage[j],
      max_vaccination_efficacy = sim_params$max_vaccine_efficacy[j],
      vaccination_time_to_max = vaccination_time_to_max,
      mean_delay_vaccination_to_infection = mean_delay_vaccination_infection,
      sd_delay_vaccination_to_infection = sd_delay_vaccination_infection,
      max_delay = max_delay,
      break_point = break_point
    )
    
    # extract summaries from the simulation:
    # number of cases, whether the outbreak was controlled, and mean values of R
    # before and after intervention
    n_case[i] <- nrow(out) 
    controlled[i] <- attr(out, "controlled")
    n_intros[i] <- attr(out, "n_intros")
    break_point_used[i] <- attr(out, "break_point")
    before <- out$date_onset <= break_point
    R_before[i] <- mean(out$R[before])
    R_after[i] <- mean(out$R[!before])
    last_force_infection[i] <- attr(out, "last_force_infection")
    R_basic[i] <- current_R_basic
    intervention_efficacy[i] <- current_intervention_efficacy
    last_intro[i] <- out %>% 
      filter(is_intro) %>% 
      pull(date_onset) %>% 
      max()
    
    gc()
  }
  
  list_results[[j]] <- tibble(
    # keep track of parameter values
    label = sim_params$label[j],
    r_intro = sim_params$r_intro[j],
    n_case,
    n_intros,
    break_point_used,
    R_before,
    R_after, 
    R_basic,
    intervention_efficacy, 
    controlled,
    last_intro,
    last_force_infection
  )
  
}
# make a tibble for this scenario
res_reactive_mass <- bind_rows(list_results) %>%
  mutate(scenario = "reactive mass") %>%
  select(scenario, everything())

#################################################################

# Reactive Targeted Vaccination
source("branching_basics_reactive_targetted.R")

## These parameters are relevant to targeted vaccination schemes
prop_hcw <- 0.06

# The amount of time between vaccination and infection
# For ring vaccination, a reactive strategy, we sample from the slope of
# the logistic curve
mean_delay_vaccination_infection <- 9
sd_delay_vaccination_infection <- 4

## Define the parameters values for different simulations
# sim_params <- expand.grid(
#   max_vaccine_efficacy = c(1, .9, .7),
#   vaccine_coverage = c(.9, .7, .5),
#   r_intro = c(0.002, 0.06)
# ) %>% mutate(label = paste(
#   "E:",
#   round(max_vaccine_efficacy * 100),
#   "/",
#   "C:",
#   round(vaccine_coverage * 100),
#   sep = ""
# ))

sim_params <- expand.grid(
  max_vaccine_efficacy = c(.9),
  vaccine_coverage = c(.9),
  r_intro = c(0.002, 0.06)
) %>% mutate(label = paste(
  "E:",
  round(max_vaccine_efficacy * 100),
  "/",
  "C:",
  round(vaccine_coverage * 100),
  sep = ""
))


# #Include the no vaccine scenario
# NV_low_intros  <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.002,
#   label = "E:0/C:0"
# )
# NV_high_intros <-  data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.06,
#   label = "E:0/C:0"
# )
# sim_params <- bind_rows(sim_params, NV_low_intros, NV_high_intros)

list_results <- list()

for (j in seq_len(nrow(sim_params))) {
  n_case <- integer(n_sims)
  #Also record number of intros
  n_intros <- integer(n_sims)
  controlled <- logical(n_sims)
  R_before <- numeric(n_sims)
  R_after <- numeric(n_sims)
  last_force_infection <- numeric(n_sims)
  R_basic <- current_R_basic <- numeric(n_sims)
  intervention_efficacy <- numeric(n_sims)
  last_intro <- integer(n_sims)
  break_point_used <- integer(n_sims)
  
  for (i in 1:n_sims) {
    # select posterior sample to use
    current_idx <- sample.int(nrow(sampled_post), 1)
    current_R_basic <- sampled_post$R_basic[current_idx]
    current_intervention_efficacy <- sampled_post$intervention_efficacy[current_idx]
    
    out <- branching_process_model_reactive_target(
      R_basic = current_R_basic,
      intervention_efficacy = current_intervention_efficacy,
      serial_interval_mean = SI_mean,
      serial_interval_sd = SI_sd,
      r_daily_intro = sim_params$r_intro[j],
      max_duration = max_duration,
      vaccination_coverage = sim_params$vaccine_coverage[j],
      max_vaccination_efficacy = sim_params$max_vaccine_efficacy[j],
      vaccination_time_to_max = vaccination_time_to_max,
      mean_delay_vaccination_to_infection = mean_delay_vaccination_infection,
      sd_delay_vaccination_to_infection = sd_delay_vaccination_infection,
      max_delay = max_delay,
      prop_hcw = prop_hcw,
      break_point = break_point
    )
    
    # extract summaries from the simulation:
    # number of cases, whether the outbreak was controlled, and mean values of R
    # before and after intervention
    n_case[i] <- nrow(out) 
    controlled[i] <- attr(out, "controlled")
    n_intros[i] <- attr(out, "n_intros")
    break_point_used[i] <- attr(out, "break_point")
    before <- out$date_onset <= break_point
    R_before[i] <- mean(out$R[before])
    R_after[i] <- mean(out$R[!before])
    last_force_infection[i] <- attr(out, "last_force_infection")
    R_basic[i] <- current_R_basic
    intervention_efficacy[i] <- current_intervention_efficacy
    last_intro[i] <- out %>% 
      filter(is_intro) %>% 
      pull(date_onset) %>% 
      max()
    
    gc()
  }
  
  list_results[[j]] <- tibble(
    # keep track of parameter values
    label = sim_params$label[j],
    r_intro = sim_params$r_intro[j],
    n_case,
    n_intros,
    break_point_used,
    R_before,
    R_after, 
    R_basic,
    intervention_efficacy, 
    controlled,
    last_intro,
    last_force_infection
  )
  
}
# make a tibble for this scenario
res_reactive_targeted <- bind_rows(list_results) %>%
  mutate(scenario = "reactive targeted") %>%
  select(scenario, everything())


##############################################################################
# Simulate Ring + Reactive Targeted Vaccination

source("branching_basics_reactive_targetted_ring.R")

## These parameters are relevant to ring/reactive

# Proportion of HCWs in all Marburg cases
prop_hcw = 0.06

# Fraction of cases that are reported after intervention
reporting <- 0.9

# The amount of time between vaccination and infection
# For reactive strategies, we sample from the slope of
# the logistic curve
mean_delay_vaccination_to_infection_target <- 9
sd_delay_vaccination_to_infection_target <- 4

mean_delay_vaccination_to_infection_ring <- 9
sd_delay_vaccination_to_infection_ring <- 4

## Define the parameters values for different simulations
# sim_params <- expand.grid(
#   max_vaccine_efficacy = c(1, .9, .7),
#   vaccine_coverage = c(.9, .7, .5),
#   r_intro = c(0.002, 0.06)
# ) %>% mutate(label = paste(
#   "E:",
#   round(max_vaccine_efficacy * 100),
#   "/",
#   "C:",
#   round(vaccine_coverage * 100),
#   sep = ""
# ))

sim_params <- expand.grid(
  max_vaccine_efficacy = c(.9),
  vaccine_coverage = c(.9),
  r_intro = c(0.002, 0.06)
) %>% mutate(label = paste(
  "E:",
  round(max_vaccine_efficacy * 100),
  "/",
  "C:",
  round(vaccine_coverage * 100),
  sep = ""
))

# # Include the no vaccine scenario
# NV_low_intros <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.002,
#   label = "E:0/C:0"
# )
# NV_high_intros <- data.frame(
#   max_vaccine_efficacy = 0,
#   vaccine_coverage = 0,
#   r_intro = 0.06,
#   label = "E:0/C:0"
# )
# sim_params <- bind_rows(sim_params, NV_low_intros, NV_high_intros)

list_results <- list()

for (j in seq_len(nrow(sim_params))) {
  n_case <- integer(n_sims)
  #Also record number of intros
  n_intros <- integer(n_sims)
  controlled <- logical(n_sims)
  R_before <- numeric(n_sims)
  R_after <- numeric(n_sims)
  last_force_infection <- numeric(n_sims)
  R_basic <- current_R_basic <- numeric(n_sims)
  intervention_efficacy <- numeric(n_sims)
  last_intro <- integer(n_sims)
  break_point_used <- integer(n_sims)
  
  for (i in 1:n_sims) {
    # select posterior sample to use
    current_idx <- sample.int(nrow(sampled_post), 1)
    current_R_basic <- sampled_post$R_basic[current_idx]
    current_intervention_efficacy <- sampled_post$intervention_efficacy[current_idx]
    
    out <- branching_process_model_reactive_target_ring(
      R_basic = current_R_basic,
      intervention_efficacy = current_intervention_efficacy,
      serial_interval_mean = SI_mean,
      serial_interval_sd = SI_sd,
      r_daily_intro = sim_params$r_intro[j],
      max_duration = max_duration,
      reporting = reporting,
      vaccination_coverage = sim_params$vaccine_coverage[j],
      max_vaccination_efficacy = sim_params$max_vaccine_efficacy[j],
      vaccination_time_to_max = vaccination_time_to_max,
      mean_delay_vaccination_to_infection_target = mean_delay_vaccination_to_infection_target,
      sd_delay_vaccination_to_infection_target = sd_delay_vaccination_to_infection_target,
      mean_delay_vaccination_to_infection_ring = mean_delay_vaccination_to_infection_ring,
      sd_delay_vaccination_to_infection_ring = sd_delay_vaccination_to_infection_ring,
      max_delay = max_delay,
      prop_hcw = prop_hcw,
      break_point = break_point
    )
    
    # extract summaries from the simulation:
    # number of cases, whether the outbreak was controlled, and mean values of R
    # before and after intervention
    n_case[i] <- nrow(out) # TODO: be consistent with plurals
    controlled[i] <- attr(out, "controlled")
    n_intros[i] <- attr(out, "n_intros")
    break_point_used[i] <- attr(out, "break_point")
    before <- out$date_onset <= break_point
    R_before[i] <- mean(out$R[before])
    R_after[i] <- mean(out$R[!before])
    last_force_infection[i] <- attr(out, "last_force_infection")
    R_basic[i] <- current_R_basic
    intervention_efficacy[i] <- current_intervention_efficacy
    last_intro[i] <- out %>% 
      filter(is_intro) %>% 
      pull(date_onset) %>% 
      max()
    
    gc()
  }
  
  list_results[[j]] <- tibble(
    # keep track of parameter values
    label = sim_params$label[j],
    r_intro = sim_params$r_intro[j],
    n_case,
    n_intros,
    break_point_used,
    R_before,
    R_after, 
    R_basic,
    intervention_efficacy, 
    controlled,
    last_intro,
    last_force_infection
  )
}

# make a tibble for this scenario
res_ring_target <- bind_rows(list_results) %>%
  mutate(scenario = "ring and targeted") %>%
  select(scenario, everything())


##############################################################################

# Bind all rows

res_all <- bind_rows(res_prophylactic_target,
                     res_prophylactic_mass,
                     res_reactive_targeted,
                     res_reactive_mass,
                     res_ring,
                     res_ring_target,
                     res_no_vaccine)

# Plot the relevant graphs after running the simulations above
library(tidyverse)
library(cowplot)
library(incidence2)
source("CI_prop_binomial.R")

# Rename the output tibble to just 'res' for easier handling and consistency
res <- res_ring_target #res_ring #res_reactive_targeted #res_reactive_mass #res_prophylactic_mass #res_prophylactic_target 

# Filter the tibble so that we group the output wrt their intro rates
res_high_intros <- res %>% filter(r_intro == 0.06)
res_low_intros <- res %>% filter(r_intro == 0.002)

# First, plot an ECDF using the stat_ecdf() function
# High rate of introductions
ECDF_high_intros <-
  ggplot(res_high_intros, aes(x = n_case, col = label)) +
  stat_ecdf() + xlab("Predicted number of cases") +
  ylab("Cumulative probability") +
  labs(col = "Vaccination Status") + theme_bw()

# Select the 'muted' colour palette from incidence2 and plot
n_colours <- length(unique(res$label))
col_pal <- muted(n_colours)
ECDF_h_intros <-
  ECDF_high_intros +  scale_colour_manual(values = col_pal)

# Low rate of introductions
ECDF_low_intros <-
  ggplot(res_low_intros, aes(x = n_case, col = label)) +
  stat_ecdf() + xlab("Predicted number of cases") +
  ylab("Cumulative probability") +
  labs(col = "Vaccination Status") + theme_bw()

# Select the 'muted' colour palette from incidence2 and plot
ECDF_l_intros <-
  ECDF_low_intros +  scale_colour_manual(values = col_pal)

# Plot the 2 ecdfs together, using plot_grid (from the cowplot package)
plot_grid(
  ECDF_l_intros,
  ECDF_h_intros,
  labels = c('A', 'B'),
  align = 'h',
  ncol = 2
)

# Second plot shows the distribution of controlled cases

# Step 1: Find all controlled simulations
controlled_out <- subset(res, controlled)

# Step 2: Split the dataframe into sims with high and low intro rates

controlled_out_low_intros <-
  subset(controlled_out, r_intro == 0.002)
controlled_out_high_intros <-
  subset(controlled_out, r_intro == 0.06)

# Step 3: Plot the distributions of cases

# Violin plot when rate of intros is low
violin_low_intros <- ggplot(controlled_out_low_intros,
                            aes(x = label, y = n_case, fill = label)) +
  geom_violin(trim = FALSE, alpha = 0.4) + xlab("Vaccination Parameters")

# Change to log10 scale and use a brewer colour palette
violin_low_intros + scale_y_continuous(trans = 'log10') +
  scale_fill_brewer(palette = "YlGnBu") +
  theme_cowplot() + ylab("Number of cases")

# Violin plot when introduction rate is high
violin_high_intros <- ggplot(controlled_out_high_intros,
                             aes(x = label, y = n_case, fill = label)) +
  geom_violin(trim = FALSE, alpha = 0.4) + xlab("Vaccination Parameters")

# Change to log10 scale and use a brewer colour palette
violin_high_intros + scale_y_continuous(trans = 'log10') +
  scale_fill_brewer(palette = "YlGnBu") +
  theme_cowplot() + ylab("Number of cases")

violin_high_intros + ylim(0,100) +
  scale_fill_brewer(palette = "YlGnBu") +
  theme_cowplot() + ylab("Number of cases")


# Run an ANOVA on these different distributions
# First, look at the low introduction rate
# Base level (intercept) is the No Vaccine control
controlled_out_low_intros$label <-
  relevel(factor(controlled_out_low_intros$label), ref = "E:0/C:0")

# Build linear model
lm_low_intros <-
  lm(n_case ~ label, data = controlled_out_low_intros)
# Test for general differences between groups
anova(lm_low_intros)
# Test individual group against no-vaccine intercept
summary(lm_low_intros)

# Then, look at the high introduction rate
# Our base level (intercept) is the No Vaccine control
controlled_out_high_intros$label <-
  relevel(factor(controlled_out_high_intros$label), ref = "E:0/C:0")

# Build linear model
lm_high_intros <-
  lm(n_case ~ label, data = controlled_out_high_intros)
# Test for general differences between groups
anova(lm_high_intros)
# Test individual group against no-vaccine intercept
summary(lm_high_intros)

summary(aov(lm_high_intros))
TukeyHSD(aov(lm_high_intros), conf.level=.95)

# Third plot shows the % of controlled outbreaks given each vaccination scheme

# Obtain the frequencies of controlled outbreaks in table format
# Then convert into a dataframe
df_controlled_out_high_intros <- 
  as.data.frame(table(controlled_out_high_intros$label))

df_controlled_out_low_intros <- 
  as.data.frame(table(controlled_out_low_intros$label))

# Find the binomial confidence intervals of % controlled outbreaks
df_controlled_out_l_intros <- df_controlled_out_low_intros %>%
  mutate(p_low  = prop_ci(Freq, n_sims, "lower"),
         p_high = prop_ci(Freq, n_sims, "upper"))

df_controlled_out_h_intros <- df_controlled_out_high_intros %>%
  mutate(p_low  = prop_ci(Freq, n_sims, "lower"),
         p_high = prop_ci(Freq, n_sims, "upper"))

# Plot the results, including binomial CIs

# Low rate of intros
controlled_out_l_intros_plot <-
  ggplot(data = df_controlled_out_l_intros, aes(x = Var1, y = Freq / n_sims)) +
  geom_point() + geom_errorbar(aes(ymin = p_low, ymax = p_high)) +
  xlab("Vaccine Parameters") + ylab("Controlled Outbreaks")  +
  theme_cowplot()

# High rate of intros
controlled_out_h_intros_plot <-
  ggplot(data = df_controlled_out_h_intros, aes(x = Var1, y = Freq / n_sims)) +
  geom_point() + geom_errorbar(aes(ymin = p_low, ymax = p_high)) +
  xlab("Vaccine Parameters") + ylab("Controlled Outbreaks")  +
  theme_cowplot()

# View plots
controlled_out_l_intros_plot
controlled_out_h_intros_plot

# Fourth plot shows the mean R values for each vaccination scenario
# (before and after intervention)

ggplot(controlled_out_low_intros, aes(x = R_before)) +
  geom_histogram(color = "red",
                 fill = "orange",
                 alpha = 0.4) +
  facet_grid(label ~ .)

ggplot(controlled_out_low_intros, aes(x = R_after)) +
  geom_histogram(color = "darkblue",
                 fill = "lightblue",
                 alpha = 0.4) +
  facet_grid(label ~ .)

# Plot all R before and also after intervention in a single plot and add facets 
# for each vaccination scenario

# Low rate of introductions first
ggplot(controlled_out_low_intros) +
  geom_histogram(aes(x = R_before),
                 color = "red",
                 fill = "orange",
                 alpha = 0.4) +
  geom_histogram(aes(x = R_after),
                 color = "darkblue",
                 fill = "lightblue",
                 alpha = 0.4) +
  facet_grid(label ~ .) + theme_cowplot() + xlab("R value")

# Repeat for simulations using a high rate of introductions
ggplot(controlled_out_high_intros, aes(x = R_before)) +
  geom_histogram(color = "red",
                 fill = "orange",
                 alpha = 0.4) +
  facet_grid(label ~ .)

ggplot(controlled_out_high_intros, aes(x = R_after)) +
  geom_histogram(color = "darkblue",
                 fill = "lightblue",
                 alpha = 0.4) +
  facet_grid(label ~ .)

# Plot both R before and after interventions in a single plot and add facets for
# each vaccination scenario

ggplot(controlled_out_high_intros) +
  geom_histogram(aes(x = R_before),
                 color = "red",
                 fill = "orange",
                 alpha = 0.4) +
  geom_histogram(aes(x = R_after),
                 color = "darkblue",
                 fill = "lightblue",
                 alpha = 0.4) +
  facet_grid(label ~ .) + theme_cowplot() + xlab("R value")

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
 
res$R_after_modified <- res$R_after
res$R_after_modified[is.nan(res$R_after_modified)] <- 0

# reorder columns

# plot R values before/after
res %>%
  pivot_longer(cols = R_before:R_after, 
               names_to = "period", 
               values_to = "R") %>%
  mutate(period = gsub("R_", "", period), 
         factor(period, levels = c("before", "after"))) %>% 
  ggplot(aes(x = R)) + theme_bw() + 
  geom_violin(aes(y = label, color = label, x = R, linetype = period)) + 
  labs(y = "") + guides(color = "none") + 
  scale_linetype_manual("Time period", values = c(before = 1, after = 2)) + 
  theme(legend.position = "top")



# Bind rows corresponding to all vaccination strategies together
res_all <- bind_rows(
  res_prophylactic_target,
  res_prophylactic_mass,
  res_reactive_targeted,
  res_reactive_mass,
  res_ring,
  res_ring_target
)
saveRDS(res_all, file = "res_all.rds")
