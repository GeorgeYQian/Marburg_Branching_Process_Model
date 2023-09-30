


# A branching process simulator that incorporates mass vaccination of 

# George Qian 2023 - MIT License


# ------------
# Requirements
# ------------

# Required packages
library(distcrete)
library(epitrix)
library(here)
library(tidyverse) # total overkill, we merely use tibble and dplyr


#------------------------------------------------
# Define the function. Its inputs are as follows:
#------------------------------------------------
#
# 1. R_basic, the basic reproduction number when a case is undetected
#
# 2. intervention_efficacy, the prop reduction of R when cases detected
#
# 3. serial_interval_mean, the mean serial interval observed for MVD
#
# 4. serial_interval_sd, the standard deviation of this serial interval 
#
# 5. r_daily_intro, the daily rate of intro from reservoir
#
# 6. max_duration, the maximum duration of the outbreak
#
# 7. vaccination_coverage, the proportion of the population we expect to 
#    vaccinate
#
# 8. vaccination_time_to_max, time between vaccination and maximum vaccine 
#    protection
#
# 9. mean_delay_vaccination_to_infection, the mean delay between vaccination 
#    and infection
#
# 10. sd_delay_vaccination_to_infection, the standard deviation in delays 
#     between vaccination and infection
#
# 11. max_delay, the maximum delay between infections in a cluster
#
# 12. break_point, the delay between the first case and the implementation of 
#     non-pharmaceutical interventions
#

branching_process_model_mass <- function(R_basic = 1.2,
                                         intervention_efficacy = 0.5,
                                         serial_interval_mean = 9,
                                         serial_interval_sd = 4,
                                         r_daily_intro = 50 / 784,
                                         max_duration = 365,
                                         vaccination_coverage = 0.75,
                                         max_vaccination_efficacy = 0.9,
                                         vaccination_time_to_max = 7,
                                         mean_delay_vaccination_to_infection = 20,
                                         sd_delay_vaccination_to_infection = 5,
                                         max_delay = 30,
                                         break_point = 90) {
  source("make_disc_gamma.R")
  source("helper_functions_marburg_mass.R")
  source("create_logistic_function.R")
  
  serial_interval <-
    make_disc_gamma(serial_interval_mean, serial_interval_sd)
  
  ## Create a vector of times for the loop
  vec_time <- seq_len(max_duration)
  
  # --------------
  # Initialization
  # --------------
  
  # This is simulation stuff, but only done once at the beginning
  
  ## Introductions: we impose the first one on day 1, and others are drawn
  ## randomly
  
  vaccination_efficacy_intros <- make_logistic(
    n_draws = 1,
    max_VE = max_vaccination_efficacy,
    time_to_max = vaccination_time_to_max,
    mean_delay = mean_delay_vaccination_to_infection,
    sd_delay = sd_delay_vaccination_to_infection,
    max_delay = max_delay
  )
  
  r_daily_intro_vaccinated <- 
    r_daily_intro * (1 - vaccination_efficacy_intros * vaccination_coverage)
  
  n_daily_intro <- rpois(max_duration, r_daily_intro_vaccinated)  
    
  
  
  # Ensure there is one introduced case on the first day of the outbreak
  n_daily_intro[1] <- 1L
  intro_onset <- rep(vec_time, n_daily_intro)
  
  ## Build a tibble - not absolutely needed, but nicer to inspect results
  # All cases will be stored in 'out'
  out <- tibble::tibble(case_id = seq_along(intro_onset),
                        date_onset = intro_onset)
  
  # Add a label signifying that these case are introductions
  out <- mutate(out,
                is_intro = TRUE)
  
  # Add vaccination efficacy; drawn from external function
  vaccination_efficacy <- make_logistic(
    n_draws = nrow(out),
    max_VE = max_vaccination_efficacy,
    time_to_max = vaccination_time_to_max,
    mean_delay = mean_delay_vaccination_to_infection,
    sd_delay = sd_delay_vaccination_to_infection,
    max_delay = max_delay
  )
  
  # Draw values of R using external function
  R <- draw_R_mass(
    n = nrow(out),
    date_onset = out$date_onset,
    R_basic = R_basic,
    intervention_efficacy = intervention_efficacy,
    vaccination_efficacy = vaccination_efficacy,
    vaccination_coverage = vaccination_coverage,
    break_point = break_point
  )
  
  out <- mutate(out,
                vaccination_efficacy = vaccination_efficacy,
                R = R)
  
  # Record the number of introductions: this is the number of rows in 'out'
  n_intros <- nrow(out)
  
  
  # -------------------------------------
  # Time iteration: what happens each day
  # -------------------------------------
  
  # Here we generate new cases from a branching process using the state of 'out'
  # above. The algorithm is:
  
  # 1. Determine the force of infection, i.e. the average number of new cases
  # generated.
  #
  # 2. Draw the number of new cases.
  #
  # 3. Draw features of the new cases.
  #
  # 4. Append new cases to the linelist of cases.
  #
  #
  # The main issue lies in determining the force of infection. For a single
  # individual 'i', it is determined as:
  #
  # lambda_i = R_i * w(t - onset_i)
  #
  # where 'w' is the PMF of the serial interval, 't' the current time, and
  # 'onset_i' the date of onset of case 'i'.
  #
  # The global force of infection at a given point in time is obtained by summing
  # all individual forces of infection.
  #
  
  
  # Time iteration
  for (t in 2:max_duration) {
    
    # Step 1
    lambda_i <- out$R * serial_interval$d(t - out$date_onset)
    force_infection <- sum(lambda_i)
    
    # Step 2
    n_new_cases <- rpois(1, lambda = force_infection)
    
    # Step 3
    
    if (n_new_cases > 0) {
      last_id <- max(out$case_id)
      new_cases <- tibble(
        case_id = seq(
          from = last_id + 1,
          length.out = n_new_cases,
          by = 1L
        ),
        date_onset = rep(t, n_new_cases)
      )
      
      new_cases <- mutate(new_cases,
                          is_intro = FALSE)
      
      # Draw features of the new cases
      # Add vaccination efficacy; drawn from external function
      vaccination_efficacy <- make_logistic(
        n_draws = n_new_cases,
        max_VE = max_vaccination_efficacy,
        time_to_max = vaccination_time_to_max,
        mean_delay = mean_delay_vaccination_to_infection,
        sd_delay = sd_delay_vaccination_to_infection,
        max_delay = max_delay
      )
      
      # Draw values of R using external function
      R <- draw_R_mass(
        n = n_new_cases,
        date_onset = new_cases$date_onset,
        R_basic = R_basic,
        intervention_efficacy = intervention_efficacy,
        vaccination_efficacy = vaccination_efficacy,
        vaccination_coverage = vaccination_coverage,
        break_point = break_point
      )
      
      new_cases <- mutate(new_cases,
                          vaccination_efficacy = vaccination_efficacy,
                          R = R)
      
      # Step 4
      
      # Now add each new case to the 'out' dataframe (which stores all cases)
      out <- bind_rows(out, new_cases)
      
    }
    
    # Insert a break if we already have 5k cases
    if (nrow(out) > 5000) {
      break
    }
  }
  
  # Arrange cases by chronological order
  out <- arrange(out, date_onset)
  
  # Redefine IDs to match chronological order
  out <- mutate(out, case_id = seq_len(n()))
  
  # Check if the outbreak is controlled; criteria is the last force of infection
  # is negligible (less than 0.05); note that force_infection is whatever the
  # last value was in the for loop
  controlled <- force_infection < 0.05
  attr(out, "controlled") <- controlled
  attr(out, "last_force_infection") <- force_infection
  attr(out, "n_intros") <- n_intros
  attr(out, "break_point") <- break_point
  
  out
}

# # Example
# set.seed(1)
# x <- branching_process_model_mass()
# attr(x, "controlled")
# attr(x, "last_force_infection")
