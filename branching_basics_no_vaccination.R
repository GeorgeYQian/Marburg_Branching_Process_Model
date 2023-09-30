


# A branching process simulator; here, we simulate an outbreak with no 
# vaccination scheme

# George Qian 2022 - MIT License


# ------------
# Requirements
# ------------

# Note: these should ultimately be better handled through a package structure
# etc.

# Required packages
library(distcrete)
library(epitrix)
library(here)
library(tidyverse) # total overkill, we merely use tibble and dplyr


#------------------------------------------------
# Define the function. Its inputs are as follows:
#------------------------------------------------
#
# 1. R_undetected, the basic reproduction number when a case is undetected
#
# 2. intervention_efficacy, the prop reduction of R when cases detected
#
# 3. serial_interval, a distcrete object representing the serial interval
#    Use the 'make_disc_gamma' function to create this object
#
# 4. r_daily_intro, the daily rate of intro from reservoir
#
# 5. max_duration, the maximum duration of the outbreak
#
# 6. p_detected, proportion of cases detected

# Note that targeted vaccination is very similar to mass vacc - the difference
# being that we only apply vaccination to cases who are introductions

branching_process_model_no_vacc <- function(R_basic = 1.2,
                                         intervention_efficacy = 0.9,
                                         serial_interval_mean = 9,
                                         serial_interval_sd = 4,
                                         r_daily_intro = 50 / 784,
                                         max_duration = 365,
                                         break_point = 90) {
  source("make_disc_gamma.R")
  
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
  
  n_daily_intro <- rpois(max_duration, r_daily_intro)
  # Ensure there is one introduced case on the first day of the outbreak
  n_daily_intro[1] <- 1L
  intro_onset <- rep(vec_time, n_daily_intro)
  
  ## Build a tibble - not absolutely needed, but nicer to inspect results
  out <- tibble::tibble(
    case_id = seq_along(intro_onset),
    date_onset = intro_onset,
    is_intro = TRUE
  )
  
  # Outside the for loop, we only deal with introduced cases, which will not be
  # contact traced and, therefore, remain unvaccinated.
  
  # For now, all intros have an R value of R_basic
  R_intros <- rep(R_basic, nrow(out))
  
  # Add the effect of intervention
  after_intervention <- out$date_onset > break_point
  
  R_intros[after_intervention] <-   R_intros[after_intervention] *
    (1 - intervention_efficacy)
  
  out <- mutate(out,
                R = R_intros)
  
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
      
      # For now, all non-introduced cases have an R value of R_basic
      R_new <- rep(R_basic, nrow(new_cases))
      
      # Add the effect of intervention
      after_intervention <- new_cases$date_onset > break_point
      
      R_new[after_intervention] <-   R_new[after_intervention] *
        (1 - intervention_efficacy)
      
      new_cases <- mutate(new_cases,
                          R = R_new,
                          is_intro = FALSE)
      
      # Step 4
      out <- bind_rows(out, new_cases)
      
    }
    
    # Insert a break if we already have 5k cases
    if (nrow(out) > 5000) {
      break
    }
  }
  
  # ------------
  # Check output
  # ------------
  
  # arrange cases by chronological order
  out <- arrange(out, date_onset)
  
  # redefine IDs to match chrono order
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
  
} #end function here


# # Example
# set.seed(1)
# x <- branching_process_model_no_vacc()
