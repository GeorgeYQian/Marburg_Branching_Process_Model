make_logistic  <-  function(n_draws = 1,
                            max_VE = 0.9,
                            time_to_max = 7,
                            mean_delay = 10,
                            sd_delay = 4,
                            max_delay = 30)
{
  if (max_VE == 0) {
    vacc_eff = 0
    return(vacc_eff)
  }
  
  # Delay in vaccination
  
  vacc_delay <- 4
  
  day_delay_det <- round(vacc_delay)
  days_delay <- seq(from = 1, to = day_delay_det)
  
  #Vaccine efficacy on these first few days is 0
  vacc_eff_start <- rep(0.000, day_delay_det)
  
  #Delay for vaccine to reach full efficacy
  delay_vacc_eff <- time_to_max
  days_delay_vacc_eff <- seq(from = 1, to = delay_vacc_eff) +
    day_delay_det
  
  
  #Maximum vaccine efficacy
  vacc_eff_max <- max_VE
  
  #Vaccine efficacy between 0 and maximum efficacy
  vac_eff_2 <-
    seq(from = 0,
        to = vacc_eff_max,
        length.out = delay_vacc_eff + 1)
  
  #days at maximum
  cluster_period <- max_delay + vacc_delay
  days_at_max <- max(0,
                     cluster_period - (day_delay_det + delay_vacc_eff))
  
  
  if (days_at_max > 0) {
    days_max_eff <- seq(from = 1, to = days_at_max) +
      day_delay_det +
      delay_vacc_eff
  }  else{
    days_max_eff <- vector()
  }
  
  #vaccine efficacy during this period
  vacc_eff_max_period <- rep(vacc_eff_max, days_at_max)
  
  # Combine into dataframe
  
  #First the days
  
  day <- c(days_delay,  days_delay_vacc_eff, days_max_eff)
  
  if (length(day) != cluster_period) {
    day <- day[1:cluster_period]
    
  }
  
  
  #Then the vaccine efficacy on each day is:
  vaccine_efficacy <-
    c(vacc_eff_start, vac_eff_2, vacc_eff_max_period)
  
  
  
  if (length(vaccine_efficacy) != cluster_period) {
    vaccine_efficacy <-  vaccine_efficacy[1:cluster_period]
    
  }
  
  vacc_eff_df <- data.frame(day, vaccine_efficacy)
  
  #find the parameters for the equation
  SS <-
    getInitial(vaccine_efficacy ~ SSlogis(day, alpha, xmid, scale), data =  vacc_eff_df)
  
  
  
  #we used a different parametrization
  K_start <- SS["alpha"]
  R_start <- 1 / SS["scale"]
  N0_start <- SS["alpha"] / (exp(SS["xmid"] / SS["scale"]) + 1)
  
  #the formula for the model
  log_formula <-
    formula(vaccine_efficacy ~ K * N0 * exp(R * day) / (K + N0 * (exp(R * day) -
                                                                    1)))
  #fit the model
  m <- nls(log_formula, start = list(K = K_start, R = R_start, N0 = N0_start))
  
  efficacies <- predict(m)
  efficacies_df <- data.frame(day, efficacies)
  
  eff_dist <- tail(efficacies_df, -vacc_delay)
  
  eff_dist[, 1] <- eff_dist[, 1] - vacc_delay
  
  
  #Draw delay from vaccination to infection from a normal distribution
  
  vacc_delay <-
    round(rnorm(n = n_draws, mean = mean_delay,  sd = sd_delay))
  day_draw_min <- pmin(vacc_delay, max_delay)
  day_draw <- pmax(day_draw_min, 1)
  vacc_eff <- pmin(1, eff_dist[day_draw, 2])
  
  return(vacc_eff)
}
