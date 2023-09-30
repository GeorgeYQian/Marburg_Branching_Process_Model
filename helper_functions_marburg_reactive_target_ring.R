# -------
# Helpers
# -------
# These are small functions to avoid repetitive code in the simulation

# Draw individual R based on group values and frequencies
# n: number of R values to draw
# date_onset: a vector of length n of dates of onset
# R_basic: a single value of R, before intervention
# intervention_efficacy: a single value, reduction of transmission after intervention
# vaccination_efficacy: a vector of length n, reduction of transmission due to vaccination
# vaccination_coverage: the probability of cases getting vaccinated
# break_point: the time after which intervention kicks in
# is_intro: a boolean vector of length n, signifying whether the case is an introduction
draw_R_reactive_target_ring <- function(n,
                                        date_onset,
                                        R_basic,
                                        intervention_efficacy,
                                        reporting,
                                        vaccination_coverage,
                                        vaccination_efficacy_target,
                                        vaccination_efficacy_ring,
                                        break_point,
                                        is_intro,
                                        prop_hcw) {
  ## Determine R for each case
  
  # Generate output: a vector of R values of size n
  out <- rep(R_basic, n)
  
  # Add the effect of intervention
  after_intervention <- date_onset > break_point
  out[after_intervention] <-
    out[after_intervention] * (1 - intervention_efficacy)
  
  # Vaccination of miners
  ## Determine who is vaccinated
  vaccinated <-  sample(
    c(TRUE, FALSE),
    size = n,
    prob = c(vaccination_coverage, (1 - vaccination_coverage)),
    replace = TRUE
  )
  
  ## Make sure there are no vaccinations in cases who are not introductions
  vaccinated[!is_intro] <- FALSE
  
  ## Also ensure there are no vaccinations before intervention
  vaccinated[!after_intervention] <- FALSE
  
  ## Add effect of vaccination for miners
  out[vaccinated] <- out[vaccinated] *
    (1 - vaccination_efficacy_target[vaccinated])
  
  # Vaccination of HCW
  ## Determine who is vaccinated among HCW
  vaccinated_HCW <-  sample(
    c(TRUE, FALSE),
    size = n,
    prob = c(
      vaccination_coverage * prop_hcw,
      (1 - vaccination_coverage * prop_hcw)
    ),
    replace = TRUE
  )
  
  ## Make sure these people are not miners (i.e. potential intros)
  vaccinated_HCW[is_intro] <- FALSE
  
  ## Also ensure there are no vaccinations before intervention
  vaccinated_HCW[!after_intervention] <- FALSE
  
  ## Add effect of HCW vaccination (at targeted vaccination efficacy)
  out[vaccinated_HCW] <- out[vaccinated_HCW] *
    (1 - vaccination_efficacy_target[vaccinated_HCW]) 
  
  
  # Add effects of ring vaccination (contacts and contacts of contacts)
  ## Determine who is vaccinated
  vaccinated_ring <- sample(
    c(TRUE, FALSE),
    size = n,
    prob = c(vaccination_coverage * reporting,
             (1 - ((vaccination_coverage) * reporting
             ))),
    replace = TRUE
  )
  
  ## Make sure these people are not potential intros
  vaccinated_ring[is_intro] <- FALSE
  
  ## Make sure there are no vaccinations before intervention
  vaccinated_ring[!after_intervention] <- FALSE
  
  # We need to ensure we don't simulate vaccination twice, if an individual is
  # both a HCW and part of the ring of contacts. Hence, we will identify these
  # individuals and only have them as part of targeted vaccination (as we assume
  # that this is quicker than ring)
  
  vaccinated_ring[vaccinated_HCW & vaccinated_ring] <- FALSE
  
  ## Add effect of vaccination
  out[vaccinated_ring] <- out[vaccinated_ring] *
    (1 - vaccination_efficacy_ring[vaccinated_ring])
  
  out
}


# Example
# x <- draw_R_reactive_target_ring(
#   n = 50,
#   date_onset = 1:50,
#   R_basic = 1.3,
#   intervention_efficacy = 0.5,
#   reporting = 0.8,
#   vaccination_coverage = 0.7,
#   vaccination_efficacy_target = rep(0.7,50),
#   vaccination_efficacy_ring = rep(0.6,50),
#   prop_hcw = 0.15,
#   break_point = 5,
#   is_intro = c(rep(TRUE, 5), rep(FALSE, 45))
# )
# x
# hist(x, nclass = 15, border = "white", col = 3,
#      xlab = "R value",
#      ylab = "frequency",
#      main = "Distribution of R")
# 
# 
# # Extreme values example
# x <- draw_R_reactive_target_ring(
#   n = 50,
#   date_onset = 1:50,
#   R_basic = 1.3,
#   intervention_efficacy = 0.3,
#   reporting = 1.0,
#   vaccination_coverage = 0.9,
#   vaccination_efficacy_target = rep(1,50),
#   vaccination_efficacy_ring = rep(0.6,50),
#   prop_hcw = 1,
#   break_point = 5,
#   is_intro = c(rep(TRUE, 5), rep(FALSE, 45))
# )
# x
# hist(x, nclass = 15, border = "white", col = 3,
#      xlab = "R value",
#      ylab = "frequency",
#      main = "Distribution of R")
# 
