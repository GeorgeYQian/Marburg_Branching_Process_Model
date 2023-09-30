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
draw_R_target <- function(n,
                          date_onset,
                          R_basic,
                          intervention_efficacy,
                          vaccination_efficacy,
                          vaccination_coverage,
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
  ## Determine who is vaccinated among potential introductions (miners)
  vaccinated <-  sample(
    c(TRUE, FALSE),
    size = n,
    prob = c(vaccination_coverage, (1 - vaccination_coverage)),
    replace = TRUE
  )
  
  ## Make sure there are no vaccinations in cases who are not introductions
  vaccinated[!is_intro] <- FALSE
  
  ## Add effect of vaccination for miners
  out[vaccinated] <- out[vaccinated] *
    (1 - vaccination_efficacy[vaccinated])
  
  
  # Vaccination of HCW
  ## Determine who is vaccinated among HCW
  vaccinated <-  sample(
    c(TRUE, FALSE),
    size = n,
    prob = c(vaccination_coverage * prop_hcw, (1 - vaccination_coverage * prop_hcw)),
    replace = TRUE
  )
  
  ## Make sure these people are not miners (i.e. potential intros)
  vaccinated[is_intro] <- FALSE
  
  ## Add effect of vaccination for HCW
  out[vaccinated] <- out[vaccinated] *
    (1 - vaccination_efficacy[vaccinated])
  
  
  out
}


# # Example
# x <- draw_R_target(
#   n = 20,
#   date_onset = 1:20,
#   R_basic = 1.3,
#   intervention_efficacy = 0.5,
#   vaccination_efficacy = rep(0.9,20),
#   vaccination_coverage = 0.7,
#   break_point = 9,
#   is_intro = c(rep(TRUE, 10), rep(FALSE, 10)),
#   prop_hcw = 1
# )
# x
# hist(x, nclass = 15, border = "white", col = 3,
#      xlab = "R value",
#      ylab = "frequency",
#      main = "Distribution of R")
