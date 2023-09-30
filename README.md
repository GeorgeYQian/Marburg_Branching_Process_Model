# MVD-Branching-Process-Model-Repository

This repository provides the source code and related data files used in:

    "Modelling Vaccination Strategies for the Control of Marburg Virus Disease Outbreaks.".
    [Full journal info and link to be provided]
 

### Data files

1. Posterior_Distrib_MVD.csv

    This is a csv file that contains the pooled posterior distribution of our estimates of the basic (pre-intervention) reproduction number and intervention efficacy of each previous MVD outbreak. In addition, there are further columns stating the number of introductions and outbreak duration. 



### Code

# Branching process model files

1. branching_basics_no_vaccination.R

   This is the main code for the branching process model, where we simulate an MVD outbreak in the absence of vaccination.

2. branching_basics_prophylactic_mass.R

      This is the main code for the branching process model, where we simulate an MVD outbreak with prophylactic mass vaccination.

3. branching_basics_prophylactic_targetted.R

      This is the main code for the branching process model, where we simulate an MVD outbreak with prophylactic targetted vaccination.

4. branching_basics_reactive_mass.R

      This is the main code for the branching process model, where we simulate an MVD outbreak with reactive mass vaccination.

5. branching_basics_reactive_targetted.R

      This is the main code for the branching process model, where we simulate an MVD outbreak with reactive targetted vaccination.

6. branching_basics_ring.R

      This is the main code for the branching process model, where we simulate an MVD outbreak with ring vaccination.

7. branching_basics_reactive_targetted_ring.R

      This is the main code for the branching process model, where we simulate an MVD outbreak with a combination of ring and reactive targetted vaccination.
      
# Helper functions

1. helper_functions_marburg_mass.R

     A helper function corresponding to the main branching process model code for prophylactic mass vaccination. The function helps avoid repetitive code in the simulation.
     
2. helper_functions_marburg_target.R

    A helper function corresponding to the main branching process model code for prophylactic targetted vaccination. The function helps avoid repetitive code in the simulation.
    
3. helper_functions_marburg_reactive_mass.R

     A helper function corresponding to the main branching process model code for reactive mass vaccination. The function helps avoid repetitive code in the simulation.
     
4. helper_functions_marburg_reactive_target.R

    A helper function corresponding to the main branching process model code for reactive targetted vaccination. The function helps avoid repetitive code in the simulation.
    
5. helper_functions_marburg_ring.R

    A helper function corresponding to the main branching process model code for ring vaccination. The function helps avoid repetitive code in the simulation.
    
6. helper_functions_marburg_reactive_target_ring.R

    A helper function corresponding to the main branching process model code for the combined ring plus reactive targetted vaccination. The function helps avoid repetitive code in the simulation.

7. make_disc_gamma.R

      A helper function that creates a discretised Gamma distribution with given
 mean and standard deviration, returning a `distcrete` object.
 
 8. create_logistic_function.R

      This file creates a logistic function that models vaccine efficacy over time. 

### RStudio and packages

The above code was written with R v4.1.0 in mind.

Furthermore, the following packages were used:

-	tidyverse 
-	stats 
-	Rmpfr 
-	distcrete
-	epitrix
-	here
-	incidence2
