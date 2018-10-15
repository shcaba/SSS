############################################################
# EXAMPLE RUN SSS
############################################################
# load the library
library(sss)
vignette("sss")


?SSS
######################################################################################################
# Inputs to the SSS function
######################################################################################################
## filepath location where the SSS will look for model files and run the executable
## file.name vector of file names for the data and control file where the expected input is c("data file", "control file") 
## reps number of mcmc draws to perform
## seed.in seed number 
## Dep.in vector defining distribution, mean, sd, and bounds for depletion prior.  Expected input is c(distribution shape, mean, sd). The distribution options are 2 = 1 - beta, 4 = uniform, 10 = truncated normal.
## M.in vector defining natural mortality distribuition, mean, and . Expected input is c(distrbution shape for females, mean for females, sd for females, distribution shape for males, mean for males, sd for males). The distibution options are 0 = normal, 3 = lognormal, and 4 = uniform.
## SR_type The shape of the stock-recruitment curve. Options are based on SS stock-recruit options. Option 3 = Beverton-holt, 8 = Shepherd 3-parameter, 9 = Ricker 3-parameter
## h.in vector defining the steepness distribution, mean, and sd. Expected input is c(distribution, mean, sd). Distribution options are 2 = truncated beta, 10 = truncated normal, 30 = truncated lognormal, 4 = uniform.
## FMSY_M.in vector defining the Fmsy/M ratio distribution, mean, and sd. Expected input is c(distribution, mean, sd). Distribution options are; negative value = ?, 2 = truncated beta, 10 = truncated normal, 30 = truncated lognormal, 4 = uniform.
## BMSY_B0.in vector defining the Bmsy/B0 ratio distribution, mean, and sd. Expected input is c(distribution, mean, sd). Distribution options are; negative value = ?, 2 = truncated beta, 10 = truncated normal, 30 = truncated lognormal, 4 = uniform.
## L1.in vector defining the minimum length. This is an optional feature.  A vector of zeros will not draw values for this value. Expected input values are c(female mean, female sd, male mean, male sd)
## Linf.in vector defining the maximum length. This is an optional feature.  A vector of zeros will not draw values for this value. Expected input values are c(female mean, female sd, male mean, male sd)
## k.in vector defining the growth coefficient k. This is an optional feature.  A vector of zeros will not draw values for this value. Expected input values are c(female mean, female sd, male mean, male sd)
## Zfrac.Beta.in
## R_start vector allowing the user to control the starting R0 value in the control file. Expected value is c( switch option, input value) where the switch optionas are 1= draw from a random draw from a truncated lognormal distribution based on the input value and 0 = start from the input value.
## doR0.loop
## sum_age summary age for total biomass calculation used by SS
## sb_ofl_yrs the years for which OFL values should be calculated
## f_yr final model year
## year0 initial model year
## genders TRUE/FALSE allows for the user to specify whether or not sexes should have the same values from drawn parameters (e.g., natural mortality, L1, Linf)
## BH_FMSY_comp
########################################################################################################

# Assume a Beverton-Holt Stock Recruit relationship
# define the directory
Dir.in <-"C:/SSS/sssexample_BH"

POP.SSS.BH<-SSS(filepath= Dir.in, # location to run the model
                file.name=c("simple_pop.dat","simple_pop.ctl"), # data and ctl file names
                reps=10, # number of model runs
                seed.in=19, # seed number
                Dep.in=c(2, 0.60, 0.1), # c(distribution shape, mean, sd)
                M.in=c(3, 0.05, 0.2, 3, 0.05, 0.2), # c(distrbution shape for females, mean for females, sd for females, distribution shape for males, mean for males, sd for males)
                h.in=c(2, 0.70, 0.152), # c(distribution, mean, sd)
                FMSY_M.in=c(-30,0.8,0.2), # c(distribution, mean, sd)
                BMSY_B0.in=c(-2,0.4,0.05), # c(distribution, mean, sd)
                R_start=c(0,9), # c( switch option, input value)
                sb_ofl_yrs=c(2017,2018,2019), # projection years 
                f_yr=2016, # final year of the model
                year0=1918, # start year of the model
                genders=T # unique male and female M
)


####################################################################################
# Alternative SRR
####################################################################################
library(sss)
# define the directory
Dir.in <-"C:/SSS/sssexample_RickPow"

# Ricker power curve example
POP.SSS.gR<-SSS(filepath=Dir.in,
                file.name=c("simple_pop.dat","simple_pop.ctl"),
                reps=10, # Iterations
                seed.in=19, # Seed number
                M.in=c(3,0.035,0.4,3,0.037,0.4),
                Dep.in=c(2,0.32,0.2),
                SR_type=9,
                h.in=c(2,0.779,0.152),
                FMSY_M.in=c(30,0.8,0.2),
                BMSY_B0.in=c(2,0.4,0.05),
                R_start=c(0,9),
                sb_ofl_yrs=c(2017,2018,2019),
                f_yr=2016,
                year0=1918,
                genders=T)




