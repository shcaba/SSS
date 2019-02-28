*Directions to making full assessment into exSSS
#STARTER
1) Turn off using par file
2) 0 # run display detail (0,1,2). Option 0 makes less output on the R terminal.
3) 2 # detailed age-structured reports in REPORT.SSO (0,1,2). Option 2 writes a smaller Report file, so should be quicker.
2) No jitter
3) Turn off all screen output, detailed reports
4) Make convergence 0.001 (don't want it too small or it will take forever)

#DATA FILE
1) Add another survey; make it for the beginning of the year (timing = 0.01) and add the area assignment
2) +2 for the number of index lines
3) Add fleet/survey line before indices
4) Add two entries (beginning year, year of depletion prior). Make sure first year matches beginning year of model and second year matches the year of the depletion estimated.
5) Make sure total removals are in catch (e.g., if discards are estimated in the
   model, take discard estiamtes and add to landings

#CONTROL FILE
1) Use the control_new file as the CTL file
2) Set first_age_at_maturity to 0
3) Turn off all estimated parameters except lnR0
4) Turn all parameters priors to negative numbers to avoid their contribution in the total likelihood
5) If parameters were estimated, set input values to those values using report file
6) No rec devs estimated; set est years to last year in model
7) Add q line for depletion. Turn on extra variance for all surveys except the depletion survey.
8) Add selectivity line for depletion (option 30)
9) Check selectivity inputs; use assessment output or set to maturity assumption
10) Turn off other data (e.g. dicards, mean weights, comps.) likelihoods (i.e., lambda=0): Add number of changes; make lambdas = 0
11) Make sure the M and h labels are "NatM_p_1_Fem_GP_1", "NatM_p_1_mAL_GP_1", AND "SR_steep" respectively

#Additional Files
1) Copy the CTL file and add "_orig" to the end of the name. The other CTL will be modified using the AIS draws.


#Calculate the selectivity parameter width from matuirty
slope<- -0.5683
intercept<-21.84
x<-seq(1,50,0.05)
mat.vec<-1/(1+exp(slope*(x-intercept)))
mat.vec[mat.vec>0.95]
width.95per<-(x[mat.vec>0.95]-intercept)[1]
width.95per

### Make sure the model is able to complete a run before running the AIS implementation


*Modifying an exSSS file to make into SSS
#DATA FILE 

#CONTROL FILE
1) Turn off lambda of all indices except depletion


#Check MLE output before running
library(r4ss)
update_r4ss_files()
mydir<-"C:/Users/Jason.Cope/Desktop/ex_SSS/"
base <- SS_output(file.path(mydir,"Petrale"),covar=FALSE)
SS_plots(base)


Sys.sleep(0.001)

#Debugging tricks
#Make sure 