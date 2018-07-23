R and FORTRAN functions related to the paper:

"A shared-parameter continuous-time hidden Markov and survival model for longitudinal data with informative drop-out"
by F.Bartolucci and A.Farcomeni (2018)

In order to use these functions, a preliminary step is to compile the
FORTRAN code to obtain executable binaries. To do so, type the
following commands in your local R version:

system("R CMD SHLIB lk_rec.f --preclean")
system("R CMD SHLIB lk_rec_surv.f --preclean")
system("R CMD SHLIB lk_comp_joint.f --preclean")
system("R CMD SHLIB lk_comp_joint_surv.f --preclean")
system("R CMD SHLIB dlk_comp_joint.f --preclean")
system("R CMD SHLIB dlk_comp_joint_surv.f --preclean")

Additionally, the following libraries must be installed:

library(expm)
library(survival)

The repository is made of .f FORTRAN functions for computing the
log-likelihood and its gradient by means of the recursions described
in the paper. Additionally, there are the following R files:

est_comp_joint_HMcontinuous.r: it contains the functions for
estimating the process with complete data. Its main function is
est_comp_joint_HMcontinuous. All comments and details are in the file
as R comments. 

est_incomp_joint_HMcontinuous.r: it contains the functions for
estimating the process with incomplete data (real data scenarios).
Its main function is
est_incomp_joint_HMcontinuous. All comments and details are in the file
as R comments. 

simulate_joint_HMcontinuous.r: it containts a function to simulate
data as in the simulation setting in the paper. The function is called
simulate_joint_HMcontinuous and is duly commented in the file. 

example_joint_HMcontinuous.r: it contains two examples based on
Gaussian and then binomial longitudinal outcomes. All details are in
the file as R comments.


