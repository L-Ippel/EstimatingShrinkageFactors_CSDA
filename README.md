# EstimatingShrinkageFactors_CSDA
This repository contains the code used in the paper "Online Estimation of Individual-Level Effects using Streaming Shrinkage Factors" published in Computational Statistics and Data Analysis. 

The "parallel_streaming_shrinkage_factors_study.R" runs the simulation study presented in the paper. This file loads the "sim_fun.R" file which contains all the functions required for the study. 

In order to run the application study, please request access to the LISS panel data. At this point, after agreed upon by Centerdata, the data file file can be obtained through me, called "Participation_panelmember_fieldworkperiod" document, containing participation from Nov. 2007 - Dec. 2011, start/stop date of participation and an id number --better solution is still in the making--.  The "liss.R" file reorganizes the data into a long format needed for the "lissStudyRevision.R", which runs the liss predict participation study. 

