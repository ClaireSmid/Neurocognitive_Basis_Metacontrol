# Read me
Repo for the data and code used in the publication "Neurocognitive underpinnings of model-based decision making and metacontrol in childhood"

Initial Update 14/07/2022: uploaded data files, starting adding analysis scripts. 

Update 24/08/2022: added data files, analysis files and linked to the surfstat page

Update 02/09/2022: The Developmental Cognitive Neuroscience Journal does not allow supplemental material, so I have included additional information on the executive function tasks here, which are not directly relevant to the main manuscript.

This is the main folder to replicate the analyses conducted in the paper: Neurocognitive underpinnings of model-based decision making and metacontrol in childhood. Preprint will be accessible on bioRxiv shortly.

The computational modelling code is nearly identical to the one used in https://github.com/ClaireSmid/Model-based_Model-free_Developmental.

Follow the steps below to recreate the analysis as reported. For any questions, please email to: claire.r.smid@gmail.com

## Task

The current task was identical to the task used in Smid et al. 2022 (except the number of trials was dropped from 140 to 102). Stimuli for this task has been developed and shared previously by the authors of 'From Creatures of Habit to goal-directed learners' by Decker et al. 2016: https://pubmed.ncbi.nlm.nih.gov/27084852/

## Data

The behavioral data (N = 69) and the structural cortical thickness data (N = 44) is included in the 'data' folder. 

## Analysis

  Step 1. Running the behavioral analyses

The data from the sequential decison making task and the executive function tasks, as well as the demographic information is included in the analyses folder. Missing data was imputed (mainly for the executive function measures) and the procedure for this is included in the R file. 

  Step 2. Running the cortical thickness analyses
  
Cortical thickness analysis was conducted using SurfStat on smoothed data, prior to this structural data was ran through FreeSurfer, visually inspected and corrected if needed. The SurfStat package can be downloaded from here: https://www.math.mcgill.ca/keith/surfstat/

