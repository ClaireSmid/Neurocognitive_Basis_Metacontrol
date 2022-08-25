# Overview of data

File created: 25/08/2022

## Behavioral measures of the task
The file 'MB_Rocket_Check_T0.csv' contains the data from the data check when participants were asked to confirm the transition structure of the rockets. This was done both directly after training, as well as at the end of the task. 

The file 'DMTask_Performance.csv' contains the performance measure, which was a participants's individual rate of reward controlled for by their potential reward. If the value is positive, it means they outperformed the average reward rate, and if it is negative they underperformed. The first performance measure reflects performance across the whole task, the last two columns reflect performance for the low and high stake trials respectively. 

To obtain behavioral measures of model-free and model-based decision making based on stay behavior, a regression approach was used (similar as in Kool et al. 2016 and Smid et al. 2022), the code for this is included in the analysis folder. 

## Model-free and model based measures and executive functions (final imputed dataset)
The file 'Imputed_Data_11Jul22.csv' contains the final behavioral measures used in the analyses reported in the paper. This file has imputed data, and can be recreated by using the R script in the analyses folder. This file was created by using the file in the section below.

## Model-free and model based measures and executive functions (pre-imputed dataset)
Data was collected as part of a larger psychological battery (e.g. alongside the data from the executive function tasks), and as part of a training paradigm (three testing sessions). Due to Covid (from March 2020 until Jul 2021) data collection of the sequential decision making task was halted, leading to a much smaller sample at the second testing session and no data for the final (third) time point. For the current paper, only data at the first timepoint was analysed, however, data from the first and second session were imputed to add one more participant to the pre-training timepoint who had completed the sequential task at the second session but was unable to complete it at the first session.

The file 'MBMF_EF_BehavioralData.csv' contains all data across time points. 

The imputation procedure is included in the analysis folder in the R script. The data from the second time point was not further used for the current paper, and only the data from the first time point analysed.

## FreeSurfer
This folder contains the final cortical thickness measures used in this paper. The folder 'fsaverage5' has the common template used for the analyses. 
