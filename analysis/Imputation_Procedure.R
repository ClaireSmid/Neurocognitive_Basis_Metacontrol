# CR Smid 25 August 2022: 
# this script creates the imputed file to be used in the final analyses, effectively 
# only one participant gets imputed at the relevant (first) time point for the model
# based model free data. However some EF measures do get imputed for the first time point.
# missing data can be inspected via the VIM plots below. 

# clears workspace
rm(list=ls())

library(dplyr)
library(ggplot2)
library(extrafont)
library(tidyr)
library(effectsize)
library(rstatix)
library(ggpubr)

require(ggiraph)
require(ggiraphExtra)
require(plyr)
library(ggeffects)
library(stats)

library(confint)
library(boot)
library(lme4)
library(lmerTest)

# for missing data and such
library(mice)
library(VIM)

library(dplyr)
library(Rmisc)
library(ggpubr)

##################################################################################

### load all data
DF_s <- read.csv('data/MBMF_EF_BehavioralData.csv')
names(DF_s)
dim(DF_s)


### these names might have shifted in the final file - to be checked before running script
# these measures are mainly used for imputation

# names(DF)[4] <- "Date"
# names(DF)[6] <- "Online"
# names(DF)[42] <- "Perc_Delay" # tot_delay_perc
# names(DF)[48] <- "Perc_HiConf_Delay" # tot_delay_perc
# names(DF)[49] <- "Perc_LoConf_Delay" # tot_delay_perc
# names(DF)[50] <- "Perc_Conflict_Del_Diff" # hi_lo_conf_del_diff
# 
# names(DF)[53] <- "TD_AvgRT" # this was capped by 3 SD
# names(DF)[64] <- "TD_Conflict_RT_Diff" #hi_lo_conf_RT_3SD_Diff this was capped by 3 SD
# 
# names(DF)[88] <- "Corr_r_Exp_Delay" # tot_delay_perc
# names(DF)[89] <- "Corr_int_Exp_Delay" # tot_delay_perc
# names(DF)[90] <- "Corr_rsq_Exp_Delay" # tot_delay_perc
# names(DF)[91] <- "Corr_r_Exp_Reward" # tot_delay_perc
# names(DF)[92] <- "Corr_int_Exp_Reward" # tot_delay_perc
# names(DF)[93] <- "Corr_rsq_Exp_Reward" # tot_delay_perc
# 
# names(DF)[123] <- "Corsi_WM_Span" # corsi_max_wm_t
# names(DF)[126] <- "AXCPT_t" # PBI_t
# names(DF)[127] <- "CogFlex_SwitchRT_t" # cf_switchrt_t
# names(DF)[128] <- "CogFlex_t" # cogflex_t
# names(DF)[129] <- "SSRT" # REzssrtnr_t
# names(DF)[130] <- "FlankerSwitch_t" # REflankerswitch_t
# names(DF)[131] <- "FlankerInhib_t" # REflankerinh_t
# names(DF)[132] <- "Stroop_t" # stroop
# names(DF)[133] <- "OneBack_WM_t" # dprimeONEBACK_t
# names(DF)[134] <- "TwoBack_WM_t" # dprimeTWOBACK_t

DF$ID_no <- as.integer(DF$ID)

# remove the unknown participant numbers
DF = DF[which (DF$ID_no < 266),]

# only take session 0 and 1 (no MB testing in Session 2)
DF = DF[which (DF$Session < 2),]

nrow(DF)

# total NAs
sum(is.na(DF$Online))

#### only add these when imputing with medium ######
DF$Online[is.na(DF$Online) & DF$ID < 135] <- 0
sum(is.na(DF$Online))

DF$Online[is.na(DF$Online) & DF$ID > 134] <- 1
sum(is.na(DF$Online))

### checl session
sum(is.na(DF$Session))


# list the participants we want (the subset who completed model-based decision making)
# there is one participant who did MB task at the second session, but not the first one,
# we are only using the first session so effectively we are imputing one participant only
MB_IDs <- array(DF %>%
                  filter(!is.na(DF$lr)) %>%
                  distinct(ID))

MB_IDs <- as.numeric(unlist(MB_IDs))
length(MB_IDs) # 69 max participants
MB_IDs

# extract the participants we want
MBs <- DF[DF$ID %in% MB_IDs, ]

# this should say 138, 67
dim(MBs)
length(unique(MBs$ID))

## check missing for MBMF
keeps <- c("it_P6","lr_P6","eg_P6","w_P6",
           "st_P6","repst_P6","w_lo","w_hi",
           "w_diff")


mis <- MBs[, (names(MBs) %in% keeps)]
names(mis)

### Plotting missings for T0 only, 30% of the data is missing (for T1 mostly)
aggr_plot<-aggr(mis,col=c('navyblue','red'),numbers=T, sortVars=T,
                labels=T,cex.axis=0.99,oma=c(10,5,5,3))


# display missing data for T0
DF_T0 <- MBs[which (MBs$Session==0),]
names(DF_T0)
dim(DF_T0)

# display T0 demographics
table(DF_T0['Gender'])
table(DF_T0['School'])
table(DF_T0['Online'])
table(DF_T0['SES'])

# report age ranges
min(DF_T0$Age_Frac,na.rm=T)
max(DF_T0$Age_Frac,na.rm=T)
mean(DF_T0$Age_Frac,na.rm=T)
sd(DF_T0$Age_Frac,na.rm=T)

### scale some vars
names(MBs)

# these names might be different in final script
MBs$Age_Frac_Imp <- MBs$Age_Frac
MBs$SES_inv = ((MBs$SES - max(MBs$SES,na.rm=T)) * -1) + min(MBs$SES,na.rm=T)
MBs$SES_inv_z <- scale(MBs$SES_inv, center = TRUE, scale = TRUE)
MBs$SES_inv_z <- as.numeric(MBs$SES_inv_z)

MBs$DG_UG_Diff_z <- scale(MBs$DG_UG_Diff, center = TRUE, scale = TRUE)
MBs$DG_UG_Diff_z <- as.numeric(MBs$DG_UG_Diff_z)
MBs$DG_Coins_Given_z <- scale(MBs$DG_Coins_Given, center = TRUE, scale = TRUE)
MBs$DG_Coins_Given_z <- as.numeric(MBs$DG_Coins_Given_z)
MBs$UG_Coins_Given_z <- scale(MBs$UG_Coins_Given, center = TRUE, scale = TRUE)
MBs$UG_Coins_Given_z <- as.numeric(MBs$UG_Coins_Given_z)

# log transform some reaction times
MBs$DG_RT_log <- log(1/MBs$DG_RT)
MBs$UG_RT_log <- log(1/MBs$UG_RT)
MBs$UG_Decision_RT_log <- log(1/MBs$UG_Decision_RT)
MBs$DG_UG_RT_Diff_log <- log(1/MBs$DG_UG_RT_Diff)
MBs$TD_AvgRT_log <- log(1/MBs$TD_AvgRT)


# UPDATE: School has been left out of the online data for anonymity 

keeps <- c("ID","Session","Online","School","Age_Frac","SES_inv_z","Age_Frac_Imp",
           "Online","Gender",
           # MBMF measures
           "it_P6","lr_P6","eg_P6","w_P6","st_P6","repst_P6","w_lo","w_hi","w_diff",
           # EF measures
           "T_Vocab","T_Matrix",
           "Corsi_WM_Span","AXCPT_t","CogFlex_t","SSRT","FlankerSwitch_t","FlankerInhib_t",
           "Stroop_t","OneBack_WM_t","TwoBack_WM_t","Flank_Cong_Cost_t","Flank_Incon_Cost_t",
           # DM measures
           "DG_Coins_Given_z","UG_Coins_Given_z","DG_UG_Diff_z","UG_Offer_Accept",
           "DG_RT_log","UG_RT_log", "UG_Decision_RT_log","DG_UG_RT_Diff_log","TD_AvgRT_log",
           "Perc_Delay","Corr_r_Exp_Delay","Corr_int_Exp_Delay",
           "Corr_r_Exp_Reward","Corr_int_Exp_Reward"
           
)

DF_for_imp <- MBs[, (names(MBs) %in% keeps)]
names(DF_for_imp)

# we also impute SES and age, maybe keep in school too.
## Take out measures not to be used in imputing
# took out School - now used in imputation also
Demdrops <- c("ID","Group","Age_Frac",
              "maths_scores","total_scores","english_scores",
              "DG_UG_RT_Diff_log"
)

# so training coefficient and online included in the imputation
testdata <- DF_for_imp[, !(names(DF_for_imp) %in% Demdrops)]
names(testdata)
dim(testdata)

library(purrr)
library(tidyr)

testdata %>% 
  keep(is.numeric) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = 'free') +
  #geom_histogram()
  geom_density()

# visualise missing data -cant be plotted, too many combinations
aggr_plot<-aggr(testdata,col=c('navyblue','red'),numbers=T, sortVars=T,
                labels=T,cex.axis=0.8,oma=c(10,5,5,3))

names(testdata)

# impute data and save it at the end
tempData <- mice(testdata, m = 5, maxit = 500, meth = 'pmm', seed = 500)
summary(tempData)

# add the columns back in
extra_cols <- DF_for_imp[, (names(DF_for_imp) %in% Demdrops)]
#extra_cols <- extra_cols[, !(names(extra_cols) %in% "School")]
names(extra_cols)
imps <- complete(tempData,1)
DF_imp <- cbind(extra_cols,complete(tempData,1))
names(DF_imp)

write.csv(DF_imp,"data/Imputed_Data.csv",row.names=FALSE)
save(tempData,file="data/Imputed_Data.RData")

