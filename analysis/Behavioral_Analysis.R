# CR Smid 25/08/2022: This script runs the main analyses in the neurocognitive underpinnings...
# paper, and uses the imputed data from the other R script. 

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

DF_imp <- read.csv('data/Imputed_Data_11Jul22.csv')
names(DF_imp)

# this should be named tempData
load("data/Imputed_Data_11Jul22.RData")

# creating factors
#DF_imp$ID <- as.factor(DF_imp$ID)
DF_imp$Online <- as.factor(DF_imp$Online)
DF_imp$SessionFactor <- as.factor(DF_imp$Session)

# name levels
levels(DF_imp$Online) <- c("Off","On")

library(psych)
library(miceadds)

names(DF_imp)

### plot new and old distributions
nrow(complete(tempData))

# filter on session
tempData_T0 <- filter(tempData, Session == 0)

# assess the distributions - first demographic variables and IQ (Dont plot train coef for this)
vnames <- c("it_P6","lr_P6","eg_P6","w_P6",
            "st_P6","repst_P6","w_lo","w_hi",
            "w_diff")

cd <- miceadds::subset_datlist(tempData, select = vnames, toclass='mids')
densityplot(cd)

# next EF measures
names(DF_imp)
vnames <- c("Corsi_WM_Span","AXCPT_t","CogFlex_t",
            "SSRT","T_Matrix","Flank_Cong_Cost_t","Flank_Incon_Cost_t",
            "FlankerSwitch_t","FlankerInhib_t","Stroop_t",
            "OneBack_WM_t","TwoBack_WM_t","T_Vocab")

cd <- mice::complete(tempData)[, vnames]

cd <- miceadds::subset_datlist(tempData, select = vnames,toclass='mids')
densityplot(cd)

# next EF measures for T0
cd <- mice::complete(tempData_T0)[, vnames]

cd <- miceadds::subset_datlist(tempData_T0, select = vnames,toclass='mids')
densityplot(cd)

##############################################################################
###### --------------- 1. Model-based decision-making data
##############################################################################

MB_imp_T0 <- DF_imp[DF_imp$Session == 0,]

names(MB_imp_T0)

# display T0 demographics
table(MB_imp_T0['Gender'])
table(MB_imp_T0['School']) # not included in online data
table(MB_imp_T0['SES_inv_z'])

### AGE
# first, some age related and sample related questions 
# does w for the 6 parameter model increase with age for the main study data (same as for the Bursted wood data?)
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$w_P6)# now at p = .035, r = .0.25 # now at p = .042, r = .25 # yes, significant at p = .033, r = .26, so just about

# stakes effect?
MB_imp_T0$w_lo <- as.numeric(MB_imp_T0$w_lo)
MB_imp_T0$w_hi <- as.numeric(MB_imp_T0$w_hi)

t.test(MB_imp_T0$w_lo,MB_imp_T0$w_hi)

# what is the mean w?
mean(MB_imp_T0$w_P6) # w = .55
mean(MB_imp_T0$w_diff) # w = .03

# do any other parameters increase with age?
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$it)   # r = .29, p = .016 # yes, p = .020, r = .28
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$lr)   # r = .24, p = .005 # yes, p = .004, r = .34
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$eg)   # r = .13, p = .310 # no,  p = .181, r = .16
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$st)   # r = -.15, p = .233 # no,  p = .221, r = -.15

# 6 parameter model
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$it_P6)    # yes, r = .30, p = .011
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$lr_P6)    # yes, r = .26, p = .002
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$eg_P6)    # no,  r = .19, p = .121
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$st_P6)    # no,  r = -.15, p = .217
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$repst_P6) # no, r = -.18, p = .137

# model fit parameters, do they change with age? - they don't, which is good, I think
#cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$bic) # no
#cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$aic) # no

# performance. does that change over age? # NOT IMPUTED
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$Avg_Pts) # yes, p = .017, r = .29 - exact same
sum(is.na(MB_imp_T0$Avg_Pts))

# same questions for the 7 parameter model:
# do the ws for the 7 parameter model increase with age for the main study data (same as for the Bursted wood data?)
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$w_lo) # no, p = .151, r = .18
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$w_hi) # yes, just about: p = .041, r = .25

# what are the means of the ws?
mean(MB_imp_T0$w_lo) # w_lo = .52
mean(MB_imp_T0$w_hi) # w_hi = .55

# correlation between the ws?
cor.test(MB_imp_T0$w_lo, MB_imp_T0$w_hi) # no? interesting.

MB_imp_T0$wdiff <- MB_imp_T0$w_hi - MB_imp_T0$w_lo

cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$wdiff) # NS
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$w_diff) # NS

# is there a stakes effect?
t.test(MB_imp_T0$w_lo, MB_imp_T0$w_hi, paired = TRUE) # no... p = .181


# performance. does that change over age? (need from other file)
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$Avg_Pts_lo) # no,  p = .081, r = .21
cor.test(MB_imp_T0$Age_Frac_Imp, MB_imp_T0$Avg_Pts_hi) # yes, p = .020, r = .28

cor.test(DF_imp_T0$w_lo, DF_imp_T0$Avg_Pts_lo) # no,  
cor.test(DF_imp_T0$w_hi, DF_imp_T0$Avg_Pts_hi) # yes, 



################################################################################
### ---------- Executive functions relations
################################################################################

### EF relationships
names(MB_imp_T0)

# recode these for clarity
MB_imp_T0$FlankerInhib_t <- MB_imp_T0$FlankerInhib_t*-1
MB_imp_T0$FlankerSwitch_t <- MB_imp_T0$FlankerSwitch_t*-1

MB_imp_T0$Age_Frac_Imp_z <- scale(MB_imp_T0$Age_Frac_Imp, center = TRUE, scale = TRUE)

# metacontrol strongest correlations
cor.test(MB_imp_T0$w_diff,MB_imp_T0$FlankerInhib_t) # SIG, r = -.37, p = .002
cor.test(MB_imp_T0$w_diff,MB_imp_T0$Stroop_t) # marginal, r = .21, p = .051
cor.test(MB_imp_T0$w_diff,MB_imp_T0$SSRT )  # NS r = .13, p = .271

# model-based strongest correlations
cor.test(MB_imp_T0$w_P6,MB_imp_T0$CogFlex_t) # SIG, r = -.26, p = .034
cor.test(MB_imp_T0$w_P6,MB_imp_T0$Corsi_WM_Span) # Sig, at p = .039
cor.test(MB_imp_T0$w_P6,MB_imp_T0$SSRT) # marginal p = .057

library(ppcor)

pcor.test(MB_imp_T0$w_diff,MB_imp_T0$FlankerInhib_t,MB_imp_T0$Age_Frac_Imp)
pcor.test(MB_imp_T0$w_diff,MB_imp_T0$Stroop_t,MB_imp_T0$Age_Frac_Imp)
pcor.test(MB_imp_T0$w_diff,MB_imp_T0$SSRT,MB_imp_T0$Age_Frac_Imp)

pcor.test(MB_imp_T0$w_P6,MB_imp_T0$CogFlex_t,MB_imp_T0$Age_Frac_Imp) # sig
pcor.test(MB_imp_T0$w_P6,MB_imp_T0$Corsi_WM_Span,MB_imp_T0$Age_Frac_Imp) # NS
pcor.test(MB_imp_T0$w_P6,MB_imp_T0$Stroop_t,MB_imp_T0$SSRT) # NS


# check whether relations hold when corrected for age
mod <- lm(w_diff ~ FlankerInhib_t + Age_Frac_Imp_z, data = MB_imp_T0) # Sig
summary(mod)

mod <- lm(w_diff ~ Stroop_t + Age_Frac_Imp_z, data = MB_imp_T0) # Sig
summary(mod)

mod <- lm(w_P6 ~ CogFlex_t + Age_Frac_Imp_z, data = MB_imp_T0) # Sig
summary(mod)

##
mod <- lm(w_P6 ~ Corsi_WM_Span + Age_Frac_Imp_z, data = MB_imp_T0) # NS
summary(mod)

mod <- lm(w_P6 ~ SSRT + Age_Frac_Imp_z, data = MB_imp_T0) # NS
summary(mod)


### split up flanker inhibition -- to be calculated separately
cor.test(MB_imp_T0$w_diff,MB_imp_T0$Flank_Cong_Cost_t)
cor.test(MB_imp_T0$w_diff,MB_imp_T0$Flank_Incon_Cost_t)
cor.test(MB_imp_T0$w_diff,MB_imp_T0$FlankerInhib_t)

#################################################################################
#################################################################################
### --- PLOTS
#################################################################################
#################################################################################

### themes for plots

claire_theme <-   theme(
  #legend.position = "none",
  plot.title = element_text(family="Arial",color="black", size=26 ,hjust = 0.5,margin=margin(0,0,10,0)),
  text = element_text(family="Arial", size=20),
  legend.title = element_text(size = 26),
  legend.text = element_text(size = 22),
  axis.title.y = element_text(color="black", size=24),
  axis.title.x = element_text(color="black", size=24, margin=margin(5,0,0,0)),
  axis.text.x = element_text(size = 22, margin=margin(5,0,0,0)),
  axis.text.y = element_text(size = 22, margin=margin(0,5,0,10))
)

fig3_theme <- theme(
  text = element_text(family="Arial", size=16),
  plot.title = element_text(color="black", size=20 ,hjust = 0.5,margin=margin(0,0,10,0)),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black"),
  #legend.title = element_blank(),
  axis.title.x = element_text(color="black", size=18, margin=margin(20,0,0,0)),
  axis.title.y = element_text(color="black", size=18),
  axis.text.x = element_text(size = 16, margin=margin(5,0,0,0)),
  axis.text.y = element_text(size = 16, margin=margin(0,5,0,10))
  # legend.text = element_text(size = 16),
  # axis.title.y = element_text(color="black", size=18),
  # axis.title.x = element_text(color="black", size=18, margin=margin(5,0,0,0)),
  
)

plotsfolder = 'analysis/'

names(MB_imp_T0)

MB_imp_T0$Age_Round <- round(MB_imp_T0$Age_Frac_Imp)
MB_imp_T0$Age_Round <- as.factor(MB_imp_T0$Age_Round)


Corrs <- c("w_diff","w_P6","T_Vocab","T_Matrix","Corsi_WM_Span","CogFlex_t",
           "SSRT","FlankerSwitch_t","FlankerInhib_t","Stroop_t",
           "AXCPT_t","OneBack_WM_t","TwoBack_WM_t") # ,"SSRT_SSD_z"

names(MB_imp_T0)
# so training coefficient and online included in the imputation
df <- MB_imp_T0[, (names(MB_imp_T0) %in% Corrs)]
#df <- MB_imp_T0[, (!names(MB_imp_T0) %in% Exclude)]
names(df)



typeof(df)

mean(df$FlankerInhib_t, na.rm=TRUE)

names(df)

names(df)[1] <- "Model-based DM"
names(df)[2] <- "Metacontrol"
names(df)[3] <- "WASI_Vocab"
names(df)[4] <- "WASI_Matrix"
names(df)[5] <- "WM_Span"
names(df)[6] <- "AX-CPT"
names(df)[7] <- "CogFlex"#"CogFlex_Switch"
names(df)[8] <- "SSRT"
names(df)[9] <- "Flanker_Switch"
names(df)[10] <- "Flanker_Inhib"
names(df)[11] <- "Stroop"
names(df)[12] <- "WM_1back"
names(df)[13] <- "WM_2back"

#names(df)[14] <- "CogFlex_Mix"

Corr_Dat_Complete <- na.omit(df)

df[] <- lapply(df, function(x) {
  if(is.integer(x)) as.numeric(as.character(x)) else x
})

df[] <- lapply(df, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

df2 <- df[complete.cases(df),]



cormat <- round(cor(Corr_Dat_Complete),2)
head(cormat)


# from this matrix we can see that:
# 'ders' mental health seems to relate to learning
#install.packages('grDevices')
library(grDevices)
library(corrplot)

jpeg(file = paste(plotsfolder,"Corrplot_4Aug2022.jpg"),
     res=300,width = 20, height = 16, units = "cm")
corrplot(cormat,type="upper",method="square",tl.srt=35,tl.cex=0.9,
         tl.col ="black",diag=FALSE,addCoef.col=TRUE)

#corrplot(cormat,type="lower",order="hclust",cl.pos="n",tl.srt=90,tl.cex=0.4)
dev.off()

#install.packages("ggcorrplot")
library(ggcorrplot)

ggcorrplot(cormat,hc.order=FALSE,type="lower",lab=TRUE)



### adjust for multiple comparisons

corrplot(cormat,order="hclust",method ='square',type='upper')

cormat <- round(cor(df),2)

# function to adjust for multiple comparisons
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

## plotting the resulting corrected significance
res1 <- cor.mtest(df, 0.95)
res2 <- cor.mtest(df, 0.99)
corrplot(cormat, p.mat = res1[[1]], sig.level=0.05, insig="blank", cl.align="r", tl.cex=0.6, order="hclust", type="lower", tl.srt=60, cl.ratio=0.1)
# adding this to this function?
pAdj <- p.adjust(c(res1[[1]]), method = "fdr")
resAdj <- matrix(pAdj,ncol=dim(res1[[1]])[1])
corrplot(cormat, p.mat = resAdj, sig.level = 0.05, insig="blank",cl.align="r", 
         tl.cex=0.6, order="hclust", type="lower", tl.srt=60, cl.ratio=0.1)

# adjusted to see p values
corrplot(cormat, p.mat = resAdj, sig.level = -1, insig="p-value",cl.align="r", 
         tl.cex=0.6, order="hclust", type="lower", tl.srt=60, cl.ratio=0.1)

# alternative?
cormat <- round(cor(df),2)
pAdj <- p.adjust(cormat, method = "fdr")
resAdj <- matrix(pAdj,ncol=dim(res1[[1]])[1])
corrplot(cormat, p.mat = resAdj, sig.level = 0.05, insig="blank",cl.align="r", tl.cex=0.6, order="hclust", type="lower", tl.srt=60, cl.ratio=0.1)




################################################################
################ ----- graphs
################################################################

##----------------------------------------------------------------------------------------------------
claire_theme <-   theme(
  plot.title = element_text(family="Arial",color="black", size=28 ,hjust = 0.5,margin=margin(0,0,10,0)),
  text = element_text(family="Arial", size=20),
  legend.title = element_text(size = 26),
  legend.text = element_text(size = 22),
  axis.title.y = element_text(color="black", size=24),
  axis.title.x = element_text(color="black", size=24, margin=margin(5,0,0,0)),
  axis.text.x = element_text(size = 22, margin=margin(5,0,0,0)),
  axis.text.y = element_text(size = 22, margin=margin(0,5,0,10))
)

fig3_theme <- theme(
  text = element_text(family="Arial", size=16),
  plot.title = element_text(color="black", size=20 ,hjust = 0.5,margin=margin(0,0,10,0)),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black"),
  #legend.title = element_blank(),
  axis.title.x = element_text(color="black", size=18, margin=margin(20,0,0,0)),
  axis.title.y = element_text(color="black", size=18),
  axis.text.x = element_text(size = 16, margin=margin(5,0,0,0)),
  axis.text.y = element_text(size = 16, margin=margin(0,5,0,10))
  # legend.text = element_text(size = 16),
  # axis.title.y = element_text(color="black", size=18),
  # axis.title.x = element_text(color="black", size=18, margin=margin(5,0,0,0)),
  
)

fig3_adlt_theme <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y=element_blank(),
  
  #axis.title.x = element_blank(),
  axis.text.x = element_text(size = 16, margin=margin(5,0,0,0)),
  axis.title.x = element_text(color="black", size=18, margin=margin(20,0,0,0)),
  axis.ticks.x = element_blank(),
  
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line.y = element_blank(),
  #axis.line = element_blank(),
)


# ok so to recap, we find some of the same age effects, although the stake effect is now gone (maybe due to the 
# latest changes, e.g. the 4.5 q value initialisation).
# we should probably also look at stay probability and stuff, but first, let's plot some distributions of the data.


#### distribution of w values
ggplot(MB_imp_T0, aes(x = w_P6)) +
  geom_histogram(aes(y=..density..), position="identity", color = 'black', fill = 'lightblue', alpha=0.5,binwidth = 0.05) +
  geom_density(alpha=0.4) +
  scale_fill_discrete(name = "Group") +
  scale_x_continuous(breaks=seq(0,1,0.25),limits=c(0,1.25)) +
  geom_vline(xintercept = mean(MB_imp_T0$w_P6),linetype = 'dashed', color = 'black', size = 1) +
  labs(title="Distribution of w (single)",x="w", y = "density") +
  guides(colour = F, fill = F) +
  theme_light() +
  fig3_theme +
  theme(strip.text.x = element_text(size = 16, colour = "white"))


# Saving plots with fixed dimensions, quality etc.
ggsave("Single_w_distribution_31May2022.png", plot = last_plot(), path = plots_folder,
       scale = 1, width = 20, height = 16, units = "cm",
       dpi = 300)

MB_imp_T0$Age_Round <- round(MB_imp_T0$Age_Frac_Imp)
MB_imp_T0$Age_Round <- as.factor(MB_imp_T0$Age_Round)

# plotting the single w over age
MBplot <- ggplot(data = MB_imp_T0) +
  geom_point(data=MB_imp_T0,aes(x = Age_Frac_Imp, y = w_P6),
             shape = 21,stroke = 1, size = 4) +
  stat_smooth(data = MB_imp_T0, inherit.aes = F, aes(x = Age_Frac_Imp, y = w_P6), 
              color = "#c5c5c5", fill ="#c5c5c5", 
              alpha = 0.5, size = 1, linetype = "solid", method = "lm", 
              formula = y ~ x, se = T) +
  stat_smooth(data = MB_imp_T0, inherit.aes = F, aes(x = Age_Frac_Imp, y = w_P6),
              color = "black", size = 1.25,  
              linetype = "solid", method = "lm", formula = y ~ x, se = F) +
  scale_x_continuous(breaks=seq(6,13,1), lim = c(6,13)) +
  scale_y_continuous(breaks=seq(0,1,0.25), lim = c(0,1)) +
  guides(color = "none", fill = "none", shape = "none") +
  scale_color_brewer(palette = "Blues") +
  scale_fill_brewer(palette = "Blues") +
  ggtitle('Model-based decision making over age') +
  xlab('Age (in years)') +
  ylab('Model-based decision making') +
  theme_light() +
  claire_theme

MBplot

# Saving plots with fixed dimensions, quality etc.
ggsave("W_overAge_BW_V2_4Aug2022.png", plot = last_plot(), path = plots_folder,
       scale = 1, width = 20, height = 16, units = "cm",
       dpi = 300)


MB_imp_T0$Age_Round <- round(MB_imp_T0$Age_Frac_Imp)
MB_imp_T0$Age_Round <- as.factor(MB_imp_T0$Age_Round)

# plotting the single w over age
MBplot <- ggplot(data = MB_imp_T0) +
  geom_point(data=MB_imp_T0,aes(x = Age_Frac_Imp, y = w_diff),
             shape = 21,stroke = 1, size = 4) +
  stat_smooth(data = MB_imp_T0, inherit.aes = F, aes(x = Age_Frac_Imp, y = w_diff), 
              color = "#c5c5c5", fill ="#c5c5c5", 
              alpha = 0.5, size = 1, linetype = "solid", method = "lm", 
              formula = y ~ x, se = T) +
  stat_smooth(data = MB_imp_T0, inherit.aes = F, aes(x = Age_Frac_Imp, y = w_diff),
              color = "black", size = 1.25,  
              linetype = "solid", method = "lm", formula = y ~ x, se = F) +
  scale_x_continuous(breaks=seq(6,13,1), lim = c(6,13)) +
  #scale_y_continuous(breaks=seq(0,1,0.25), lim = c(0,1)) +
  guides(color = "none", fill = "none", shape = "none") +
  #scale_color_brewer(palette = "Blues") +
  #scale_fill_brewer(palette = "Blues") +
  ggtitle('Metacontrol over age') +
  xlab('Age (in years)') +
  ylab('Metacontrol') +
  theme_light() +
  claire_theme

MBplot

# Saving plots with fixed dimensions, quality etc.
ggsave("Metacontrol_overAge_BW_4Aug2022.png", plot = last_plot(), path = plots_folder,
       scale = 1, width = 20, height = 16, units = "cm",
       dpi = 300)


### COMBINED WITH FACET


w_long <- gather(MB_imp_T0, stake, w, w_lo:w_hi, factor_key=TRUE)

Combo_MB <- ggplot(data = w_long, aes(x = stake, y = w, fill = stake)) +
  geom_violin(aes(x = stake, fill = stake),alpha = 0.7) +
  geom_line(aes(group=ID),position=position_dodge(0.25), lwd=0.5,
            linetype="dotted") +
  geom_point(aes(group=ID), color = "black",shape=21,
             stroke = 1, size = 3, alpha = 0.8, position=position_dodge(0.25)) +
  scale_fill_manual(name = "Stakes", values = c("#c5c5c5","#fc8d39"),
                    labels=c("Low","High"), guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values=c(21, 23),labels = c("Children", "Adults")) +
  scale_y_continuous(name="Model-based decision making",breaks=seq(0,1,0.2),lim=c(0,1)) +
  scale_x_discrete(name="",labels=c("Low","High")) +
  ggtitle("Model-based decision making over stakes") +
  guides(fill="none",color="none",shape="none") +
  #facet_grid(. ~ Group) +
  theme_light() +
  claire_theme 


Combo_MB  


# Saving plots with fixed dimensions, quality etc.
ggsave("Combo_Stakes_paired_31May2022.png", plot = last_plot(), path = plots_folder,
       scale = 1, width = 20, height = 16, units = "cm",
       dpi = 300)



###############################################################
### -------------------- MEDIATION ANALYSIS -------------------
###############################################################

Atlas <- read.csv("data/MeanDLPFC_Thickness_Wdiff_Inhib.csv")
names(Atlas)

Atlas$lROI_noage_z <- scale(Atlas$lROI_noage,scale=TRUE,center=TRUE)
Atlas$rROI_noage_z <- scale(Atlas$rROI_noage,scale=TRUE,center=TRUE)


library(lavaan)

model <- '
w_diff ~ c*lROI

flanker ~ a*lROI
w_diff ~ b*flanker

ab := a*b

total := c + (a*b)
'
fit <- sem(model,data = Atlas)
summary(fit)


model <- '
w_diff ~ c*rROI_noage_z + Age_frac_scaled

flanker ~ a*rROI_noage_z + Age_frac_scaled
w_diff ~ b*flanker

ab := a*b

total := c + (a*b)
'
fit <- sem(model,data = Atlas,se = "bootstrap")
summary(fit)


### the following piece of code investigates whether a model with or without covariate fits the data better.

#specify variables in your model:

X<-(Atlas$lROI)
Y<-(Atlas$w_diff)
M<-(Atlas$flanker)
Cov<-(Atlas$Age_frac_scaled) #demeaned age can be used as a covariate if wanted, if age is correlated to the outcome of the EQ/AQ score
#group<-(dataAQ_EQ$Sex) #sex differences will be investigated in following code, first we looked at women and men combines

#mediationdata<-data.frame(X,Y,M,Cov,group)
mediationdata<-data.frame(X,Y,M,Cov) #group

#model without covariate
MediationModel<-
  ' Y~b*M +c*X
M~a*X
ab:=a*b
total:=c+(a*b)'

fitMediationModel=sem(MediationModel, data=mediationdata, se='bootstrap')

# get stats:
summary(fitMediationModel,fit.measures=T, standardized=T, modindices = TRUE)
fitMeasures (fitMediationModel, (c("df","chisq","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","wrmr","aic","bic")))


X<-(Atlas$lROI_noage_z)

#model with covariate
CovMediationModel<-
  ' Y~b*M +c*X
M~a*X
ab:=a*b
total:=c+(a*b)
Y~Cov'

fitCovMediationModel=sem(CovMediationModel, data=mediationdata, se='bootstrap')

# get stats:
summary(fitCovMediationModel,fit.measures=T, standardized=T, modindices = TRUE)
fitMeasures (fitCovMediationModel, (c("df","chisq","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","wrmr","aic","bic")))


# get stats:
summary(fitMediationModel,fit.measures=T, standardized=T, modindices = TRUE)
fitMeasures (fitMediationModel, (c("df","chisq","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","wrmr","aic","bic")))



#fill in est for a,b and c from summary (std.all) above
#variance explained:
coef(fitMediationModel)#only if variables are centered
a=  -0.093 #Estimate for a from summary
b=  -1.08 #Estimate for b from summary
c=  0.058 #Estimate for c from summary
varexplained=(a*b)/((a*b)+c)
varexplained

# get AIC and/or BIC value to test significance (the smaller the value, the better the fit of the model)
#AIC(fitMediationModel)
#BIC(fitMediationModel) 

# when interpreting the output, c is read as c' and the total path is c+a*b
#indirect path is ab in the output
#total is c is direct path, and c' is association between X and Y.



################################################################################
################################################################################
### ---------- GLM analysis
################################################################################
################################################################################

Stay_kids <- read.csv('data/Stay_Prob_Kids_P7.csv')

# initial values:
# same: -1 is different
# same; 1 is same
names(Stay_kids)
range(Stay_kids$prevpoints) # goes from 0 to 9
range(Stay_kids$stay) # from 0 to 1
range(Stay_kids$same) # from -1 to 1
range(Stay_kids$stake) # from 0 to 9
range(Stay_kids$prevrewdiff) # from 0 to 1

library(dplyr)
DF_imp_T0$ID <- as.numeric(MB_imp_T0$ID)
# MERGE IN AGE FOR KIDS to use as predictor (based on the Data_P dataframes from other R script)
Age_kids <- MB_imp_T0 %>%
  select(ID,Age_Frac_Imp)

names(Age_kids)[1] <- "subnr"

# creating combos of children and adults
Stay_kids_age <- merge(Stay_kids,Age_kids)

# this will create z scored variable
Stay_kids_age$Age_Frac_Imp_z <- scale(Stay_kids_age$Age_Frac_Imp, center = TRUE, scale = TRUE)
Stay_kids_age$prevpoints_z <- scale(Stay_kids_age$prevpoints, center = TRUE, scale = TRUE)

#install.packages('boot')
library(boot)
library(lme4)

#install.packages('stargazer')
library(stargazer)
library(ggeffects)
library(AICcmodavg)

# model comparison with AICcmodavg
Cand.mod1 <- list()

# first model - basic measures
Cand.mod1[[1]] <- glmer(stay ~ prevpoints_z * same + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))


# plot this
mydf <- ggpredict(Cand.mod1[[1]], terms = c("prevpoints","same"))
plot(mydf)

# first mo

# model comparison with AICcmodavg
Cand.mod1 <- list()

# first model - basic measures
Cand.mod1[[1]] <- glmer(stay ~ prevpoints_z * same + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# first model - adding in age
Cand.mod1[[2]] <- glmer(stay ~ prevpoints_z * same + Age_Frac_Imp_z + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# first model - adding in age
Cand.mod1[[3]] <- glmer(stay ~ prevpoints_z * same * Age_Frac_Imp_z + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# adding in stake
Cand.mod1[[4]] <- glmer(stay ~ prevpoints_z * same * stake + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# adding in stake and age
Cand.mod1[[5]] <- glmer(stay ~ prevpoints_z * same * stake * Age_Frac_Imp_z + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# adding in stake and age
Cand.mod1[[6]] <- glmer(stay ~ prevpoints_z * same * stake + Age_Frac_Imp_z + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# adding in stake and age
Cand.mod1[[7]] <- glmer(stay ~ prevpoints_z * same * prevrewdiff + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# adding in stake and age
Cand.mod1[[8]] <- glmer(stay ~ prevpoints_z * same * prevrewdiff * stake + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# adding in stake and age
Cand.mod1[[9]] <- glmer(stay ~ prevpoints_z * same * prevrewdiff * stake * Age_Frac_Imp_z + (1| subnr), data = Stay_kids_age, 
                        family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# nulll model
Cand.mod1[[10]] <- glmer(stay ~ 1 + (1| subnr), data = Stay_kids_age, 
                         family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))


# assign names to each model
# Modnames <- c("1prevp*same*prevd*stake","2prevp*same*prevd+stake","3prevp*same+prevd+stake",
#               "4prevp+same+prevd+stake")

# model selection table based on AICc
#aictab(cand.set = Cand.mod1,modnames=Modnames)
aictab(cand.set = Cand.mod1)

# compute evidence ratio - based on raw method (what is raw?)
# confset(cand.set = Cand.mod9, modnames = Modnames, second.ord = TRUE,
#         method = "ordinal")
confset(cand.set = Cand.mod1, second.ord = TRUE,
        method = "ordinal")

# for the kids, model 5 was the best (prevp * same * stake * age)
summary(Cand.mod1[[5]])


# plot this
mydf <- ggpredict(Cand.mod1[[5]], terms = c("prevpoints","same","stake","Age_Frac_Imp_z"))
plot(mydf)



############ ############ ############ ############ ############ ############ 
############ PLOTTING THIS
############ ############ ############ ############ ############ ############ 


## plotting
#install.packages("ggiraphExtra")
#install.packages("ggiraph")
require(ggiraph)
require(ggiraphExtra)
require(plyr)
library(ggeffects)
library(ggplot2)
library(cowplot)


claire_theme <-   theme(
  #legend.position = "none",
  plot.title = element_text(family="Arial",color="black", size=18 ,hjust = 0.5,margin=margin(0,0,10,0)),
  text = element_text(family="Arial", size=16),
  legend.title = element_text(size = 20),
  legend.text = element_text(size = 16),
  axis.title.y = element_text(color="black", size=18,margin=margin(r = 20)),
  axis.title.x = element_text(color="black", size=18, margin=margin(5,0,0,0)),
  axis.text.x = element_text(size = 16, margin=margin(5,0,0,0)),
  axis.text.y = element_text(size = 16, margin=margin(5,5,0,10))
)


### -------------------- win model
# winning model:
# adding in stake and age
# Cand.mod1[[5]] <- glmer(stay ~ prevpoints * same * stake * Age_Frac_Imp_z + (1| subnr), data = Stay_kids_age, 
#                         family = binomial(link="logit"), control = glmerControl(optimizer = "bobyqa"))

# plot with the full model --------------------- x -------group - facet -- panel
mydf <- ggpredict(Cand.mod1[[5]], terms = c("prevpoints_z","same","stake","Age_Frac_Imp_z"))
plot(mydf)
names(mydf)

levels(mydf$facet) <- c("Low Stakes","High Stakes")
levels(mydf$panel) <- c("Youngest","Middle","Oldest")

p1df <- mydf[which (mydf$panel == -1),]
p2df <- mydf[which (mydf$panel == 0),]
p3df <- mydf[which (mydf$panel == 1),]

names(p1df)
plot(p1df)

p1 <- ggplot(p1df, aes(x = x, y = predicted, color= group)) +
  geom_line(lwd=2) +
  facet_wrap(~facet) +
  scale_y_continuous(name = "Predicted Stay Probability",breaks=seq(0.25,0.75,0.25),lim = c(0.25,0.85)) +
  scale_x_continuous(name = "Previous reward", breaks=seq(-2,2,1), lim = c(-2,2)) +
  scale_color_discrete(name = "Starting State",labels = c("Different","Same"),
                       guide = guide_legend(reverse=TRUE)) +
  ggtitle("Youngest Children") +
  guides(color="none") +
  theme_light() +
  claire_theme +
  theme(
    axis.title.y = element_blank(),
    plot.margin=margin(0,0,0,2,"cm")
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank()
  )

p1

p2 <- ggplot(p2df, aes(x = x, y = predicted, color= group)) +
  geom_line(lwd=2) +
  facet_wrap(~facet) +
  scale_y_continuous(name = "Predicted stay probability",breaks=seq(0.25,0.75,0.25),lim = c(0.25,0.85)) +
  scale_x_continuous(name = "Previous reward", breaks=seq(-2,2,1), lim = c(-2,2)) +
  scale_color_discrete(name = "Starting State",labels = c("Different","Same"),
                       guide = guide_legend(reverse=TRUE)) +
  ggtitle("Middle Children") +
  guides(color="none") +
  theme_light() +
  claire_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p2

p3 <- ggplot(p3df, aes(x = x, y = predicted, color= group)) +
  geom_line(lwd=2) +
  facet_wrap(~facet) +
  scale_y_continuous(name = "Predicted Stay Probability",breaks=seq(0.25,0.75,0.25),lim = c(0.25,0.85)) +
  scale_x_continuous(name = "Previous reward", breaks=seq(-2,2,1), lim = c(-2,2)) +
  scale_color_discrete(name = "Starting State",labels = c("Different","Same"),
                       guide = guide_legend(reverse=TRUE)) +
  ggtitle("Oldest Children") +
  guides(color="none") +
  theme_light() +
  claire_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p3

# # Saving plots with fixed dimensions, quality etc.
# ggsave("Regression_P3_Legend.png", plot = last_plot(), path = plots_folder,
#        scale = 1, width = 26, height = 16, units = "cm",
#        dpi = 300)

###
title <- ggdraw() + 
  draw_label(
    "Behavioral Markers of Model-based Decision Making across Age and Stake",
    #"Model-based decision-making for children over\nage and adults plotted separately",
    fontfamily = "Arial",
    fontface = "plain",
    size = 22,
    hjust = 0.5
  )

gg <- plot_grid(p1,p2,p3,nrow=1, align = "h", rel_widths = c(4,3,3))

gg

g_title1 <- plot_grid(title,gg, nrow=2, rel_heights = c(1,9))

g_title1




# Saving plots with fixed dimensions, quality etc.
ggsave("Kids_Regression_6Jul22.png", plot = last_plot(), path = plots_folder,
       scale = 1, width = 26, height = 16, units = "cm",
       dpi = 300)



#### calculate better performance
names(Stay_kids_age)

unique(Stay_kids_age$stake)

Stay_kids_age$stake_val <- Stay_kids_age$stake

Stay_kids_age$stake_val[Stay_kids_age$stake_val == -1] <- 0
Stay_kids_age$stake_val[Stay_kids_age$stake_val == 1] <- 5

unique(Stay_kids_age$stake_val)

Stay_kids_age$real_points <- Stay_kids_age$prevpoints * Stay_kids_age$stake_val

real_rews <- summarySE(data = Stay_kids_age, measurevar = "real_points", groupvars = c("subnr"))
names(real_rews)[1] <- "ID"

# merge real rews into the data
Data_rews <- left_join(DF_imp_T0,real_rews,by=c("ID"))

Data_rews$Age_Round <- round(Data_rews$Age_Frac_Imp)
Data_rews$Age_Round <- as.factor(Data_rews$Age_Round)

cor.test(Data_rews$real_points,Data_rews$Age_Frac_Imp)

cor.test(Data_rews$Avg_Pts,Data_rews$Age_Frac_Imp)

Data_rews$rews_2 <- Data_rews$Avg_Pts_hi * 5

Data_rews$Mean_Rews = rowMeans(Data_rews[,c("rews_2","Avg_Pts_lo")])

t.test(Data_rews$Avg_Pts, mu = 0, alternative ="two.sided")
cohens_d(Data_rews$Avg_Pts, mu = 0)

Mean = mean(Data_rews$Avg_Pts)
Mu = 0
Diff = Mean - Mu
SD = sd(Data_rews$Avg_Pts)
CohenD = (Mean-Mu) / SD
CohenD



