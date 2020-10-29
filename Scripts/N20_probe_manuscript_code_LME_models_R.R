### Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses ###
# Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

#### Linear-mixed-effects models on relation between pre-stimulus alpha and early SEP amplitudes ####
  
# load packages
library(R.matlab)
library(lme4)
library(lmerTest)

# load data
data_all <- readMat("/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Data/LME_prestim_peak_amp_MANUSCRIPT.mat")

subj_id <- c(data_all$subj.id)
subj_id <- as.factor(subj_id)
prestim_emp <- c(data_all$alpha.prestim.tang.concat)
peak_emp <- c(data_all$peak.amp.N20.concat)

# random-slope LME
mod <- lmer(peak_emp ~ prestim_emp + (prestim_emp|subj_id))
summary(mod)





