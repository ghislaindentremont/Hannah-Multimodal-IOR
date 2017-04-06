library(ggplot2)
library(plyr)
library(reshape)
library(scales)
library(ggthemes)

setwd("/Users/ghislaindentremont/Documents/Experiments/Multimodal_IOR/Hannah/P_analysis")


##################################################################
####                  Grand Averages                          ####
##################################################################

b_list2 = list.files(
  path = "."
  , pattern = "grand_averages"
  , recursive = T
)
b_list = b_list2[grep(".csv", b_list2)]

b_data = b_list
b = ldply(
  .data = b_data
  , .fun = function(file) {
    df = read.csv(
      file
      , header = TRUE
    )
    P_col = data.frame(id = substr(file,1,3))
    to_return = cbind(P_col, df)
    return(to_return)
  }
)

# copy
b2 = b 

# melt 
b2 = melt(b2[,-c(2:3)])

# add time variable
b2$X = -200:500

# average across participants
b2_agg = aggregate(value ~ variable + X, data = b2, FUN = mean)

# add cueing/cue/target/laterality columns
b2_agg$ctoa = ifelse(
  substr(b2_agg$variable,4,4) == "0"
  , "1000"
  , "700"  
)

# convert to factor
b2_agg$ctoa = factor(b2_agg$ctoa)

# change names
levels(b2_agg$ctoa) = c('1000 ms', '700 ms')

gg = ggplot(
  b2_agg
  , aes(x = X, y = value*1e6)
) +
  geom_line(size = 1) +
  facet_grid(. ~ ctoa, scales = "free")+
  geom_vline(
    xintercept =  0
    , linetype = 2
  )+ 
  geom_hline(
    yintercept = 0
  )+
  scale_x_continuous("Time (ms)")+
  scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
  theme_gray(base_size = 20)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)

##################################################################
####                      For Presentation                    ####
##################################################################

b2_agg_avg = aggregate(value ~ X, data = b2_agg, FUN = mean)

plot_grands = function(dat, vline1, vline2, vline3, vline4) {
  gg = ggplot(
    dat
    , aes(x = X, y = value*1e6)
  ) +
    geom_line(size = 1) +
    geom_vline(
      xintercept =  0
      , linetype = 2
    )+
    geom_hline(
      yintercept = 0
    )+
    geom_vline(xintercept = vline1, color = "red", size = 1)+
    geom_vline(xintercept = vline2, color = "red", size = 1)+
    geom_vline(xintercept = vline3, color = "red", size = 1)+
    geom_vline(xintercept = vline4, color = "red", size = 1)+
    scale_x_continuous("Time (ms)", limits = c(-100, 400))+
    scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
    theme_gray(base_size = 20)+
    theme(panel.grid.major = element_line(size = 1.5)
          ,panel.grid.minor = element_line(size = 1))

  return(gg)
}

# set bounds
P1lo = 70
P1hi = 170
N1lo = 170
N1hi = 210
plot_grands(b2_agg_avg ,P1lo, P1hi, N1lo, N1hi)



##################################################################
####                 Condition Averages                       ####
##################################################################

a_list2 = list.files(
path = "."
, pattern = "condition_averages"
, recursive = T
)

a_list = a_list2[grep(".csv", a_list2)]

a_data = a_list

a = ldply(
.data = a_data
, .fun = function(file) {
    df = read.csv(
    file
    , header = TRUE
    )
    P_col = data.frame(id = substr(file,1,3), stringsAsFactors=FALSE)
    to_return = cbind(P_col, df)
  return(to_return)
}
)

# batch 1 df
a2 = a

# melt
a2 = melt(a2[,-c(2:3)])

# add time variable
a2$time = -200:500

# add cueing/laterality/cue modality/CTOA columns
a2$cue_modality = ifelse(
  substr(a2$variable,3,3) == "T"
  , "tactile"
  , "visual"
)
a2$ctoa = ifelse(
  substr(a2$variable,4,4) == "7"
  , "700"
  , "1000"
)
a2$laterality = ifelse(
  substr(a2$variable,2,2) == "C"
  , "contra"
  , "ipsi"
)
a2$cueing = ifelse(
  substr(a2$variable,1,1) == "C"
  , "cued"
  , "uncued"
)

# average across participants
a2_agg = aggregate(value ~ variable + time + cue_modality + ctoa + laterality + cueing, data = a2, FUN = mean)

# change factor order (note: this does not change DF)
a2_agg$cueing = factor(a2_agg$cueing, c("uncued", "cued"))

# convert to factor
a2_agg$cue_modality = factor(a2_agg$cue_modality)
a2_agg$ctoa = factor(a2_agg$ctoa)
a2_agg$laterality = factor(a2_agg$laterality)

# change names
levels(a2_agg$cue_modality) = c('Tactile Cue', 'Visual Cue')
levels(a2_agg$ctoa) = c('1000 ms', '700 ms')
levels(a2_agg$laterality) = c("Contralateral", "Ipsilateral")
levels(a2_agg$cueing) = c("Uncued", "Cued")

gg = ggplot(
  a2_agg
  , aes(x = time, y = value*1e6, group = cueing, color = cueing)
) +
  geom_line(size = 1) +
  facet_grid(laterality ~ cue_modality + ctoa)+
  geom_vline(
    xintercept =  0
    , linetype = 2
  )+
  geom_hline(
    yintercept = 0
  )+
  scale_x_continuous("Time (ms)", limits = c(-100, 400))+
  scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
  labs(color = "Cueing")+
  theme_gray(base_size = 16)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))

print(gg)





##################################################################
####                     Analysis                             ####
##################################################################

lower_bound = P1lo
upper_bound = P1hi
component = "P1"
  
eeg_comp = a2[a2$time >= lower_bound & a2$time <= upper_bound, ] 

# get average voltage in window for each condition (8) x cueing (2) x participant (n)
eeg_comp_agg = aggregate(value ~ ctoa + cue_modality + laterality + cueing + id, data = eeg_comp, FUN = mean)

# get means
eeg_comp_agg_means = aggregate(value ~ ctoa + cue_modality + laterality + cueing, data = eeg_comp_agg, FUN = mean)

# get sds
eeg_comp_agg_sds = aggregate(value ~ ctoa + cue_modality + laterality + cueing, data = eeg_comp_agg, FUN = sd)

# get cueing effects
cueing_effects = aggregate(value ~ ctoa + cue_modality + laterality + id, data = eeg_comp_agg, FUN = diff)
avg_cueing_effects = aggregate (value ~ ctoa + cue_modality + laterality, data = cueing_effects, mean)

# data frame of summary stats
eeg_comp_agg_sum = cbind(eeg_comp_agg_means, eeg_comp_agg_sds$value)
names(eeg_comp_agg_sum)[5:6] = c("M", "SD")
print(eeg_comp_agg_sum)

# ANOVA
m = aov(value ~
          cue_modality*laterality*cueing*ctoa
        + Error(id/(cue_modality*laterality*cueing*ctoa))
        , data = eeg_comp_agg)
m_summary = summary(m)
print(m_summary)



#----------------------------- Two-way interaction ------------------------------------#
eeg_comp_agg_cueing2 = aggregate(value ~ laterality + cue_modality + ctoa + id, data = eeg_comp_agg, FUN = diff)
eeg_comp_agg_cueing = aggregate(value ~ cue_modality + id, data = eeg_comp_agg_cueing2, FUN = mean)

# calculate paired differences (from zero!)
mods = unique(eeg_comp_agg_cueing$cue_modality)
MEs = NULL
alpha = 0.05
for (mod in mods) {
  dat = eeg_comp_agg_cueing[eeg_comp_agg_cueing$cue_modality == mod,]$value

  MEAN = mean(dat)
  SD = sd(dat)
  N = length(dat)
  DF = N - 1
  T_CRIT = qt(1-alpha/2, DF)

  SE = SD/sqrt(N)

  ME = T_CRIT * SE
  MEs = c(MEs, ME)
}

eeg_comp_agg_cueing$mix_factor = paste(eeg_comp_agg_cueing$cue_modality)
# model
m3 = aov(value ~ mix_factor
         + Error(id/mix_factor)
         , data = eeg_comp_agg_cueing)
m3_summary = summary(m3)

MSE3 = m3_summary$`Error: id:mix_factor`[1][[1]][[3]][2]
df3 = m3_summary$`Error: id:mix_factor`[1][[1]][[1]][2]

eeg_comp_agg_2way = aggregate(value ~ id, data = eeg_comp_agg_cueing, FUN = diff)
twoway_var = var(eeg_comp_agg_2way$value)
(MSE3*2 - twoway_var) < 1.0e-8

SEM_LM3 = sqrt(MSE3/8)

# get critical t value and LSD
t_crit3 = qt(1-alpha/2, df3)
ME_LM3 = abs(t_crit3 * SEM_LM3)

EEG_CIs2 = aggregate(value ~ cue_modality + cueing, data = eeg_comp_agg_means, FUN = mean)
EEG_CIs = aggregate(value ~ cue_modality, data = EEG_CIs2, FUN = diff)
EEG_CIs$CI = ME_LM3  # use margin of error !

names(EEG_CIs)[2] = "M"

EEG_CIs$cue_modality = as.factor(EEG_CIs$cue_modality)
levels(EEG_CIs$cue_modality) = c("Tactile Cue", "Visual Cue")

# add pair-wise
EEG_CIs$from_zero = MEs

# plot!
ggplot(EEG_CIs, aes(x = cue_modality, y = M*1e6))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = (M - CI)*1e6, ymax = (M + CI)*1e6)
                , width = 0.1
                , size = 1
  )+
  geom_errorbar(aes(ymin = (M - from_zero)*1e6, ymax = (M + from_zero)*1e6)
                , width = 0.1
                , size = 1
                , color = "red"
  )+
  labs(x = "Cue Modality",y = paste(component, " Reduction (", expression(u),"V)", sep = ""), color = "Laterality")+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))
#----------------------------- Two-way interaction ------------------------------------#



#----------------------------- Four-way interaction ----------------------------------#
twoways = aggregate(value~laterality + ctoa + id, data=eeg_comp_agg_cueing2, FUN=diff)
twoways_meanz = aggregate(value~laterality + ctoa, data=twoways, FUN=mean)
twoways_varz = aggregate(value~laterality + ctoa, data=twoways, FUN=var)

ciw = sqrt(twoways_varz$value/8) * qt(1-alpha/2, 8-1)

twoways_meanz$CI = ciw

names(twoways_meanz)[3:4] = c("M", "CI")

twoways_meanz$laterality = as.factor(twoways_meanz$laterality)
levels(twoways_meanz$laterality) = c("Contralateral", "Ipsilateral")

twoways_meanz$ctoa = as.factor(twoways_meanz$ctoa)
levels(twoways_meanz$ctoa) = c("1000 ms", "700 ms")

ggplot(data=twoways_meanz, aes(x=ctoa, y=M*1e6, color = laterality))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = (M - CI)*1e6, ymax = (M + CI)*1e6)
                , width = 0.1
                , size = 1
                # , color = "red"
  )+
  labs(x = "Cue-Target Onset Asynchrony",y = paste(component, ": Visual - Tactile Reduction (", expression(u),"V)", sep = ""), color = "Laterality")+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))
#----------------------------- Four-way interaction ----------------------------------#
