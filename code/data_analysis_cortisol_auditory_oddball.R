# ---- READS preprocessed auditory oddball data and merges to cortisol data

# PRELIMINARY FINDINGS (May2024)
# - higher rpd_low associated with attenuated saliva cortisol reduction
# - higher saliva cortisol before is associated with lower BPS (rpd_low)
# - higher hair cortisol is associated with higher rpd responses to all trials before manipulation
# - higher saliva cortisol before experiment differentiates pupillary responses

# - self reported internalizing with lower hair cortisol
# - self reported externalizing with higher hair cortisol
# - self reported externalizing with lower saliva cortisol change (pre-püost change)


# PRELIMINARY FINDINGS (September 2024)
# - higher hair cortisol predicts the average salivary cortisol level (b = 0.98)
# - higher saliva cortisol before is associated with lower BPS
# - higher saliva cortisol differentiates pupillary response to oddballs
# - higher BPS assoicated with attenuated cortisol decrease

### TODO: control for BMI and time in models

# SETUP ####

# PATHS
if (Sys.info()["sysname"] == "Linux") {
  home_path <- "~"
  project_path <- "/PowerFolders/project_sega"
  data_path <- "/PowerFolders/project_sega/data/AuditoryOddball"
  data_path_eeg <- "/PowerFolders/project_sega/data/AuditoryOddball_EEG"
  datapath <- paste0(home_path, data_path) # .csv + .hdf5 input files
  datapath_eeg <- paste0(home_path, data_path_eeg) # .txt input files (eeg)
  # List all .hdf and .csv files
  data_files <- list.files(path = datapath, full.names = TRUE)
}

if (Sys.info()["sysname"] == "Windows") {
  home_path <- "C:/Users/Nico"
  project_path <- "/PowerFolders/project_sega"
  data_path <- "/PowerFolders/project_sega/data/et auditory oddball"
  data_path_task<-"/PowerFolders/project_sega/data/task auditory oddball"
  data_path_eeg <- "/PowerFolders/project_sega/data/tests/eeg preprocessed auditory oddball"
  datapath <- paste0(home_path, data_path) # hdf5 input files
  datapath_task <- paste0(home_path, data_path_task) # hdf5 input files
  datapath_eeg <- paste0(home_path, data_path_eeg) # .txt input files (eeg)
  # List all .hdf and .csv files
  data_files <- c(list.files(path = datapath, full.names = TRUE),
                  list.files(path = datapath_task, full.names = TRUE))
}

if (Sys.info()["sysname"] == "Darwin") {
  home_path <- "~"
  project_path <- "/code"
  data_path <- "/code/input/AuditoryOddball/eyetracking"
  data_path_task <- "/code/input/AuditoryOddball/taskdata"
  data_path_eeg <- "/code/input/AuditoryOddball/eeg"
  datapath <- paste0(home_path, data_path) # .hdf5 input files
  datapath_task <- paste0(home_path, data_path_task) # .csv input files
  datapath_eeg <- paste0(home_path, data_path_eeg) # .txt input files (eeg)
  # List all .hdf and .csv files
  data_files <- c(
    list.files(path = datapath, full.names = TRUE),
    list.files(path = datapath_task, full.names = TRUE))
}


# Can be used to skip preprocessing and directly read proprocessed data from .rds file:
df_trial <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata_24092024.rds'))

#changed random intercept to a factor
df_trial$id<-as.factor(df_trial$id) #change ID to factor
df_trial<-df_trial[df_trial$manipulation %in% c('before','after'),]
df_trial$block<-as.factor(df_trial$block)
df_trial$manipulation<-droplevels(df_trial$manipulation)


require(performance)
require(ggplot2)
require(lme4)
require(lmerTest)
require(emmeans)


# CORISOL PREPROCESSING ####

### - load hair and saliva cortisol data from xlsx ####

require(readxl)
#hair cortisol
df_hair<-read_xlsx(path="C:/Users/nico/PowerFolders/project_sega/data/hair_cortisol_1439_23082023.xlsx",
                   skip=11) #starting reading in line 12

df_hair2<-read_xlsx(path="C:/Users/nico/PowerFolders/project_sega/data/hair_cortisol_1548_08052024.xlsx",
                    skip=14) #starting reading in line 15


# #in file "hair_cortisol_1548_08052024" use label on foil as reference, is same as in our sheets:
# df_hair2_labels<-read_xlsx(path="C:/Users/nico/PowerFolders/project_sega/data/hair_labelling_second_batch.xlsx")
# cbind(df_hair2$`ID on the foil`,df_hair2_labels$`Patient site and ID`)

#remove wrong labels
df_hair2<-df_hair2[,-c(1,2)]

#change variable to appropriate length
df_hair$`externe ID`<-substr(df_hair$`externe ID`,nchar(df_hair$`externe ID`)-7,nchar(df_hair$`externe ID`))

#change names
names(df_hair)<-c('id','hair_mass','hair_comments','hair_cortisol','hair_cortisone')
names(df_hair2)<-c('id','hair_mass','hair_comments','hair_cortisol','hair_cortisone')

df_hair<-rbind(df_hair,df_hair2)  


#saliva cortisol
df_saliva<-read_xlsx(path="C:/Users/nico/PowerFolders/project_sega/data/saliva_cortisol_20230912_Results_SEGA_Salimetrics.xlsx",
                     skip=9) #starting reading in line 10

df_saliva2<-read_xlsx(path="C:/Users/nico/PowerFolders/project_sega/data/saliva_cortisol_20240521_Results2_SEGA_Salimetrics.xlsx",
                      skip=12) #starting reading in line 13

names(df_saliva)<-c('saliva_sample_id','id','saliva_date','saliva_time',
                    'saliva_cortisol_1','saliva_cortisol_2','saliva_cortisol_avg','saliva_cv','saliva_comments')
names(df_saliva2)<-c('saliva_sample_id','id','saliva_date','saliva_time',
                     'saliva_cortisol_1','saliva_cortisol_2','saliva_cortisol_avg','saliva_cv','saliva_comments')

df_saliva<-rbind(df_saliva,df_saliva2)

#remove missing
df_saliva<-df_saliva[!(grepl('65',df_saliva$saliva_sample_id)),]

#categorize sample
df_saliva$saliva_when<-rep(c('before','after'),nrow(df_saliva)/2)

#get time as numeric variable
df_saliva$saliva_time_num<-as.numeric(format(df_saliva$saliva_time, "%H"))


### - create SEGA id variable to match df_hair and df_saliva ####

sega_id<-ifelse(nchar(levels(df_trial$id))==2,
                paste0('0',levels(df_trial$id)),
                ifelse(nchar(levels(df_trial$id))==1,
                       paste0('00',levels(df_trial$id)),levels(df_trial$id)))

sega_id<-paste0('SEGA_',sega_id)

levels(df_trial$id)<-sega_id


### - df_agg: estimate pupillometry measures from trial data #####
df_agg<-aggregate(df_trial[,names(df_trial) %in% c('rpd_auc','rpd_high','rpd_low','rpd')],
                  by=with(df_trial,list(id,trial,manipulation,block)),
                  mean,na.rm=T)

names(df_agg)<-c('id','oddball','manipulation','reverse','rpd_auc','rpd_high','rpd_low','rpd')


### - saliva change ####

saliva_cortisol_change<-with(df_saliva,saliva_cortisol_avg[saliva_when=='after']-
                               saliva_cortisol_avg[saliva_when=='before'])

id<-unique(df_saliva$id)

df_saliva_change<-data.frame(id,saliva_cortisol_change)


### - merge data ####
df_agg_hair<-merge(df_agg,df_hair,by='id')
df_agg_hair<-merge(df_agg_hair,df_saliva_change,by='id',all.x=T)


df_agg_saliva<-merge(df_agg,df_saliva,by='id')
df_agg<-merge(df_agg_saliva,df_hair,by='id')
df_agg<-merge(df_agg,df_saliva_change,by='id')


### - CHECK data distribution ####

hist(df_hair$hair_cortisol)

  plot(df_hair$hair_cortisol,df_hair$hair_mass)
  ###--> remove participants with very low hair mass

ggplot(df_hair,aes(hair_cortisone,hair_cortisol,color=hair_mass))+geom_point()+
  labs(x='hair cortisone (pg/mg, log-scaled)', y='hair cortisol (pg/mg, log-scaled)', color='sample mass (mg)')+
  geom_smooth(color='grey',method='lm')+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_bw()

hist(df_hair$hair_cortisol[df_hair$hair_cortisol<5])

hist(df_saliva$saliva_cortisol_avg)
hist(df_saliva$saliva_cortisol_avg[df_saliva$saliva_when=='before'])
hist(df_saliva$saliva_cortisol_avg[df_saliva$saliva_when=='after'])
hist(df_saliva$saliva_cv)

ggplot(df_saliva,aes(saliva_cortisol_1,saliva_cortisol_2,color=scale(saliva_cv)))+geom_point()+
  labs(x='saliva sample A (nmol/l, log-scaled)', y='saliva sample B (nmol/l, log-scaled)', color='variation index (z)')+
  geom_smooth(color='grey',method='lm')+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_bw()

# DATA ANALYSIS ####

length(unique(df_agg$id))
names(df_agg)
require(performance)

## - association of hair variables ####
summary(lm(scale(hair_cortisol)~scale(hair_mass)+scale(hair_cortisone),df_agg[!duplicated(df_agg$id),]))
###--> hair cortisone and cortisol highly correlated (0.98)

#association of saliva variables
lmm<-lmer(scale(saliva_cortisol_avg)~saliva_when*saliva_time_num+(1|id),df_agg)
summary(lmm)
r2_nakagawa(lmm)


confint(contrast(emmeans(lmm,~saliva_when),'pairwise'))
###--> saliva cortisol reduction after experiment
contrast(emtrends(lmm,~saliva_when,var = 'saliva_time_num'),'pairwise')
###--> time of day is more strongly associated with cortisol after experiment

##effect of experiment
ggplot(df_agg,aes(x=factor(saliva_when,level=c('before','after')),log(saliva_cortisol_avg),fill=saliva_when))+
  geom_violin(adjust=2,alpha=0.5)+geom_boxplot(notch=T,width=0.4)+
  labs(x='saliva sample',y='cortisol (nmol/l, log-scaled)')+
  theme_bw()

ggplot(df_agg,aes(x=factor(saliva_when,level=c('before','after')),log(saliva_cortisol_avg),group=id))+
  geom_point()+geom_line()+
  labs(x='saliva sample',y='cortisol (nmol/l), log-scaled')+
  theme_bw()


## - association of hair and saliva ####
lmm<-lmer(scale(saliva_cortisol_avg)~scale(hair_cortisol)+(1|id),df_agg)
anova(lmm)

lmm<-lmer(scale(saliva_cortisol_avg)~saliva_when*saliva_time_num*scale(hair_cortisol)+(1|id),df_agg)
anova(lmm)
fixef(lmm)
r2_nakagawa(lmm)
#---> higher hair cortisol predicts the average salivary cortisol level

hist(df_agg$saliva_time_num)
ggplot(df_agg[df_agg$saliva_time_num<12,],aes(hair_cortisol,saliva_cortisol_avg))+geom_point()+geom_smooth(method='lm')
ggplot(df_agg[df_agg$saliva_time_num>16,],aes(hair_cortisol,saliva_cortisol_avg))+geom_point()+geom_smooth(method='lm')
###--> to few data points to draw strong conclusion

# ## - asocaition of pupillometric and saliva and hair ###
# lm<-lm(scale(rpd_low)~scale(hair_cortisol),
#        df_agg_hair[!(duplicated(df_agg_hair$id)),])
# summary(lm)
# 
# lm<-lm(scale(rpd_low)~scale(saliva_cortisol_avg)*saliva_when,
#        df_agg[!(duplicated(interaction(df_agg$id,df_agg$saliva_when))),])
# summary(lm)
# ###--> pd_baseline associated with saliva cortisol

### - SEPR and hair and saliva cortisol ####

df_agg$hair_cortisol_z<-scale(df_agg$hair_cortisol) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd)~(oddball*manipulation*reverse)*hair_cortisol_z+(1|id),df_agg)
anova(lmm)
r2_nakagawa(lmm)

# contrast(emmeans(lmm,~oddball),'pairwise')
# ###--> substantial oddball effect
emtrends(lmm,~manipulation,var = 'hair_cortisol_z')
confint(contrast(emtrends(lmm,~manipulation,var = 'hair_cortisol_z'),'pairwise'))
###--> higher hair cortisol is associated with higher rpd responses to all trials before manipulation

ggplot(df_agg[df_agg$manipulation=='before' &
                df_agg$hair_cortisol<15,],aes(hair_cortisol,rpd))+geom_point()+geom_smooth(method='lm')+
  labs(x='hair corisol (pg/mg)',y='pupillary response (z)')+
  theme_bw()
###--> rather driven by outliers

df_agg$saliva_cortisol_z<-scale(df_agg$saliva_cortisol_avg) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd)~(oddball*manipulation)*saliva_cortisol_z+saliva_cv+(1|id),
          df_agg[df_agg$reverse=='forward',])
anova(lmm)
r2_nakagawa(lmm)

plot(emtrends(lmm,~oddball|manipulation,var='saliva_cortisol_z'))
confint(emtrends(lmm,~oddball|manipulation,var='saliva_cortisol_z'))

ggplot(df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward' &
                df_agg$saliva_cortisol_avg<10,],
       aes(x=saliva_cortisol_avg,y=rpd,color=oddball))+geom_point()+geom_smooth(method='lm')+
  labs(x='salivary cortisol (nmol/l)',y='pupillary repsonse (mm)')+
  theme_bw()
#--> higher saliva cortisol differentiates pupillary responses

ggplot(df_agg[df_agg$manipulation=='after' & df_agg$reverse=='forward' &
                df_agg$saliva_cortisol_avg<10,],
       aes(x=saliva_cortisol_avg,y=rpd,color=oddball))+geom_point()+geom_smooth(method='lm')+
  labs(x='saliva cortisol (nmol/l)',y='pupillary repsonse (mm)')+
  coord_cartesian(ylim=c(0,0.08))+
  theme_bw()
###--> higher saliva cortisol might attenuate the oddball effect after the manipulation

### - BPS and hair and saliva cortisol ####  

lmm<-lmer(scale(rpd_low)~(oddball*manipulation*reverse)*hair_cortisol_z+(1|id),df_agg)
anova(lmm)
r2_nakagawa(lmm)

contrast(emmeans(lmm,~reverse),'pairwise')
emtrends(lmm,~reverse,var='hair_cortisol_z')
confint(contrast(emtrends(lmm,~reverse,var='hair_cortisol_z'),'pairwise'))
###--> higher hair cortisol leads to higher BPS in reverse condition


lmm<-lmer(scale(rpd_low)~(oddball*manipulation*reverse)*saliva_cortisol_z+(1|id),
          df_agg[df_agg$saliva_when=='before',])
anova(lmm)
r2_nakagawa(lmm)

cbind(round(fixef(lmm)['saliva_cortisol_z'],2),
      round(confint(lmm,parm = 'saliva_cortisol_z'),2))
#--> higher saliva cortisol before is associated with lower BPS

lmm<-lmer(scale(rpd_low)~(oddball*manipulation*reverse)*saliva_cortisol_z+(1|id),
          df_agg[df_agg$saliva_when=='after',])
anova(lmm)
r2_nakagawa(lmm)

cbind(round(fixef(lmm)['saliva_cortisol_z'],2),
      round(confint(lmm,parm = 'saliva_cortisol_z'),2))
#--> higher saliva cortisol after is associated with lower BPS


ggplot(df_agg[df_agg$oddball=='oddball' & df_agg$saliva_cortisol_z<3,],
       aes(x=saliva_cortisol_z,y=rpd_low))+geom_point()+geom_smooth(method='lm')+
  labs(x='saliva cortisol (z)',y='pupil size (mm)')+
  theme_bw()


### - cortisol change ####

hist(df_agg$saliva_cortisol_change)
lm<-lm(scale(saliva_cortisol_change)~(manipulation*reverse)*scale(rpd_low),df_agg[!duplicated(df_agg$id),])
anova(lm)
###---> higher baseline PD is marginally associated with less reduction in saliva cortisol 

summary(lm(scale(saliva_cortisol_change)~scale(rpd_low)+scale(hair_cortisol),
           df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward',]))
summary(lm(scale(saliva_cortisol_change)~scale(rpd_low)+scale(hair_cortisol),
           df_agg[df_agg$manipulation=='after' & df_agg$reverse=='forward',]))
summary(lm(scale(saliva_cortisol_change)~scale(rpd_low)+scale(hair_cortisol),
           df_agg[df_agg$manipulation=='after' & df_agg$reverse=='reverse',]))
##--> higher BPS associated with attenuated cortisol decrease

    ggplot(df_agg[!duplicated(df_agg$id),],
           aes(x=rpd_low,y=saliva_cortisol_change))+geom_point()+geom_smooth(method='lm')+
      labs(x='pupil size (mm)',y='cortisol change (nmol/l)')+
      theme_bw()

summary(lmer(scale(hair_cortisol)~scale(rpd)*oddball*manipulation*reverse+(1|id),df_agg))

## --> analysis in ERC Grant 2024 ####

##assocaition of hair cortisol and pupillary response
df_agg$hair_cortisol_z<-scale(df_agg$hair_cortisol) #z standardize before to use emtrends later

lmm<-lmer(scale(rpd)~(oddball*manipulation*reverse)*hair_cortisol_z+(1|id),df_agg)

anova(lmm)
r2_nakagawa(lmm)

contrast(emmeans(lmm,~oddball),'pairwise')
###--> substantial oddball effect
emtrends(lmm,~manipulation,var = 'hair_cortisol_z')
confint(contrast(emtrends(lmm,~manipulation,var = 'hair_cortisol_z'),'pairwise'))

#association of baseline pupil size and saliva cortisol change 
lm<-lm(scale(saliva_cortisol_change)~scale(rpd_low),df_agg[!duplicated(df_agg$id),])
summary(lm)

#association of saliva cortisol before the task and baseline pupil size
df_agg$saliva_cortisol_z<-scale(df_agg$saliva_cortisol_avg) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd_low)~(oddball*manipulation*reverse)*saliva_cortisol_z+(1|id),
          df_agg[df_agg$saliva_when=='before',])
anova(lmm)
r2_nakagawa(lmm)

cbind(round(fixef(lmm)['saliva_cortisol_z'],2),
      round(confint(lmm,parm = 'saliva_cortisol_z'),2))
#--> higher saliva cortisol before is associated with lower BPS

#association of saliva cortisol and pupillary response
df_agg$saliva_cortisol_z<-scale(df_agg$saliva_cortisol_avg) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd)~(oddball*manipulation*reverse)*saliva_cortisol_z+(1|id),
          df_agg)
anova(lmm)
r2_nakagawa(lmm)

plot(emtrends(lmm,~oddball|manipulation+reverse,var='saliva_cortisol_z'))
##--> higher saliva cortisol is associated with higher oddball response before manipulation
##--> higher saliva cortisol is associated with lower oddball response after manipulation

require(wesanderson) #custom color palettes
custom_condition_colors <- wes_palette('Darjeeling1',2,type='discrete') #reverse custom colors to match color coding in other figures


ggplot(df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward' &
                df_agg$saliva_cortisol_avg<10,],
       aes(x=saliva_cortisol_avg,y=rpd,color=oddball))+geom_point()+geom_smooth(method='lm')+
  scale_color_manual(values = custom_condition_colors)+
  labs(x='salivary cortisol (nmol/l)',y='pupillary repsonse (mm)')+
  theme_bw()
#--> higher saliva cortisol differentiates pupillary responses

# ADD demographic data (database export) ####

#data is generated with: data_preprocessing_database.Rmd
df_dem<-readRDS("C:/Users/nico/PowerFolders/project_sega/data/demographics_25092024.rds")
#df_dem<-readRDS("C:/Users/nico/PowerFolders/project_sega/data/demographics_28052024.rds")

df_agg<-merge(df_agg,df_dem,by='id',all.x=T) # only data with saliva and hair data
df_trial<-merge(df_trial,df_dem,by='id',all.x=T) #all data 

## - association of cortisol with BMI ####
df_agg$BMI<-with(df_agg,Gewicht/((Groeße/100)*(Groeße/100)))
hist(df_agg$BMI)


lm<-lm(scale(saliva_cortisol_avg)~scale(BMI),
       df_agg[df_agg$saliva_when=='before' & df_agg$manipulation=='before' & df_agg$oddball=='oddball' & df_agg$reverse=='forward',])
summary(lm)
###--> higher BMI is associated with higher cortisol before assessments

lm<-lm(scale(saliva_cortisol_change)~scale(BMI),
       df_agg[!duplicated(df_agg$id),])
summary(lm)
###--> higher BMI also associated with stronger cortisol dropoff

lm<-lm(scale(hair_cortisol)~scale(BMI),
       df_agg[!duplicated(df_agg$id),])
anova(lm)
###--> but not associated with hair cortisol

lm<-lm(scale(rpd_low)~scale(BMI),
       df_agg[!duplicated(df_agg$id),])
anova(lm)
###--> but not associated with baseline pupil size


## - analysis Internalizing and Externalizing with core measures ####

### explore in df_trial
lmm<-lmer(scale(rpd)~trial*(SDQ_int_s+SDQ_ext_s)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$block=='forward',])
summary(lmm)
anova(lmm)
r2_nakagawa(lmm)

lmm<-lmer(scale(rpd_low)~trial*(SDQ_int_s+SDQ_ext_s)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$block=='forward',])
anova(lmm)
###--> no association

lmm<-lmer(scale(rpd)~trial*(CBCL_T_INT+CBCL_T_EXT)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$block=='forward',])
anova(lmm)
###--> no association

lmm<-lmer(scale(rpd_low)~trial*(CBCL_T_INT+CBCL_T_EXT)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$block=='forward',])
anova(lmm)
###--> no association

df_trial$YSR_T_INT_z<-scale(df_trial$YSR_T_INT)
df_trial$YSR_T_EXT_z<-scale(df_trial$YSR_T_EXT)
df_trial$YSR_T_GES_z<-scale(df_trial$YSR_T_GES)
lmm<-lmer(scale(rpd)~trial*(YSR_T_INT_z+YSR_T_EXT_z+YSR_T_GES_z)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$block=='forward',])
summary(lmm)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm,~trial,var="YSR_T_INT_z")
## higher self-rated Internalizing is associated with lower pupillary response to oddballs

###--> explore for an per-participant aggregated measure
df_agg_firstblock<-aggregate(df_trial[df_trial$manipulation=='before' & df_trial$block=='forward',
                                      names(df_trial) %in% c('rpd_auc','rpd_high','rpd_low','rpd')],
                             by=with(df_trial[df_trial$manipulation=='before' & df_trial$block=='forward',],
                                     list(id,trial)),
                             mean,na.rm=T)
names(df_agg_firstblock)<-c('id','oddball','rpd_auc','rpd_high','rpd_low','rpd')


df_agg_firstblock<-merge(df_agg_firstblock,df_dem,by='id')

df_agg_firstblock$YSR_T_INT_z<-scale(df_agg_firstblock$YSR_T_INT)
df_agg_firstblock$YSR_T_EXT_z<-scale(df_agg_firstblock$YSR_T_EXT)
df_agg_firstblock$YSR_T_GES_z<-scale(df_agg_firstblock$YSR_T_GES)
df_agg_firstblock$YSR_T_split<-ifelse(df_agg_firstblock$YSR_T_GES>69,'above clinical cutoff','below clinical cutoff')
df_agg_firstblock<-df_agg_firstblock[!is.na(df_agg_firstblock$YSR_T_INT),]

lmm<-lmer(scale(rpd)~oddball*(YSR_T_INT_z+YSR_T_EXT_z+YSR_T_GES_z)+(1|id),
          df_agg_firstblock)
summary(lmm)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm,~oddball,var='YSR_T_INT_z')    
emtrends(lmm,~oddball,var='YSR_T_GES_z')    

lmm<-lmer(rpd~oddball*(YSR_T_INT+YSR_T_EXT+YSR_T_GES)+(1|id),
          df_agg_firstblock)
df_agg_firstblock$predicted_rpd<-predict(lmm)
ggplot(df_agg_firstblock,aes(YSR_T_GES,predicted_rpd,color=oddball))+
  geom_point()+
  labs(x='p-factor (YSR total t-score)',y='pupillary response (mm)')+
  geom_smooth(method='lm')

ggplot(df_agg_firstblock,aes(YSR_T_GES,predicted_rpd,group=interaction(YSR_T_split,oddball),color=oddball))+
  geom_point()+
  geom_vline(xintercept=65,linetype=2)+
  labs(x='p-factor (YSR total t-score)',y='pupillary response (mm)')+
  geom_smooth(method='lm')+theme_bw()


require(wesanderson) #custom color palettes
df_agg_firstblock$custom_group<-with(df_agg_firstblock,
                                     ifelse(YSR_T_GES<65 & YSR_T_INT<65,'non-clinical',
                                            ifelse(YSR_T_GES>YSR_T_INT,'high p-factor',
                                                   ifelse(YSR_T_GES<=YSR_T_INT,'high int','else'))))


custom_condition_colors <- wes_palette('Cavalcanti1',3,type='discrete') #reverse custom colors to match color coding in other figures

ggplot(df_agg_firstblock[df_agg_firstblock$oddball=='oddball',],
       aes(YSR_T_INT,YSR_T_GES,color=custom_group,
           size=predicted_rpd,alpha=predicted_rpd))+
  geom_point()+
  scale_color_manual(values = custom_condition_colors)+
  labs(x='internalizing (t-score)',y='p-factor (t-score)',
       color = "cluster group",size= 'pupillary resp. (mm)',alpha= 'pupillary resp. (mm)')+ theme_bw()




ggplot(df_agg_firstblock,aes(YSR_T_INT_z,predicted_rpd,color=oddball))+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x='internalizing psychopathology (YSR, t-score)',y='pupillary response (mm)')+ 
  coord_cartesian(ylim=c(0,0.08))+
  theme_bw()





lmm<-lmer(scale(rpd_low)~trial*(YSR_T_INT+YSR_T_EXT+YSR_T_GES)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$block=='forward',])
summary(lmm)
anova(lmm)


lm<-lm(scale(rpd_low)~oddball*(YSR_T_INT+YSR_T_EXT+YSR_T_GES),
       df_agg_firstblock)
summary(lm)
anova(lm)


## - cortisol and psychopathology ####

lm<-lm(scale(hair_cortisol)~scale(YSR_T_INT)+scale(YSR_T_EXT),
       df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward' & df_agg$oddball=='oddball',])

summary(lm)

lm<-lm(hair_cortisol~SDQ_int_s+SDQ_ext_s,
       df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward' & df_agg$oddball=='oddball',])

summary(lm)


lm<-lm(saliva_cortisol_change~YSR_T_INT+YSR_T_EXT,
       df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward'  & df_agg$oddball=='oddball',])

summary(lm)


lm<-lm(saliva_cortisol_change~SDQ_int_s+SDQ_ext_s,
       df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward'  & df_agg$oddball=='oddball',])

summary(lm)


lm<-lm(saliva_cortisol_z~YSR_T_INT+YSR_T_EXT,
       df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward'  & df_agg$oddball=='oddball',])

anova(lm)


lm<-lm(saliva_cortisol_avg~YSR_T_INT+YSR_T_EXT,
       df_agg[df_agg$manipulation=='before' & df_agg$reverse=='reverse'  & df_agg$saliva_when=='before',])

summary(lm)




