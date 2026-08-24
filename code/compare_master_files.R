
root_path<-'C:/Users/nico/Nextcloud/project_sega'
require(openxlsx) #read xls

df<-read.xlsx(paste0(root_path,'/data/master_file_SEGA_data.xlsx'),na.strings = c("NA","-"),detectDates = T)
df_sep<-read.xlsx(paste0(root_path,'/data/master_file_SEGA_data_Sept2025.xlsx'),startRow = 3)
df_oct<-read.xlsx(paste0(root_path,'/data/master_file_SEGA_Oct2025.xlsx'),na.strings = c("NA","-"),detectDates = T)


#SEP: df + SDQ + handdyn + cortisol additional + exp data
names(df_sep)[!names(df_sep) %in% names(df)]

# all but SDQ and sometimes other names
names(df_oct)[!names(df_oct) %in% names(df)]
