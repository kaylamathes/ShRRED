############# Analysis for Short Term Rhizosphere Response to Disturbance (ShRRED) ##################
####### Script for manuscript analysis 
####### Created on July 5th, 2023 
##### Kayla Mathes 


##### Library 

library(googledrive)
library(dplyr)
library(plotrix)
library(car)
library(ggpubr)
library(gridExtra)
library(AICcmodavg)
library(agricolae)
library(gridExtra)
library(ggplot2)
library(rstatix)


##Upload All Data from googledrive 
# stem map
as_id("https://drive.google.com/drive/folders/1i6H86C8OSXAQVhMJgeEKMtbEIMk7M20U") %>% 
  drive_ls ->
  gdfiles

# soil respiration
as_id("https://drive.google.com/drive/folders/1DTSowyhfDoIr6MhZlUacvtg_C-YonOWF") %>% 
  drive_ls ->
  gdfiles

### NSCs
as_id("https://drive.google.com/drive/folders/1r00ZJCCPAAQIrhQH3YT0b9xNETzBDDmM") %>% 
  drive_ls ->
  gdfiles

# Nitrogen
as_id("https://drive.google.com/drive/folders/1m1oaYbm4vtJtOoi9uvTl9HTnUF02GaDD") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "googledrive_data/"
if(!dir.exists(data_dir)) dir.create(data_dir)

#Download data from google drive folders 
###Stem map folder 
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}

#soil resipiration folder
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}

##NSC folder 
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}

##Nitrogen 
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}

## Import downloaded date from new data directory "googledrive_data"

##Stem Map 
stem_map <- read.csv("googledrive_data/ShRRED_Stem_map_data.csv", na.strings = c("NA", "na"))%>%
  select(1:18)%>%
  na.omit()%>%
  mutate(DBH_cm = Tree_Dia/10)%>%
  mutate(Longitude_backwards= Longitude*-1)%>%
  rename(species = Tree_Spc)

##Soil Respiration
Rs <- read.csv("googledrive_data/Soil_Respiration_ShRRED.csv", na.strings = c("NA", "na"))

##NSC
NSC <- read.csv("googledrive_data/NSC_summary_mathes.csv", na.strings = c("NA", "na"))%>%
  select(Sample_Name, Date, Core_type, Plot_ID, Treatment_type, Tree_type, NSC)%>%
  filter(!is.na(NSC))

#Nitrogen
nitrogen <- read.csv("googledrive_data/CalculatedN_ShRRED_REU2022.csv", na.strings = c("NA", "na"))%>%filter(!is.na(ammonium_calculated))

#######Stem map ########## 

stem_map <- stem_map %>% 
  mutate(a = case_when(
    species == "ACRU" ~ 0.03177,
    species == "FAGR" ~ 0.1892,
    species == "PIRE" ~ 0.0526,
    species == "PIST" ~ 0.0408,
    species == "POGR" ~ 0.1387,
    species == "QURU" ~ 0.0398
  )) %>% 
  mutate(b = case_when(
    species == "ACRU" ~ 2.7780,
    species == "FAGR" ~ 2.3097,
    species == "PIRE" ~ 2.5258,
    species == "PIST" ~ 2.5735,
    species == "POGR" ~ 2.3498,
    species == "QURU" ~ 2.7734,
    
  ))%>%
  mutate(biomass_kg = a*DBH_cm^b)

stem_map_summary <- stem_map%>%
  group_by(Plot_ID)%>%
  summarize(biomass_kg_sum = sum(biomass_kg))%>%
  mutate(Type = case_when(Plot_ID == "C1" |Plot_ID == "C2" |Plot_ID == "C3" |Plot_ID == "C4" ~ "Control", 
                          Plot_ID == "T1" |Plot_ID == "T2"|Plot_ID == "T3"|Plot_ID == "T4" ~ "Treatment"))

#####################################Soil Respiration #############################

##Clean Data
Rs <- Rs%>%
  select(!Comments)%>%
  filter(!is.na(Efflux))%>%
  rename(Plot = Plot_.)

Rs_stats <- Rs%>%
  group_by(Type, Date, Plot_ID, Plot)%>%
  summarize(ave_Efflux = mean(Efflux), std_efflux = std.error(Efflux),
            ave_soil_T = mean(Soil_T), std_soil_T = std.error(Soil_T),
            ave_VWC = mean(Soil_VWC), std_VWC = std.error(Soil_VWC)) %>%
  ungroup()


Rs_stats <- Rs_stats%>%
  mutate(Days_after_girdle = case_when(Date == "2022-05-26" ~ "-4", 
                                       Date == "2022-06-22" ~ "22", 
                                       Date == "2022-06-28" ~ "28", 
                                       Date == "2022-07-18" ~ "48", 
                                       Date == "2022-07-22" ~ "52", 
                                       Date == "2022-07-29" ~ "59", 
                                       Date == "2022-09-09" ~ "99"))

########################################### Rs Statistics ######################################
########Convert to factors
Rs_stats$Plot_ID <- as.factor(Rs_stats$Plot_ID)
Rs_stats$Type <- as.factor(Rs_stats$Type)
Rs_stats$Date <- as.factor(Rs_stats$Date)

##########Testing Assumptions 
###Test for outliers test: no extreme outliers
outliers_Rs <- Rs_stats %>%
  group_by(Type, Date) %>%
  identify_outliers(ave_Efflux)


##Equality of variance test for Type: Data are equal
leveneTest(ave_Efflux ~ Type, data = Rs_stats)

##Normality: Data not normal 
Rs_stats %>%
  group_by(Type)%>%
  shapiro_test(ave_Efflux)%>%
  ungroup()

##Log transform data: normal 
##transform data
Rs_stats <- Rs_stats%>%
  mutate(ave_Efflux_transformed = log(ave_Efflux))

##Normality: Data normal 
Rs_stats%>%
  group_by(Type)%>%
  shapiro_test(ave_Efflux_transformed)%>%
  ungroup()

#####Trying different models with soil temperature and moisture as covariates
Rs_anova_VWC_Temp <- aov(ave_Efflux_transformed ~ Type*Date + ave_VWC + ave_soil_T+ Error(Plot), data = Rs_stats)
summary(Rs_anova_VWC_Temp)

Rs_anova_Temp <- aov(ave_Efflux_transformed ~ Type * Date +  ave_soil_T+ Error(Plot), data = Rs_stats)
summary(Rs_anova_Temp)

Rs_anova_VWC <- aov(ave_Efflux_transformed ~ Type * Date + ave_VWC + Error(Plot), data = Rs_stats)
summary(Rs_anova_VWC)

Rs_anova <- aov(ave_Efflux_transformed ~ Type * Date + Error(Plot), data = Rs_stats)
summary(Rs_anova)


##Turn the models into lm 
Rs_lm_VWC_Temp <- lm(ave_Efflux_transformed ~ Type * Date + ave_VWC + ave_soil_T, data = Rs_stats)
Rs_lm_Temp <- lm(ave_Efflux_transformed ~ Type * Date +  ave_soil_T, data = Rs_stats)
Rs_lm_VWC <- lm(ave_Efflux_transformed ~ Type * Date + ave_VWC, data = Rs_stats)
Rs_lm <- lm(ave_Efflux_transformed ~ Type * Date, data = Rs_stats)

#Define list of models 
rs_models <- (list(Rs_lm,Rs_lm_VWC,Rs_lm_Temp,Rs_lm_VWC_Temp))

#Specify model names
rs_mod.names <- c("Rs_lm","Rs_lm_VWC","Rs_lm_Temp","Rs_lm_VWC_Temp")

#Calculate AIC of each model: Model without temp and moisture is best fit model   
aictab(cand.set = rs_models, modnames = rs_mod.names)

##################################### NSCs ###################################

#######Clean the dataframe and make labels consistent 
NSC_summary <- NSC%>%
  mutate(Type_Location = case_when(Treatment_type == "Control" ~ "C_Focal", 
                                   Treatment_type == "Treatment" & Tree_type == "F" ~ "T_focal", 
                                   Treatment_type == "Treatment" & Tree_type == "N" ~ "T_adjacent"))%>%
  mutate(Core_type_2 = case_when(Core_type == "B" ~ "Bole", 
                                 Core_type == "R" ~ "Root")) %>%
  rename(Type = Treatment_type)%>%
  rename(Plot = Plot_ID)%>%
  mutate(Plot_ID = case_when(Type == "Treatment" & Plot == "1" ~ "T1", 
                             Type == "Treatment" & Plot == "2" ~ "T2",
                             Type == "Treatment" & Plot == "3" ~ "T3",
                             Type == "Treatment" & Plot == "4" ~ "T4",
                             Type == "Control" & Plot == "1" ~ "C1",
                             Type == "Control" & Plot == "2" ~ "C2",
                             Type == "Control" & Plot == "3" ~ "C3",
                             Type == "Control" & Plot == "4" ~ "C4"))%>%
  mutate(Days_after_girdle = case_when(Date == "2022-07-13" ~ "44", 
                                       Date == "2022-07-19" ~ "50",
                                       Date == "2022-07-25" ~ "56"))


##### Make a root core girdled and control data frame only
NSC_summary_root_girdled = NSC_summary%>%
  filter(Core_type_2 == "Root")%>%
  filter(Type_Location == "C_Focal" | Type_Location == "T_adjacent")

##### Make a root core non-girdled and control data frame only
NSC_summary_root_non_girdled = NSC_summary%>%
  filter(Core_type_2 == "Root")%>%
  filter(Type_Location == "C_Focal" | Type_Location == "T_focal")


################################################# NSC Statistics #######################################
###ANOVA Repeated Measures for girdled/control NSC roots
##Convert to factors

NSC_summary_root_girdled$Type <- as.factor(NSC_summary_root_girdled$Type)
NSC_summary_root_girdled$Date <- as.factor(NSC_summary_root_girdled$Date)


####Testing Assumptions
##Test for outliers test: no outliers
outliers_NSC_R <- NSC_summary_root_girdled %>%
  group_by(Type, Date) %>%
  identify_outliers(NSC)%>%
  ungroup()


##Equality of variance test for Type: Data are equal
leveneTest(NSC ~ Type*Date, data = NSC_summary_root_girdled)

##Normality: Not normal

NSC_summary_root_girdled %>%
  shapiro_test(NSC)

NSC_summary_root_girdled <- NSC_summary_root_girdled%>%
  mutate(NSC_log = log(NSC))

#####Data are normal with a log transformation
NSC_summary_root_girdled %>%
  shapiro_test(NSC_log)


NSC_root_anova <- aov(NSC_log ~  Type +Date + Error(Plot), data = NSC_summary_root_girdled)
summary(NSC_root_anova)



###ANOVA Repeated Measures for nongirdled/control NSC roots
##Convert to factors

NSC_summary_root_non_girdled$Type <- as.factor(NSC_summary_root_non_girdled$Type)
NSC_summary_root_non_girdled$Date <- as.factor(NSC_summary_root_non_girdled$Date)


####Testing Assumptions
##Test for outliers test: no outliers
outliers_NSC_R <- NSC_summary_root_non_girdled %>%
  group_by(Type, Date) %>%
  identify_outliers(NSC)%>%
  ungroup()


##Equality of variance test for Type: Data are equal
leveneTest(NSC ~ Type*Date, data = NSC_summary_root_non_girdled)

##Normality: Not normal

NSC_summary_root_non_girdled %>%
  shapiro_test(NSC)

NSC_summary_root_non_girdled <- NSC_summary_root_non_girdled%>%
  mutate(NSC_log = log(NSC))

#####Data are normal with a log transformation
NSC_summary_root_non_girdled %>%
  shapiro_test(NSC_log)


NSC_root_anova_non_girdle <- aov(NSC_log~ Type +Date + Error(Plot), data = NSC_summary_root_non_girdled)
summary(NSC_root_anova_non_girdle)



