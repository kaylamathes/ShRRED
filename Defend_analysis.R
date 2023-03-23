###########Analysis for Defense #############

##library
library(googledrive)
library(ggplot2)
library(dplyr)
library(rstatix)
library(plotrix)
library(car)
library(ggpubr)
library(gridExtra)
library(AICcmodavg)
library(agricolae)
library(gridExtra)
library(ggbeeswarm)

#######Stem map for ShRRED project########## 


##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data/FoRTE_prelim_data"
as_id("https://drive.google.com/drive/folders/1i6H86C8OSXAQVhMJgeEKMtbEIMk7M20U") %>% 
  drive_ls ->
  gdfiles


# Create a new data directory for files, if necessary
data_dir <- "googledrive_data/"
if(!dir.exists(data_dir)) dir.create(data_dir)

#Download date
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}


## Import downloaded date from new data directory "googledrive_data"
stem_map <- read.csv("googledrive_data/ShRRED_Stem_map_data.csv", na.strings = c("NA", "na"))%>%
  select(1:18)%>%
  na.omit()%>%
  mutate(DBH_cm = Tree_Dia/10)%>%
  mutate(Longitude_backwards= Longitude*-1)%>%
  rename(species = Tree_Spc)

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

##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data"
as_id("https://drive.google.com/drive/folders/1DTSowyhfDoIr6MhZlUacvtg_C-YonOWF") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "googledrive_data/"
if(!dir.exists(data_dir)) dir.create(data_dir)

#Download data
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}

## Import downloaded date from new data directory "googledrive_data"
Rs <- read.csv("googledrive_data/Soil_Respiration_ShRRED.csv", na.strings = c("NA", "na"))

##Clean Data
Rs <- Rs%>%
  select(!Comments)%>%
  filter(!is.na(Efflux))%>%
  rename(Plot = Plot_.)


Rs_stats <- Rs%>%
  group_by(Type,Date,Plot_ID, Plot)%>%
  summarize(ave_Efflux = mean(Efflux), std_efflux = std.error(Efflux),
            ave_soil_T = mean(Soil_T), std_soil_T = std.error(Soil_T),
            ave_VWC = mean(Soil_VWC), std_VWC = std.error(Soil_VWC)) %>%
  ungroup()

Rs_stats$Type <- as.factor(Rs_stats$Type)
Rs_stats$Date <- as.factor(Rs_stats$Date)

Rs_stats <- Rs_stats%>%
  mutate(Days_after_girdle = case_when(Date == "2022-05-26" ~ "-4", 
                                       Date == "2022-06-22" ~ "22", 
                                       Date == "2022-06-28" ~ "28", 
                                       Date == "2022-07-18" ~ "48", 
                                       Date == "2022-07-22" ~ "52", 
                                       Date == "2022-07-29" ~ "59", 
                                       Date == "2022-09-09" ~ "99"))



#################################### Rs statistics###########################


##Convert to factors
Rs_stats$Plot_ID <- as.factor(Rs_stats$Plot_ID)
Rs_stats$Type <- as.factor(Rs_stats$Type)
Rs_stats$Date <- as.factor(Rs_stats$Date)

####Testing Assumptions 
##Test for outliers test: no extreme outliers
outliers <- Rs_stats %>%
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

##turn the models into lm 

Rs_lm_VWC_Temp <- lm(ave_Efflux_transformed ~ Type * Date + ave_VWC + ave_soil_T, data = Rs_stats)
Rs_lm_Temp <- lm(ave_Efflux_transformed ~ Type * Date +  ave_soil_T, data = Rs_stats)
Rs_lm_VWC <- lm(ave_Efflux_transformed ~ Type * Date + ave_VWC, data = Rs_stats)
Rs_lm <- lm(ave_Efflux_transformed ~ Type * Date, data = Rs_stats)

#define list of models 
rs_models <- (list(Rs_lm,Rs_lm_VWC,Rs_lm_Temp,Rs_lm_VWC_Temp))

#specify model names
rs_mod.names <- c("Rs_lm","Rs_lm_VWC","Rs_lm_Temp","Rs_lm_VWC_Temp")

#calculate AIC of each model  
aictab(cand.set = rs_models, modnames = rs_mod.names)



###Create a pre-disturbance soil met dataframe 

Soil_met <- Rs_stats%>%
  filter(Date == "2022-05-26")%>%
  select(Type, Plot_ID, Plot, ave_soil_T, ave_VWC)


###Soil respiration Data frame for multivariates 

Soil_summary <- Rs_stats%>%
  filter(Date == "2022-07-18" | Date == "2022-07-22" | Date == "2022-07-29")%>%
  mutate(week = case_when(Date == "2022-07-18" ~ "1", 
                           Date == "2022-07-22" ~ "2",
                           Date == "2022-07-29" ~ "3"))%>%
  select(ave_Efflux, Plot_ID, Type, week)

##################################NSC###############################################



##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data/"
as_id("https://drive.google.com/drive/folders/1r00ZJCCPAAQIrhQH3YT0b9xNETzBDDmM") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "googledrive_data/"
if(!dir.exists(data_dir)) dir.create(data_dir)

#Download date
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}


## Import downloaded date from new data directory "googledrive_data"
NSC <- read.csv("googledrive_data/NSC_summary_mathes.csv", na.strings = c("NA", "na"))%>%
  select(Sample_Name, Date, Core_type, Plot_ID, Treatment_type, Tree_type, NSC)%>%
  filter(!is.na(NSC))


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


##### Make a girdled and control data frame only
NSC_summary_root_girdled = NSC_summary%>%
  filter(Core_type_2 == "Root")%>%
  filter(Type_Location == "C_Focal" | Type_Location == "T_adjacent")

##### Make a non-girdled and control data frame only
NSC_summary_root_non_girdled = NSC_summary%>%
  filter(Core_type_2 == "Root")%>%
  filter(Type_Location == "C_Focal" | Type_Location == "T_focal")


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


NSC_root_anova <- aov(NSC_log ~  Type*Date + Error(Plot), data = NSC_summary_root_girdled)
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


NSC_root_anova_non_girdle <- aov(NSC_log~ Type*Date + Error(Plot), data = NSC_summary_root_non_girdled)
summary(NSC_root_anova_non_girdle)


###################################### Nitrogen ####################################


##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data/FoRTE_prelim_data"
as_id("https://drive.google.com/drive/folders/1m1oaYbm4vtJtOoi9uvTl9HTnUF02GaDD") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "googledrive_data/"
if(!dir.exists(data_dir)) dir.create(data_dir)

#Download date
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}


## Import downloaded date from new data directory "googledrive_data"
nitrogen <- read.csv("googledrive_data/CalculatedN_ShRRED_REU2022.csv", na.strings = c("NA", "na"))%>%filter(!is.na(ammonium_calculated))



#################Summarize Data##########

nitrogen_summary <- nitrogen%>%
  group_by(Treatment, Date_Collected, Plot, Plot_ID)%>%
  summarize(ave_ammonium = mean(ammonium_calculated), ave_nitrate = mean(nitrate_calculated), 
            Average_Soil_Temp = mean(Average_Soil_Temp))%>%
  ungroup()

nitrogen$Treatment <- as.factor(nitrogen$Treatment)
nitrogen$Date_Collected <- as.factor(nitrogen$Date_Collected)

nitrogen_summary <- nitrogen_summary%>%
rename(Type = Treatment)%>%
  rename(Date = Date_Collected)

nitrogen_summary <- nitrogen_summary%>%
  mutate(Days_after_girdle = case_when(Date == "2022-06-28" ~ "28",
                                       Date == "2022-07-13" ~ "44", 
                                       Date == "2022-07-19" ~ "50",
                                       Date == "2022-07-25" ~ "56"))
  


nitrogen_summary_total <- nitrogen_summary%>%
  group_by(Type, Plot_ID)%>%
  summarize(ave_ammonium = mean(ave_ammonium))

###Means Test

##Convert to factors
nitrogen_summary$Plot_ID <- as.factor(nitrogen_summary$Plot_ID)
nitrogen_summary$Type <- as.factor(nitrogen_summary$Type)
nitrogen_summary$Date <- as.factor(nitrogen_summary$Date)

nitrogen_summary <- nitrogen_summary%>%
  mutate(id = 1:32)

####Testing Assumptions 
##Test for outliers test: no outliers
outliers <- nitrogen_summary %>%
  group_by(Type) %>%
  identify_outliers(ave_ammonium)


##Equality of variance test for Type: Data are equal
leveneTest(ave_ammonium ~ Type, data = nitrogen_summary)

##Normality: Not normal 
nitrogen_summary %>%
  group_by(Type)%>%
  shapiro_test(ave_ammonium)%>%
  ungroup()


##Log transform data 

nitrogen_summary <- nitrogen_summary%>%
  mutate(ave_ammonium_transformed = log(ave_ammonium))

nitrogen_summary %>%
  group_by(Type)%>%
  shapiro_test(ave_ammonium_transformed)%>%
  ungroup()



Ammonium_test <- aov(ave_ammonium_transformed ~ Type +Date +Error(Plot), data = nitrogen_summary)
summary(Ammonium_test)

nitrate_test <- aov(ave_nitrate ~ Type + Date +Error(Plot), data = nitrogen_summary)
summary(nitrate_test)


nitrogen_post_hoc <- with(nitrogen_summary, LSD.test(ave_ammonium_transformed,Type,26, 0.5312, console = TRUE, alpha = 0.1))


nitrogen_summary1 <- nitrogen_summary%>%
  filter(Type == "Treatment")

## Import downloaded date from new data directory "googledrive_data"
leaf_C_N <- read.csv("googledrive_data/Leaf_N_and_C.csv", na.strings = c("NA", "na"))






####################################################Fine-Root Structure################



##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data/FoRTE_prelim_data"
as_id("https://drive.google.com/drive/folders/1kdQ_jNnmR45p9U2Cnqo_taW5XvDPgs3o") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "googledrive_data/"
if(!dir.exists(data_dir)) dir.create(data_dir)

#Download date
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}

## Import downloaded date from new data directory "googledrive_data"
FR_6_26 <- read.csv("googledrive_data/ShRRED_root_structure_6_26_22.csv", na.strings = c("NA", "na"))
FR_7_13 <- read.csv("googledrive_data/ShRRED_root_structure_7_13_22.csv", na.strings = c("NA","na"))
FR_7_25 <- read.csv("googledrive_data/ShRRED_root_structure_7_25_22..csv", na.strings = c("NA","na"))

###Clean Data frames 

FR_6_26_clean <- FR_6_26%>%
  select(Label, Type, Plot, Dry_root_weight_g_measured, Dry_root_adjusted_g,  SoilVol.m3., Length.cm., ProjArea.cm2., SurfArea.cm2., AvgDiam.mm., RootVolume.cm3., Tips, Forks, Crossings)%>%
  filter(!is.na(Dry_root_weight_g_measured))%>%
  cbind(Date = "2022-06-26")


FR_7_13_clean <- FR_7_13%>%
  select(RHIZO.2017a, Dry_root_weight_g_measured, percent_root, SoilVol.m3., Length.cm., ProjArea.cm2., SurfArea.cm2., AvgDiam.mm., RootVolume.cm3., Tips, Forks, Crossings)%>%
  filter(!is.na(Dry_root_weight_g_measured))%>%
  cbind(Date = "2022-07-13")

FR_7_25_clean <- FR_7_25%>%
  select(RHIZO.2017a, Dry_root_weight_g_measured, SoilVol.m3., Length.cm., ProjArea.cm2., SurfArea.cm2., AvgDiam.mm., RootVolume.cm3., Tips, Forks, Crossings)%>%
  filter(!is.na(Dry_root_weight_g_measured))%>%
  cbind(Date = "2022-07-25")


###Add a type and plot column for 7-13 and 7-25 data frames 

FR_7_25_clean <- FR_7_25_clean%>%
  mutate(Type = case_when(RHIZO.2017a == "C1_1" | RHIZO.2017a == "C1_2" |RHIZO.2017a == "C1_3" |RHIZO.2017a == "C1_4" |RHIZO.2017a == "C2_1" |RHIZO.2017a == "C2_2"
                          |RHIZO.2017a == "C2_3" |RHIZO.2017a == "C2_4" |RHIZO.2017a == "C3_1" |RHIZO.2017a == "C3_2" |RHIZO.2017a == "C3_3" | RHIZO.2017a == "C3_4"
                          |RHIZO.2017a == "C4_1"
                          |RHIZO.2017a == "C4_2" |RHIZO.2017a == "C4_3" |RHIZO.2017a == "C4_4" ~ "Control",
                          RHIZO.2017a == "T1_1" |  RHIZO.2017a == "T1_2" |  RHIZO.2017a == "T1_3" |  RHIZO.2017a == "T1_4" |  RHIZO.2017a == "T2_1" | 
                            RHIZO.2017a == "T2_2" |  RHIZO.2017a == "T2_3" |  RHIZO.2017a == "T2_4" |  RHIZO.2017a == "T3_1" |  RHIZO.2017a == "T3_2" |  
                            RHIZO.2017a == "T3_3" |  RHIZO.2017a == "T3_4" |  RHIZO.2017a == "T3_4" |  RHIZO.2017a == "T4_1" |  RHIZO.2017a == "T4_2" | 
                            RHIZO.2017a == "T4_3" |  RHIZO.2017a == "T4_4" ~ "Treatment"))%>%
  mutate(Plot = case_when(RHIZO.2017a == "C1_1" | RHIZO.2017a == "C1_2" |RHIZO.2017a == "C1_3" |RHIZO.2017a == "C1_4" | RHIZO.2017a == "T1_1" |  
                            RHIZO.2017a ==  "T1_2" |  RHIZO.2017a == "T1_3" |  RHIZO.2017a == "T1_4" ~ "1",
                          RHIZO.2017a == "C2_1" |RHIZO.2017a == "C2_2"
                          |RHIZO.2017a == "C2_3" |RHIZO.2017a == "C2_4"|RHIZO.2017a == "T2_1" | 
                            RHIZO.2017a == "T2_2" |  RHIZO.2017a == "T2_3" |  RHIZO.2017a == "T2_4" ~ "2", 
                          RHIZO.2017a == "C3_1" |RHIZO.2017a == "C3_2" |RHIZO.2017a == "C3_3" | RHIZO.2017a == "C3_4" |RHIZO.2017a == "T3_1" | 
                            RHIZO.2017a == "T3_2" |  RHIZO.2017a == "T3_3" |  RHIZO.2017a == "T3_4" ~ "3", 
                          RHIZO.2017a == "C4_1"
                          |RHIZO.2017a == "C4_2" |RHIZO.2017a == "C4_3" |RHIZO.2017a == "C4_4"|  RHIZO.2017a == "T4_1" |  RHIZO.2017a == "T4_2" | 
                            RHIZO.2017a == "T4_3" |  RHIZO.2017a == "T4_4" ~ "4" ))


FR_7_13_clean <- FR_7_13_clean%>%
  mutate(Type = case_when(RHIZO.2017a == "C1_1" | RHIZO.2017a == "C1_2" |RHIZO.2017a == "C1_3" |RHIZO.2017a == "C1_4" |RHIZO.2017a == "C2_1" |RHIZO.2017a == "C2_2"
                          |RHIZO.2017a == "C2_3" |RHIZO.2017a == "C2_4" |RHIZO.2017a == "C3_1" |RHIZO.2017a == "C3_2" |RHIZO.2017a == "C3_3" | RHIZO.2017a == "C3_4"
                          |RHIZO.2017a == "C4_1"
                          |RHIZO.2017a == "C4_2" |RHIZO.2017a == "C4_3" |RHIZO.2017a == "C4_4" ~ "Control",
                          RHIZO.2017a == "T1_1" |  RHIZO.2017a == "T1_2" |  RHIZO.2017a == "T1_3" |  RHIZO.2017a == "T1_4" |  RHIZO.2017a == "T2_1" | 
                            RHIZO.2017a == "T2_2" |  RHIZO.2017a == "T2_3" |  RHIZO.2017a == "T2_4" |  RHIZO.2017a == "T3_1" |  RHIZO.2017a == "T3_2" |  
                            RHIZO.2017a == "T3_3" |  RHIZO.2017a == "T3_4" |  RHIZO.2017a == "T3_4" |  RHIZO.2017a == "T4_1" |  RHIZO.2017a == "T4_2" | 
                            RHIZO.2017a == "T4_3" |  RHIZO.2017a == "T4_4" ~ "Treatment"))%>%
  mutate(Plot = case_when(RHIZO.2017a == "C1_1" | RHIZO.2017a == "C1_2" |RHIZO.2017a == "C1_3" |RHIZO.2017a == "C1_4" | RHIZO.2017a == "T1_1" |  
                            RHIZO.2017a ==  "T1_2" |  RHIZO.2017a == "T1_3" |  RHIZO.2017a == "T1_4" ~ "1",
                          RHIZO.2017a == "C2_1" |RHIZO.2017a == "C2_2"
                          |RHIZO.2017a == "C2_3" |RHIZO.2017a == "C2_4"|RHIZO.2017a == "T2_1" | 
                            RHIZO.2017a == "T2_2" |  RHIZO.2017a == "T2_3" |  RHIZO.2017a == "T2_4" ~ "2", 
                          RHIZO.2017a == "C3_1" |RHIZO.2017a == "C3_2" |RHIZO.2017a == "C3_3" | RHIZO.2017a == "C3_4" |RHIZO.2017a == "T3_1" | 
                            RHIZO.2017a == "T3_2" |  RHIZO.2017a == "T3_3" |  RHIZO.2017a == "T3_4" ~ "3", 
                          RHIZO.2017a == "C4_1"
                          |RHIZO.2017a == "C4_2" |RHIZO.2017a == "C4_3" |RHIZO.2017a == "C4_4"|  RHIZO.2017a == "T4_1" |  RHIZO.2017a == "T4_2" | 
                            RHIZO.2017a == "T4_3" |  RHIZO.2017a == "T4_4" ~ "4" ))


###Add an adjusted root mass (per ash free weight) But 7-13 and 7-25 will use 95.4% root (from the 7-13 procedure)

FR_7_13_clean <-FR_7_13_clean%>%
  mutate(Dry_root_adjusted_g = Dry_root_weight_g_measured*0.954)%>%
  select(!percent_root)


FR_7_25_clean <-FR_7_25_clean%>%
  mutate(Dry_root_adjusted_g = Dry_root_weight_g_measured*0.954)

###Change the label name in 7-13 and 7-25 files 
FR_7_13_clean <- rename(FR_7_13_clean, Label = RHIZO.2017a)
FR_7_25_clean <- rename(FR_7_25_clean, Label = RHIZO.2017a)

##Reordering columns in 7-13 and 7-25 to match 
FR_7_13_clean <- FR_7_13_clean[, c(1,13,14,2,15,3,4,5,6,7,8,9,10,11,12)]

FR_7_25_clean <- FR_7_25_clean[, c(1,13,14,2,15,3,4,5,6,7,8,9,10,11,12)]


##Combine all three data frames together 
FR_total <- rbind(FR_6_26_clean,FR_7_13_clean,FR_7_25_clean)

#######Calculate Specific Root Length (m/g) as the ratio between Fine-root length and dry root mass
FR_total$Length.cm. <- as.numeric(FR_total$Length.cm.)
FR_total$SurfArea.cm2. <- as.numeric(FR_total$SurfArea.cm2.)
FR_total$AvgDiam.mm. <- as.numeric(FR_total$AvgDiam.mm.)

FR_total <- FR_total%>%
  mutate(SRL_mg = (Length.cm./Dry_root_adjusted_g)/100)

###Calculate Fine-root tissue density (g/cm^3): Divide fine-root dry wt by fine-root volume

FR_total <- FR_total%>%
  mutate(Tissue_density_gcm3 = (Dry_root_adjusted_g/RootVolume.cm3.))

###Calculate Fine-root Biomass as kg_m3 from g/cm3

FR_total <- FR_total%>%
  mutate(biomass_kg_m3 = ((Dry_root_adjusted_g/314.16)*100000)/1000)

####Calculate surface area to volume ratio 

FR_total <- FR_total%>%
  mutate(SA_V = SurfArea.cm2./RootVolume.cm3.)

##Calculate

###Summarize Root Traits 

FR_total_summary <- FR_total%>%
  group_by(Type, Date, Plot)%>%
  summarize(ave_SA_V = mean(SA_V), ave_SRL = mean(SRL_mg), ave_root_mass_g = mean(Dry_root_adjusted_g), ave_root_length = mean(Length.cm.), ave_root_volume = mean(RootVolume.cm3.), ave_root_surf_area = mean(SurfArea.cm2.), ave_root_diam = mean(AvgDiam.mm.), ave_tissue_density = mean(Tissue_density_gcm3), ave_biomass = mean(biomass_kg_m3))%>%
  mutate(Plot_ID = case_when(Type == "Treatment" & Plot == "1" ~ "T1", 
                             Type == "Treatment" & Plot == "2" ~ "T2",
                             Type == "Treatment" & Plot == "3" ~ "T3",
                             Type == "Treatment" & Plot == "4" ~ "T4",
                             Type == "Control" & Plot == "1" ~ "C1",
                             Type == "Control" & Plot == "2" ~ "C2",
                             Type == "Control" & Plot == "3" ~ "C3",
                             Type == "Control" & Plot == "4" ~ "C4"))%>%
  mutate(week = case_when(Date == "2022-06-26" ~ "1", 
                          Date == "2022-07-13" ~ "2",
                          Date == "2022-07-25" ~ "3"))%>%
  ungroup()


FR_srl_summary_total <- FR_total_summary%>%
  mutate(Plot_ID = case_when(Type == "Treatment" & Plot == "1" ~ "T1", 
                             Type == "Treatment" & Plot == "2" ~ "T2",
                             Type == "Treatment" & Plot == "3" ~ "T3",
                             Type == "Treatment" & Plot == "4" ~ "T4",
                             Type == "Control" & Plot == "1" ~ "C1",
                             Type == "Control" & Plot == "2" ~ "C2",
                             Type == "Control" & Plot == "3" ~ "C3",
                             Type == "Control" & Plot == "4" ~ "C4"))%>%
  select(Plot_ID, Type, ave_SA_V, ave_root_diam, ave_biomass)%>%
  group_by(Plot_ID, Type)%>%
  summarize(ave_SA_V = mean(ave_SA_V), ave_root_diam = mean(ave_root_diam),ave_biomass = mean(ave_biomass))%>%
  ungroup()

##SAV ANOVA
##Equality of variance test for Type:
leveneTest(ave_SA_V ~ Type, data = FR_total_summary)


##Normality
FR_total_summary %>%
  group_by(Type)%>%
  shapiro_test(ave_SA_V)%>%
  ungroup()

FR_total_summary <- FR_total_summary%>%
  mutate(log_SAV = log(ave_SA_V))

SAV_model <- aov(log_SAV ~ Type + Date +Error(Plot), data = FR_total_summary)
summary(SAV_model)


##root diameter ANOVA
##Equality of variance test for Type:
leveneTest(ave_root_diam ~ Type, data = FR_total_summary)


##Normality
FR_total_summary %>%
  group_by(Type)%>%
  shapiro_test(ave_root_diam)%>%
  ungroup()

FR_total_summary <- FR_total_summary%>%
  mutate(log_root_diam = log(ave_root_diam))

root_diam_model <- aov(log_root_diam ~ Type + Date +Error(Plot), data = FR_total_summary)
summary(root_diam_model)





###############CV Test############

library(cvequality)

FR_total_summary_CV <- FR_total_summary%>%
  group_by(Type, week)%>%
  summarize(mean_SAV = mean(ave_SA_V), sd_SAV = sd(ave_SA_V), mean_root_diam = mean(ave_root_diam), sd_root_diam = sd(ave_root_diam))%>%
mutate(cv_SAV = sd_SAV / mean_SAV * 100)%>%
  mutate(cv_root_diam = sd_root_diam / mean_root_diam * 100)


SAV_CV_test <- 
  with(FR_total_summary, 
       asymptotic_test(ave_SA_V, Type))

root_diamter_CV_test <- 
  with(FR_total_summary, 
       asymptotic_test(ave_root_diam, Type))


 ##########################################Root Exudation#############################


##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data/FoRTE_prelim_data"
as_id("https://drive.google.com/drive/folders/1uSUJSXYMNQP3nkiFGQ4BgwwjS9_01Bxr") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "googledrive_data/"
if(!dir.exists(data_dir)) dir.create(data_dir)

#Download date
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}

## Import downloaded date from new data directory "googledrive_data"
Exudate_structure <- read.csv("googledrive_data/Root_exudate_structure.csv", na.strings = c("NA", "na"))
Exudates <- read.csv("googledrive_data/Root_exudates.csv", na.strings = c("NA","na"))


##Cleaning Root exudate Structure dataframe 
##Get rid of first four rows 
Exudate_structure <- Exudate_structure[-c(1, 2, 3, 4),]

##Select just the Plot ID columns and adjusted root mass column
Exudate_structure <- Exudate_structure%>%
  select(Plot, Type, Location, Dry_root_adjusted_g)

##Get an average adjusted root mass by Plot 
Exudate_Structure_mass_summary <- Exudate_structure%>%
  group_by(Plot)%>%
  summarize(ave_dry_root_mass_g = mean(Dry_root_adjusted_g))


##Clean Exudate Dataframe 

Exudates <- Exudates%>%
  select(Date, Plot, X.NPOC..mg.L)
colnames(Exudates)[3] = "NPOC_mgL"

##Adjusted NPOC for blank value 
Exudates <- Exudates%>%
  mutate(NPOC_mgL_adusted = NPOC_mgL-2.223)%>%
  filter(Plot != 0)

##Combine Exudate structure and Exudate values dataframes 

Exudates_combined <- merge(Exudates, Exudate_Structure_mass_summary, by = "Plot")

##Convert from mg/L to mg C/g root day

Exudates_combined <- Exudates_combined%>%
  mutate(NPOC_mg_g = (NPOC_mgL_adusted/1000)/ave_dry_root_mass_g)

##Add Treatment or Control Column 

Exudates_combined <- Exudates_combined%>%
  mutate(Type = case_when(Plot == "C1" |Plot == "C2" |Plot == "C3" |Plot == "C4" ~ "Control", 
                          Plot == "T1" |Plot == "T2"|Plot == "T3"|Plot == "T4" ~ "Treatment"))%>%
  rename(Plot_ID = Plot)

ggplot(Exudates_combined, aes(x = Type, y = NPOC_mg_g, fill = Type)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#8F6B2B","#4967D7"))+
  theme_classic() +
  xlab("Type") +ylab("Roote Exudates (NPOC mg*g-1*day-1)")+
  geom_quasirandom()

##CV Test

library(cvequality)

Exudate_CV_test <- 
  with(Exudates_combined, 
       asymptotic_test(NPOC_mg_g, Type))


Exudates_summary <- Exudates_combined%>%
  select(Plot_ID, Type, NPOC_mg_g)

Exudates1 <- Exudates_summary%>%
  filter(Type == "Treatment")


##root diameter ANOVA
##Equality of variance test for Type:
leveneTest(ave_root_diam ~ Type, data = FR_total_summary)


##Normality
FR_total_summary %>%
  group_by(Type)%>%
  shapiro_test(ave_root_diam)%>%
  ungroup()

FR_total_summary <- FR_total_summary%>%
  mutate(log_root_diam = log(ave_root_diam))

root_diam_model <- aov(log_root_diam ~ Type + Date +Error(Plot), data = FR_total_summary)
summary(root_diam_model)


###################################### Big Data Frame combined#########################

###Merge all factors into a dataframe 
NSC_summary_total_girdle <- NSC_summary%>%
  filter(Core_type_2 == "Root")%>%
  filter(Type_Location == "C_Focal" | Type_Location =="T_adjacent")%>%
  group_by(Plot_ID, Type)%>%
  summarize(ave_NSC_girdled = mean(NSC))

NSC_summary_girdle <- NSC_summary%>%
  filter(Core_type_2 == "Root")%>%
  filter(Type_Location == "C_Focal" | Type_Location =="T_adjacent")%>%
  rename(NSC_girdle = NSC)%>%
  mutate(week = case_when(Date == "2022-07-13" ~ "1", 
                          Date == "2022-07-19" ~ "2",
                          Date == "2022-07-25" ~ "3"))


NSC_summary_total_nongirdle <- NSC_summary%>%
  filter(Core_type_2 == "Root")%>%
  filter(Type_Location == "C_Focal" | Type_Location =="T_focal")%>%
  group_by(Plot_ID, Type)%>%
  summarize(ave_NSC_non_girdled = mean(NSC))

NSC_summary_nongirdle <- NSC_summary%>%
  filter(Core_type_2 == "Root")%>%
  filter(Type_Location == "C_Focal" | Type_Location =="T_focal")%>%
  rename(NSC_nongirdle = NSC)%>%
  mutate(week = case_when(Date == "2022-07-13" ~ "1", 
                          Date == "2022-07-19" ~ "2",
                          Date == "2022-07-25" ~ "3"))


nitrogen_summary_total = nitrogen_summary%>%
  group_by(Plot_ID, Type)%>%
  summarize(ave_ammonium = mean(ave_ammonium),ave_nitrate = mean(ave_nitrate))

nitrogen_summary <- nitrogen_summary%>%
  mutate(week = case_when(Date == "2022-06-28" ~ "0",
                          Date == "2022-07-13" ~ "1", 
                          Date == "2022-07-19" ~ "2",
                          Date == "2022-07-25" ~ "3"))

rhizosphere_dataframe <- merge(nitrogen_summary_total,NSC_summary_total_nongirdle,by = c("Plot_ID", "Type"))
rhizosphere_dataframe <- merge(NSC_summary_total_girdle,rhizosphere_dataframe, by = c("Plot_ID", "Type"))
rhizosphere_dataframe <- merge(Exudates_summary,rhizosphere_dataframe, by = c("Plot_ID", "Type"))
rhizosphere_dataframe <- merge(FR_srl_summary_total,rhizosphere_dataframe, by = c("Plot_ID", "Type"))
rhizosphere_dataframe <- merge(stem_map_summary,rhizosphere_dataframe, by = c("Plot_ID", "Type"))
rhizosphere_dataframe <- merge(Soil_met, rhizosphere_dataframe, by = c("Plot_ID", "Type"))


###Paired dataframes (With multiple dates included)

nitrogen_NSC <- merge(nitrogen_summary, NSC_summary_nongirdle, by = c("Plot_ID", "Type", "week"))

nitrogen_NSC <- merge(nitrogen_NSC, NSC_summary_girdle, by = c("Plot_ID", "Type", "week"))%>%
select(Plot_ID, Plot.x, Type,week, ave_ammonium, ave_nitrate, NSC_girdle, NSC_nongirdle)

nitrogen_NSC_root <- merge(nitrogen_NSC, FR_total_summary, by = c("Plot_ID", "Type", "week"))

nitrogen_NSC_root_Rs <- merge(nitrogen_NSC_root, Soil_summary,by = c("Plot_ID", "Type", "week"))

####Relationship between ammonium and non girdled focal tree (Significant!!!)
library(sjPlot)
library(sjmisc)
library(emmeans)



###Make a log NSC girdled value (normality improves)
nitrogen_NSC <- nitrogen_NSC%>%
  mutate(log_NSC_nongirdle = log(NSC_nongirdle))
##histogram 
hist(nitrogen_NSC$log_NSC_nongirdle)
  
#####Regression model      
lm_nitrogen_non_girdled_2 <- lm(log_NSC_nongirdle ~ave_ammonium + ave_ammonium*Type + week, data = nitrogen_NSC )
summary(lm_nitrogen_non_girdled_2)



###Asumptions: Assumptions met 
par(mfrow = c(2, 2))
plot(lm_nitrogen_non_girdled_2)
gvlma::gvlma(lm_nitrogen_non_girdled_2)

##Interaction plot ###
#plot_model(lm_nitrogen_non_girdled_2, type = "int")



####Relationship between nitrate and non girdled focal trees (not significant) 

lm_nitrogen_non_girdled_3 <- lm(NSC_nongirdle ~ ave_nitrate*Type +week , data = nitrogen_NSC )
summary(lm_nitrogen_non_girdled_3)

###Relationship between ammoumium and girdled trees (Not significant)
hist(nitrogen_NSC$NSC_girdle)

nitrogen_NSC <- nitrogen_NSC%>%
  mutate(log_NSC_girdle = log(NSC_girdle))

lm_nitrogen_girdle <- lm(log_NSC_girdle ~ ave_ammonium*Type + week, data = nitrogen_NSC )
summary(lm_nitrogen_girdle)

########Relationship between nitrogen and root traits (not significant)

lm_nitrogen_root_1 <- lm(ave_root_diam ~ ave_ammonium*Type, data = nitrogen_NSC_root)
summary(lm_nitrogen_root_1)


ggplot(nitrogen_NSC_root)



###relationship between nitrogen and Exudation (not significant)

lm_nitrogen_exudation <- lm(NPOC_mg_g ~ ave_ammonium*Type, data = rhizosphere_dataframe)
summary(lm_nitrogen_exudation)


###Relationship between Rs and NSC

nitrogen_NSC_root_Rs <- nitrogen_NSC_root_Rs%>%
  mutate(log_NSC_nongirdle = log(NSC_nongirdle))%>%
  mutate(log_ammonium = log(ave_ammonium))

lm_NSC_Rs <- lm(log_NSC_nongirdle ~ave_Efflux + Date, data = nitrogen_NSC_root_Rs)
summary(lm_NSC_Rs)

lm_nitrogen_Rs <- lm(log_NSC_nongirdle ~ ave_ammonium*Type + week, data = nitrogen_NSC_root_Rs)
summary(lm_nitrogen_Rs)

par(mfrow = c(2, 2))
plot(lm_nitrogen_Rs)

plot_model(lm_nitrogen_Rs, type = "int")


#########Figures ###############


######## Rs Figures (Figure 1) #############
soil_figure <- Rs_stats%>%
  group_by(Type, Days_after_girdle)%>%
  summarize(ave_Rs = mean(ave_Efflux), std_err_Rs = std.error(ave_Efflux))

soil_figure$Days_after_girdle <- as.character(soil_figure$Days_after_girdle)
soil_figure$Days_after_girdle <- as.numeric(soil_figure$Days_after_girdle)


Plot1 <- ggplot(soil_figure, aes(x = Days_after_girdle, y = ave_Rs, color = Type, group = Type))+
  geom_path(size = 2, alpha = 0.9, linetype = 2) +
  geom_point(size = 10) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 35),  axis.title.y = element_text(size = 40),axis.title.x = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = 35), legend.title = element_text(size = 35), legend.position = c(0.2,0.8), plot.margin = margin(0,0,-0.2,1.1, "cm")) +
  geom_errorbar(mapping=aes(x=Days_after_girdle, ymin=ave_Rs - std_err_Rs, ymax=ave_Rs + std_err_Rs), size = 1, width = 1)+
  scale_color_manual(values = c("#8F6B2B","#4967D7"), labels = c("Control", "Disturbance"))+
  scale_y_continuous(position = "left",breaks = seq(from = 1, to = 20, by = 1), sec.axis = dup_axis(name = NULL, labels = NULL), labels = scales::number_format(accuracy = 0.1))+
  labs(x = "Days after disturbance", y=expression(paste(" ",Soil," ",Respiration," (",mu*molCO[2]," ",m^-2," ",sec^-1,")"))) +
  geom_vline(xintercept= 0, color = "red4", linetype = 2, size = 2) +
  scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
  annotate("text", label = "A", x = -5, y = 12, size = 17) +
  annotate("text", label = "n.s.", x = 100, y = 12, size = 15)

#ggsave(path = "Figures", filename = "Rs_line.png", height = 10, width = 15, units = "in")

##Percent Change figure

##Select necessary columns
Rs_stats_percent_change <- Rs_stats%>%
  select(Type, Days_after_girdle, Plot, ave_Efflux)

##Change from long to wide and calculate log ratio
Rs_stats_percent_change <- spread(Rs_stats_percent_change, Type, ave_Efflux)%>%
  mutate(log_ratio = log(Treatment/Control))

Rs_stats_log_ratio_summary <- Rs_stats_percent_change%>%
  group_by(Days_after_girdle)%>%
  mutate(ave_log_ratio = mean(log_ratio), std_error_ratio = std.error(log_ratio))

Rs_stats_log_ratio_summary$Days_after_girdle <- as.numeric(Rs_stats_log_ratio_summary$Days_after_girdle )



Plot2 <- ggplot(Rs_stats_log_ratio_summary, aes(x = Days_after_girdle, y = ave_log_ratio)) +
  geom_path(size = 2, alpha = 0.9, color = "#4967D7", linetype = 2)+
  theme_classic()+
  theme(axis.text = element_text(size = 35), axis.title = element_text(size = 40), plot.margin = margin(-0.20,0,0,0, "cm")) +
  geom_point(size = 10, color = "#4967D7")+
  geom_errorbar(mapping=aes(x=Days_after_girdle, ymin=ave_log_ratio - std_error_ratio, ymax=ave_log_ratio + std_error_ratio), color = "#4967D7", width = 3, size = 1) +
  labs(x = "Days After Disturbance", y=expression(paste( Resistance: ln,over(R[s]*Disturbance, R[s]*Control)))) +
  geom_hline(yintercept=0, size = 1 ) +
  scale_y_continuous(position = "left", sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
  geom_vline(xintercept= 0, color = "red4", linetype = 2, size = 2) +
  annotate("text", label = "B", x = -5, y = 0.4, size = 15) 

#ggsave(path = "Figures", filename = "Soil_Respiration_Resistance.png", height = 10, width = 15, units = "in")

## Make a multipaneled figure for Rs data 
p1_grob <- ggplotGrob(Plot1)
p2_grob <- ggplotGrob(Plot2)  

layout <- rbind(c(1),
                c(2))

Rs_plots <- grid.arrange(p1_grob, p2_grob, layout_matrix=layout)

ggsave(path = "Figures", filename = "Rs_combined.png", height = 20, width = 20, units = "in",Rs_plots )


########Roots NSC (girdled/control): Figure 2 #########
NSC_summary_root_girdled$Days_after_girdle <- as.character(NSC_summary_root_girdled$Days_after_girdle)

NSC_summary_root_girdled$Days_after_girdle <- as.numeric(NSC_summary_root_girdled$Days_after_girdle)

NSC_summary_root_girdled <- NSC_summary_root_girdled%>%
  mutate(Type2 = case_when(Type == "Control" ~ "Control", 
                           Type == "Treatment" ~ "Disturbance"))

ggplot(NSC_summary_root_girdled, aes(x = Type2 , y = NSC, fill = Type))+
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B", "#4967D7" )) +
  theme(axis.title.y = element_text(size = 35), axis.text = element_text(size = 30),axis.title.x = element_blank(), legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, size = 0.7)) +
  ylab("Non-structural Carbohydates (%)") +
 annotate("text", label = "n.s.", x = 2.4, y = 31, size = 14)

ggsave(path = "Figures", filename = "NSC_roots_girdled_control.png", height = 10, width = 10, units = "in" )


NSC_summary_root_girdled1 <- NSC_summary_root_girdled%>%
filter(Type == "Treatment")


#######nitrogen Figure 3 ###########


rhizosphere_dataframe_N <- rhizosphere_dataframe%>%
  filter(ave_ammonium <10000)

rhizosphere_dataframe_N  <-rhizosphere_dataframe_N %>%
  mutate(Type2 = case_when(Type == "Control" ~ "Control", 
                           Type == "Treatment" ~ "Disturbance"))

rhizosphere_dataframe  <-rhizosphere_dataframe %>%
  mutate(Type2 = case_when(Type == "Control" ~ "Control", 
                           Type == "Treatment" ~ "Disturbance"))


Plot_Nh4 <- ggplot(rhizosphere_dataframe_N , aes(x = Type2, y = ave_ammonium, fill = Type)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B","#4967D7")) +
  theme(axis.title.y = element_text(size = 40), axis.text.y = element_text(size = 35), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.9), legend.position = "none") +
  labs(y=expression(paste(" ",Soil," ",NH[4]^"+"," ",(mg~NH[4]^"+"/kg~soil),""))) +
  annotate("text", label = "A", x = 0.5, y = 4655, size = 15) +
  annotate("text", label = "a", x = 1, y = 3400, size = 15) +
annotate("text", label = "b", x = 2, y = 4700, size = 15) +
  scale_y_continuous(breaks = seq(from = 2500, to = 4500, by = 500))
  

Plot_no3 <- ggplot(rhizosphere_dataframe , aes(x = Type2, y = ave_nitrate, fill = Type)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B","#4967D7")) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, size = 0.9), legend.position = "none",plot.margin = margin(-0.40,0.2,0,0.85, "cm")) +
  labs(x = "Plot Type", y=expression(paste(" ",Soil," ",NO[3]^"-"," ",(mg~NO[3]^"-"/kg~soil),""))) +
  annotate("text", label = "B", x = 0.5, y = 200, size = 15) +
  annotate("text", label = "n.s.", x = 2.5, y = 200, size = 15) 

## Make a multipaneled figure for Rs data 
p3_grob <- ggplotGrob(Plot_Nh4)
p4_grob <- ggplotGrob(Plot_no3)  

layout <- rbind(c(1),
                c(2))

Nitrogen_plots <- grid.arrange(p3_grob, p4_grob, layout_matrix=layout)

ggsave(path = "Figures", filename = "Nitrogen_combined.png", height = 20, width = 15, units = "in",Nitrogen_plots )



####### Ammonium and time 
ggplot(nitrogen_summary, aes(x = Days_after_girdle, y = ave_ammonium, fill = Type)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B","#4967D7")) +
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 25), legend.title = element_text(size = 30), legend.text = element_text(size = 25)) +
  labs(x = "Days after girdle", y=expression(paste(" ",Soil," ",NH[4]^"+"," ",(mg~NH[4]^"+"/kg~soil),"")))


ggsave(path = "Figures", filename = "ammonium.png", height = 10, width = 15, units = "in" )




##########################Nitrogen and root traits figure (Figure 4) ####################

Plot_FR <- ggplot(rhizosphere_dataframe, aes(x = Type2 , y = ave_SA_V, fill = Type2))+
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B", "#4967D7" )) +
  theme(axis.title.y = element_text(size = 30), axis.text.y = element_text(size = 30),axis.title.x = element_blank(),axis.text.x = element_blank(), legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, size = 0.7),plot.margin = margin(0,0,-.1,1.2, "cm")) +
  ylab("Fine-root \n surface area (cm2):Volume (cm3)") +
  annotate("text", label = "n.s.", x = 2.4, y = 88, size = 14) +
  annotate("text", label = "A", x = 0.5, y = 88, size = 14)

Plot_FR2 <- ggplot(rhizosphere_dataframe, aes(x = Type2 , y = ave_root_diam, fill = Type2))+
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B", "#4967D7" )) +
  theme(axis.title.y = element_text(size = 35), axis.text = element_text(size = 30),axis.title.x = element_blank(),axis.text.x = element_blank(), legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, size = 0.7),plot.margin = margin(-.1,0,-.10,0, "cm")) +
  ylab("Fine-root \n Average diameter (mm)") +
  annotate("text", label = "n.s.", x = 2.4, y =0.61 , size = 14) +
  annotate("text", label = "B", x = 0.5, y =0.61 , size = 14)

Plot_FR3 <- ggplot(rhizosphere_dataframe, aes(x = Type2 , y = NPOC_mg_g, fill = Type2))+
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B", "#4967D7" )) +
  theme(axis.title.y = element_text(size = 35), axis.text = element_text(size = 30),axis.title.x = element_blank(), legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, size = 0.7),plot.margin = margin(-.10,0,0,0, "cm")) +
  ylab("Fine-root Exudation \n (mg NPOC / g soil)") +
  annotate("text", label = "n.s.", x = 2.4, y =0.07 , size = 14) +
  annotate("text", label = "C", x = 0.5, y =0.07 , size = 14)

FR_total_summary_CV <- FR_total_summary_CV%>%
  mutate(Type2 = case_when(Type == "Control" ~ "Control", 
                            Type == "Treatment" ~ "Disturbance"))

FR_total_summary_CV1 <- FR_total_summary_CV%>%
  filter(Type == "Control")

SAV <- ggplot(FR_total_summary_CV, aes(x = Type2, y = cv_SAV, fill = Type2)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#8F6B2B","#4967D7"), labels = c("Control", "Non-girdled trees \n in disturbance"))+
  theme_classic() +
  xlab("Plot Type") +ylab("CV Fine-root Surface area:Volume") +
  theme(axis.title.y = element_text(size = 30), axis.text.y = element_text(size = 25), legend.title = element_text(size = 30), legend.position = "none",panel.border = element_rect(colour = "black", fill = NA, size = 0.7), axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0,0.7,-0.1,0.80, "cm")) +
  annotate("text", label = "D", x = 0.5, y = 17.5, size = 14) +
  annotate("text", label = "a", x = 1, y = 8.5, size = 10)+
  annotate("text", label = "b", x = 2, y = 17.8, size = 10) +
  scale_y_continuous(position = "right")


RD <- ggplot(FR_total_summary_CV, aes(x = Type2, y = cv_root_diam, fill = Type2)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#8F6B2B","#4967D7"), labels = c("Control", "Non-girdled trees \n in disturbance"))+
  theme_classic() +
  ylab("CV Average Fine-root Diameter") +
  theme(axis.title.y = element_text(size = 30), axis.text.y = element_text(size = 25), axis.title.x = element_blank(), axis.text.x = element_text(size = 30), legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, size = 0.7),plot.margin = margin(-0.1,0,-1,0.8, "cm")) +
  annotate("text", label = "E", x = 0.5, y = 16.8, size = 14) +
  annotate("text", label = "a", x = 1, y = 8.5, size = 10)+
  annotate("text", label = "b", x = 2, y = 15.5, size = 10) +
  scale_y_continuous(position = "right")




p5_grob <- ggplotGrob(Plot_FR)
p6_grob <- ggplotGrob(Plot_FR2)
p7_grob <- ggplotGrob(Plot_FR3) 
p8_grob <- ggplotGrob(SAV)
p9_grob <- ggplotGrob(RD)



layout <- rbind(c(1,4),
                c(2,5), 
                c(3,NA))

FR_plots <- grid.arrange(p5_grob, p6_grob, p7_grob, p8_grob, p9_grob,layout_matrix=layout)

ggsave(path = "Figures", filename = "Fr_CV_combined.png", height = 20, width = 20, units = "in",FR_plots)


################################ nitrogen and nongirdled NSC (Figure 5)  #############################


ggplot(nitrogen_NSC, aes(x = ave_ammonium, y = NSC_nongirdle, color = Type)) +
  geom_point(size = 4)+
  geom_smooth(method = "lm", se = FALSE, size = 2)+
  theme_classic() +
  scale_color_manual(values = c("#8F6B2B","#4967D7"), labels = c("Control", "Non-girdled trees \n in disturbance")) +
  labs(x=expression(paste(" ",Soil," ",NH[4]^"+"," ",(mg~NH[4]^"+"/kg~soil),"")), y = "NSC concentration (%)")+
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 25), legend.title = element_text(size = 30), legend.text = element_text(size = 25), legend.position = c(.2, .85)) +
  scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL))

ggsave(path = "Figures", filename = "ammonium_nongirdle_NSC.png", height = 10, width = 15, units = "in" )



######################## Root trait Figure (Figure 6) ###############################################

FR_total_summary <- FR_total_summary%>%
  mutate(Days_after_girdle = case_when(Date == "2022-06-26" ~ "26",
                                       Date == "2022-07-13" ~ "44", 
                                       Date == "2022-07-25" ~ "56"))

SAV <- ggplot(FR_total_summary_CV, aes(x = Type, y = cv_SAV, fill = Type)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#8F6B2B","#4967D7"), labels = c("Control", "Non-girdled trees \n in disturbance"))+
  theme_classic() +
  xlab("Plot Type") +ylab("CV of Fine-root \n Surface area:Volume") +
  theme(axis.title.y = element_text(size = 30), axis.text.y = element_text(size = 25), legend.title = element_text(size = 30), legend.position = "none",panel.border = element_rect(colour = "black", fill = NA, size = 0.7), axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0,0.4,-0.2,0, "cm")) +
  annotate("text", label = "A", x = 0.5, y = 17, size = 10) +
  annotate("text", label = "a", x = 1, y = 8.5, size = 10)+
  annotate("text", label = "b", x = 2, y = 17.5, size = 10) 


RD <- ggplot(FR_total_summary_CV, aes(x = Type, y = cv_root_diam, fill = Type)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#8F6B2B","#4967D7"), labels = c("Control", "Non-girdled trees \n in disturbance"))+
  theme_classic() +
  xlab("Plot Type") +ylab("CV of \n Average Fine-root Diameter (mm)") +
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 25), legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, size = 0.7),plot.margin = margin(-0.2,0,0,0, "cm")) +
  annotate("text", label = "B", x = 0.5, y = 17, size = 10) +
  annotate("text", label = "a", x = 1, y = 8.5, size = 10)+
  annotate("text", label = "b", x = 2, y = 15.5, size = 10) 



## Make a multipaneled figure for Rs data 
p8_grob <- ggplotGrob(SAV)
p9_grob <- ggplotGrob(RD)


layout <- rbind(c(1),
                c(2))

CV_plots <- grid.arrange(p8_grob, p9_grob, layout_matrix=layout)

ggsave(path = "Figures", filename = "CV.png", height = 20, width = 10, units = "in",CV_plots)




