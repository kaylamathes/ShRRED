###############################################
##ShRRED Fine-root Analysis 
### Kayla C. Mathes
################################################

##library
library(googledrive)
library(ggplot2)
library(dplyr)
library(rstatix)
library(plotrix)
library(car)
library(ggpubr)
library(gridExtra)


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
  select(Label, Type, Plot, Dry_root_weight_g_measured, Dry_root_adjusted_g, percent_root, SoilVol.m3., Length.cm., ProjArea.cm2., SurfArea.cm2., AvgDiam.mm., RootVolume.cm3., Tips, Forks, Crossings)%>%
  filter(!is.na(Dry_root_weight_g_measured))%>%
  cbind(Date = "2022-06-26")
 
 
 FR_7_13_clean <- FR_7_13%>%
   select(RHIZO.2017a, Dry_root_weight_g_measured, percent_root, SoilVol.m3., Length.cm., ProjArea.cm2., SurfArea.cm2., AvgDiam.mm., RootVolume.cm3., Tips, Forks, Crossings)%>%
   filter(!is.na(Dry_root_weight_g_measured))%>%
   cbind(Date = "2022-07-13")
 
 FR_7_25_clean <- FR_7_25%>%
   select(RHIZO.2017a, Dry_root_weight_g_measured, SoilVol.m3., Length.cm., ProjArea.cm2., SurfArea.cm2., AvgDiam.mm., RootVolume.cm3., Tips, Forks, Crossings)%>%
   filter(!is.na(Dry_root_weight_g_measured))%>%
   cbind(Date = "2022-07-")
 
 
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
                                 
   
   
  


