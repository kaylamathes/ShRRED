#####ShRRED Data Analysis#######


##library
library(googledrive)
library(ggplot2)
library(dplyr)
library(rstatix)
library(plotrix)
library(car)
library(ggpubr)

#####Run prelim analysis on percent aboveground biomass by vertical stratum from FoRTE D replicate (similar ecosystem type as ShRRED)#### 

##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data/FoRTE_prelim_data"
as_id("https://drive.google.com/drive/folders/1TM992s2f8tcRM6AwaFkaNRM2kfeO4lQV") %>% 
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
canopy_2021 <- read.csv("googledrive_data/canopy_biomass_2021.csv", na.strings = c("NA", "na"))
subcanopy_2021 <- read.csv("googledrive_data/subcanopy_biomass_2021.csv", na.strings = c("NA","na"))
seedling <- read.csv("googledrive_data/three_yr_seedling_bm.csv", na.strings = c("NA","na"))

##Cleaning and combining datafiles 
seedling_2021 <- seedling%>%
  filter(year == 2021)%>%
  select(subplot_id, total_biomass, stratum)

subcanopy_2021 <- subcanopy_2021%>%
  select(subplot_id, total_biomass)
subcanopy_2021$stratum <- c("subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy","subcanopy" )

canopy_2021 <- canopy_2021%>%
  select(subplot_id, total_biomass)
canopy_2021$stratum <- c("canopy", "canopy", "canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy","canopy")

##Combine dataframe
all_strata <- rbind(seedling_2021, subcanopy_2021, canopy_2021)

##Include only D replicate of 0 severity, sum both subplots together for each stratum
all_strata_D <- all_strata%>%
  filter(subplot_id == "D01E" | subplot_id == "D01W")%>%
  group_by(stratum)%>%
  summarize(total_biomass = sum(total_biomass))

##Find the sum of biomass across all strata
sum(all_strata_D$total_biomass)

#get rid of scientific notation
options(scipen=999)

#Cacluate percent biomass per stratum
percent_biomass_strata <- all_strata_D%>%
  group_by(stratum)%>%
  summarize(percent_of_biomass = (total_biomass/36848.91)*100)


###Creating a Stem map of ShRRED Plots

##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data/FoRTE_prelim_data"
as_id("https://drive.google.com/drive/folders/1lWx-Ggzz1f49S52UkDw85WQ0gUX5m-KN") %>% 
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
  mutate(Longitude_backwards= Longitude*-1)

##Make a stem map in ggplot 

ggplot(stem_map, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(size = DBH_cm, color = Type, shape = Strata),alpha = 0.5)+
  scale_color_manual(values = c("#8F6B2B","#4967D7"))+
  theme_classic()+
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.title = element_text(size =15), legend.text = element_text(size = 10))+
  scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL))

ggsave(path = "Figures", filename = "Stem_map.png", height = 10, width = 15, units = "in")


###Soil Respiration 

##Upload Data from googledrive 
# Direct Google Drive link to "ShRRED_data"
as_id("https://drive.google.com/drive/folders/1lWx-Ggzz1f49S52UkDw85WQ0gUX5m-KN") %>% 
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
Rs <- read.csv("googledrive_data/Soil_Respiration_ShRRED.csv", na.strings = c("NA", "na"))

##Clean Data
Rs <- Rs%>%
  select(!Comments)%>%
  filter(!is.na(Efflux))

Rs_summary <- Rs%>%
  group_by(Plot_ID, Type)%>%
  summarize(ave_Efflux = mean(Efflux), std_efflux = std.error(Efflux),
            ave_soil_T = mean(Soil_T), std_soil_T = std.error(Soil_T),
            ave_VWC = mean(Soil_VWC), std_VWC = std.error(Soil_VWC))


####Testing Assumptions 
##Test for outliers test: no outliers
Rs_summary %>% 
  group_by(Type) %>%
  identify_outliers(ave_Efflux)

##Equality of variance test for Type
leveneTest(ave_Efflux ~ Type, data = Rs_summary)

##Normality (Data are normal)
# Build the linear model
normality_test  <- lm(ave_Efflux ~ Type,
                      data = Rs_summary)

# Create a QQ plot of residuals
ggqqplot(residuals(normality_test))
# Shapiro test of normality 
shapiro_test(residuals(normality_test))

##T-test for Rs between Control and treatment 

t.test(ave_Efflux ~ Type, var.equal = TRUE, data = Rs_summary)


