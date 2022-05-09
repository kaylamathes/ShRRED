#####ShRRED Data Analysis#######


##library
library(googledrive)
library(ggplot2)
library(dplyr)

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


