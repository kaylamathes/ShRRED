#####ShRRED Data Analysis#######


##library
library(googledrive)
library(ggplot2)
library(dplyr)
library(rstatix)
library(plotrix)
library(car)
library(ggpubr)
library(gridExtra)

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
  filter(!is.na(Efflux))

Rs_stats <- Rs%>%
  group_by(Plot_ID,Type, Date)%>%
  summarize(ave_Efflux = mean(Efflux), std_efflux = std.error(Efflux),
            ave_soil_T = mean(Soil_T), std_soil_T = std.error(Soil_T),
            ave_VWC = mean(Soil_VWC), std_VWC = std.error(Soil_VWC))

Rs_summary <- Rs%>%
  group_by(Type, Date)%>%
  summarize(ave_Efflux = mean(Efflux), std_efflux = std.error(Efflux),
            ave_soil_T = mean(Soil_T), std_soil_T = std.error(Soil_T),
            ave_VWC = mean(Soil_VWC), std_VWC = std.error(Soil_VWC))

###Graphing 

Rs_summary$Date <- as.Date(Rs_summary$Date)

##Make summary dataframe for line graph of Rs, temperature and moisture 
Rs_figure <- Rs_summary%>%
  group_by(Date, Type)%>%
  summarize(mean_Efflux = mean(ave_Efflux), std_efflux = std.error(ave_Efflux),
mean_soil_T = mean(ave_soil_T), std_soil_T = std.error(ave_soil_T),
mean_VWC = mean(ave_VWC), std_VWC = std.error(ave_VWC))








##Make a plot of Rs, temperature and moisture by type 
Rs <- ggplot(Rs_summary, aes(x = Date, y = ave_Efflux, color = Type))+
  geom_path()+
  geom_point()+
  scale_color_manual(values = c("#8F6B2B","#4967D7"))+
  theme_classic() +
  geom_errorbar(mapping=aes(x=Date, ymin=ave_Efflux - std_efflux, ymax=ave_Efflux + std_efflux), width = 0.1) +
  theme(legend.position = c(0.2,.9))+
  scale_x_date( sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_y_continuous(position = "left",breaks = seq(from = 1, to = 10, by = 1), sec.axis = dup_axis(name = NULL, labels = NULL))+
  labs(x = "Date", y=expression(paste(" ",R[s]," (",mu*molCO[2]," ",m^-2," ",sec^-1,")"))) 

Temp <- ggplot(Rs_figure, aes(x = Date, y = mean_soil_T, color = Type))+
  geom_path()+
  geom_point()+
  scale_color_manual(values = c("#8F6B2B","#4967D7"))+
  theme_classic() +
  theme(legend.position = "none")+
  geom_errorbar(mapping=aes(x=Date, ymin=mean_soil_T - std_soil_T, ymax=mean_soil_T + std_soil_T), width = 0.1) +
  scale_x_date( sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_y_continuous(position = "left", sec.axis = dup_axis(name = NULL, labels = NULL))+
  labs(x = "Date", y= expression('Temperature ('*~degree*C*')'))

VWC <- ggplot(Rs_figure, aes(x = Date, y = mean_VWC, color = Type))+
  geom_path()+
  geom_point()+
  scale_color_manual(values = c("#8F6B2B","#4967D7"))+
  theme_classic() +
  theme(legend.position = "none")+
  scale_x_date( sec.axis = dup_axis(name = NULL, labels = NULL))+
  geom_errorbar(mapping=aes(x=Date, ymin=mean_VWC - std_VWC, ymax=mean_VWC + std_VWC), width = 0.1) +
  scale_y_continuous(position = "left", sec.axis = dup_axis(name = NULL, labels = NULL))+
  labs(x = "Date", y= "Soil Moisture (%)") 

##Rs and micrometerology Figure 
p1_grob <- ggplotGrob(Rs)
p2_grob <- ggplotGrob(Temp)
p3_grob <- ggplotGrob(VWC)

layout <- rbind(c(1),
                c(2), 
                c(3))
g_timeseries <- grid.arrange(p1_grob, p2_grob, p3_grob, layout_matrix=layout)


###Stats 

####Testing Assumptions 
##Test for outliers test: no outliers
outliers <- Rs %>%
  group_by(Type) %>%
  identify_outliers(Efflux)

Rs_stats_2 <- Rs%>%
  filter(Efflux < 20.7)%>%
  group_by(Plot_ID,Type, Date, Plot_.)%>%
  summarize(ave_Efflux = mean(Efflux), std_efflux = std.error(Efflux),
            ave_soil_T = mean(Soil_T), std_soil_T = std.error(Soil_T),
            ave_VWC = mean(Soil_VWC), std_VWC = std.error(Soil_VWC))%>%
ungroup()

##Equality of variance test for Type
leveneTest(ave_Efflux ~ Type, data = Rs_stats_2)

##Normality (Data are normal)
# Build the linear model
normality_test  <- lm(ave_Efflux ~ Type,
                      data = Rs_stats_2)

# Create a QQ plot of residuals
ggqqplot(residuals(normality_test))
# Shapiro test of normality 
shapiro_test(residuals(normality_test))

Rs_t_test <- Rs_stats_2%>%
  select(Type, Date, ave_Efflux, Plot_.)

##T-test for Rs between Control and treatment 
t.test(ave_Efflux ~ Type, var.equal = TRUE, data = Rs_t_test, paired=TRUE)

anova <- anova_test(data = Rs_t_test, dv = ave_Efflux, wid = Plot_., within = c(Date, Type))

get_anova_table(anova)

