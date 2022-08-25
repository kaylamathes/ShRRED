########################
###Nitrogen Data
########################

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

nitrogen$Treatment <- as.factor(nitrogen$Treatment)
nitrogen$Date_Collected <- as.factor(nitrogen$Date_Collected)


#################Summarize Data##########

nitrogen_summary <- nitrogen%>%
  group_by(Treatment, Date_Collected, Plot)%>%
  summarize(ave_ammonium = mean(ammonium_calculated), ave_nitrate = mean(nitrate_calculated))


################figures###################

ggplot(nitrogen_summary, aes(x = Date_Collected, y = ave_ammonium, fill = Treatment)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B","#4967D7")) +
  xlab("Date") +ylab("Soil ammonium (mg ammonium/kg soil")

ggplot(nitrogen_summary, aes(x = Date_Collected, y = ave_nitrate, fill = Treatment)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B","#4967D7")) +
  xlab("Date") +ylab("Soil nitrate (mg nitrate/kg soil")





#################Stats#####################
#Equality of variance test for Type
leveneTest(ammonium_calculated_log ~ Treatment, data = nitrogen)

##Normality (Data are normal)
# Build the linear model
normality_test  <- lm(ammonium_calculated_log ~ Treatment,
                      data = nitrogen)

# Create a QQ plot of residuals
ggqqplot(residuals(normality_test))
# Shapiro test of normality 
shapiro_test(residuals(normality_test))


nitrogen <- nitrogen%>%
  mutate(ammonium_calculated_log = log(ammonium_calculated))

nitrogen_anova <- aov(ammonium_calculated_log ~ Treatment, data = nitrogen)
summary(nitrogen_anova)


##########Leaf C:N#################

## Import downloaded date from new data directory "googledrive_data"
leaf_C_N <- read.csv("googledrive_data/Leaf_N_and_C.csv", na.strings = c("NA", "na"))


#############Figures##############

ggplot(leaf_C_N, aes(x = Type_Location, y = N_percent, fill = Type_Location)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B", "#ADD8E6", "#4967D7" )) +
  xlab("Location") +ylab("% Leaf Nitrogen")

