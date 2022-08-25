########################
###Non-Structural Carbohydrate Data
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


##Add a "Type_Location" Column that combines the Treatment_type with the Tree_type data for vizualization 
NSC_summary <- NSC%>%
  mutate(Type_Location = case_when(Treatment_type == "C" ~ "C_Focal", 
                                   Treatment_type == "T" & Tree_type == "F" ~ "T_focal", 
                                   Treatment_type == "T" & Tree_type == "N" ~ "T_adjacent"))
NSC_summary <- NSC_summary%>%
  mutate(Core_type_2 = case_when(Core_type == "B" ~ "Bole", 
                                   Core_type == "R" ~ "Root")) 
                              


##############Figures################

ggplot(NSC_summary, aes(x = Date, y = NSC, fill = Type_Location))+
  geom_boxplot()+
  facet_wrap(~ Core_type_2) +
  scale_y_continuous(position = "left", sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme_classic() +
  scale_fill_manual(values = c("#8F6B2B", "#ADD8E6", "#4967D7" )) +
  ylab("Non-structural Carbohydates (%)")
  
  

         
