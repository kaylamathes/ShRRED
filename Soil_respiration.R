##############################
####ShRRED Soil Respiration 
###Kayla C. Mathes
##############################

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
  filter(!is.na(Efflux))

Rs$Date <- as.Date(Rs$Date)

Rs_summary <- Rs%>%
  group_by(Type, Date)%>%
  summarize(ave_Efflux = mean(Efflux), std_efflux = std.error(Efflux),
            ave_soil_T = mean(Soil_T), std_soil_T = std.error(Soil_T),
            ave_VWC = mean(Soil_VWC), std_VWC = std.error(Soil_VWC))


Rs_stats <- Rs%>%
  group_by(Type,Date,Plot_ID,)%>%
  summarize(ave_Efflux = mean(Efflux), std_efflux = std.error(Efflux),
            ave_soil_T = mean(Soil_T), std_soil_T = std.error(Soil_T),
            ave_VWC = mean(Soil_VWC), std_VWC = std.error(Soil_VWC))


###Graphing 

Rs_summary$Date <- as.Date(Rs_summary$Date)
Rs_stats$Date <- as.factor(Rs_stats$Date)
Rs_stats$Type <- as.factor(Rs_stats$Type)


##Make a plot of Rs, temperature and moisture by type 
Rs <- ggplot(Rs_stats, aes(x = Date, y = ave_Efflux, fill = Type))+
  geom_boxplot() +
  scale_fill_manual(values = c("#8F6B2B","#4967D7"))+
  theme_classic() +
  theme(legend.position = c(0.2,.9))+
  scale_y_continuous(position = "left",breaks = seq(from = 1, to = 10, by = 1), sec.axis = dup_axis(name = NULL, labels = NULL))+
  labs(x = "Date", y=expression(paste(" ",R[s]," (",mu*molCO[2]," ",m^-2," ",sec^-1,")"))) 

Temp <- ggplot(Rs_stats, aes(x = Date, y = ave_soil_T, fill = Type))+
  geom_boxplot() +
  scale_fill_manual(values = c("#8F6B2B","#4967D7"))+
  theme_classic() +
  theme(legend.position = "none")+
  scale_y_continuous(position = "left", sec.axis = dup_axis(name = NULL, labels = NULL))+
  labs(x = "Date", y= expression('Temperature ('*~degree*C*')'))

VWC <- ggplot(Rs_stats, aes(x = Date, y = ave_VWC, fill = Type))+
  geom_boxplot() +
  scale_fill_manual(values = c("#8F6B2B","#4967D7"))+
  theme_classic() +
  theme(legend.position = "none")+
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
outliers <- Rs_stats %>%
  group_by(Type) %>%
  identify_outliers(ave_Efflux)


##Equality of variance test for Type: Data are equal
leveneTest(ave_Efflux ~ Type, data = Rs_stats)

##Normality 
# Build the linear model
normality_test  <- lm(ave_Efflux ~ Type,
                      data = Rs_stats)

# Create a QQ plot of residuals
ggqqplot(residuals(normality_test))
# Shapiro test of normality = DATA ARE NOT NORMAL
shapiro_test(residuals(normality_test))

Rs_stats_transformed <- Rs_stats%>%
  mutate(ave_Efflux_transformed = log(ave_Efflux))

##Normality Test with log transformed
normality_test_log  <- lm(ave_Efflux_transformed ~ Type,
                      data = Rs_stats_transformed)
##QQ plot
ggqqplot(residuals(normality_test_log))
# Shapiro test of normality = DATA ARE NORMAL
shapiro_test(residuals(normality_test_log))


##ANOVA for Rs between Control and treatment 

Rs_ANOVA <- aov(ave_Efflux_transformed ~ Type*Date,  data = Rs_stats_transformed)
summary(Rs_ANOVA)


