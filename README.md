SHort-term Rhizosphere Response to Experimental Disturbance (ShRRED) 
DOI: 10.5281/zenodo.22921444
Authors 
Kayla C. Mathes, Tatum Delaney, Lisa Nalliah and Christopher M. Gough 

Description 
All code and data files in preparation for manuscript: Dynamic rhizosphere processes stabilize soil respiration immediately following a phloem-disrupting disturbance. 

Location of study: 
University of Michigan Biological Station, Pellston, MI

All updated analysis can be found in the script:  
manuscript_script.R

R packages used: 
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
library(cvequality)
library(sjPlot)
library(sjmisc)
library(emmeans)

Data sources and metadata: 
Soil_respiration_ShRRED.csv: Tabular data file that contains all plot-level soil respiration, soil temperature and soil moisture data. 
-Unique ID: unique identifier for each row 
-Plot_ID: Plot within the experimental design (C1,C2,C3,C4,T1,T2,T3,T4) 
-Type: Disturbance treatment (Treatment) or Control 
-Plot_#: Plot nested within a treatment or control (1,2,3,4) 
-Location: Soil collar location nested within a plot (1,2,3,4)
-Efflux: soil respiration measurement (umols CO2 m^-2 sec^-1) 
-Soil_VWC: Soil Volumetric water content to 20 cm (%) 
-Soil_T: Soil temperature to 7cm (Celcius) 
-Date: Date of measurement

CalculatedN_ShRRED_REU2022.csv: Plot-level soil nitrogen data including ammonium and nitrate concentrations. 
-	Unique ID: unique identifier for each row 
-	Date_Collected: Date of field collection
-	Treatment: Disturbance (Treatment) or Control 
-	Plot: Plot within the experimental design (C1,C2,C3,C4,T1,T2,T3,T4) 
-	Average_Soil_Temp: 
-	Wet_sample_weight_analysis: weight of soil (g) used in analysis for nitrogen content 
-	Wet_sample_weight_drying: weight of soil (g) to be used to determine gravimetric water 
content 
-	Dry_soil_weight: Dry soil weight (g) after drying for 24 hours 
-	Nitrate: concentration (ug/L) 
-	Ammonium: concentration (ug/L)
-	grav_water_content: gravimetric water content (unitless) 
-	nitrate_calculated: Inorganic soil nitrate concentration mg NO3- per kg of dry soil 
-	ammonium_calculated: Inorganic soil nitrate concentration mg NH4+ per kg of dry soil 

NSC_summary_mathes.csv: Tree-level non-structural carbohydrate concentrations. 
-	Sample_Name: Unique indentifier for tree 
-	Date: Date of field collection 
-	Core_Type: Root (R) or Bole (B) sample collection 
-	Treatment_Type: Disturbance (Treatment) or Control
-	Plot_ID: Plot within the experimental design (C1,C2,C3,C4,T1,T2,T3,T4)
-	Tree_type: Focal (F) or Neighboring (N) 
-	NSC: Non-structural carbohydrate concentrations (proportion % of NSC concentration to sample dry mass) 

ShRRED_root_structure_6_26_22.csv:
ShRRED_root_structure_7_13_22.csv
ShRRED_root_structure_7_25_22.csv: Three different tabular datasets separated by date of collection (6-26, 7-13, 7-25) for fine-root structural data. 
 (selected columns used for analysis) 
-	Label: Unique ID
-	Type: Disturbance (Treatment) or Control
-	Plot: Plot within the experimental design (C1,C2,C3,C4,T1,T2,T3,T4)
-	Dry_root_weight_g_measured: root weight 
-	Dry_root_adjusted_g: root weight after ash-free adjustments 
-	SoilVol.m3: Soil sample volume m3
-	Length.cm: Total sample length (cm)
-	ProjArea.cm2: Total area of sample (cm2) 
-	SurfArea.cm2: Surface area of sample (cm2)
-	AvgDiam.mm: Average diameter of roots in sample (mm)
-	RootVolume.cm3: Total root volume of sample (cm3) 
-	Tips: Number of root tips in the sample
-	Forks: Number of root forks in the sample 
-	Crossings: number of crossings in sample 

ShRRED_Stem_map_data.csv: stem map and tree species information.  
-	Plot_ID: Plot within the experimental design (C1,C2,C3,C4,T1,T2,T3,T4)
-	Type: Disturbance (Treatment) or Control
-	Plot_Haglof: Plot ID in Haglof Instrument 
-	Plot_radius: Plot radius from center (m) 
-	Tree_Nr: Unique tree number within each plot
-	Tree_Spc: Tree species code: 
o	QURU: Quercus rubra 
o	PIST: Pinus strobus
o	PIRE: Pinus resinosa
o	ACRU: Acer rubrum 
o	POGR: Populus grandidentata
o	FAGR: Fagus grandifolia 
-	Strata: Canopy position: Canopy (can) or Subcanopy (Sub) 
-	Tree_Dia: Diameter at breast Height (cm) 
-	Tree_PosTex1, Tree_PosTex2, Tree_PosTex3, Tree_Local_x, Tree_Local_y, Tree_Local_Dist, Tree_Local_Angle, Tree_Angle_ToPlotCenter: Haglof measurements to spatially locate tree within the plot 
-	Latitude
-	Longitude

<img width="468" height="636" alt="image" src="https://github.com/user-attachments/assets/1a8a4be2-83a3-4920-97cb-caf5a03e020d" />


    








