## This script is uuuuugggggglllllyyy and needs some TLC before being uploaded. 
## General order of operations is read in the metadata, do a ton of manual 
## cleaning, then merge with respiration data and explore basic relationships
## between environmental context (meta) and respiration. Hopefully of some use
## to the WHONDRS Topic 3 team.... 
##
## 2021-11-02 (updated 2022-02-03)
## Peter Regier 
## 
# #############
# #############

# 1. Setup environment ---------------------------------------------------------

## Clear environment
rm(list = ls())

## Load packages
require(pacman)
p_load(cowplot, tidyverse, parsedate, lubridate, maps)

## Using local versions of the files, but setting a filepath so easier to generalize
local_filepath <- "data/"

## Set ggplot theme and some basic dimensions for graphs
theme_set(theme_bw())
plot_width = 5
plot_height = 3

## Hi! I'm the only function in the whole code cause everything is manual for now
give.n <- function(x){
  return(c(y = median(x) / 100, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


# 2. Import and format data ----------------------------------------------------
## Read in metadata
meta <- read_csv(paste0(local_filepath, "WHONDRS_S19S_Metadata_v3.csv"), skip = 1) %>% 
  ## rename all the interesting columns
  rename("id" = "Sample ID",
         "date" = "Sampling date",
         "stream_name" = "Stream name", 
         "stream_order" = "Stream order", 
         "veg_type" = "General vegetation type (select up to 2 if mixed)", 
         "river_morphology" = "General hydrogeomorphology of river", 
         "intermittent" = "Is the river intermittent or perennial?", 
         "river_gradient" = "General gradient of river", 
         "sed_type" = "What is the dominant sediment type of the streambed?", 
         "latitude" = "Latitude of midstream site (decimal degrees)"  , 
         "longitude" = "Longitude of midstream site (decimal degrees)", 
         "water_height" = "Height of water column at midstream sampling location (cm)", 
         "sed_deposition" = "Midstream sediment location depositional zone type", 
         "algal_mats" = "At the midstream site what is the algal mat coverage?", 
         "macrophytes" = "At the midstream site what is the macrophyte coverage?",
         "ph" = "pH value from colorimetric pH strip",
         "do_mgl" = "Dissolved oxygen at 50% depth at downstream sampling location (mg/L) if available", 
         "flow_variation" = "What are the primary sources of variation in river flow?",
         "contamination" = "Are there any upstream contamination sources?") %>% 
  ## Down-select to interesting columns
  dplyr::select(id, date, stream_name, stream_order, veg_type, river_morphology, 
         intermittent, river_gradient, sed_type, latitude, longitude, 
         water_height, sed_deposition, algal_mats, macrophytes, ph, do_mgl, 
         flow_variation, contamination) %>% 
  ## Clean up NAs
  mutate_all(funs(str_replace_all(., "Not_Provided", "NA"))) %>% 
  ## Clean up macrophyte levels
  mutate_all(funs(str_replace_all(., "Submerged aquatic macrophytes abundant without full coverage", "High"))) %>% 
  ## Clean up sediment type levels
  mutate(sed_type_recode = recode(as.factor(sed_type), 
                                  "Silt/mud (<0.0625mm)" = "Silt", 
                                  "Gravel/cobble (>2mm)" = "Gravel", 
                                  "Bedrock (primarily)" = "Bedrock"),
         ## Clean up river gradient levels
         river_gradient_recode = recode(as.factor(river_gradient), 
                                        "Relatively flat/gentle gradient (e.g. valleys)" = "Gentle gradient", 
                                        "Relatively flat/gentle gradient (e.g.. valleys)" = "Gentle gradient", 
                                        "Relatively steep gradient (e.g. mountainous or hilly terrain)" = "Steep gradient", 
                                        "Relatively steep gradient (e.g.. mountainous or hilly terrain)" = "Steep gradient"), 
         ## Clean up veg type levels
         veg_type_recode = recode(as.factor(veg_type), 
                                  "Not vegetated" = "None",
                                  "Shrub. Grass" = "Shrub", 
                                  "Broadleaf deciduous tree. Shrub" = "Broadleaf deciduous",
                                  "Needleleaf evergreen tree. Shrub" = "Needleleaf evergreen", 
                                  "Needleleaf evergreen tree" = "Needleleaf evergreen", 
                                  "Broadleaf deciduous tree" = "Broadleaf deciduous", 
                                  "Needleleaf deciduous tree. Broadleaf deciduous tree" = "Needleleaf deciduous",
                                  "Broadleaf deciduous tree. Grass" = "Broadleaf deciduous", 
                                  "Broadleaf deciduous tree. Crop" = "Broadleaf deciduous", 
                                  "Needleleaf evergreen tree. Broadleaf deciduous tree" = "Needleleaf evergreen", 
                                  "Broadleaf evergreen tree. Broadleaf deciduous tree" = "Broadleaf evergreen", 
                                  "Needleleaf evergreen tree. Grass" = "Needleleaf evergreen", 
                                  "Grass. Not vegetated" = "Grass", 
                                  "Broadleaf evergreen tree" = "Broadleaf evergreen", 
                                  "Broadleaf evergreen tree. Shrub" = "Broadleaf evergreen", 
                                  "Needleleaf evergreen tree. Broadleaf evergreen tree" = "Needleleaf evergreen", 
                                  "Needleleaf deciduous tree. Grass" = "Needleleaf deciduous"), 
         ## Second round of veg type cleanup: going for simpler yet
         veg_type_recode2 = recode(as.factor(veg_type), 
                                  "Not vegetated" = "None",
                                  "Shrub. Grass" = "Shrub", 
                                  "Broadleaf deciduous tree. Shrub" = "Broadleaf",
                                  "Needleleaf evergreen tree. Shrub" = "Needleleaf", 
                                  "Needleleaf evergreen tree" = "Needleleaf", 
                                  "Broadleaf deciduous tree" = "Broadleaf", 
                                  "Needleleaf deciduous tree. Broadleaf deciduous tree" = "Needleleaf",
                                  "Broadleaf deciduous tree. Grass" = "Broadleaf", 
                                  "Broadleaf deciduous tree. Crop" = "Broadleaf", 
                                  "Needleleaf evergreen tree. Broadleaf deciduous tree" = "Needleleaf", 
                                  "Broadleaf evergreen tree. Broadleaf deciduous tree" = "Broadleaf", 
                                  "Needleleaf evergreen tree. Grass" = "Needleleaf", 
                                  "Grass. Not vegetated" = "Grass", 
                                  "Broadleaf evergreen tree" = "Broadleaf", 
                                  "Broadleaf evergreen tree. Shrub" = "Broadleaf", 
                                  "Needleleaf evergreen tree. Broadleaf evergreen tree" = "Needleleaf", 
                                  "Needleleaf deciduous tree. Grass" = "Needleleaf"), 
         ## Third round of veg type cleanup: going for simplest
         veg_type_recode3 = recode(as.factor(veg_type), 
                                   "Not vegetated" = "None",
                                   "Shrub. Grass" = "Shrub", 
                                   "Broadleaf deciduous tree. Shrub" = "deciduous",
                                   "Needleleaf evergreen tree. Shrub" = "evergreen", 
                                   "Needleleaf evergreen tree" = "evergreen", 
                                   "Broadleaf deciduous tree" = "deciduous", 
                                   "Needleleaf deciduous tree. Broadleaf deciduous tree" = "deciduous",
                                   "Broadleaf deciduous tree. Grass" = "deciduous", 
                                   "Broadleaf deciduous tree. Crop" = "deciduous", 
                                   "Needleleaf evergreen tree. Broadleaf deciduous tree" = "evergreen", 
                                   "Broadleaf evergreen tree. Broadleaf deciduous tree" = "evergreen", 
                                   "Needleleaf evergreen tree. Grass" = "evergreen", 
                                   "Grass. Not vegetated" = "Grass", 
                                   "Broadleaf evergreen tree" = "evergreen", 
                                   "Broadleaf evergreen tree. Shrub" = "evergreen", 
                                   "Needleleaf evergreen tree. Broadleaf evergreen tree" = "evergreen", 
                                   "Needleleaf deciduous tree. Grass" = "deciduous"), 
         ## Clean up algal mats levels
         algal_mats_recode = recode(as.factor(algal_mats), 
                                   "Low. No" = "No",
                                   "Low. The algal mat coverage is not on sediment but nearby" = "Low", 
                                   "the water was very turbid. and no visual inspection was possible. I assume no algal mats present" = "No"))

## Create a dataset for the respiration data
respiration <- read_csv(paste0(local_filepath, "WHONDRS_S19S_Sediment_Incubations_Respiration_Rates.csv"))  %>% 
  dplyr::mutate(id = substr(Sample_ID, 1, 9)) %>% 
  drop_na() %>% 
  group_by(id) %>% 
  summarize("slope" = mean(slope_of_the_regression), 
         "mgL_per_hr" = mean(rate_mg_per_L_per_h))

## Join respiration and metadata into a single dataset
df <- left_join(meta, respiration, by = "id")

## Write out this dataset in case it's of use later
write_csv(df, paste0(local_filepath, "joined_metadata_respiration.csv"))

# 3. Make a whole truckload of plots -------------------------------------------

p0 <- ggplot(df %>% drop_na(), aes(stream_order, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75))
ggsave(paste0(local_filepath, "graphs/respiration_v_stream_order.png"), width = plot_width, height = plot_height)

p1 <- ggplot(df%>% drop_na(), aes(sed_type_recode, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75))
ggsave(paste0(local_filepath, "graphs/respiration_v_sed_type.png"), width = plot_width, height = plot_height)

plot_grid(p0, p1, nrow = 1) 
ggsave(paste0(local_filepath, "graphs/good_categories.png"), width = plot_width * 1.5, height = plot_height)


p2.1 <- ggplot(df %>% drop_na(), aes(macrophytes, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75))
ggsave(paste0(local_filepath, "graphs/respiration_v_macrophytes.png"), width = plot_width, height = plot_height)

p2.2 <- ggplot(df %>% drop_na(), aes(algal_mats_recode, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75))
ggsave(paste0(local_filepath, "graphs/respiration_v_algal.png"), width = plot_width, height = plot_height)

plot_grid(p2.1, p2.2, nrow = 1) 
ggsave(paste0(local_filepath, "graphs/good_categories_bio.png"), width = plot_width * 1.5, height = plot_height)



p3 <- ggplot(df %>% drop_na(), aes(as.numeric(longitude), mgL_per_hr)) + geom_point() + geom_smooth(se = F)
ggsave(paste0(local_filepath, "graphs/respiration_v_longitude.png"), width = plot_width, height = plot_height)

p4 <- ggplot(df %>% drop_na(), aes(as.numeric(ph), mgL_per_hr)) + geom_point() + geom_smooth(se = F)
ggsave(paste0(local_filepath, "graphs/respiration_v_ph.png"), width = plot_width, height = plot_height)

plot_grid(p3, p4, nrow = 1)
ggsave(paste0(local_filepath, "graphs/good_correlations.png"), width = plot_width * 1.5, height = plot_height)


q0 <- ggplot(df %>% drop_na(), aes(river_gradient_recode, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + xlab("") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

q1 <- ggplot(df %>% drop_na(), aes(river_morphology, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + xlab("") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

q2 <- ggplot(df %>% drop_na(), aes(veg_type_recode, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + xlab("") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

plot_grid(q0, q1, q2, nrow = 1, align = "hv")
ggsave(paste0(local_filepath, "graphs/bad_categories.png"), width = plot_width * 2, height = plot_height * 1.2)




