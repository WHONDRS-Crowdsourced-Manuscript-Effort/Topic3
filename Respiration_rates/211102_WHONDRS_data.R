## This script does some metadata cleaning, merges with respiration data, then
## plots respiration rates versus a number of ecosystem variables of interest 
##
## 2021-11-02 (updated 2022-05-12)
## Peter Regier 
## 
# #############
# #############

# 1. Setup environment ---------------------------------------------------------

## Load packages
require(pacman)
p_load(cowplot, # plot_grid
       tidyverse, # keep things tidy
       parsedate, # handle POSIX easier
       lubridate,  # tidy date work
       ggpubr) # stat_compare_means

## Set ggplot theme
theme_set(theme_bw())


# 2. Import and format data ----------------------------------------------------
## Read in metadata
meta <- read_csv("data/WHONDRS_S19S_Metadata_v3.csv", skip = 1) %>% 
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
respiration <- read_csv("data/WHONDRS_S19S_Sediment_Incubations_Respiration_Rates.csv") %>% 
  dplyr::mutate(id = substr(Sample_ID, 1, 9)) %>% 
  drop_na() %>% 
  group_by(id) %>% 
  summarize("slope" = mean(slope_of_the_regression), 
         "mgL_per_hr" = mean(rate_mg_per_L_per_h))

## Join respiration and metadata into a single dataset
df <- left_join(meta, respiration, by = "id")

## Write out this dataset in case it's of use later
write_csv(df, "data/joined_metadata_respiration.csv")


# 3. Make some plots -------------------------------------------

## I'm a function that prints counts on plots for each box
give.n <- function(x){
  return(c(y = median(x) / 100, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

## I'm a function that counts the number of NAs for a given variable
count_nas <- function(var){
  nrow(df %>% filter(is.na({{var}}) | {{var}} == "NA"))
}

## Plot repiration rates v stream order
p_stream_order <- ggplot(df %>% drop_na() %>% filter(stream_order != "NA"), 
             aes(stream_order, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + 
  labs(x = "Stream order", y = "Respiration rate (mg/L/hr)") + 
  stat_compare_means(label = "p.format", label.x = 7, label.y = 25) +
  annotate("text", x = 7, y = 23, label = paste0(count_nas(stream_order), " NAs"))

## Plot repiration rates v sediment type
p_sediment_type <- ggplot(df%>% drop_na(), aes(sed_type_recode, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + 
  labs(x = "Sediment type", y = "Respiration rate (mg/L/hr)") + 
  stat_compare_means(label = "p.format", label.x = 1, label.y = 25) +
  annotate("text", x = 1, y = 23, label = paste0(count_nas(sed_type_recode), " NAs"))

## Plot repiration rates v macrophytes
p_macrophytes <- ggplot(df %>% drop_na(), aes(macrophytes, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + 
  labs(x = "Macrophyte coverage", y = "Respiration rate (mg/L/hr)") + 
  stat_compare_means(label = "p.format", label.x = 1, label.y = 25) +
  annotate("text", x = 2, y = 23, label = paste0(count_nas(macrophytes), " NA"))

## Plot repiration rates v algal mats
p_algal_mats <- ggplot(df %>% drop_na(), aes(algal_mats_recode, mgL_per_hr)) + geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + 
  labs(x = "Algal mat coverage", y = "Respiration rate (mg/L/hr)") + 
  stat_compare_means(label = "p.format", label.x = 1, label.y = 25) +
  annotate("text", x = 1, y = 23, label = paste0(count_nas(algal_mats_recode), " NA"))

## Combine into a multi-plot object
plot_grid(p_stream_order, 
          plot_grid(p_sediment_type, 
                    p_macrophytes, 
                    p_algal_mats, nrow = 1), ncol = 1)

## Save that plot
ggsave("plots/220512_respiration_plots.png", 
       width = 6, height = 6)
