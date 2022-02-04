## Prelim observations by JZB for whondrs 2021-2022
## Code from Danczakre github with modifications

## NMDS of all putative biogeochemical transformations of FTICR-MS data

library(dplyr)
library(tidyverse)
library(vegan)
setwd("~/Documents/whondrs-prelim-analyses/put-trans-nmds")

############################################################################################
# Data Prep 
############################################################################################

# read in putative biogeo transformations
trans = read.csv("S19S_Sed-Water_Trans_Profiles.csv", row.names = 1)

#remove the mass column
trans = trans[,-which(colnames(trans) %in% "Mass")]
#remove any rows/cols with only 0s
trans <- trans[, colSums(trans !=0)>0]
trans <- trans[rowSums(trans[])>0,]

#make a factors file
factors = data.frame(SampleID = colnames(trans),
                     Type = "Surface Water", stringsAsFactors = F)
factors$Type[grep("Sed_Field", factors$SampleID)] = "Sediment"
# Change samples names to match metadata, match metadata
factors$Samples = str_extract(factors$SampleID, "S19S_[0-9]{4}")

#transpose the rows and columns of abundance data
rownames(trans) <- factor(rownames(trans), levels = rownames(trans)) # fix row order
trans <- as.data.frame(t(trans[rowSums(trans)!=0, ]))

############################################################################################
# Metadata prep
############################################################################################

# read in npoc and metadata
npoc_sed = read.csv("WHONDRS_S19S_Sediment_NPOC.csv", 
                    stringsAsFactors = F, row.names = 1)
npoc_sw = read.csv("WHONDRS_S19S_SW_NPOC.csv", 
                   stringsAsFactors = F, row.names = 1)
metadata = read.csv("WHONDRS_S19S_Metadata_v3.csv", 
                    stringsAsFactors = F, row.names = 1)

#remove rows called S19S_num_Sed_Inc...
npoc_sed <- npoc_sed[!grepl("INC",npoc_sed$Sample_ID),]
# and there's ONE ROW with "Inc". you don't want to know how long this took me
npoc_sed <- npoc_sed[!grepl("Inc",npoc_sed$Sample_ID),]

# Change npoc samples names to match metadata, match metadata
npoc_sed$Samples = str_extract(npoc_sed$Sample_ID, "S19S_[0-9]{4}")
metadata$Samples = as.character(rownames(metadata))
npoc_sw$Samples = str_extract(npoc_sw$Sample_ID, "S19S_[0-9]{4}")
metadata$Samples = as.character(rownames(metadata))

#match npoc sample IDs to factors file
npoc_sed$SampleID <- str_replace(npoc_sed$Sample_ID, "-", "\\.")
npoc_sed$SampleID <- sapply(npoc_sed$SampleID,
                            function(x) paste0("Sample_", x, "_P2"), USE.NAMES=FALSE)
#match npoc sample IDs to factors file
npoc_sw$SampleID <- str_replace(npoc_sw$Sample_ID, "-", "\\.")
npoc_sw$SampleID <- sapply(npoc_sw$SampleID,
                           function(x) paste0("Sample_", x, "_p05"), USE.NAMES=FALSE)

#join sed and sw npoc datasets
npoc <- full_join(npoc_sw, npoc_sed)

#################################
# metadata and factors organize #
#################################

#observations are different between factors file and npoc data\
#need to find out which are different
diff <- subset(npoc,!(SampleID%in%factors$SampleID))

#drop extra sample 0099_ICR-2 in npoc_sw
factors_npoc <- semi_join(factors, npoc)

#add npoc data to factors data
factors_npoc <- npoc %>%
  select(SampleID, X00681_NPOC_mg_per_L_as_C) %>%
  left_join(factors_npoc, npoc, by = "SampleID",type="left", match="first") %>%
  distinct()

#change "Below range.." to NA
factors_npoc$X00681_NPOC_mg_per_L_as_C[factors_npoc$X00681_NPOC_mg_per_L_as_C 
                                           == "Above_Range_Greater_Than_22"] <- NA

#subset all other relevant metadata
#these were just chosen by hand/interest.
rel.meta <- metadata[,c(8:11,14,15,17,20,21,41:44,53,57,65)]

#merge relevant metadata
factors_npoc <- inner_join(factors_npoc,rel.meta, by = "Samples")

#weird but needs to happen
#clean up later
#ignore NA warning
factors_npoc$X00681_NPOC_mg_per_L_as_C <- as.numeric(factors_npoc$X00681_NPOC_mg_per_L_as_C)
factors_npoc$SW_pH <- as.numeric(factors_npoc$SW_pH)
factors_npoc$SW_Temp_degC <- as.numeric(factors_npoc$SW_Temp_degC)
factors_npoc$DO_perc.sat <- as.numeric(factors_npoc$DO_perc.sat)
factors_npoc$DO_mg.per.L <- as.numeric(factors_npoc$DO_mg.per.L)

#############################
# put biogeo trans organize #
#############################

#filter the matrix trans_sed to match the factors_npoc_sed file
include_list <- factors_npoc$SampleID
trans_filtered <- trans[include_list, ]

# Converting to rel. abund.
trans_filtered = as.data.frame(apply(trans_filtered, 2, function(x) x/sum(x)))
trans_filtered = trans_filtered*100

############################################################################################
# SW NMDS Analyses
############################################################################################
# relative abundance as B-C dist
dist <- vegdist(trans_filtered, "bray", na.rm = TRUE)

## perform ordination using a 'sample x sample' distance matrix, visualize w/ shepard plot
set.seed(2)
ord = metaMDS(dist, k = 10, trymax = 20, trace = FALSE)
ord
stressplot(ord)
plot(ord)

# perform envfit to find correlations between your data shape and \
# environmental variables that you measured
env_fit <- envfit(ord$points, factors_npoc, perm=999, na.rm=TRUE)
env_fit_df = as.data.frame(env_fit$vectors$arrows*sqrt(env_fit$vectors$r))
env_fit_df$vector = rownames(env_fit_df)
## View which biogeocheimcal variables are significant in explaining \
## the variability of your data shape
env_fit

## extract the scaling scores from your NMDS and convert to a data.frame
data.scores = as.data.frame(scores(ord)) 
## create a column of sample names, from the rownames of data.scores
data.scores$sample = rownames(data.scores)

## add group info to the score dataframe from your metadata table
data.scores <- data.scores[order(data.scores$sample),]
factors <- factors[order(factors$Samples),]
data.scores$type = factors_npoc$Type
data.scores$streamorder = factors_npoc$Stream_Order
data.scores$city = factors_npoc$City
data.scores$country = factors_npoc$Country
data.scores$hgm = factors_npoc$Hydrogeomorphology
data.scores$sediment = factors_npoc$Sediment
data.scores$contam = factors_npoc$Contamination.Source.Upstream

## clean up the dataframe strings
data.scores = as.data.frame(data.scores, stringsAsFactors = TRUE)
data.scores$sample <- factor(data.scores$sample, levels = data.scores$sample)

## In addition to env data, let's find dispersion of your data (samples)
# using same distance matrix that was used for NMDS
heatmap(as.matrix(dist))
dist.clust <- hclust(dist,method="single")
plot(dist.clust) #super uggo don't recommend

# Image ord with SW as pretty
# Fill color is arbitrary at this point
p.nmds <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill=sediment), colour="black", pch=21, size = 3) +
  geom_segment(data=env_fit_df, aes(x=0, xend=env_fit_df$MDS1, y=0, yend=env_fit_df$MDS2), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data=env_fit_df, aes(x=env_fit_df$MDS1, y=env_fit_df$MDS2, label=vector),
            inherit.aes = FALSE, size=3) +
  theme_bw()
print(p.nmds)

## Done for now ##