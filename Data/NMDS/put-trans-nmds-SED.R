## Prelim observations by JZB for whondrs 2021-2022
## Code from Danczakre github with modifications

## NMDS of Sediment putative biogeochemical transformations of FTICR-MS data

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
# SED Only 
############################################################################################

# read in npoc and metadata
npoc_sed = read.csv("WHONDRS_S19S_Sediment_NPOC.csv", 
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

#match npoc sample IDs to factors file
npoc_sed$SampleID <- str_replace(npoc_sed$Sample_ID, "-", "\\.")
npoc_sed$SampleID <- sapply(npoc_sed$SampleID,
                           function(x) paste0("Sample_", x, "_P2"), USE.NAMES=FALSE)

#####################################
# SED metadata and factors organize #
#####################################
#subset sw from factors file 216 obs
factors_sed <- factors %>%
  filter(Type == "Sediment")

#observations are different between factors file (216 obs) and npoc data (282 obs)\
#need to find out which are different (78 obs)
diff <- subset(npoc_sed,!(SampleID%in%factors_sed$SampleID))

#drop obs which are different
factors_npoc_sed <- semi_join(factors_sed, npoc_sed)

#add npoc data to factors data
factors_npoc_sed <- npoc_sed %>%
  select(SampleID, X00681_NPOC_mg_per_L_as_C) %>%
  left_join(factors_npoc_sed, npoc_sed, by = "SampleID",type="left", match="first") %>%
  distinct()

#change "Below range.." to NA
factors_npoc_sed$X00681_NPOC_mg_per_L_as_C[factors_npoc_sed$X00681_NPOC_mg_per_L_as_C 
                                           == "Above_Range_Greater_Than_22"] <- NA

#subset all other relevant metadata
rel.meta.sed <- metadata[,c(8:11,14,15,17,20,21,41:44,53,57,65)]

#merge relevant metadata
factors_npoc_sed <- inner_join(factors_npoc_sed,rel.meta.sed, by = "Samples")

#weird but needs to happen
#clean up later
#ignore NA warning
factors_npoc_sed$X00681_NPOC_mg_per_L_as_C <- as.numeric(factors_npoc_sed$X00681_NPOC_mg_per_L_as_C)
factors_npoc_sed$SW_pH <- as.numeric(factors_npoc_sed$SW_pH)
factors_npoc_sed$SW_Temp_degC <- as.numeric(factors_npoc_sed$SW_Temp_degC)
factors_npoc_sed$DO_perc.sat <- as.numeric(factors_npoc_sed$DO_perc.sat)
factors_npoc_sed$DO_mg.per.L <- as.numeric(factors_npoc_sed$DO_mg.per.L)

################################
# SW put biogeo trans organize #
################################
#subset sw from trans
sw.loc = which(factors$Type == "Surface Water")
trans_sed = trans[-sw.loc,]
rm(sw.loc)

#filter the matrix trans_sed to match the factors_npoc_sed file
include_list <- factors_npoc_sed$SampleID
trans_sed <- trans_sed[include_list, ]

# Converting to rel. abund.
trans_sed = as.data.frame(apply(trans_sed, 2, function(x) x/sum(x)))
trans_sed = trans_sed*100

############################################################################################
# SW NMDS Analyses
############################################################################################
# relative abundance as B-C dist
dist <- vegdist(trans_sed, "bray", na.rm = TRUE)

## perform ordination using a 'sample x sample' distance matrix, visualize w/ shepard plot
set.seed(2)
ord = metaMDS(dist, k = 10, trymax = 20, trace = FALSE)
ord
stressplot(ord) ## Looks kinda weird--need to address (from subsetting?)
plot(ord) ## Looks kinda weird--need to address (from subsetting?)

# perform envfit to find correlations between your data shape and \
# environmental variables that you measured
env_fit <- envfit(ord$points, factors_npoc_sed, perm=999, na.rm=TRUE)
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
factors <- factors[order(factors_sed$Samples),]
data.scores$type = factors_npoc_sed$Type
data.scores$streamorder = factors_npoc_sed$Stream_Order
data.scores$city = factors_npoc_sed$City
data.scores$country = factors_npoc_sed$Country
data.scores$hgm = factors_npoc_sed$Hydrogeomorphology
data.scores$sediment = factors_npoc_sed$Sediment
data.scores$contam = factors_npoc_sed$Contamination.Source.Upstream

## clean up the dataframe strings
data.scores = as.data.frame(data.scores, stringsAsFactors = TRUE)
data.scores$sample <- factor(data.scores$sample, levels = data.scores$sample)

## In addition to env data, let's find dispersion of your data (samples)
# using same distance matrix that was used for NMDS
heatmap(as.matrix(dist))
dist.clust <- hclust(dist,method="single")
plot(dist.clust) #super uggo don't recommend

# Image ord with SW as pretty
# fill color is arbitrary at this point
p.nmds <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill=country), colour="black", pch=21, size = 3) +
  geom_segment(data=env_fit_df, aes(x=0, xend=env_fit_df$MDS1, y=0, yend=env_fit_df$MDS2), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data=env_fit_df, aes(x=env_fit_df$MDS1, y=env_fit_df$MDS2, label=vector),
            inherit.aes = FALSE, size=3) +
  theme_bw()
print(p.nmds)

## Done for now ##