## Prelim observations by JZB for whondrs 2021-2022
## Code from Danczakre github with modifications

## NMDS of Surface Water putative biogeochemical transformations of FTICR-MS data

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
# SW Only 
############################################################################################

# read in npoc and metadata
npoc_sw = read.csv("WHONDRS_S19S_SW_NPOC.csv", 
                   stringsAsFactors = F, row.names = 1)
metadata = read.csv("WHONDRS_S19S_Metadata_v3.csv", 
                    stringsAsFactors = F, row.names = 1)

# Change npoc samples names to match metadata, match metadata
npoc_sw$Samples = str_extract(npoc_sw$Sample_ID, "S19S_[0-9]{4}")
metadata$Samples = as.character(rownames(metadata))

#match npoc sample IDs to factors file
npoc_sw$SampleID <- str_replace(npoc_sw$Sample_ID, "-", "\\.")
npoc_sw$SampleID <- sapply(npoc_sw$SampleID,
       function(x) paste0("Sample_", x, "_p05"), USE.NAMES=FALSE)

####################################
# SW metadata and factors organize #
####################################
#subset sw from factors file 
factors_sw <- factors %>%
  filter(Type == "Surface Water")

#drop extra sample 0099_ICR-2 in npoc_sw
factors_npoc_sw <- semi_join(factors_sw, npoc_sw)

#add npoc data to factors data
factors_npoc_sw <- npoc_sw %>%
  select(SampleID, X00681_NPOC_mg_per_L_as_C) %>%
  left_join(factors_npoc_sw, npoc_sw, by = "SampleID",type="left", match="first") %>%
  distinct()

#change "Below range.." to NA
factors_npoc_sw$X00681_NPOC_mg_per_L_as_C[factors_npoc_sw$X00681_NPOC_mg_per_L_as_C 
                                           == "Below_Range_Less_Than_0.45"] <- NA

#subset all other relevant metadata
rel.meta.sw <- metadata[,c(8:11,14,15,17,20,21,41:44,53,57,65)]

#merge relevant metadata
factors_npoc_sw <- inner_join(factors_npoc_sw,rel.meta.sw, by = "Samples")

#weird but needs to happen
#clean up later
#ignore NA warning
factors_npoc_sw$X00681_NPOC_mg_per_L_as_C <- as.numeric(factors_npoc_sw$X00681_NPOC_mg_per_L_as_C)
factors_npoc_sw$SW_pH <- as.numeric(factors_npoc_sw$SW_pH)
factors_npoc_sw$SW_Temp_degC <- as.numeric(factors_npoc_sw$SW_Temp_degC)
factors_npoc_sw$DO_perc.sat <- as.numeric(factors_npoc_sw$DO_perc.sat)
factors_npoc_sw$DO_mg.per.L <- as.numeric(factors_npoc_sw$DO_mg.per.L)

################################
# SW put biogeo trans organize #
################################
#subset sw from trans
sed.loc = which(factors$Type == "Sediment")
trans_sw = trans[-sed.loc,]
rm(sed.loc)

# Converting to rel. abund.
trans_sw = as.data.frame(apply(trans_sw, 2, function(x) x/sum(x)))
trans_sw = trans_sw*100

############################################################################################
# SW NMDS Analyses
############################################################################################
# relative abundance as B-C dist
dist <- vegdist(trans_sw, "bray", na.rm = TRUE)

## perform ordination using a 'sample x sample' distance matrix, visualize w/ shepard plot
set.seed(2)
ord = metaMDS(dist, k = 10, trymax = 20, trace = FALSE)
ord
stressplot(ord)
plot(ord)

# perform envfit to find correlations between your data shape and \
# environmental variables that you measured
env_fit <- envfit(ord$points, factors_npoc_sw, perm=999, na.rm=TRUE)
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
factors <- factors[order(factors_sw$Samples),]
data.scores$type = factors_npoc_sw$Type
data.scores$streamorder = factors_npoc_sw$Stream_Order
data.scores$city = factors_npoc_sw$City
data.scores$country = factors_npoc_sw$Country
data.scores$hgm = factors_npoc_sw$Hydrogeomorphology
data.scores$sediment = factors_npoc_sw$Sediment
data.scores$contam = factors_npoc_sw$Contamination.Source.Upstream

## clean up the dataframe strings
data.scores = as.data.frame(data.scores, stringsAsFactors = TRUE)
data.scores$sample <- factor(data.scores$sample, levels = data.scores$sample)

## In addition to env data, let's find dispersion of your data (samples)
# using same distance matrix that was used for NMDS
heatmap(as.matrix(dist))
dist.clust <- hclust(dist,method="single")
plot(dist.clust) #super uggo don't recommend

# Image ord with SW as pretty
# fill color is pretty arbitrary at this point
p.nmds <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill=country), colour="black", pch=21, size = 3) +
  geom_segment(data=env_fit_df, aes(x=0, xend=env_fit_df$MDS1, y=0, yend=env_fit_df$MDS2), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data=env_fit_df, aes(x=env_fit_df$MDS1, y=env_fit_df$MDS2, label=vector),
            inherit.aes = FALSE, size=3) +
  theme_bw()
print(p.nmds)

## Done for now ##