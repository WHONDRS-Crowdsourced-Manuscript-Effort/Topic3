### Data analysis ####
### This script is for cleaning the data anions and NPOC
### Topic 3####
### 04-Feb-22 Patricia Garcia
####################################

#libraries
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(ggpubr)
## Set the WD
##user shoud change this

setwd("data/")

### data import

Isotop <- read.csv("WHONDRS_S19S_SW_Isotopes.csv") 
NPOC <- read.csv("WHONDRS_S19S_SW_NPOC.csv")
Anions <- read.csv("WHONDRS_S19S_SWS_SpC_Anions_TN.csv")
NPOC_sed <- read.csv("WHONDRS_S19S_Sediment_NPOC.csv")


NPOC$NPOC <- as.numeric(NPOC$X000681_NPOC_mg_per_L_as_C)

### split the sample ID

res <- data.frame(str_split_fixed(Anions$Sample_ID, "_", 3))## data Anions

res$X2 <- as.numeric(res$X2)
res$SampleID <- paste0(res$X1,res$X2)

Anions$Sample_n <- res$X2
Anions$Sample_n <- res$SampleID
remove(res)

res <- data.frame(str_split_fixed(NPOC$Sample_ID, "_", 3)) ##data NPOC

res$SampleID <- paste0(res$X1,res$X2)
res$X2 <- as.numeric(res$X2)
NPOC$Sample_n <- res$X2
NPOC$Sample_n <- res$SampleID


remove(res)
NPOC_m <- NPOC [,c(2,4:5)]
NPOC_m$Sample_ID <- NULL

### data merge manually

Anion_SWDOC <- read_excel("Data_Anions_SWDOC.xlsx", 
                          +     col_types = c("text", "text", "numeric", 
                                              +         "numeric", "numeric", "numeric", 
                                              +         "numeric", "numeric", "numeric", 
                                              +         "text", "text", "numeric"))
summary(Anion_SWDOC)
res.cor <- rcorr(as.matrix((Anion_SWDOC [,c(3:9,12)]))) # correlacion como matriz
corrplot(res.cor$r, type="upper", order="hclust", 
         p.mat = res.cor$P, sig.level = 0.05, insig = "blank", tl.col="black", tl.srt=45)## grafico
res.cor$r

data.stan <- scale(Anion_SWDOC [,c(3:9,12)], center=TRUE)

res.pca <- PCA(data.stan, graph=FALSE)
eigenvalues <- res.pca$eig# eigenvalues
head(eigenvalues[, 1:3])

fviz_eig(res.pca, addlabels = TRUE)
fviz_pca_biplot(res.pca,  repel=TRUE, invisible = "quali")+theme_classic()


data.stan_i <- imputePCA(data.stan, ncp=2)

res.pca <- PCA(data.stan_i$completeObs, graph=FALSE)
eigenvalues <- res.pca$eig# eigenvalues
head(eigenvalues[, 1:3])
fviz_pca_biplot(res.pca,  repel=TRUE, invisible = "quali")+theme_classic()
fviz_pca_biplot(res.pca,  repel=TRUE, invisible = "quali", habillage=as.factor(Anion_SWDOC$Sample_ID...2))+theme_classic()


### use dataset anion_average

res.cor <- rcorr(as.matrix((anion_doc_SW[1:8])))
corrplot(res.cor$r, type="upper", order="hclust", 
         p.mat = res.cor$P, sig.level = 0.05, insig = "blank", tl.col="black", tl.srt=45)## grafico
res.cor$r

ggplot(anion_doc_SW, aes(x=US_Latitude_dec.deg, y= NPOC))+ geom_point()

ggplot(anion_doc_SW, aes(x=US_Latitude_dec.deg, y= NPOC, color=Stream_Order ))+ geom_point() + labs(x="Latitude" , y="DOC (mg/L)")


ggplot(anion_doc_SW, aes(x=US_Longitude_dec.deg, y= NPOC, color=Stream_Order ))+ geom_point() + labs(x="Longitude" , y="DOC (mg/L)")


ggplot(anion_doc_SW, aes(x=NPOC, y= TN, color=Stream_Order ))+ geom_point() + labs(x="DOC (mg/L)" , y="TN")

ggplot(anion_doc_SW, aes(x=SpC, y= NPOC, color=Stream_Order ))+ geom_point() + labs(x="Specific Conductivity" , y="DOC (mg/L)")

new <- anion_doc_SW %>%  filter(SpC<1000)

ggplot(new, aes(x=SpC, y= NPOC, color=Stream_Order ))+ geom_point() + labs(x="Specific Conductivity" , y="DOC (mg/L)")

ggplot(new, aes(x=SpC, y= NPOC, color=Stream_Order ))+ geom_point() + labs(x="Specific Conductivity" , y="DOC (mg/L)") + scale_color_brewer(palette="Blues")


ggplot(new, aes(x=NO3, y= TN, color=Stream_Order ))+ geom_point() + labs(x="NO3" , y="TN") 


ggbarplot(data=new, x=Stream_Order, y=NPOC)
ggbarplot(data=new, x="Stream_Order", y="NPOC", add = c("mean_sd"), color="red", fill="darkred", xlab=c("Stream order"), ylab=c("DOC mg /L"))