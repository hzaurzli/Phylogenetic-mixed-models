#===============================================================================
#-------------------------------------------------------------------------------
#       Assessing heritability in susceptibility to ISP using spot tests
#-------------------------------------------------------------------------------
#===============================================================================

# 0.1 Purpose-------------------------------------------------------------------

'This script analyses data using plaque assays (refered to as spot tests in the
 data and script), optical density assays, and qPCR to assess the suseptibility 
 of 64 Staphylococcus isolates to a virus, ISP, using the phylogenetic mixed 
 model approach discussed in the following papers:
  
  - Hadfield JD, Nakagawa S. General quantitative genetic methods for 
    comparative biology: phylogenies, taxonomies and multi-trait models for 
    continuous and categorical characters. Journal of evolutionary biology. 
    2010 Mar;23(3):494-508.
  
  - Housworth EA, Martins EP, Lynch M. The phylogenetic mixed model. The 
    American Naturalist. 2004 Jan;163(1):84-96.
  
  - Lynch M. Methods for the analysis of comparative data in evolutionary 
    biology. Evolution. 1991 Aug;45(5):1065-80.
  
  - Losos JB. Seeing the Forest for the Trees: The Limitations of Phylogenies 
    in Comparative Biology: (American Society of Naturalists Address). The 
    American Naturalist. 2011 Jun 1;177(6):709-27.'

# ------------------------------------------------------------------------------
# -------------------------- [1] Initialisation --------------------------------
# ------------------------------------------------------------------------------
# 1.1 Load packages ---------------------------------------------------------

library(tidyverse); library(MCMCglmm); library(ape); library(ggtree); 
library(geiger); library(lme4); library(nortest); library(lmtest); library(gt);
library(tidyr); library(phytools); library(plotrix)

# 1.2 Set directory ------------------------------------------------------------

"Chose your working directory"

setwd("")

# 1.2 Load data ----------------------------------------------------------------

Data <- read.csv("PhageSusceptibilityData.csv", header = TRUE)

head(Data)

"Columns are as follows:
  - Sample: the name of the bacterial strain as needed to join with the 
            individual data sheets. 
  - Block: experimental block. 
  - ISP: measure of virus presence, i.e., PFU/uL for plaque assay, proportion 
         change in OD for the OD assay, and log10 fold change in viral load 
         for qPCR. 
  - trait: the method being used to assess ISP presence. 
  - Family: the distribution of the data, either threshold (binary) or normal
            here. 
  - PhotoshopName: the name to be displayed as a tiplabel on the figure.
  - animal: allows the model to calculate effect of phylogeny, must be the same
            as the tiplabels in the tree.
  - TreeOrder: order that the samples appear in on the ladderised tree, this is
               needed for plotting the figures. 
  - Sample_Full: allows the model to calculate the species specific effect."

Data$trait <- as.factor(Data$trait)

Data$trait <- factor(Data$trait, levels = c("OD", "qPCR", "ST_Cont", "ST_Bin"))

# 1.3 Load host phylogeny ------------------------------------------------------

Tree <- read.tree("Staph_Phylogeny.nwk")

ggtree(Tree, size = 0.25, right = T)+ 
  coord_cartesian(xlim = c(0, 0.6)) + 
  geom_treescale() +
  geom_tiplab(size = 2.5)

# 1.6 Tree lengths -------------------------------------------------------------

mean(diag(vcv(Tree)))

"Check - length of the whole phylogeny should be 0.477006"

# ------------------------------------------------------------------------------
# ---------------------- [2] Phylogenetic Mixed Models  ------------------------
# ------------------------------------------------------------------------------

"If you have previously run the models, skip to 3.1 Load models"

# 2.1 Define priors ------------------------------------------------------------

Prior.With.Species <- list(G = list(
  G1 = list(V = diag(4), nu = 4, alpha.mu = rep(0,4), alpha.V = diag(4) * 1000),
  G2 = list(V = diag(4), nu = 4, alpha.mu = rep(0,4), alpha.V = diag(4) * 1000)),
  R = list(V = diag(4), nu = 0.002, fix = 4))

Prior.Without.Species <- list(G = list(
  G1 = list(V = diag(4), nu = 4, alpha.mu = rep(0,4), alpha.V = diag(4) * 1000)),
  R = list(V = diag(4), nu = 0.002, fix = 4))

# 2.2 Run main model -----------------------------------------------------------

Ainv<-inverseA(Tree, scale=FALSE)$Ainv

"As the whole phylogeny (64 strains), within-aureus phylogeny, and 
 among-species phylogenies have different lengths, Ainv allows us to run the 
 model with the length of the tree (rather than the root-to-tip = 1). When 
 calcualting repeatability and phylogenetic heritability at the end, we then 
 scale for tree length in these equations to make the resulting numbers 
 comparable between trees of differing lengths."

WithSpeciesModel <- 
  MCMCglmm(ISP ~ trait, random = ~us(trait) : animal + us(trait) : Sample_Full,
           rcov = ~idh(trait) : units, ginv = list(animal=Ainv), prior = Prior.With.Species, data = Data,
           family = NULL, nitt = 13000000, thin = 5000, burnin = 3000000, pr = TRUE)

save(WithSpeciesModel, file = "WithSpeciesModel.Rdata")

WithoutSpeciesModel <- 
  MCMCglmm(ISP ~ trait, random = ~us(trait) : animal,
           rcov = ~idh(trait) : units, ginv = list(animal=Ainv), prior = Prior.Without.Species, data = Data,
           family = NULL, nitt = 13000000, thin = 5000, burnin = 3000000, pr = TRUE)

save(WithoutSpeciesModel, file = "WithoutSpeciesModel.Rdata")

# ------------------------------------------------------------------------------
# ---------------------------- [3] Model Exploration ---------------------------
# ------------------------------------------------------------------------------

"From this point onwards, the code will be shown as using the models run for 
 main dataframe. To look at any other models, the name of the model should be
 substituted and parts of the VCV matrix called checked. For more information 
 see the end of the script, Section 6 - Alternative Models."

# 3.1 Load models --------------------------------------------------------------

load("WithSpeciesModel.Rdata")

load("WithoutSpeciesModel.Rdata")

# 3.2 Model summaries ----------------------------------------------------------

options(max.print=100000)

summary(WithSpeciesModel) # Structure
View(summary(WithSpeciesModel$Sol)) # Fixed Effects (too big for Rstudio console, write to file TBC)
summary(WithSpeciesModel$VCV) # Variance/Covariance

summary(WithoutSpeciesModel) # Structure
View(summary(WithoutSpeciesModel$Sol)) # Fixed Effects (too big for Rstudio console, write to file TBC)
summary(WithoutSpeciesModel$VCV) # Variance/Covariance

# 3.3 Model plots --------------------------------------------------------------

plot(WithSpeciesModel)

plot(WithoutSpeciesModel)

# 3.4 Fixed Effects Plots ------------------------------------------------------

# Note: These are very long

plot(WithSpeciesModel$Sol)

plot(WithoutSpeciesModel$Sol)

# 3.5 Variance/Covariance Plots ------------------------------------------------

plot(WithSpeciesModel$VCV)

plot(WithoutSpeciesModel$VCV)

# 3.6 Autocorrelation ----------------------------------------------------------

"With species"

WithSpeciesModel.Autocorr <- data.frame(t(autocorr.diag(WithSpeciesModel$VCV, lags = c(0:10))))

print(WithSpeciesModel.Autocorr)

WithSpeciesModel.Autocorr$Comparison <- rownames(WithSpeciesModel.Autocorr)
WithSpeciesModel.Autocorr <- gather(WithSpeciesModel.Autocorr, key = "Lag", value = "Autocorrelation", -Comparison)
WithSpeciesModel.Autocorr$Lag <- as.numeric(str_split(WithSpeciesModel.Autocorr$Lag, pattern = "\\.", simplify = TRUE)[,2])
WithSpeciesModel.Autocorr$Comparison <- factor(WithSpeciesModel.Autocorr$Comparison, levels = unique(WithSpeciesModel.Autocorr$Comparison))

ggplot(data = WithSpeciesModel.Autocorr) +
  ggtitle("Autocorrelation of Parameter Estimates") +
  geom_line(mapping = aes(x = Lag, y = Autocorrelation)) +
  geom_point(mapping = aes(x = Lag, y = Autocorrelation)) +
  geom_hline(yintercept = 0.1, color = "red") +
  facet_wrap(~Comparison, ncol = 4)

"Without species"

WithoutSpeciesModel.Autocorr <- data.frame(t(autocorr.diag(WithoutSpeciesModel$VCV, lags = c(0:10))))

print(WithoutSpeciesModel.Autocorr)

WithoutSpeciesModel.Autocorr$Comparison <- rownames(WithoutSpeciesModel.Autocorr)
WithoutSpeciesModel.Autocorr <- gather(WithoutSpeciesModel.Autocorr, key = "Lag", value = "Autocorrelation", -Comparison)
WithoutSpeciesModel.Autocorr$Lag <- as.numeric(str_split(WithoutSpeciesModel.Autocorr$Lag, pattern = "\\.", simplify = TRUE)[,2])
WithoutSpeciesModel.Autocorr$Comparison <- factor(WithoutSpeciesModel.Autocorr$Comparison, levels = unique(WithoutSpeciesModel.Autocorr$Comparison))

ggplot(data = WithoutSpeciesModel.Autocorr) +
  ggtitle("Autocorrelation of Parameter Estimates") +
  geom_line(mapping = aes(x = Lag, y = Autocorrelation)) +
  geom_point(mapping = aes(x = Lag, y = Autocorrelation)) +
  geom_hline(yintercept = 0.1, color = "red") +
  facet_wrap(~Comparison, ncol = 4)

# ------------------------------------------------------------------------------
# --------------- [4] Repeatability and Phylogenetic Heritability --------------
# ------------------------------------------------------------------------------
# 4.1 Repeatability (without species effects) ----------------------------------
# 4.1.1. Collect animal variances ----------------------------------------------

"As previously mentioned, we are using *mean(diag(vcv(Tree))) here to scale the
 estimates of animal varaince by the length of the phylognetic tree used in the
 model. This allows us to compare across models with different tree lengths."

dimnames(WithoutSpeciesModel$VCV)

Without.Var.Animal.STCont <- WithoutSpeciesModel$VCV[,11] * mean(diag(vcv(Tree)))

Without.Var.Animal.STBin <- WithoutSpeciesModel$VCV[,16] * mean(diag(vcv(Tree)))

Without.Var.Animal.OD <- WithoutSpeciesModel$VCV[,1] * mean(diag(vcv(Tree)))

Without.Var.Animal.qPCR <- WithoutSpeciesModel$VCV[,6] * mean(diag(vcv(Tree)))

# 4.1.2 Collecting residuals ---------------------------------------------------

dimnames(WithoutSpeciesModel$VCV)

Without.Var.units.STCont <- WithoutSpeciesModel$VCV[,19]

Without.Var.units.STBin <- WithoutSpeciesModel$VCV[,20]

Without.Var.units.OD <- WithoutSpeciesModel$VCV[,17]

Without.Var.units.qPCR <- WithoutSpeciesModel$VCV[,18]

# 4.1.3 Calculating repeatabilities --------------------------------------------

"ST Continuous"

repeatability.STCont <- Without.Var.Animal.STCont / 
  (Without.Var.Animal.STCont  + Without.Var.units.STCont)

mean(repeatability.STCont)
HPDinterval(as.mcmc(repeatability.STCont))

"ST Binary"

repeatability.STBin <- Without.Var.Animal.STBin / 
  (Without.Var.Animal.STBin  + Without.Var.units.STBin)

mean(repeatability.STBin)
HPDinterval(as.mcmc(repeatability.STBin))

"Optical density"

repeatability.OD <- Without.Var.Animal.OD / 
  (Without.Var.Animal.OD  + Without.Var.units.OD)

mean(repeatability.OD)
HPDinterval(as.mcmc(repeatability.OD))

"qPCR"

repeatability.qPCR <- Without.Var.Animal.qPCR / 
  (Without.Var.Animal.qPCR  + Without.Var.units.qPCR)

mean(repeatability.qPCR)
HPDinterval(as.mcmc(repeatability.qPCR))

# 4.2 Heritability (with species effects) --------------------------------------
# 4.2.1. Collect animal variances ----------------------------------------------

dimnames(WithSpeciesModel$VCV)

With.Var.Animal.STCont <- WithSpeciesModel$VCV[,11] * mean(diag(vcv(Tree)))

With.Var.Animal.STBin <- WithSpeciesModel$VCV[,16] * mean(diag(vcv(Tree)))

With.Var.Animal.OD <- WithSpeciesModel$VCV[,1] * mean(diag(vcv(Tree)))

With.Var.Animal.qPCR <- WithSpeciesModel$VCV[,6] * mean(diag(vcv(Tree)))

# 4.2.2. Collecting Species Effects --------------------------------------------

dimnames(WithSpeciesModel$VCV)

With.Var.Spc.STCont <- WithSpeciesModel$VCV[,27]

With.Var.Spc.STBin <- WithSpeciesModel$VCV[,32]

With.Var.Spc.OD <- WithSpeciesModel$VCV[,17]

With.Var.Spc.qPCR <- WithSpeciesModel$VCV[,22]

# 4.2.3. Calculating heritabilities --------------------------------------------

"ST continuous"

heritability.STCont <- With.Var.Animal.STCont /
  (With.Var.Animal.STCont + With.Var.Spc.STCont)

mean(heritability.STCont)
HPDinterval(as.mcmc(heritability.STCont))

"ST binary"

heritability.STBin <- With.Var.Animal.STBin /
  (With.Var.Animal.STBin + With.Var.Spc.STBin)

mean(heritability.STBin)
HPDinterval(as.mcmc(heritability.STBin))

"Optical density"

heritability.OD <- With.Var.Animal.OD /
  (With.Var.Animal.OD + With.Var.Spc.OD)

mean(heritability.OD)
HPDinterval(as.mcmc(heritability.OD))

"qPCR"

heritability.qPCR <- With.Var.Animal.qPCR /
  (With.Var.Animal.qPCR + With.Var.Spc.qPCR)

mean(heritability.qPCR)
HPDinterval(as.mcmc(heritability.qPCR))

# 4.3 Heritability as a product of total variance ------------------------------
# 4.3.1. Collect animal variances ----------------------------------------------

dimnames(WithSpeciesModel$VCV)

With.Var.Animal.STCont <- WithSpeciesModel$VCV[,11] * mean(diag(vcv(Tree)))

With.Var.Animal.STBin <- WithSpeciesModel$VCV[,16] * mean(diag(vcv(Tree)))

With.Var.Animal.OD <- WithSpeciesModel$VCV[,1] * mean(diag(vcv(Tree)))

With.Var.Animal.qPCR <- WithSpeciesModel$VCV[,6] * mean(diag(vcv(Tree)))

# 4.3.2. Collecting Species Effects --------------------------------------------

dimnames(WithSpeciesModel$VCV)

With.Var.Spc.STCont <- WithSpeciesModel$VCV[,27]

With.Var.Spc.STBin <- WithSpeciesModel$VCV[,32]

With.Var.Spc.OD <- WithSpeciesModel$VCV[,17]

With.Var.Spc.qPCR <- WithSpeciesModel$VCV[,22]

# 4.3.3. Collecting the within-strain variance ---------------------------------

With.Var.Units.STCont <- WithSpeciesModel$VCV[,35]

With.Var.Units.STBin <- WithSpeciesModel$VCV[,36]

With.Var.Units.OD <- WithSpeciesModel$VCV[,33]

With.Var.Units.qPCR <- WithSpeciesModel$VCV[,34]

# 4.3.3. Calculating heritability as a product of total variation --------------

"ST continuous"

TotalVarHeritability.STCont <- With.Var.Animal.STCont /
  (With.Var.Animal.STCont + With.Var.Spc.STCont + With.Var.Units.STCont)

mean(TotalVarHeritability.STCont)
HPDinterval(as.mcmc(TotalVarHeritability.STCont))

"ST binary"

TotalVarHeritability.STBin <- With.Var.Animal.STBin /
  (With.Var.Animal.STBin + With.Var.Spc.STBin + With.Var.Units.STBin)

mean(TotalVarHeritability.STBin)
HPDinterval(as.mcmc(TotalVarHeritability.STBin))

"Optical density"

TotalVarHeritability.OD <- With.Var.Animal.OD /
  (With.Var.Animal.OD + With.Var.Spc.OD + With.Var.Units.OD)

mean(TotalVarHeritability.OD)
HPDinterval(as.mcmc(TotalVarHeritability.OD))

"qPCR"

TotalVarHeritability.qPCR <- With.Var.Animal.qPCR /
  (With.Var.Animal.qPCR + With.Var.Spc.qPCR + With.Var.Units.qPCR)

mean(TotalVarHeritability.qPCR)
HPDinterval(as.mcmc(TotalVarHeritability.qPCR))

# ------------------------------------------------------------------------------
# ------------------------- [5] Correlations -----------------------------------
# ------------------------------------------------------------------------------

# 5.1. Calculating the correlation coefficient ---------------------------------

"The correlation coefficient is calculated as follows: 

                        R = cov.x.y / sqrt(var.x + var.y)"

dimnames(WithoutSpeciesModel$VCV)

"STcon:OD" 
mean(WithoutSpeciesModel$VCV[,3] / sqrt(WithoutSpeciesModel$VCV[,11] * WithoutSpeciesModel$VCV[,1]))

HPDinterval(WithoutSpeciesModel$VCV[,3] / sqrt(WithoutSpeciesModel$VCV[,11] * WithoutSpeciesModel$VCV[,1]))

"OD:STcon" 
mean(WithoutSpeciesModel$VCV[,9] / sqrt(WithoutSpeciesModel$VCV[,11] * WithoutSpeciesModel$VCV[,1]))

HPDinterval(WithoutSpeciesModel$VCV[,9] / sqrt(WithoutSpeciesModel$VCV[,11] * WithoutSpeciesModel$VCV[,1]))

"STcon:qPCR" 
mean(WithoutSpeciesModel$VCV[,7] / sqrt(WithoutSpeciesModel$VCV[,11] * WithoutSpeciesModel$VCV[,6]))

HPDinterval(WithoutSpeciesModel$VCV[,7] / sqrt(WithoutSpeciesModel$VCV[,11] * WithoutSpeciesModel$VCV[,6]))

"qPCR:STcon" 
mean(WithoutSpeciesModel$VCV[,10] / sqrt(WithoutSpeciesModel$VCV[,11] * WithoutSpeciesModel$VCV[,6]))

HPDinterval(WithoutSpeciesModel$VCV[,10] / sqrt(WithoutSpeciesModel$VCV[,11] * WithoutSpeciesModel$VCV[,6]))

"OD:qPCR" 
mean(WithoutSpeciesModel$VCV[,5] / sqrt(WithoutSpeciesModel$VCV[,1] * WithoutSpeciesModel$VCV[,6]))

HPDinterval(WithoutSpeciesModel$VCV[,5] / sqrt(WithoutSpeciesModel$VCV[,1] * WithoutSpeciesModel$VCV[,6]))

"qPCR:OD" 
mean(WithoutSpeciesModel$VCV[,2] / sqrt(WithoutSpeciesModel$VCV[,1] * WithoutSpeciesModel$VCV[,6]))

HPDinterval(WithoutSpeciesModel$VCV[,2] / sqrt(WithoutSpeciesModel$VCV[,1] * WithoutSpeciesModel$VCV[,6]))

"binary ST:OD"
mean(WithoutSpeciesModel$VCV[,4] / sqrt(WithoutSpeciesModel$VCV[,16] * WithoutSpeciesModel$VCV[,1]))

HPDinterval(WithoutSpeciesModel$VCV[,4] / sqrt(WithoutSpeciesModel$VCV[,16] * WithoutSpeciesModel$VCV[,1]))

"OD:binary ST"
mean(WithoutSpeciesModel$VCV[,13] / sqrt(WithoutSpeciesModel$VCV[,16] * WithoutSpeciesModel$VCV[,1]))

HPDinterval(WithoutSpeciesModel$VCV[,13] / sqrt(WithoutSpeciesModel$VCV[,16] * WithoutSpeciesModel$VCV[,1]))

"binary ST:qPCR"
mean(WithoutSpeciesModel$VCV[,8] / sqrt(WithoutSpeciesModel$VCV[,16] * WithoutSpeciesModel$VCV[,6]))

HPDinterval(WithoutSpeciesModel$VCV[,8] / sqrt(WithoutSpeciesModel$VCV[,16] * WithoutSpeciesModel$VCV[,6]))

"qPCR:binary ST"
mean(WithoutSpeciesModel$VCV[,14] / sqrt(WithoutSpeciesModel$VCV[,16] * WithoutSpeciesModel$VCV[,6]))

HPDinterval(WithoutSpeciesModel$VCV[,14] / sqrt(WithoutSpeciesModel$VCV[,16] * WithoutSpeciesModel$VCV[,6]))

# 5.2. Calculating the slope ---------------------------------------------------

"The slope of the correlation is calculated as follows: 

                            slope = cov.x.y / var.x "

dimnames(WithoutSpeciesModel$VCV)

"STcon:OD" 
mean(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,3] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,11])

HPDinterval(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,3] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,11])

"OD:STcon" 
mean(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,9] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,1])

HPDinterval(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,9] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,1])

"STcon:qPCR" 
mean(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,7] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,11])

HPDinterval(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,7] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,11])

"qPCR:STcon" 
mean(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,10] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,6])

HPDinterval(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,10] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,6])

"OD:qPCR" 
mean(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,5] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,1])

HPDinterval(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,5] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,1])

"qPCR:OD" 
mean(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,2] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,6])

HPDinterval(V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,2] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,6])

"binary ST:OD"
mean(pnorm((V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,4] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,16]), 0, sqrt(1)))

HPDinterval(pnorm((V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,4] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,16]), 0, sqrt(1)))

"OD:binary ST"
mean(pnorm((V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,13] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,1]), 0, sqrt(1)))

HPDinterval(pnorm((V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,13] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,1]), 0, sqrt(1)))

"binary ST:qPCR"
mean(pnorm((V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,8] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,16]), 0, sqrt(1)))

HPDinterval(pnorm((V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,8] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,16]), 0, sqrt(1)))

"qPCR:binary ST"
mean(pnorm((V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,14] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,6]), 0, sqrt(1)))

HPDinterval(pnorm((V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,14] / V003_MCMCglmm_WholePhylogeny_WithoutSpecies$VCV[,6]), 0, sqrt(1)))

# ------------------------------------------------------------------------------
# ------------------------ [6] Alternative Models ------------------------------
# ------------------------------------------------------------------------------

"The following section shows you how to run variation of the model that were
 mentioned in the manuscript."


# 6.1. Data for the OD assay outlier analysis ----------------------------------

"If you want to run the outlier analysis instead of the main model, load in 
 either one of the following dataframes and continue the analysis from 1.2. The  
 remainder of code can be run as presented but with model, tree, and dataframe 
 names changed to match the model you are trying to run."

Data <- read.csv("PhageSusceptibilityDataNoOutlierRemoval.csv", header = TRUE) 

Data <- read.csv("PhageSusceptibilityDataMajorOutliersRemoved.csv", header = TRUE) 

# 6.2. Data for the drop tip among-species and within-aureus trees -------------

"If you want to run the analysis for the among-species and within-aureus trees, 
 run these additional lines of code after 1.3 and before 1.4. The remainder of 
 code can be run as presented but with model, tree, and dataframe names changed
 to match the model you are trying to run."

Key <- read_csv("WholePhylogenyToSmallPhylogeniesKey.csv")

colnames(Key)[1] <- "Sample"

Data <- left_join(Data, Key, by = "Sample")

Data.Aureus <- Data[!is.na(Data$Aureus), ]
Data.Aureus <- Data.Aureus[-c(10, 12)]

Data.Species <- Data[!is.na(Data$Species), ]
Data.Species <- Data.Species[-c(11, 12)]

Tree.Aureus <- drop.tip(Tree, sort(Tree$tip.label[
  as.character(unique(Tree$tip.label)) %in% Data.Aureus$animal == FALSE]))

plot(ladderize(Tree.Aureus), cex = 0.75)

Tree.Species <- drop.tip(Tree, sort(Tree$tip.label[
  as.character(unique(Tree$tip.label)) %in% Data.Species$animal == FALSE]))

plot(ladderize(Tree.Species), cex = 0.75)

# 6.3. Compare the drop-tip phylogenies to the individually constructed ---------
#      phylogenies -------------------------------------------------------------

"6.2 shows you how to drop tips from the 64 strain phylogeny to create the 
 individual trees used to run the within-aureus and among-species models. The 
 topology of these trees was compared with the topology of the within-aureus 
 and among-species trees constructed originally. To load the original trees and
 compare them to the drop-tip trees, run the following code."

Tree.Species <- read.tree("Staph_Phylogeny_Species.nwk")

Tree.aureus <- read.tree("Staph_Phylongey_Aureus.nwk")

# 6.4. Alternative priors ------------------------------------------------------

"If you want to run the analysis with the alternative (Inverse Wishart) priors, 
 load the following and re-run the analysis."

Prior.With.Species <- list(G = list( 
  G1 = list(V = diag(4), n = 4.004),
  G2 = list(V = diag(4), n = 4.004)),
  R = list(V = diag(4), nu = 0.004, fix = 4))

Prior.Without.Species <- list(G = list( 
  G1 = list(V = diag(4), n = 4.004)),
  R = list(V = diag(4), nu = 0.004, fix = 4))

# 6.5 Caluclating distance effects ---------------------------------------------
# 6.5.1. Distance from amplification host --------------------------------------

"Gets the distances between tips in the tree in a matrix."

Distances <- cophenetic.phylo(Tree)

"Check which number the amplification host (13S44S9) is, for me it was the 
 38th column."

colnames(Distances) 

"Collect the distances from the amplification host"

DistanceFromAmp <- data.frame(Distance = Distances[,38])
DistanceFromAmp$animal <- rownames(DistanceFromAmp)
rownames(DistanceFromAmp) <- c(1:64)

"Calculate relative distance from the amplification host"

DistanceFromAmp$RelativeDistance = ((DistanceFromAmp$Distance/max(DistanceFromAmp$Distance))*2)

"Rejoin the data"

Data <- left_join(Data, DistanceFromAmp, by = "animal")

# 6.5.2. Model structure for distance effects ----------------------------------

WithSpeciesModel_Distance <- 
  MCMCglmm(ISP ~ trait * RelativeDistance, random = ~us(trait) : animal + us(trait) : Sample_Full,
           rcov = ~idh(trait) : units, ginv = list(animal=Ainv), prior = Prior.With.Species, data = Data,
           family = NULL, nitt = 13000000, thin = 5000, burnin = 3000000, pr = TRUE)







