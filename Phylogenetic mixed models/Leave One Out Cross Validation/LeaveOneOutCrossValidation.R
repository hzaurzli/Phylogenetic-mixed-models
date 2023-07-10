#===============================================================================
#-------------------------------------------------------------------------------
#       Leave-one-out cross validation of phylogenetic GLMMs
#-------------------------------------------------------------------------------
#===============================================================================

# 0.1 Purpose-------------------------------------------------------------------

'This script validates the findings of phylogenetic GLMMs in 002 Susceptibility,
which all appear to have heritabilities and repeatabilities of 1 (1, 1), by
asking models to predict the trait values of each Staph strain based on its
evolutionary relationships and the trait values of its relatives.

These values are plotted and a root mean square error calculated for qPCR, OD,
and both components of ST.

The majority of the explanation for the phylogenetic models and data structures
used below can be found in the file Phage Susceptibility Analysis.R'

# ------------------------------------------------------------------------------
# -------------------------- [1] Initialisation --------------------------------
# ------------------------------------------------------------------------------
# 1.1 Install packages ---------------------------------------------------------

library(tidyverse); library(MCMCglmm); library(ape); library(ggtree); 
library(geiger); library(lme4); library(nortest); library(lmtest); library(gt);
library(tidyr); library(phytools); library(plotrix)

# 1.2 Set directory ------------------------------------------------------------

"Chose your working directory"

setwd("")

# 1.3 Load data ----------------------------------------------------------------

Data <- read.csv("PhageSusceptibilityData.csv", header = TRUE)

colnames(Data)[4] <- "method"

Data$family <- NULL

# 1.4 Load host phylogeny ------------------------------------------------------

Tree <- read.tree("Staph_Phylogeny.nwk")

ggtree(Tree, size = 0.25, right = T)+ 
  coord_cartesian(xlim = c(0, 0.6)) + 
  geom_treescale() +
  geom_tiplab(size = 2.5)


# 1.5 Tree lengths -------------------------------------------------------------

mean(diag(vcv(Tree))) # 0.477006


# ------------------------------------------------------------------------------
# ---------------------- [2] Drop One Models  ----------------------------------
# ------------------------------------------------------------------------------

"If you have previously run the models, skip to 3.1 Load models"

# 2.1 Define priors ------------------------------------------------------------

Prior.qPCR <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))

Prior.OD <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))

Prior.STcon <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))

Prior.STbin <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002, fix = 1))


# 2.2 Run models ---------------------------------------------------------------

Ainv<-inverseA(Tree, scale=FALSE)$Ainv

strains <- unique(Data$animal)


# ----- 2.2.1 qPCR -----

"The code block below does the following:
1. Create an empty dataframe to store predicted values for each strain.
2. Create a dataframe containing only qPCR data from experiment.
3. For each strain in the qPCR data:
    a) Filter out the qPCR data for this strain.
    b) Run a phylogenetic mixed model with this strain present in the phylogeny
       but absent from the qPCR data. This causes the model to predict the
       missing value as best it can from the evolutionary relationships of the
       missing strain and the qPCR values of its relatives.
    c) Collect the mean predicted qPCR value for the missing strain and store it
       in the dataframe from (1.)
4. Calculate the root-mean-squared-error of the predicted values compared to the
   measured values for each strain.
5. Plot predicted vs measured values for each strain.
"

data.drop.one.qPCR <- data.frame(strains = "",
                            ISP = 0)

data.model.qPCR <- filter(Data, method == "qPCR")


for (i in c(1:length(strains))){
  
  current.data <- filter(data.model.qPCR, animal != strains[i])
  
  current.model <- MCMCglmm(ISP ~ 1, random = ~animal,
                            rcov = ~units, ginv = list(animal=Ainv), prior = Prior.qPCR, data = current.data,
                            family = "gaussian", nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)
  
  x <- 0
  found <- F
  
  while(found == F){
    x <- x + 1
    if (grepl(paste("animal.", strains[i], sep = ""), dimnames(current.model$Sol)[[2]][x])){
      found <- T
    }
  }
  
  data.drop.one.qPCR[i,] <- c(strains[i],
                         mean(current.model$Sol[,x] + current.model$Sol[,1]))
  
}

data.qPCR <- data.model.qPCR %>% group_by(animal) %>% summarise(ISP.actual = mean(ISP, na.rm = T))

colnames(data.drop.one.qPCR) <- c("animal", "ISP.predicted")

data.qPCR <- left_join(data.qPCR, data.drop.one.qPCR)

data.qPCR$ISP.predicted <- as.numeric(data.qPCR$ISP.predicted)


RMSE.qPCR <- data.qPCR$ISP.predicted -
  data.qPCR$ISP.actual

RMSE.qPCR <- RMSE.qPCR^2

RMSE.qPCR <- mean(RMSE.qPCR)

RMSE.qPCR <- sqrt(RMSE.qPCR)

plot.qPCR <- ggplot(data.qPCR) +
  geom_abline(color = "#abadb2") +
  geom_vline(xintercept = 4.47, color = "#fc9b63") +
  geom_abline(intercept = RMSE.qPCR, linetype = "dotted", color = "#abadb2") +
  geom_abline(intercept = -RMSE.qPCR, linetype = "dotted", color = "#abadb2") +
  geom_point(aes(x = ISP.predicted, y = ISP.actual), size = 2, alpha = 0.35, color = "#2E3440") +
  coord_fixed() +
  geom_text(aes(x = ISP.predicted, y = ISP.actual, label = animal)) +
  scale_y_continuous(limits = c(0,6), expand = c(0,0),
                     name = "Measured change in viral load",
                     breaks = c(0,2,4,6),
                     labels = c(expression(10^0),
                                expression(10^2),
                                expression(10^4),
                                expression(10^6))) +
  scale_x_continuous(limits = c(0,6), expand = c(0,0),
                     name = "Predicted change in viral load",
                     breaks = c(0,2,4,6),
                     labels = c(expression(10^0),
                                expression(10^2),
                                expression(10^4),
                                expression(10^6))) +
  theme_bw() +
  theme(text = element_text(size = 8.5, color = "#2E3440"),
        panel.border = element_blank(),
        axis.line = element_line(color = "#2E3440"),
        panel.grid.minor = element_line(size = 0.5))

plot.qPCR

# ----- 2.2.2 OD -----

"For a description of the code below see section 2.2.1 above."

data.drop.one.OD <- data.frame(strains = "",
                               ISP = 0)

data.model.OD <- filter(Data, method == "OD")


for (i in c(1:length(strains))){
  
  current.data <- filter(data.model.OD, animal != strains[i])
  
  current.model <- MCMCglmm(ISP ~ 1, random = ~animal,
                            rcov = ~units, ginv = list(animal=Ainv), prior = Prior.OD, data = current.data,
                            family = "gaussian", nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)
  
  x <- 0
  found <- F
  
  while(found == F){
    x <- x + 1
    if (grepl(paste("animal.", strains[i], sep = ""), dimnames(current.model$Sol)[[2]][x])){
      found <- T
    }
  }
  
  data.drop.one.OD[i,] <- c(strains[i],
                            mean(current.model$Sol[,x] + current.model$Sol[,1]))
  
}

data.OD <- data.model.OD %>% group_by(animal) %>% summarise(ISP.actual = mean(ISP, na.rm = T))

colnames(data.drop.one.OD) <- c("animal", "ISP.predicted")

data.OD <- left_join(data.OD, data.drop.one.OD)

data.OD$ISP.predicted <- as.numeric(data.OD$ISP.predicted)


RMSE.OD <- data.OD$ISP.predicted -
  data.OD$ISP.actual

RMSE.OD <- RMSE.OD^2

RMSE.OD <- mean(RMSE.OD)

RMSE.OD <- sqrt(RMSE.OD)

plot.OD <- ggplot(data.OD) +
  geom_abline(color = "#abadb2") +
  geom_vline(xintercept = 0.675, color = "#fc9b63") +
  geom_abline(intercept = RMSE.OD, linetype = "dotted", color = "#abadb2") +
  geom_abline(intercept = -RMSE.OD, linetype = "dotted", color = "#abadb2") +
  geom_point(aes(x = ISP.predicted, y = ISP.actual), size = 2, alpha = 0.35, color = "#2E3440") +
  coord_fixed() +
  geom_text(aes(x = ISP.predicted, y = ISP.actual, label = animal)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0),
                     name = "Measured change in OD",
                     breaks = c(0,0.25,0.5,0.75, 1)) +
  scale_x_continuous(limits = c(0,1), expand = c(0,0),
                     name = "Predicted change in OD",
                     breaks = c(0,0.25,0.5,0.75, 1)) +
  theme_bw() +
  theme(text = element_text(size = 8.5, color = "#2E3440"),
        panel.border = element_blank(),
        axis.line = element_line(color = "#2E3440"),
        panel.grid.minor = element_line(size = 0.5))

plot.OD

# ----- 2.2.3 ST continuous -----

"For a description of the code below see section 2.2.1 above."

data.drop.one.STcon <- data.frame(strains = "",
                                  ISP = 0)

data.model.STcon <- filter(Data, method == "ST_Cont")

strains.STcon <- unique(na.omit(data.model.STcon)$animal)

data.model.STcon <- filter(data.model.STcon, animal %in% strains.STcon)

Tree.STcon <- drop.tip(Tree, sort(Tree$tip.label[
  as.character(unique(Tree$tip.label)) %in% strains.STcon == FALSE]))


Ainv.STcon<-inverseA(Tree.STcon, scale=FALSE)$Ainv


for (i in c(1:length(strains.STcon))){
  
  current.data <- filter(data.model.STcon, animal != strains.STcon[i])
  
  current.model <- MCMCglmm(ISP ~ 1, random = ~animal,
                            rcov = ~units, ginv = list(animal=Ainv.STcon), prior = Prior.STcon, data = current.data,
                            family = "gaussian", nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)
  
  x <- 0
  found <- F
  
  while(found == F){
    x <- x + 1
    if (grepl(paste("animal.", strains.STcon[i], sep = ""), dimnames(current.model$Sol)[[2]][x])){
      found <- T
    }
  }
  
  data.drop.one.STcon[i,] <- c(strains.STcon[i],
                               mean(current.model$Sol[,x] + current.model$Sol[,1]))
  
}

data.STcon <- data.model.STcon %>% group_by(animal) %>% summarise(ISP.actual = mean(ISP, na.rm = T))

colnames(data.drop.one.STcon) <- c("animal", "ISP.predicted")

data.STcon <- left_join(data.STcon, data.drop.one.STcon)

data.STcon$ISP.predicted <- as.numeric(data.STcon$ISP.predicted)


RMSE.STcon <- data.STcon$ISP.predicted -
  data.STcon$ISP.actual

RMSE.STcon <- RMSE.STcon^2

RMSE.STcon <- mean(RMSE.STcon)

RMSE.STcon <- sqrt(RMSE.STcon)

plot.STcon <- ggplot(data.STcon) +
  geom_abline(color = "#abadb2") +
  geom_vline(xintercept = 3.92e+05, color = "#fc9b63") +
  geom_point(aes(x = ISP.predicted, y = ISP.actual), size = 2, alpha = 0.35, color = "#2E3440") +
  coord_fixed() +
  scale_y_continuous(limits = c(0,900000), expand = c(0,0),
                     name = "Measured PFU/uL",
                     breaks = c(0, 300000, 600000, 900000),
                     labels = c(expression(0),
                                expression(3*"x"*10^5),
                                expression(6*"x"*10^5),
                                expression(9*"x"*10^5))) +
  scale_x_continuous(limits = c(0,900000), expand = c(0,0),
                     name = "Predicted PFU/uL",
                     breaks = c(0, 300000, 600000, 900000),
                     labels = c(expression(0),
                                expression(3*"x"*10^5),
                                expression(6*"x"*10^5),
                                expression(9*"x"*10^5))) +
  theme_bw() +
  theme(text = element_text(size = 8.5, color = "#2E3440"),
        panel.border = element_blank(),
        axis.line = element_line(color = "#2E3440"),
        panel.grid.minor = element_line(size = 0.5))

plot.STcon

# ----- 2.2.4 ST binary -----

"For a description of the code below see section 2.2.1 above. To get binary
predictions (0/1) from these models, predicted values have been rounded to the
nearest 0 or 1 and then classified depending on whether they were predicted to
be permissive (1) or non-permissive (0) and whether the prediction was correct."

data.drop.one.STbin <- data.frame(strains = "",
                               ISP = 0)

data.model.STbin <- filter(Data, method == "ST_Bin")


for (i in c(1:length(strains))){
  
  current.data <- filter(data.model.STbin, animal != strains[i])
  
  current.model <- MCMCglmm(ISP ~ 1, random = ~animal,
                            rcov = ~units, ginv = list(animal=Ainv), prior = Prior.STbin, data = current.data,
                            family = "threshold", nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)
  
  x <- 0
  found <- F
  
  while(found == F){
    x <- x + 1
    if (grepl(paste("animal.", strains[i], sep = ""), dimnames(current.model$Sol)[[2]][x])){
      found <- T
    }
  }
  
  data.drop.one.STbin[i,] <- c(strains[i],
                            mean(current.model$Sol[,x] + current.model$Sol[,1]))
  
}

data.drop.one.STbin$ISP <- pnorm(as.numeric(data.drop.one.STbin$ISP), 0, sqrt(1))


data.STbin <- data.model.STbin %>% group_by(animal) %>% summarise(ISP.actual = median(ISP, na.rm = T))

colnames(data.drop.one.STbin) <- c("animal", "ISP.predicted")

data.STbin <- left_join(data.STbin, data.drop.one.STbin)

data.STbin$ISP.predicted <- ifelse(data.STbin$ISP.predicted - 0.5 > 0, 1, 0)

data.STbin$result <- ifelse(data.STbin$ISP.actual == 0 & data.STbin$ISP.predicted == 0, "Non-permissive, prediction correct",
                             ifelse(data.STbin$ISP.actual == 0 & data.STbin$ISP.predicted == 1, "Non-permissive, prediction incorrect",
                                    ifelse(data.STbin$ISP.actual == 1 & data.STbin$ISP.predicted == 1, "Permissive, prediction correct", "Permissive, prediction incorrect")))

data.STbin$prediction <- data.STbin$ISP.actual == data.STbin$ISP.predicted

data.STbin.bars <- data.STbin %>% group_by(ISP.actual) %>%
  summarise(Success = sum(prediction),
            Failure = sum(prediction == F))

data.STbin.bars <- gather(data.STbin.bars, key = "Prediction", value = "count", 2:3)

data.STbin.bars$ISP.actual <- factor(data.STbin.bars$ISP.actual, levels = c("0", "1"))

plot.STbin <- ggplot(data.STbin.bars) +
  geom_bar(aes(x = ISP.actual, y = count, fill = Prediction), stat = "identity", position = "dodge") +
  theme_bw() +
  theme(text = element_text(size = 8.5, color = "#2E3440"),
        panel.border = element_rect(color = "#ebebeb"),
        axis.line = element_line(color = "#2E3440"),
        axis.line.y.right  = element_line(color = 'red'),
        panel.grid.minor = element_line(size = 0.5),
        legend.position = c(0.3, 0.775)) +
  scale_x_discrete(labels = c("Non-permissive", "Permissive"), name = "Measured plaque assay permissiveness") +
  scale_y_continuous(name = "Number of strains", expand = c(0,0), limits = c(0,35)) +
  scale_fill_manual(name = "Prediction outcome:", values = c("#353b47", "#b6b8bc"))

plot.STbin


# ------------------------------------------------------------------------------
# ---------------------- [3] Null Models  --------------------------------------
# ------------------------------------------------------------------------------

# 3.1 Priors -------------------------------------------------------------------

Prior.qPCR.null <- list(R = list(V = diag(1), nu = 0.002))

Prior.OD.null <- list(R = list(V = diag(1), nu = 0.002))

Prior.STcon.null <- list(R = list(V = diag(1), nu = 0.002))

Prior.STbin.null <- list(R = list(V = diag(1), nu = 0.002, fix = 1))


# 3.2 Run Models ---------------------------------------------------------------

# ----- 3.2.1 qPCR -----

model.qPCR.null <- MCMCglmm(ISP ~ 1,
                            rcov = ~units, prior = Prior.qPCR.null, data = data.model.qPCR,
                            family = "gaussian", nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)

"Model will have an intercept value."

summary(model.qPCR.null)

"Post mean is 4.470"

"Create dataframe with measured, predicted, and null predicted values for each
 strain:"

data.qPCR.comparison <- data.model.qPCR %>% group_by(Sample) %>% summarise(ISP.actual = mean(ISP))

data.qPCR.comparison$ISP.null <- 4.470

key <- unique(dplyr::select(data.model.qPCR, Sample, animal))

data.qPCR.comparison <- left_join(data.qPCR.comparison, key)

data.qPCR.comparison <- left_join(data.qPCR.comparison, data.drop.one.qPCR)


# ----- 3.2.2 OD -----

model.OD.null <- MCMCglmm(ISP ~ 1,
                          rcov = ~units, prior = Prior.OD.null, data = data.model.OD,
                          family = "gaussian", nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)

"Model will have an intercept value."

summary(model.OD.null)

"Post mean is 0.6744"

"Create dataframe with measured, predicted, and null predicted values for each
 strain:"

data.OD.comparison <- data.model.OD %>% group_by(Sample) %>% summarise(ISP.actual = mean(ISP))

data.OD.comparison$ISP.null <- 0.6744

data.OD.comparison <- left_join(data.OD.comparison, key)

data.OD.comparison <- left_join(data.OD.comparison, data.drop.one.OD)


# ----- 3.2.3 STcon -----

model.STcon.null <- MCMCglmm(ISP ~ 1,
                             rcov = ~units, prior = Prior.STcon.null, data = data.model.STcon,
                             family = "gaussian", nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)

"Model will have an intercept value."

summary(model.STcon.null)

"Post mean is 392165"

"Create dataframe with measured, predicted, and null predicted values for each
 strain:"

data.STcon.comparison <- data.model.STcon %>% group_by(Sample) %>% summarise(ISP.actual = mean(ISP))

data.STcon.comparison$ISP.null <- 392165

data.STcon.comparison <- left_join(data.STcon.comparison, key)

data.STcon.comparison <- left_join(data.STcon.comparison, data.drop.one.STcon)


# ----- 3.2.4 STbin -----

model.STbin.null <- MCMCglmm(ISP ~ 1,
                             rcov = ~units, prior = Prior.STbin.null, data = data.model.STbin,
                             family = "gaussian", nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)

"Model will have an intercept value."

summary(model.STbin.null)

"Post mean is 0.6094"

pnorm(0.6094, 0, sqrt(1))

"Converted to probability is 0.73"

"=> most probable prediction from null model is permissive (1)"

"Create dataframe with measured, predicted, and null predicted values for each
 strain:"

data.STbin.comparison <- data.model.STbin %>% group_by(Sample) %>% summarise(ISP.actual = median(ISP))

data.STbin.comparison$ISP.null <- 1

data.STbin.comparison <- left_join(data.STbin.comparison, key)

data.STbin.comparison <- left_join(data.STbin.comparison, data.drop.one.STbin)

data.STbin.comparison$ISP.predicted <- ifelse(data.STbin.comparison$ISP.predicted > 0.5, 1, 0)


# ------------------------------------------------------------------------------
# ---------------------- [4] Comparison of errors ------------------------------
# ------------------------------------------------------------------------------

# 4.1 qPCR ---------------------------------------------------------------------

data.qPCR.comparison$ISP.predicted <- as.numeric(data.qPCR.comparison$ISP.predicted)

data.qPCR.comparison$Errors.null <- abs(data.qPCR.comparison$ISP.null - data.qPCR.comparison$ISP.actual)

data.qPCR.comparison$Errors.predicted <- abs(data.qPCR.comparison$ISP.predicted - data.qPCR.comparison$ISP.actual)

ggplot(data.qPCR.comparison) +
  geom_boxplot(aes(x = "Null", y = Errors.null)) +
  geom_boxplot(aes(x = "Predicted", y = Errors.predicted)) +
  scale_y_continuous(name = "Error")

wilcox.test(data.qPCR.comparison$Errors.null,
            data.qPCR.comparison$Errors.predicted)


# 4.2 OD -----------------------------------------------------------------------

data.OD.comparison$ISP.predicted <- as.numeric(data.OD.comparison$ISP.predicted)

data.OD.comparison$Errors.null <- abs(data.OD.comparison$ISP.null - data.OD.comparison$ISP.actual)

data.OD.comparison$Errors.predicted <- abs(data.OD.comparison$ISP.predicted - data.OD.comparison$ISP.actual)

ggplot(data.OD.comparison) +
  geom_boxplot(aes(x = "Null", y = Errors.null)) +
  geom_boxplot(aes(x = "Predicted", y = Errors.predicted)) +
  scale_y_continuous(name = "Error")

wilcox.test(data.OD.comparison$Errors.null,
            data.OD.comparison$Errors.predicted)


# 4.3 STcon --------------------------------------------------------------------

data.STcon.comparison$ISP.predicted <- as.numeric(data.STcon.comparison$ISP.predicted)

data.STcon.comparison$Errors.null <- abs(data.STcon.comparison$ISP.null - data.STcon.comparison$ISP.actual)

data.STcon.comparison$Errors.predicted <- abs(data.STcon.comparison$ISP.predicted - data.STcon.comparison$ISP.actual)

ggplot(data.STcon.comparison) +
  geom_boxplot(aes(x = "Null", y = Errors.null)) +
  geom_boxplot(aes(x = "Predicted", y = Errors.predicted)) +
  scale_y_continuous(name = "Error")

wilcox.test(data.STcon.comparison$Errors.null,
            data.STcon.comparison$Errors.predicted)

