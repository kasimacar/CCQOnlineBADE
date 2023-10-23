# Title: COVID-19 conspiracy ideation is associated with delusion proneness and resistance to belief update
# Author: Kasim Acar, Karolinska Institutet, Department of Clinical Neuroscience

# Load libraries

library(openxlsx)
library(dplyr)
library(reghelper)
library(ggplot2)
library(ggstance)
library(radarchart)
library(coefplot)
library(effsize)
library(lm.beta)
library(jtools)
library(xlsx)
library(stats)
library(ggpubr)
library(data.table)
library(jmv)
library(lavaan)
library(semPlot)
library(broom)
library(tidyverse)
library(semTable)
library(broom.mixed)
library(DiagrammeR)
library(MASS)
library(rcompanion)
library(bestNormalize)
library(GGally)
library(psych)


# Load data
load("CCQDataGitHub_Final.RData")


# Find best-performing normalization method for variables of interest
(BNobjectCon <- bestNormalize(databadeall$CCQCon))
(BNobjectFeaAct <- bestNormalize(databadeall$CCQFeaAct))
(BNobjectDis <- bestNormalize(databadeall$CCQDis))
(BNobjectCCQTot <- bestNormalize(databadeall$CCQTot))
(BNobjectPDI <- bestNormalize(databadeall$PDI_total))
(BNobjectEII <- bestNormalize(databadeall$EII1score))
(BNobjectPDIAll <- bestNormalize(AllCovid$PDI))
(BNobjectCCQTotAll <- bestNormalize(AllCovid$CCQTot))
(BNobjectPDICMQ <- bestNormalize(AllCMQ$PDI_total))



# add normalized CCQ conspiracy variable to BADE data 
databadeall$CCQConNorm <- (BNobjectCon$chosen_transform$x.t)


# add normalized CCQ Fear/Acttion variable to BADE data 
databadeall$CCQFeaActNorm <- (BNobjectFeaAct$chosen_transform$x.t)


# add normalized CCQ Distrust variable to BADE data 
databadeall$CCQDisNorm <- (BNobjectDis$chosen_transform$x.t)


# add normalized CCQ Total variable to BADE data 
databadeall$CCQTotNorm <- (BNobjectCCQTot$chosen_transform$x.t)


# add normalized EII variable to BADE data 
databadeall$EIINorm <- (BNobjectEII$chosen_transform$x.t)


# add normalized PDI variable to BADE data 
databadeall$PDINorm <- (BNobjectPDI$chosen_transform$x.t)


# add normalized PDI variable All Covid  to data 
AllCovid$PDINorm <- (BNobjectPDIAll$chosen_transform$x.t)


# add normalized PDI variable to CMQ data 
AllCMQ$PDINorm <- (BNobjectPDICMQ$chosen_transform$x.t)


# add normalized PDI (All Covid) variable to data 
AllCovid$PDINorm <- (BNobjectPDIAll$chosen_transform$x.t)


# add normalized PDI (All Covid) variable to data 
AllCovid$CCQTotNorm <- (BNobjectCCQTotAll$chosen_transform$x.t)



# plot distributions
plotNormalHistogram(databadeall$CCQFeaActNorm)
plotNormalHistogram(databadeall$CCQDisNorm)
plotNormalHistogram(databadeall$CCQConNorm)
plotNormalHistogram(databadeall$CCQTotNorm)
plotNormalHistogram(databadeall$CCQTot)
plotNormalHistogram(databadeall$EIINorm)
plotNormalHistogram(databadeall$PDINorm)
plotNormalHistogram(AllCMQ$PDINorm)
plotNormalHistogram(AllCovid$PDINorm)
plotNormalHistogram(AllCovid$CCQTotNorm)

# Pearson's correlation CCQ Total and PDI
cor.test(AllCovid$PDINorm, AllCovid$CCQTot)

# Pearson's correlation CCQ total and CMQ
cor.test(AllCMQ$CONS5, AllCMQ$CCQTot)

# Pearson's correlation CCQ total and EII
cor.test(databadeall$EIINorm, databadeall$CCQTotNorm)


# Pearson's correlation PDI-Non-paranoia

cor.test(AllCovidNP$PDI, AllCovidNP$CCQTot)
cor.test(AllCovidNP$PDI, AllCovidNP$CCQFeaAct)
cor.test(AllCovidNP$PDI, AllCovidNP$CCQCon)
cor.test(AllCovidNP$PDI, AllCovidNP$CCQDis)

jmv::corrMatrix(
  data = AllCovid,
  vars = vars(CCQFeaAct, CCQDis, CCQCon, CCQTot, PDI),
  plotDens = TRUE,
  plotStats = TRUE)


cor.test(databadeall$PDI_total, databadeall$CCQTot)
cor.test(databadeall$PDI_total, databadeall$CCQFeaAct)
cor.test(databadeall$PDI_total, databadeall$CCQCon)
cor.test(databadeall$PDI_total, databadeall$CCQDis)

jmv::corrMatrix(
  data = databadeall,
  vars = vars(CCQFeaAct, CCQDis, CCQCon, CCQTot, PDI_total),
  plotDens = TRUE,
  plotStats = TRUE)


#create data frame for correlogram
myDataCorr <- data.frame(PDI = databadeall$PDINorm, EII = databadeall$EIINorm, CCQTot = databadeall$CCQTotNorm)


# Create correlation and scatter plot with densities for PDI, EII & CCQ Total
pairs.panels(myDataCorr,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             col = "black",
             pch = 16,           # pch symbol
             lm = TRUE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,          # Jittering factor
             smoother = FALSE,   # If TRUE,then smooth.scatter the data points
             hist.col = 4,       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = TRUE)         # If TRUE, adds confidence intervals




## Principal Component Analysis of CCQ ##

jmv::pca(
  data = SCRCVP,
  vars = vars(CV1, CV2, CV3, CV4, CV5, CV6, CV7, CV8, CV9, CV10, CV11, CV12, CV13, CV14, CV15, CV16, CV17, CV18, CV19, CV20, CV21, CV22),
  hideLoadings = 0.4,
  screePlot = TRUE,
  factorSummary = TRUE)

# Remove items CV16 due to significant cross-loadings and CV22 due to no significant loading

jmv::pca(
  data = SCRCVP,
  vars = vars(CV1, CV2, CV3, CV4, CV5, CV6, CV7, CV8, CV9, CV10, CV11, CV12, CV13, CV14, CV15, CV17, CV18, CV19, CV20, CV21),
  hideLoadings = 0.4,
  screePlot = TRUE,
  factorSummary = TRUE)


# Mutliple Regression: PDI + CCQ. N = 617
## PDI + control for demographics + diagnosis + CCQTot ##

ModelSCRCVPt <- lm(PDI ~ Age+Sex+Education+Diagnosis+CCQTot+ASRS+RAADSY+RAADSN+RAADSB,
                   data= AllCovid)

# Excluding variables of interest for residual extraction
ModelSCRCVPt <- lm(PDI ~ Age+Sex+Education+Diagnosis,
                   data= AllCovid)

# Extract residuals
ModelSCRCVPt.res <- residuals.lm(ModelSCRCVPt)


# plot residuals
plotNormalHistogram(ModelSCRCVPt.res)


# Summarize results
summary(ModelSUMt)
summary(ModelSCRCVPt)
summ(ModelSCRCVPt)
beta(ModelSCRCVPt)

# PDI + control for demographics + diagnosis
ModelSCRCVPp <- lm(PDI ~ Age+Sex+Education+Diagnosis+CCQFeaAct+CCQDis+CCQCon+ASRS+RAADSY+RAADSN+RAADSB,
                   data= AllCovid)

# excluding variables of interest
ModelSCRCVPp <- lm(PDI ~ Age+Sex+Education+Diagnosis,
                   data= AllCovid)

# Extract residuals
ModelSCRCVPp.res <- residuals.lm(ModelSCRCVPp)

# plot residuals
plotNormalHistogram(ModelSCRCVPp.res)


# Summarize results
summary(ModelSCRCVPp)
summ(ModelSCRCVPp)

# Multiple Regression: EII and CCQ. N = 313

# Create data frame with wanted variables 
myDataAll <- data.frame(Age = databadeall$age.x, Sex = databadeall$genderm1, CCQFeaAct = databadeall$CCQFeaAct, CCQDis = databadeall$CCQDis,
                        CCQCon = databadeall$CCQCon, PDI = databadeall$PDI_total, Schizotypy = databadeall$DP, OLIFE = databadeall$OLIFE_tot,
                        EII = databadeall$EII1score, EII2 = databadeall$EII2score, Education = databadeall$education, CCQTot = databadeall$CCQTot,
                        Diagnosis = databadeall$PsychDiagAny, ASRS = databadeall$ASRSLog, RAADSY = databadeall$raads_young, RAADSN = databadeall$raads_now, RAADSB = databadeall$raads_both)


# PDI and EII + control for demographics + diagnosis
ModelPDI <- lm(EII ~ Age+Sex+Education+Diagnosis+PDI+ASRS+RAADSY+RAADSN+RAADSB,
               data= myDataAll)


# Excluding variables of interest
ModelPDI <- lm(EII ~ Age+Sex+Education+Diagnosis,
               data= myDataAll)

# Extract residuals
ModelPDI.res <- residuals.lm(ModelPDI)

# plot residuals
plotNormalHistogram(ModelPDI.res)

# Summarize results
summary(ModelPDI)


# EII + Total CCQ
ModelFull1 <- lm(EII ~ Age+Sex+Education+Diagnosis+PDI+CCQTot+ASRS+RAADSY+RAADSN+RAADSB,
                 data= myDataAll)


# Extract residuals
ModelFull1.res <- residuals.lm(ModelFull1)

# plot residuals
plotNormalHistogram(ModelFull1.res)


# Summarize results
summary(ModelFull1)
beta(ModelFull1)

# Full model with Covid questions as separate components
ModelFull2 <- lm(EII ~ Age+Sex+Education+Diagnosis+PDI+CCQFeaAct+CCQDis+CCQCon+ASRS+RAADSY+RAADSN+RAADSB,
                 data= myDataAll)




# Summarize results
summary(ModelFull2)
#summ(ModelFull2)



# Create plot for all  models
plot_summs(ModelPDI, ModelFull1, ModelFull2, scale = TRUE, plot.distributions = FALSE, inner_ci_level = .9,
           model.names = c("PDI Only", "PDI + CCQ Total", "PDI + CCQ Components"))




# Path Models: PDI -> EII -> CCQ. N = 313

# Create data frame for CCQ Path Models
MeMo.Data <- data.frame(Age = databadeall$age.x, Sex = databadeall$sex, Education = databadeall$education,
                        CCQFeaActNorm = databadeall$CCQFeaActNorm, CCQDisNorm = databadeall$CCQDisNorm,
                        CCQConNorm = databadeall$CCQConNorm, PDINorm = databadeall$PDINorm,
                        EIINorm = databadeall$EIINorm, CCQTotNorm = databadeall$CCQTotNorm,
                        CCQFeaAct = databadeall$CCQFeaAct, CCQDis = databadeall$CCQDis,
                        CCQCon = databadeall$CCQCon, PDI = databadeall$PDI_total,
                        EII = databadeall$EII1score, CCQTot = databadeall$CCQTot)


# CCQTot path-model - non-normalized
model1 <- ' # direct effect
             CCQTot ~ c * PDI
           # mediator 1
             EII ~ a1*PDI
             CCQTot ~ b1*EII
           # indirect effect
             a1b1 := a1*b1
           # total effect
             total := c + (a1*b1)
             '


# Fit data to model
fit1 <- sem(model1, data = MeMo.Data, estimator = "WLS") # add estimator = "wls" for WLS estimates


# Extract estimates
summary(fit1, fit.measures=T, rsq=T, ci = TRUE, standardized = T)

# Create plot
semPaths(fit1, "model", "std", residuals = F, exoCov = F)


# Create plot with diagrammeR for non-normalized path-model
MeMoCCQTot <- grViz('
      digraph {
      
      #graph statements
      graph [rankdir = LR,
      nodesep = .5]
      
      #node statements
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica]
      "PDI";
      
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica,
      style = filled,
      fillcolor = "#50b8fa"]
      "EII";
      
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica,
      style = filled,
      fillcolor = white]
      "CCQ \\n Total";
      
      
      #edge statements
      edge [color = grey,
      fontname = "helvetica-bold",
      fontstyle = bold]
      "PDI"->"EII" [label ="β = 0.18***", fontsize = 10, color = "#50b8fa", penwidth = 2, fontcolor= "#50b8fa"];
      "EII"-> "CCQ \\n Total" [label ="   β = 0.15**", fontsize = 10, color = "#50b8fa", fontcolor = "#50b8fa"];
      "PDI"->"CCQ \\n Total" [label ="β = 0.32***", fontsize = 10, penwidth = 2, fontcolor = grey];
      }
      ')

# Show plot
MeMoCCQTot



# Non-transformed variables
# CCQ factors Path-model
model2 <- ' # direct effect
             CCQDis ~ c1 * PDI
             CCQFeaAct ~ c2 * PDI
             CCQCon ~ c3 * PDI
           # mediator
             EII ~ a1*PDI
             CCQDis ~ b1*EII
             CCQFeaAct ~ b2*EII
             CCQCon ~ b3*EII
           # indirect effect
             a1b1 := a1*b1
             a1b2 := a1*b2
             a1b3 := a1*b3
           # total effect
             total := c1 + c2 + c3 + (a1*b1) + (a1*b2) + (a1*b3)
             '

# Fit data to model 
fit2 <- sem(model2, data = MeMo.Data, estimator = "wls") # add estimator = "WLS" for wls estimates


# Extract estimates
summary(fit2, fit.measures=T, rsq=T, standardized = T)

# Create plot
semPaths(fit2, "model", "std", residuals = F, exoCov = F)





# Create plot with diagrammeR
MeMoCCQComp <- grViz('
      digraph {
      
      #graph statements
      graph [rankdir = LR,
      nodesep = .5]
      
      #node statements
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica]
      "PDI";
      
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica,
      style = filled,
      fillcolor = "#50b8fa"]
      "EII";
      
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica,
      style = filled,
      fillcolor = white]
      "CCQ \\n Conspiracy";
      
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica,
      style = filled,
      fillcolor = white]
      "CCQ \\n Fear/Action";
      
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica,
      style = filled,
      fillcolor = white]
      "CCQ \\n Distrust";
      
      #edge statements
      edge [color = grey,
      fontname = "helvetica-bold",
      fontstyle = bold]
      "PDI"->"EII" [label ="β = 0.18***", fontsize = 10, penwidth = 2 color = "#50b8fa", fontcolor= "#50b8fa"];
      "EII"-> "CCQ \\n Distrust" [label ="  β = 0.15**", fontsize = 10, color = "#50b8fa", fontcolor = "#50b8fa"];
      "EII"-> "CCQ \\n Conspiracy" [label ="β = 0.24***", fontsize = 10, penwidth = 2, color = "#50b8fa", fontcolor = "#50b8fa"];
      "EII"-> "CCQ \\n Fear/Action" [label ="    β = -0.11", fontsize = 10, color = "#50b8fa",  style = dotted, fontcolor = "#50b8fa"];
      "PDI"->"CCQ \\n Distrust" [label ="    β = 0.33***", fontsize = 10, penwidth = 2, fontcolor = grey];
      "PDI"->"CCQ \\n Conspiracy" [label ="β = 0.28***", fontsize = 10, penwidth = 2, fontcolor = grey];
      "PDI"->"CCQ \\n Fear/Action" [label ="β = -0.04", fontsize = 10, style = dotted, fontcolor = grey];
      }
      ')


# show plot
MeMoCCQComp


