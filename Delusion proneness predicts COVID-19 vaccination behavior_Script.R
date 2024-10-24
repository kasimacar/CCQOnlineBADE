## Title: Delusion proneness predicts COVID-19 vaccination behavior
## DOI: 10.3389/fpsyt.2024.1450429 
## Script Authors: Kasim Acar & Ariadni Karagiannidou

# Clear workspace
rm(list = ls())

# Load data

load("Delusion proneness predicts COVID-19 vaccination behavior.RData")


### subset only subject under the age of 36 --- conduct the same analyses after running the line below ###
#SCRCVP <- subset(SCRCVP, (SCRCVP$age.x<36))


##################
# load libraries #
##################

pacman::p_load(
  lmSupport, 
  lubridate, #time and date manipulations
  dplyr,
  lavaan, #SEM/Path model fitting
  tidytext, # text mining using dplyr, ggplot2 
  tm, # load library for the text data
  jtools,
  rstatix,
  ggplot2,
  gtsummary,
  interactions,
  DiagrammeR,
  wordcloud,
  wordcloud2,
  semPlot,
  reshape2,
  rms,
  ggrepel,
  flextable,
  here,
  stringr
)


# Raw Data overview
skimr::skim(SCRCVP)
head(SCRCVP)
names(SCRCVP)
str(SCRCVP) # check what classes columns are


##############
# Clean data #
##############


# Make variables numeric

SCRCVP$raads_now <- as.numeric(SCRCVP$raads_now)
SCRCVP$genderm1 <- as.numeric(SCRCVP$genderm1)


# Vaccine

SCRCVP$vaccination <- replace(SCRCVP$vaccinated, SCRCVP$vaccinated==2, 0)    # where 0 indicates the non vaccination and 1 indicates the vaccinated people
SCRCVP$ccq_sum <- varScore(SCRCVP, Forward = c('reason_vaccinate_1' , 'reason_vaccinate_2' , 'reason_vaccinate_3'))


# Summarize CCQ Components

SCRCVP$CCQDis <- varScore(SCRCVP, Forward =  c('CV3', 'CV8' , 'CV12' , 'CV13'))
SCRCVP$CCQCon <- varScore(SCRCVP , Forward = c( 'CV9' , 'CV10' , 'CV11' , 'CV14' , 'CV15' , 'CV17' ,'CV18' , 'CV19' , 'CV20' , 'CV21'))
SCRCVP$CCQFeaAct <- varScore(SCRCVP,  Forward = c( 'CV1' , 'CV2' , 'CV4' , 'CV5' , 'CV6' ) , Reverse = 'CV7', Range = c(1,11))
SCRCVP$CCQTot <- varScore(SCRCVP, Forward = c( 'CV1' , 'CV2' , 'CV3' , 'CV4' , 'CV5' , 'CV6' , 'CV8' , 'CV9' , 'CV10' , 'CV11' , 'CV12', 'CV13' , 'CV14' , 'CV15' , 'CV16' , 'CV17', 'CV18' , 'CV19','CV20' , 'CV21' , 'CV22'), Reverse = 'CV7' , Range = c(1,11)) 


# Summarize STAI - trait score

SCRCVP$STAI_Tscore <- varScore(SCRCVP,  Forward = c('STAI_T_2', 'STAI_T_4', 'STAI_T_5', 'STAI_T_8', 'STAI_T_9', 'STAI_T_11', 'STAI_T_12', 'STAI_T_15', 'STAI_T_17', 'STAI_T_18', 'STAI_T_20'),
                               Reverse = c('STAI_T_1', 'STAI_T_3', 'STAI_T_6', 'STAI_T_7', 'STAI_T_10', 'STAI_T_13', 'STAI_T_14', 'STAI_T_16', 'STAI_T_19'), Range=c(1,4)) -20


# Create Delta-T for Vaccination

SCRCVP$VaccinationStartDate <- as.Date("2020-12-27") # start date for vaccination in Sweden
SCRCVP$subjectVaccDate <- as.Date(SCRCVP$first_dose, origin = "1900-01-01")
format(SCRCVP$subjectVaccDate, "%m/%Y")
SCRCVP$Delta_time <- interval(SCRCVP$VaccinationStartDate, SCRCVP$subjectVaccDate) %/% months(1)


# change time to vaccinate for some subjects (who reported getting vaccinated before it was possible to get vaccinated)
SCRCVP["Delta_time"][SCRCVP["Delta_time"] < 0] <- NA


# group comparison: PDI and vaccination
wilcox.test(SCRCVP$PDI_total[SCRCVP$vaccination==1], SCRCVP$PDI_total[SCRCVP$vaccination==0])


#effect size
wilcox_effsize(SCRCVP, PDI_total ~ vaccination, paired = F)



### PDI and time to vaccinate correlation ###

#check correlation
cor(SCRCVP$Delta_time, SCRCVP$PDI_total, method = "spearman", use = "complete.obs")


# plot correlation
smoothScatter(SCRCVP$Delta_time ~ SCRCVP$PDI_total,
              pch = 16, col = "black",
              nrpoints = 50, bandwidth = 0.8,
              transformation = function(x) x^0.5,
              xlab = "Delusion Proneness (PDI)", ylab = "Time to get vaccinated (in months)")
abline(glm(SCRCVP$Delta_time~SCRCVP$PDI_total), col = 1, lwd = 3)


# Creating logistic regression models for the prediction of vaccination

# model1 -  DI total and vaccination
logmodel1=glm(SCRCVP$vaccination~PDI_total, data = SCRCVP, family='binomial')

# summarize logmodel results
summary(logmodel1, standardized = TRUE)
summ(logmodel1, standardized = TRUE)

# show odds ratio for logmodel 
exp(coef(logmodel1))

# calculcate 95% CI interval
confint(logmodel1)

# model2 -PDI total, ASRS, RAADS, STAI and Vaccination
logmodel2=glm(SCRCVP$vaccination~PDI_total + ASRS + STAI_Tscore + raads_now, data = SCRCVP, family='binomial')

# summarize logmodel results
summary(logmodel2, standardized = TRUE)
summ(logmodel2, standardized = TRUE)

# show odds ratio for logmodel 
exp(coef(logmodel2))

# calculcate 95% CI interval
confint(logmodel2)

# model 3 - pdi total, ASRS, RAADS, STAI, gender 0 psychiatric diagnosis and Vaccination
logmodel3=glm(SCRCVP$vaccination~PDI_total + ASRS + STAI_Tscore + raads_now + PsychDiagAny, data = SCRCVP, family='binomial')

summary(logmodel3)
summ(logmodel3, standardized = TRUE)

# show odds ratio for logmodel 
exp(coef(logmodel3))

# calculcate 95% CI interval
confint(logmodel3)

# model 4 - 
logmodel4=glm(SCRCVP$vaccination~PDI_total + ASRS + STAI_Tscore + raads_now + PsychDiagAny + genderm1 + education, data = SCRCVP, family='binomial')
summary(logmodel4)
summ(logmodel4, standardized = TRUE)

# show odds ratio for logmodel 
exp(coef(logmodel4))

# calculcate 95% CI interval
confint(logmodel4)

# model 5 (add age as covariate)
logmodel5=glm(SCRCVP$vaccination~PDI_total + age.x + ASRS + STAI_Tscore + raads_now + PsychDiagAny + genderm1 + education, data = SCRCVP, family='binomial')
summary(logmodel5)
summ(logmodel5, standardized = TRUE)

# show odds ratio for logmodel 
exp(coef(logmodel5))

# calculate 95% CI interval
confint(logmodel5)

# display variance inflation factor (check for multicollinearity)
vif(logmodel1)
vif(logmodel2)
vif(logmodel3)
vif(logmodel4)
vif(logmodel5)

# test interaction PDI * STAI
logmodelint=glm(SCRCVP$vaccination~PDI_total * STAI_Tscore ,data = SCRCVP, family='binomial')

summary(logmodelint, standardized = TRUE)

# create table for interaction model
logtbl_int <- tbl_regression(logmodelint, exponentiate= T) |> bold_p(t=0.05)

# create interaction plot
interact_plot(logmodelint, pred = PDI_total, modx = STAI_Tscore, plot.points = FALSE, interval = TRUE,
              x.label = "PDI", y.label = "Vaccinated (yes or  no)", legend.main = "STAI-T",
              int.width = 0.9, robust = TRUE, jitter = F, data = SCRCVP)


# create regression table for log models
logtbl1 <- tbl_regression(logmodel1, exponentiate = T, list("PDI_total" = "PDI")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(AIC))

logtbl2 <- tbl_regression(logmodel2, exponentiate=T, list("PDI_total" = "PDI", "ASRS" = "ASRS",
                                                        "STAI_Tscore" = "STAI-T", "raads_now" = "RAADS-N")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(AIC))

logtbl3 <- tbl_regression(logmodel3, exponentiate=T, list("PDI_total" = "PDI", "ASRS" = "ASRS",
                                                        "STAI_Tscore" = "STAI-T", "raads_now" = "RAADS-N",
                                                        "PsychDiagAny" = "Psychiatric Diagnosis")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(AIC))

logtbl4 <- tbl_regression(logmodel4, exponentiate=T, list("PDI_total" = "PDI", "ASRS" = "ASRS",
                                                        "STAI_Tscore" = "STAI-T", "raads_now" = "RAADS-N",
                                                        "PsychDiagAny" = "Psychiatric Diagnosis", "genderm1" = "Sex (M1F2)",
                                                        "education" = "Education")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(AIC))


logtbl5 <- tbl_regression(logmodel5, exponentiate=T, list("PDI_total" = "PDI", "age.x" = "Age", "ASRS" = "ASRS",
                                                          "STAI_Tscore" = "STAI-T", "raads_now" = "RAADS-N",
                                                          "PsychDiagAny" = "Psychiatric Diagnosis", "genderm1" = "Sex (M1F2)",
                                                          "education" = "Education")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(AIC))


# merge tables
logtbls_merged <- tbl_merge(list(logtbl1, logtbl2, logtbl3, logtbl4),
          tab_spanner = c("Model 1", "Model 2", "Model 3", "Model 4"))


# merge tables for controlling age
logtbls_merged_age <- tbl_merge(list(logtbl1, logtbl2, logtbl3, logtbl5),
                            tab_spanner = c("Model 1", "Model 2", "Model 3", "Model 4"))



#########
# g l m #
#########

# Creating glm to predict the time of vaccination

model_1 <- lm(SCRCVP$Delta_time ~ SCRCVP$PDI_total, data = SCRCVP)
summary(model_1)
summ(model_1)

model_2 <- lm(SCRCVP$Delta_time ~ SCRCVP$PDI_total + SCRCVP$ASRS + SCRCVP$STAI_Tscore + SCRCVP$raads_now, data = SCRCVP)
summary(model_2)
summ(model_2)

model_3 <- lm(SCRCVP$Delta_time ~ SCRCVP$PDI_total + SCRCVP$ASRS + SCRCVP$STAI_Tscore + SCRCVP$raads_now + SCRCVP$PsychDiagAny, data = SCRCVP)
summary(model_3)
summ(model_3)

model_4 <- lm(SCRCVP$Delta_time ~ SCRCVP$PDI_total + SCRCVP$ASRS + SCRCVP$STAI_Tscore +  SCRCVP$raads_now + SCRCVP$PsychDiagAny + SCRCVP$genderm1 + SCRCVP$education, data = SCRCVP)
summary(model_4)
summ(model_4)

model_5 <- lm(SCRCVP$Delta_time ~ SCRCVP$PDI_total + SCRCVP$age.x + SCRCVP$ASRS + SCRCVP$STAI_Tscore +  SCRCVP$raads_now + SCRCVP$PsychDiagAny + SCRCVP$genderm1 + SCRCVP$education, data = SCRCVP)
summary(model_5)
summ(model_5)

## create tables for linear regression models

lintbl1 <- tbl_regression(model_1, exponentiate=F, list("SCRCVP$PDI_total" = "PDI")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(adj.r.squared, AIC))

lintbl2 <- tbl_regression(model_2, exponentiate=F, list("SCRCVP$PDI_total" = "PDI", "SCRCVP$ASRS" = "ASRS",
                          "SCRCVP$STAI_Tscore" = "STAI-T", "SCRCVP$raads_now" = "RAADS-N")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(adj.r.squared, AIC))


lintbl3 <- tbl_regression(model_3, exponentiate=F, list("SCRCVP$PDI_total" = "PDI", "SCRCVP$ASRS" = "ASRS",
                          "SCRCVP$STAI_Tscore" = "STAI-T", "SCRCVP$raads_now" = "RAADS-N",
                          "SCRCVP$PsychDiagAny" = "Psychiatric Diagnosis")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(adj.r.squared, AIC))


lintbl4 <- tbl_regression(model_4, exponentiate=F, list("SCRCVP$PDI_total" = "PDI", "SCRCVP$ASRS" = "ASRS",
                          "SCRCVP$STAI_Tscore" = "STAI-T", "SCRCVP$raads_now" = "RAADS-N",
                          "SCRCVP$PsychDiagAny" = "Psychiatric Diagnosis", "SCRCVP$genderm1" = "Sex (M1F2)",
                          "SCRCVP$education" = "Education")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(adj.r.squared, AIC))


lintbl5 <- tbl_regression(model_5, exponentiate=F, list("SCRCVP$PDI_total" = "PDI", "SCRCVP$age.x" = "Age", "SCRCVP$ASRS" = "ASRS",
                                                        "SCRCVP$STAI_Tscore" = "STAI-T", "SCRCVP$raads_now" = "RAADS-N",
                                                        "SCRCVP$PsychDiagAny" = "Psychiatric Diagnosis", "SCRCVP$genderm1" = "Sex (M1F2)",
                                                        "SCRCVP$education" = "Education")) |> bold_p(t=0.05) |>
  add_glance_table(
    include = c(adj.r.squared, AIC))

# merge linear tables
lintbls_merged <- tbl_merge(list(lintbl1, lintbl2, lintbl3, lintbl4),
          tab_spanner = c("Model 1", "Model 2", "Model 3", "Model 4"))


#merge linear tables -  age correction
lintbls_merged_age <- tbl_merge(list(lintbl1, lintbl2, lintbl3, lintbl5),
                            tab_spanner = c("Model 1", "Model 2", "Model 3", "Model 4"))

#save table as editable word

#lintbls_merged %>% as_flex_table() %>%
# flextable::save_as_docx(path = here(paste0('linear_tbls_merged', ".docx")))


# test interaction PDI * STAI
model_int <- glm(Delta_time ~ PDI_total * STAI_Tscore, data = SCRCVP)
summary(model_int)

# create table
logtbl_int <- tbl_regression(model_int, exponentiate= F) |> bold_p(t = 0.05)

# display variance inflation factor (check for multicollinearity)
vif(model_1)
vif(model_2)
vif(model_3)
vif(model_4)
vif(model_5)

# create plot
interact_plot(model_int, pred = PDI_total, modx = STAI_Tscore, plot.points = FALSE, interval = TRUE,
              x.label = "PDI", y.label = "Time to get vaccinated (months)", legend.main = "STAI-T",
              int.width = 0.9, robust = TRUE, jitter = F, data = SCRCVP)


# Create data frame for Path/mediation Model
MeMo.Data <- data.frame(PDI = scale(SCRCVP$PDI_total), Vaccinated = (SCRCVP$vaccination), CCQTot = (SCRCVP$CCQTot), DeltaT = (SCRCVP$Delta_time))


# CCQTotal model
model1 <- ' # direct effect
             Vaccinated ~ c*PDI
           # mediator
             CCQTot ~ a1*PDI
             Vaccinated ~ b1*CCQTot
           # indirect effect
             a1b1 := a1*b1
           # total effect
             total := c + (a1*b1)
             '


# Fit data to model
fit1 <- sem(model1, data = MeMo.Data)

# Extract estimates

summary(fit1, fit.measures=T, rsq=T, ci = TRUE, standardize = TRUE)


# Create plot
semPaths(fit1, "model", "std", residuals = F, exoCov = F)


# Create Path Model plotwith diagrammeR
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
      "CCQ Total";
      
      node [shape = box,
      color = black,
      height = 0.7,
      fontsize = 12,
      fontname = helvetica,
      style = filled,
      fillcolor = white]
      "Vaccinated";
      
      
      #edge statements
      edge [color = grey,
      fontname = "helvetica-bold",
      fontstyle = bold]
      "PDI"->"CCQ Total" [label ="β = 0.29***", fontsize = 10, penwidth = 2, color = "#50b8fa", fontcolor= "#50b8fa"];
      "CCQ Total"-> "Vaccinated" [label ="      β = -0.34***", fontsize = 10, penwidth = 2, color = "#50b8fa", fontcolor = "#50b8fa"];
      "PDI"->"Vaccinated" [label ="β = -0.18**", fontsize = 10, fontcolor = grey];
      }
      ')

# show mediation figure
MeMoCCQTot


#### WORDCLOUD ####

# how many NAs/empty strings?
# line 222 must be removed
length(vaccinated_english[vaccinated_english == ''])/length(vaccinated_english) #14.5% NAs
length(unvaccinated_english[unvaccinated_english == '']) # 0% NAs
# + the ones originally not answered


# Remove the NAs from the dataframe
unvaccinated_english <- unvaccinated_english %>% na.omit()
vaccinated_english <- vaccinated_english %>% na.omit()


# Do not use "and", "the", because it cuts down sentences
keyword1 <- c("@\\S*", 
              "[\r\n]",  
              "[[:punct:]]",
              "https\\S*")


# convert all the strings of the dataframe vaccinated_english to lowercase with
vaccinated_english$text <- casefold(vaccinated_english$text, upper = FALSE)
unvaccinated_english$text <- casefold(unvaccinated_english$text, upper = FALSE)


# Subset
# create a data frame from $text
# a matrix can be more easily manipulated
vaccinated_engelska <- vaccinated_english[,2]
unvaccinated_engelska <- unvaccinated_english[,2]

# turn data frame into matrix
vaccinated_engelska <- as.matrix(vaccinated_engelska)
unvaccinated_engelska <- as.matrix(unvaccinated_engelska) # now it's an array

# remove punctuation and different characters
# for vaccinated
for (j in (keyword1)){ 
  
  vaccinated_engelska <- gsub(j, "", vaccinated_engelska)}


#Unvaccinated
for (j in (keyword1)){ 
  
  unvaccinated_engelska <- gsub(j, "", unvaccinated_engelska)}


########################
# create the wordcloud #
########################

# Create a vector for the loop
Q_un <- c()
Q_vac <- c()

for (i in 1:length(unvaccinated_engelska)) {
  Q_un <- c(Q_un, unvaccinated_engelska[i,1])
}

for (i in 1:length(vaccinated_engelska)) {
  Q_vac <- c(Q_vac, vaccinated_engelska[i,1])
}

# Create the doc
Doc_unvaccinated = Corpus(VectorSource(Q_un))
Doc_vaccinated = Corpus(VectorSource(Q_vac))

# Create a term document matrix
tdm_unvac=TermDocumentMatrix(Doc_unvaccinated)
tdm_vac=TermDocumentMatrix(Doc_vaccinated)

# Creating a word matrix
mat_unvac=as.matrix(tdm_unvac)
mat_vac=as.matrix(tdm_vac)

# Finding the words
words_unvac = sort(rowSums(mat_unvac),decreasing=TRUE) 
words_vac = sort(rowSums(mat_vac),decreasing=TRUE) 

# Actually, creating the data frame
df_unvac = data.frame(word = names(words_unvac),count=words_unvac)
df_vac = data.frame(word = names(words_vac),count=words_vac)

# plot the wordclouds
### bing sentiments
bing<-get_sentiments('bing')
df_unvac_bingsent<-inner_join(df_unvac, bing)
df_vac_bingsent<-inner_join(df_vac, bing)


matrix_unvac_bingsent<-acast(df_unvac_bingsent, word~sentiment, value.var='count', fill=0)
matrix_vac_bingsent<-acast(df_vac_bingsent, word~sentiment, value.var='count', fill=0)


## unvaccinated wordcloud
comparison.cloud(matrix_unvac_bingsent, colors=c('black', 'black'), title.size=1.4)

## vaccinated wordcloud
comparison.cloud(matrix_vac_bingsent, colors=c('black', 'black'), title.size=1.4)


## count number of words
SCRCVP$nwords <- str_count(SCRCVP$will_you_vaccinateC, "\\w+")


#### compare groups average number of words
t.test(SCRCVP$nwords[SCRCVP$vaccination==1], SCRCVP$nwords[SCRCVP$vaccination==0])


##############################
## create timeline of study ##
##############################

  # create data about timeline
  timeline.tb <-
    data.frame(what = c("Measurement of\nDelusion Proneness\n(N = 1032)", "WHO declare COVID-19\na pandemic",
                        "COVID-19 Conspiracy Questionnaire &\nBias Against Disconfirmatory Evidence\n(N = 577)",
                        "Vaccination start in Sweden", "Vaccines available\n for all age groups", "Vaccination survey (N= 273)"),
               when = ymd(c("2018-11-18", "2020-03-11", "2020-11-28", "2020-12-27",
                            "2021-07-14", "2021-10-21")),
               series = "Timeline")


# Now the labels would overlap, so we let R find a place for them using geom_text_repel() instead of geom_text().
  
  
  ggplot(timeline.tb, aes(x = when, y = series, label = what)) +
    geom_line() +
    geom_point() +
    geom_text_repel(direction = "y",
                    point.padding = 0.5,
                    hjust = 0,
                    box.padding = 1,
                    seed = 123) +
    scale_x_date(name = "", date_breaks = "5 months", date_labels = "%d %B %Y",
                 expand = expansion(mult = c(0.12, 0.12))) +
    scale_y_discrete(name = "") +
    theme_minimal()
  
  
# create second data frame 
  
  timeline_periods.tb <-
    data.frame(Study = c("Acar et al., 2022",
                           "Present study"),
               start = ymd(c("2018-11-18", "2020-12-27")),
               end = ymd(c("2020-12-15", "2021-11-10")),
               series = "Timeline")
  
  
# We highlight two periods using colours, and move the corresponding key to the top.
  
  ggplot(timeline.tb, aes(x = when, y = series)) +
    geom_line() +
    geom_segment(data = timeline_periods.tb,
                 mapping = aes(x = start, xend = end,
                               y = series, yend = series,
                               colour = Study),
                 linewidth = 2) +
    geom_point() +
    geom_text_repel(aes(label = what),
                    direction = "y",
                    point.padding = 0.5,
                    hjust = 0,
                    box.padding = 1,
                    seed = 123) +
    scale_x_date(name = "", date_breaks = "5 months", date_labels = "%d %B %Y",
                 expand = expansion(mult = c(0.12, 0.12))) +
    scale_y_discrete(name = "") +
    theme_minimal() +
    theme(legend.position = "top")



###############################################
## Descriptive table using gtsummary package ##
###############################################
  
# first, replace gender code to male/female
SCRCVP$genderm1[SCRCVP$genderm1 == 1] <- 'Male'
SCRCVP$genderm1[SCRCVP$genderm1 == 2] <- 'Female'



# Create table
descriptive_tbl <-
    SCRCVP |> 
    tbl_summary(include = c(age.x, genderm1, PDI_total, vaccinated, education, PsychDiagAny),
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c("{min} - {max}", 
                                     "{mean} ({sd})",
                                     "{median} [{p25} - {p75}]"),
    missing = "no",
    list(age.x = "Age", genderm1 = "Sex", PDI_total = "PDI", education = "Education", PsychDiagAny = "Psychiatric Diagnosis"),
    by = vaccinated) |>
    add_overall(last = TRUE, col_label = "**Total** \nN = {style_number(N)}") |>
    modify_header(
      stat_1 ~ "**Vaccinated**, \nN = {n}",
      stat_2 ~ "**Unvaccinated**, \nN = {n}")


