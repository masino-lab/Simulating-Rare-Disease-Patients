# load packages
# install.packages("COMPoissonReg")
library(COMPoissonReg)

# read in data
patinet_dataframe <- read.csv("patients_dataframe.csv")
# fit the model 
# model <- glm.cmp (y ~ x1+x2, data = my_data); y = outcome; x1/x2 = predictors
model <- glm.cmp(positive_phenotypes_count ~ age, data = patinet_dataframe)
# nrow(patinet_dataframe)
# mean(patinet_dataframe$positive_phenotypes_count)
# var(patinet_dataframe$positive_phenotypes_count)
# range(patinet_dataframe$positive_phenotypes_count)
# var(patinet_dataframe$age)


# stratified dataframe
stratified_patients <- read.csv("stratified_patients_dataframe.csv")
# fit the model 
stratified_model <- glm.cmp(positive_phenotypes_count ~ age, data = stratified_patients)
# nrow(stratified_patients)
# mean(stratified_patients$positive_phenotypes_count)
# var(stratified_patients$positive_phenotypes_count)
# range(stratified_patients$positive_phenotypes_count)
# var(stratified_patients$age)


# File name
file_name <- "COMP_output.txt"

# Write unstratified results (overwrite file)
cat("===== Unstratified Data =====\n\n", file = "COMP_output.txt")
capture.output(summary(model), file = "COMP_output.txt", append = TRUE)
cat ("\n\n\n", file = "COMP_output.txt", append = TRUE)

# Append the statsified results
cat("===== Stratsfied Data =====\n\n", file = "COMP_output.txt", append = TRUE)
capture.output(summary(stratified_model), file = "COMP_output.txt", append = TRUE)

