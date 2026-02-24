# load packages
library(COMPoissonReg)

# read in data
patinet_dataframe <- read.csv("patients_dataframe.csv")
# fit the model 
model <- glm.cmp(positive_phenotypes_count ~ age, data = patinet_dataframe)

# stratified dataframe
stratified_patients <- read.csv("stratified_patients_dataframe.csv")
# fit the model 
stratified_model <- glm.cmp(positive_phenotypes_count ~ age, data = stratified_patients)

sink ("COMP_output.txt")
# print a summary to output file
print(summary(model))
print(summary(stratified_model))

