import json
import pandas as pd
import numpy as np 
import statsmodels.api as sm
import jsonlines
import sys
import config
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

from statsmodels.discrete.discrete_model import Poisson

def create_dataframe(filename):
    # create dataframe from csv of patients
    patients_dataframe = pd.read_csv(filename)
    
    # drop the columns we don't need
    patients_dataframe = patients_dataframe.drop(columns=['true_genes', 'n_distractor_genes', 
        'distractor_genes', 'dropout_phenotypes', 'corruption_phenotypes', 'id'])
    # add columns: counts for +/- HPO 
    patients_dataframe.insert(2, 'positive_phenotypes_count', patients_dataframe['positive_phenotypes'].apply(len))
    patients_dataframe.insert(4,'negative_phenotypes_count', patients_dataframe['negative_phenotypes'].apply(len))

    # write to csv for R script
    output_file = "patients_dataframe.csv"
    try:
        patients_dataframe.to_csv(output_file)
    except:
        print (f"Could not open file {output_file}")

    # return dataframe
    return patients_dataframe

def validate_strata(patients_dataframe):
    ''' 
    long term: want to find all disease/age stratum <50
    short term: for this case, let's use 5
    every disease group, want to check the age. if count age >=5, use in strata 
    '''

    # create dataframe (counts) for each disease/age combo and then count(n_patients)
    counts = (patients_dataframe.groupby(['disease_id', 'age']).size().reset_index(name='n_patients'))
    # filtering data to only include n_patients >= 5 and storing in dataframe (valid_strata)
    valid_strata = counts[counts['n_patients'] >= 5]
    # merge -> inner merge to only keep patients_dataframe patients whose rows appear in valid_strata
    patients_dataframe_filtered = patients_dataframe.merge(valid_strata[['disease_id', 'age']], 
        on=['disease_id', 'age'], how='inner')

    # write to csv for R script
    output_file = "stratified_patients_dataframe.csv"
    try:
        patients_dataframe_filtered.to_csv(output_file)
    except:
        print (f"Could not open file {output_file}")

    # return this merged dataframe
    return patients_dataframe_filtered

def regression(patients_dataframe, stratified = False): 
    # dictionary for the regression model and outfile
    models = {
        sm.Poisson: "poisson.txt",
        sm.NegativeBinomial: "negative_binomial_output.txt", 
        sm.ZeroInflatedPoisson: "zero_inflated_poisson_output.txt",
        # sm.ZeroInflatedNegativeBinomialP: "zero_inflated_negative_binomial_output.txt"
    }

    # iterate through every item in dictionary
    for key, value in models.items():
        regression_model = key
        filename = value
        # generate regression models
        pos_results = model(patients_dataframe, regression_model, "positive")
        neg_results = model(patients_dataframe, regression_model, "negative")
            
        # stratified data called second and therefore appended to output file
        if stratified:
            access = "a"
            description = "Stratified"
        else:
            # non-stratified data called first and write to output file
            access = "w"
            description = ""
        
        # open file and write positive and negative summaries
        with open(filename, access) as file:
            try:
                print(f"{description} Positive HPO Count {pos_results.summary()} \n\n", file=file)
                print(f"{description} Negative HPO Count {neg_results.summary()} \n\n", file=file)
            except:
                print (f"Could not open file {file}")

        # goodness of fit test
        goodness_of_fit(patients_dataframe, regression_model, pos_results, neg_results, filename)

def model(patients_dataframe, regression_model, category):
    # defining ind and dep variables for regression
    ind_var = patients_dataframe[['age']]
    if category == "positive":
        dep_var = patients_dataframe['positive_phenotypes_count']
    else:
        dep_var = patients_dataframe['negative_phenotypes_count']

    # HPO regression model
    # Zero Inflated Poisson Model
    if regression_model == sm.ZeroInflatedPoisson:
        results_regression_model = regression_model(dep_var, ind_var, exog_infl=ind_var, inflation='logit')
        results = results_regression_model.fit()
    else:
        # need to add independent variable constant for regression
        ind_var = sm.add_constant(ind_var)
        # Zero Inflated Negative Binomial Model
        if regression_model == sm.ZeroInflatedNegativeBinomialP:
            results_regression_model = regression_model(dep_var, ind_var, exog_infl=ind_var, p=2)
            results = results_regression_model.fit()
        else:
            # Poisson and Negative Binomial models
            results_regression_model = regression_model(dep_var, ind_var)
            results = results_regression_model.fit()
    return results

def goodness_of_fit(patients_dataframe, regression_model, pos_results, neg_results, filename):
    with open (filename, "a") as file:
        try:
            for results in [pos_results, neg_results]:
                print (f"{'POSITIVE' if results == pos_results  else 'NEGATIVE'} HPO COUNTS:", file=file)
                # Poisson model requires mean and variance measures; doing the same of NB
                if (regression_model == sm.Poisson):
                    print (f"   Mean: {np.mean(patients_dataframe['positive_phenotypes_count' if results == pos_results  else 'negative_phenotypes_count' ])}", file=file)
                    print (f"   Variance: {np.var(patients_dataframe['positive_phenotypes_count'])}", file=file)
                # all models get AIC and BIC to compare
                print (f"   AIC: {results.aic}", file=file)
                print (f"   BIC: {results.aic} \n\n", file=file)
        except:
            print (f"Could not open file {file}")
    

def main():
    print ('Reading in patients...')
    # create dataframe
    patients_dataframe = create_dataframe('abbrev_mock_data.csv')
    # before stratifying the data, want to do the regressions
    regression(patients_dataframe)

    # then stratify 
    patients_dataframe = validate_strata(patients_dataframe)
    # regression again  
    regression(patients_dataframe, stratified=True)


if __name__ == "__main__":
    main()
