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
    # we have an array of patients
    patients_dataframe = pd.read_json(filename, lines=True)
    
    # drop the columns we don't need
    patients_dataframe = patients_dataframe.drop(columns=['true_genes', 'n_distractor_genes', 
        'distractor_genes', 'dropout_phenotypes', 'corruption_phenotypes', 'id'])
    # add columns: counts for +/- HPO 
    patients_dataframe.insert(2, 'positive_phenotypes_count', patients_dataframe['positive_phenotypes'].apply(len))
    patients_dataframe.insert(4,'negative_phenotypes_count', patients_dataframe['negative_phenotypes'].apply(len))

    return patients_dataframe

def write_output(patients_dataframe):
    # write file function
    filename = "output.jsonl"
    try:
        patients_dataframe.to_json(filename, orient="records", lines=True)
    except:
        print (f"Could not open file {filename}")

def validate_strata(patients_dataframe):
    # want to find all disease/age stratum <50
    # for this case, let's use 5
   
    # every disease group, want to check the age. if count age >=5, use in strata

    # create dataframe (counts) for each disease/age combo and then count(n_patients)
    counts = (patients_dataframe.groupby(['disease_id', 'age']).size().reset_index(name='n_patients'))

    # filtering data to only include n_patients >= 5 and storing in dataframe (valid_strata)
    valid_strata = counts[counts['n_patients'] >= 5]

    # merge -> inner merge to only keep patients_dataframe patients whose rows appear in valid_strata
    patients_dataframe_filtered = patients_dataframe.merge(valid_strata[['disease_id', 'age']], 
        on=['disease_id', 'age'], how='inner')

    # return this merged dataframe
    return patients_dataframe_filtered

def regression(patients_dataframe, stratified = False): 
    models = {
        sm.Poisson: "poisson.txt",
        sm.NegativeBinomial: "negative_binomial_output.txt", 
        sm.ZeroInflatedPoisson: "zero_inflated_poisson_output.txt",
        # sm.ZeroInflatedNegativeBinomialP: "zero_inflated_negative_binomial_output.txt"
    }

    for key, value in models.items():
        regression_model = key
        filename = value
        
        # positive HPO
        pos_ind_var = patients_dataframe[['age']]
        pos_dep_var = patients_dataframe['positive_phenotypes_count']

        # positive HPO regression model
        if regression_model == sm.ZeroInflatedPoisson:
            pos_regression_model = regression_model(pos_dep_var, pos_ind_var, exog_infl=pos_ind_var, inflation='logit')
            pos_results = pos_regression_model.fit()
        else:
            pos_ind_var = sm.add_constant(pos_ind_var)
            if regression_model == sm.ZeroInflatedNegativeBinomialP:
                pos_regression_model = regression_model(pos_dep_var, pos_ind_var, exog_infl=pos_ind_var, p=2)
                pos_results = pos_regression_model.fit()
            else:
                pos_regression_model = regression_model(pos_dep_var, pos_ind_var)
                pos_results = pos_regression_model.fit()

        # negative HPO regression model
        neg_ind_var = patients_dataframe[['age']]
        neg_dep_var = patients_dataframe['negative_phenotypes_count']

        if regression_model == sm.ZeroInflatedPoisson:
            neg_regression_model = regression_model(neg_dep_var, neg_ind_var, exog_infl=neg_ind_var, inflation='logit')
            neg_results = neg_regression_model.fit()
        else:
            neg_ind_var = sm.add_constant(neg_ind_var)
            if regression_model == sm.ZeroInflatedNegativeBinomialP:
                neg_regression_model = regression_model(pos_dep_var, neg_ind_var, exog_infl=neg_ind_var, p=2)
                neg_results = neg_regression_model.fit()
            else:
                neg_regression_model = regression_model(pos_dep_var, neg_ind_var)
                neg_results = neg_regression_model.fit()
            
        if stratified:
            access = "a"
            description = "Stratified"
        else:
            access = "w"
            description = ""
        
        with open(filename, access) as file:
            try:
                print(f"{description} Positive HPO Count {pos_results.summary()} \n\n", file=file)
                print(f"{description} Negative HPO Count {neg_results.summary()} \n\n", file=file)
            except:
                print (f"Could not open file {file}")

        # goodness of fit test
        goodness_of_fit(patients_dataframe, regression_model, pos_results, neg_results, filename)

def goodness_of_fit(patients_dataframe, regression_model, pos_results, neg_results, filename):

    # poisson first
    # need the mean and var comparison 
    # also need to do AIC and BIC
    with open (filename, "a") as file:
        try:
            for results in [pos_results, neg_results]:
                print (f"{'POSITIVE' if results == pos_results  else 'NEGATIVE'} HPO COUNTS:", file=file)
                if (regression_model == sm.Poisson):
                    print (f"   Mean: {np.mean(patients_dataframe['positive_phenotypes_count' if results == pos_results  else 'negative_phenotypes_count' ])}", file=file)
                    print (f"   Variance: {np.var(patients_dataframe['positive_phenotypes_count'])}", file=file)
                print (f"   AIC: {results.aic}", file=file)
                print (f"   BIC: {results.aic} \n\n", file=file)
            
        except:
            print (f"Could not open file {file}")
    

def main():
    # read in patients from file
    print ('Reading in patients...')
    patients_dataframe = create_dataframe('abbrev_mock_data.jsonl')

    # before stratfing the data, want to do the regressions
    # poisson(patients_dataframe)
    regression(patients_dataframe)

    ### then stratify ###
    patients_dataframe = validate_strata(patients_dataframe)

    # regression again  
    regression(patients_dataframe, stratified=True)


    # write_output(patients_dataframe) 


if __name__ == "__main__":
    main()