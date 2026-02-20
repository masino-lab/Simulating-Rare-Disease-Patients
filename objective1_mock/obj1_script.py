'''
Script to achieve objective 1:

Objective 1: model the prob distribution for the number of disease-positive and disease-negative HPO 
terms; assign to simulated patient conditioned on age and disease

Steps:
    1. find mock data: run pipeline and use first 100 patients for mock_data
        - simulated patients path: Simulating-Rare-Disease-Patients/data_release/simulated_patients/simulated_patients.jsonl
        - mock data path: Simulating-Rare-Disease-Patients/objective1_mock/mock_data.jsonl

### Create a script for the rest of the objectives [this] ###
### Simulating-Rare-Disease-Patients/objective1_mock/obj1_script.py ###

    2. with this extracted data, make a table
        - rows: patients
        - columns: number HPO positive terms, number HPO negative terms, age, disease???
    3. stratify this data by age and disease
        - if n<50, create a group
    4. each stratum, apply Python stats models
    5. simulate patients and compare distribution to test set
    
'''

'''
2. make a table:
    read in data
    make a table
        the rows are the rows in mock_data.jsonl
        store count of # +/- terms patient
        store disease, age

3. stratify the data by disease and age
    group into smaller and non-overlapping groups
    runtime huge to iterate a bunch
    prob a way to do this pandas
    i have the data sorted-- i need to figure out how to stratify data

4. apply statsmodels
    for each stratum, apply
    poisson, negative binomial, zero inflated poisson, zero inflated negative binomial, R mpcmp
    goal: fit the data to obtain the parameters
'''

import json
import pandas as pd
import numpy as np 
import statsmodels.api as sm
import jsonlines
import sys
import config

from statsmodels.discrete.discrete_model import Poisson


# function to read in jsonl of mock data
def read_jsonl(filename):
    # read in the mock dat patients
    print ('filename: ', filename)
    patients = []
    with jsonlines.open(filename) as reader:
        for patient in reader:
            patients.append(patient)
    
    return patients

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

# def stratify_data(patients_dataframe):
#     # stratify by disease and age
#     # patients_dataframe = patients_dataframe.sort_values(['disease_id', 'age']).reset_index(drop=True)
#     patients_dataframe = patients_dataframe.sort_values(['disease_id', 'age'])
#     return patients_dataframe

def statsmodels(patients_dataframe):
    # goal: fit the data to obtain the parameters
    # for each stratum in the dataframe, apply various stats models
    # maybe then call a goodness of fit test function

    # ok, this function: want to call a bunch of different stats functions
    poisson(patients_dataframe)

def poisson(patients_dataframe, stratified = False):
    # want to be able to make this work with the raw data and the stratified data
        # could always just make another function later on
    # first, find the mean of the dataframe of choice

    # start with just the positive
    mean_positive_hpo = patients_dataframe['positive_phenotypes_count'].mean()
    var_pos = patients_dataframe['positive_phenotypes_count'].var()
    # print (f"mean positive = {mean_positive_hpo}")
    # print (f"var positive = {var_pos}")

    # ind var is the age
    # dep is the HPO #
    pos_ind_var = patients_dataframe['age']
    pos_dep_var = patients_dataframe['positive_phenotypes_count']

    # # fit the Poisson model
    poisson_model = Poisson(pos_ind_var, pos_dep_var)
    pos_poisson_results = poisson_model.fit()
    # print summary
    # print (pos_poisson_results.summary())   

    # then the negative
    mean_negative_hpo = patients_dataframe['negative_phenotypes_count'].mean()
    var_neg = patients_dataframe['negative_phenotypes_count'].var()
    # print(f"mean neg = {mean_negative_hpo}")
    # print(f"var neg = {var_neg}")

    # ind var is the age
    # dep is the HPO #
    neg_ind_var = patients_dataframe['age']
    neg_dep_var = patients_dataframe['negative_phenotypes_count']

    # # fit the Poisson model
    poisson_model = Poisson(neg_ind_var, neg_dep_var)
    neg_poisson_results = poisson_model.fit()
    # print summary
    # print (neg_poisson_results.summary())

    filename = "poisson_output.txt"
    if stratified:
        with open("poisson_output.txt", "a") as file:
            try:
                print ("Stratified Positive HPO Terms Poisson Results", pos_poisson_results.summary(), "\n\n", file=file)
                print ("Stratified Negative HPO Terms Poisson Results", neg_poisson_results.summary(), file=file)
            except:    
                print (f"Could not open file {file}")
    else:
        with open("poisson_output.txt", "w") as file:
            try:    
                print ("Positive HPO Terms Poisson Results", pos_poisson_results.summary(), "\n\n", file=file)
                print ("Negative HPO Terms Poisson Results", neg_poisson_results.summary(), "\n\n", file=file)
            except:    
                print (f"Could not open file {file}")

    print('AIC pos:', pos_poisson_results.aic)
    print('AIC neg:', neg_poisson_results.aic)
    print('BIC pos:', pos_poisson_results.aic)
    print('BIC neg:', neg_poisson_results.aic)

def validate_strata(patients_dataframe):
    # want to find all disease/age stratum <50
    # for this case, let's use 2
   
    # every disease group, want to check the age. if count age >2, use in strata

    counts = (
        patients_dataframe.groupby(['disease_id', 'age']).size().reset_index(name='n_patients')
    )

    valid_strata = counts[counts['n_patients'] >= 5]
    patients_dataframe_filtered = patients_dataframe.merge(valid_strata[['disease_id', 'age']], 
    on=['disease_id', 'age'], how='inner')

    return patients_dataframe_filtered




def main():
    # read in patients from file
    print ('Reading in patients...')
    patients_dataframe = create_dataframe('abbrev_mock_data.jsonl')
    # before stratfing the data, want to do the poisson distb
    poisson(patients_dataframe)

    # then stratify
    patients_dataframe = validate_strata(patients_dataframe)
    # poisson again
    poisson(patients_dataframe, stratified=True)   


    # print(type(patients_dataframe))

    # print(type(patients_dataframe))
    write_output(patients_dataframe) 


if __name__ == "__main__":
    main()