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
'''

import json
import pandas as pd
import jsonlinesgit 
import sys
import config



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
    # we have an array of patinets
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


def main():
    # read in patients from file
    print ('Reading in patients...')
    patients_dataframe = create_dataframe('mock_data.jsonl')
    write_output(patients_dataframe) 


if __name__ == "__main__":
    main()