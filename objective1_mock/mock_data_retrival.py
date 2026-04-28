# script to get all 20 patients of 10 diseases from the simulated_patients_formatted.jsonl

import json
import pandas as pd
import jsonlines
import sys
import config


def main():
# make array of diseases of interest
# open file to read
# if disease_id == disease in array, 
    # add to list
# send all of list to output file

    # first 10 diseases
    diseases = ["404493", "1496", "329971", "87", "417", "242", "254704", "564", "7", "370927"]
    
    input_file = "simulated_patients_formatted.jsonl"
    input_dataframe = pd.read_json(input_file, lines=True)
    
    output_dataframe = input_dataframe[input_dataframe["disease_id"].isin(diseases)].copy()

    # mapping age category to number
    mapping = {'Onset_Infant': 0, 'Onset_Child': 1, 'Onset_Adolescent': 2, 'Onset_Adult': 3, 'Onset_Elderly': 4}
    output_dataframe['age'] = output_dataframe['age'].map(mapping)

    # map to csv
    output_file = "abbrev_mock_data.csv"
    try:
        output_dataframe.to_csv(output_file)
    except:
        print (f"Could not open file {output_file}")
    

if __name__ == "__main__":
    main()