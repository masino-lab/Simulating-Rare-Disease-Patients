# want to first look at the distribution of the OG file
# "simulated_patients_formatted.jsonl"
# then, want to change the parameters and compare

import json
import pandas as pd
import jsonlines
import sys
import config


def main():
    
    files = ["simulated_patients_formatted.jsonl", "NB_simulated_patients_formatted.jsonl"]

    for file in files:
        input_file = file
        patients_dataframe = pd.read_json(input_file, lines=True)
        
        # add columns: counts for +/- HPO 
        patients_dataframe.insert(3, 'positive_phenotypes_count', patients_dataframe['positive_phenotypes'].apply(len))
        patients_dataframe.insert(5,'negative_phenotypes_count', patients_dataframe['negative_phenotypes'].apply(len))

        average_pos = patients_dataframe['positive_phenotypes_count'].mean()
        var_pos = patients_dataframe['positive_phenotypes_count'].var()
        average_neg = patients_dataframe['negative_phenotypes_count'].mean()
        var_neg = patients_dataframe['negative_phenotypes_count'].var()

        print (average_pos, var_pos, average_neg, var_neg)


        filename = "test_data.txt"
        with open (filename, "a") as outfile:
            if file == "simulated_patients_formatted.jsonl":
                print (f"original test set", file=outfile)
            else:
                print (f"NB set", file=outfile)
            try:

                print (f"average pos: {average_pos} var pos: {var_pos} \naverage neg: {average_neg} var neg: {var_neg}", file=outfile)
            except:
                print (f"Could not open file {outfile}")

    

if __name__ == "__main__":
    main()