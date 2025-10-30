#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

df = pd.read_csv("mock_phenotypes.csv")

# all patients present in the dataset
all_patients = pd.Index(df["patient_id"].unique(), name="patient_id")

# first want to find the number noisy HPOs per pateints (is_true == 0)
noisy_counts = df.query("is_true == 0").groupby("patient_id").size().reindex(all_patients, fill_value=0)

# number of true phenotypes per disease
true_counts = (
    pd.Series({d: len(hpos) for d, hpos in {
        #Multiple epiphyseal dysplasia, Al-Gazali type
        "D001": ["HP:0000256", "HP:0000272", "HP:0000316", "HP:0000369", "HP:0000470", "HP:0000767",
            "HP:0001274", "HP:0001373", "HP:0001513", "HP:0002007", "HP:0002758", "HP:0002857", 
            "HP:0005930", "HP:0006101", "HP:0012444", "HP:0030084"],
        #Beta-mannosidosis
        "D002": ["HP:0000365", "HP:0001249", "HP:0001250", "HP:0001999", "HP:0002205", "HP:0005247"],
        #Glycogen storage disease due to muscle phosphofructokinase deficiency
        "D003": ["HP:0001324", "HP:0001903", "HP:0002149", "HP:0002486"],
        #Isolated osteopoikilosis
        "D004": ["HP:0000086", "HP:0000252", "HP:0001482", "HP:0002652", "HP:0004322", "HP:0005789"],
        #Hyperekplexia-epilepsy syndrome
        "D005": ["HP:0001276", "HP:0002267", "HP:0002376", "HP:0002384", "HP:0007333", 
            "HP:0010818", "HP:0012018", "HP:0200134", "HP:0000243"],
    }.items()})
)

# estimate how many phenotypes are missing per patient based on 
# their known disease
obs_true_counts = df.query("is_true == 1").groupby(["patient_id", "disease_id"]).size()
obs_true_counts = obs_true_counts.reset_index(name="n_true_obs")
obs_true_counts["expected_true"] = obs_true_counts ["disease_id"].map(true_counts)
obs_true_counts["n_missing"] = obs_true_counts["expected_true"] - obs_true_counts["n_true_obs"]

# compute parametes
prob_noisy_positive = (noisy_counts > 0).mean()
noisy_pos_lambda = noisy_counts.mean()
noisy_neg_lambda = obs_true_counts["n_missing"].mean()

# output the variables
print ("PROB_NOISY_POSITIVE =", round(prob_noisy_positive, 3))
print ("NOISY_POS_PHEN_SAMPLES_LAMBDA =", round(noisy_pos_lambda, 3))
print ("NOISY_NEG_PHEN_SAMPLES_LAMBDA =", round (noisy_neg_lambda, 3))

