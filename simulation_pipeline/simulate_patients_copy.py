
import sys
import os
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import numpy as np
import random
import networkx
import obonet
import pickle
import json
from random import shuffle
import argparse
import logging
from collections import defaultdict

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns

import config
from simulation_pipeline.modules.patient_simulator_copy import PatientSimulator
from simulation_pipeline.modules.disease import Disease
from simulation_pipeline.modules.patient_copy import Patient
from simulation_pipeline.utils.util import create_disease_dict
pd.options.mode.chained_assignment = None

SEED = 42

'''
This is the entry point of the simulation. This program simulates rare disease patient profiles by combining Orphanet 
disease data with probabilistic models adding noise, dropout, and distractor genes. 

This file coordinates the components of the simulation pipeline by reading in data, running simulation logic, and 
exporting results. 
'''


# dict to keep track of the the # of genes added from each gene module
simulation_outcome_counts = {
    'non_syndromic_phenotype' : defaultdict(int),
    'common_false_positive' : defaultdict(int),
    'tissue_distractor' : defaultdict(int),
    'pathogenic_phenotype_irrelevant' : defaultdict(int),
    'insufficient_explainer' : defaultdict(int),
    'universal_distractor' : defaultdict(int),
    'phenotype_distractor' : defaultdict(int)
}

###################################################################################################
# HELPER UTILITIES (I/O AND NAMING)

def read_orphanet_data(args):
    '''
    Read in the phenotype, gene, and metadata from Orphanet data
    '''
    base_path = config.ORPHANET_PATH
    phenotypes = pd.read_csv(base_path / 'orphanet_final_disease_hpo_normalized_2015.tsv', sep='\t', dtype=str)
    genes = pd.read_csv(base_path / 'orphanet_final_disease_genes_normalized_2015.tsv', sep='\t', dtype=str)
    metadata = pd.read_csv(base_path / 'orphanet_final_disease_metadata_normalized_2015.tsv', sep='\t', dtype=str)

    print(f'There are {len(metadata.index)} diseases to simulate from Orphanet.')
    return phenotypes, genes, metadata

def get_dataset_statistics(patients):
    '''
    Calculate basic statistics about the simulated data:
        Average # of distractor genes
        Average # of positive and negative phenotypes
    '''
    
    n_distractors = [(len(patient.get_distractor_genes())) for patient in patients]
    n_positive_phenotypes = [(len(patient.get_hpo_set(is_positive = True))) for patient in patients]
    n_negative_phenotypes = [(len(patient.get_hpo_set(is_positive = False))) for patient in patients]

    logging.info('Total number of patients simulated: {}'.format(len(patients)))
    logging.info('Average number of distractor genes: {}'.format(sum(n_distractors)/len(n_distractors)))
    logging.info('Average number of positive phenotypes: {}'.format(sum(n_positive_phenotypes)/len(n_positive_phenotypes)))
    logging.info('Average number of negative phenotypes: {}'.format(sum(n_negative_phenotypes)/len(n_negative_phenotypes)))

def get_output_filename(args):
    '''
    Returns the output filename associated with the simulation run
    '''

    descriptors = []
    if args.random_genes: descriptors.append('_rand_genes')
    if args.no_phen_corruption: descriptors.append('_no_phencorrupt')
    if args.no_phen_dropout: descriptors.append('_no_phendrop')
    if args.no_phen_noise: descriptors.append('_no_phennoise')
    if args.no_gene_module_phen: descriptors.append('_no_phengenemod')
    if args.equal_probs: descriptors.append('_equal_probs')
    if args.sim_many_genes: descriptors.append('_many_genes')

    suffix = "".join(descriptors) if descriptors else ""

    filename = f"simulated_patients{suffix}.jsonl"

    if descriptors:
        filename = f"ablations/{filename}"

    return filename

###################################################################################################
# CORE SIMULATION LOGIC

def simulate_patient(args, simulator, patient, disease):
    '''
    simulates a single patient:
        1. initialize phenotypes
        2. for sampled n_distractor_genes, run distractor gene module 
        3. add noisy phenotypes
    '''
    logging.info('--- Initialize phenotypes --- ')
    
    # whether to perform phenotype corruption/dropout
    perform_phenotype_corruption = not args.no_phen_corruption
    perform_phenotype_dropout = not args.no_phen_dropout

    # initialize phenotypes for patient with disease
    simulator.initialize_phenotypes(patient, disease, \
        config.PROB_DROPOUT_POS, config.PROB_DROPOUT_NEG, config.PROB_CORRUPTION_POS, config.PROB_CORRUPTION_NEG, \
            perform_phenotype_dropout=perform_phenotype_dropout, perform_phenotype_corruption=perform_phenotype_corruption)

    # sample n_distractor_genes. Distractor genes are realistic "noise" in the gene list
    logging.info('--- Sampling Distractor Genes ---')

    # First, determine the number of distractor genes for the patient. This can be a random number based on N_DISTRACTORS_LAMBDA, or 
    # a hard coded number.

    n_distractor_genes = 1 + np.random.poisson(config.N_DISTRACTORS_LAMBDA - 1) 
    patient.n_distractor_genes = n_distractor_genes
    if args.sim_many_genes: # for the ablation analysis, we want to simulate patients with large number of genes and then down-sample
        n_distractor_genes = 100

    # initialize the number of distractor genes added, the distractor gene sampler, and the finish state of the simulation
    n_distractor_genes_added = 0
    gene_sampler_names = list(simulator.gene_samplers.keys())
    unfinished = False
    
   # Determine how to sample the distractor genes. 
    if args.random_genes: #randomly add genes
        while n_distractor_genes_added < n_distractor_genes: 
            did_add_gene = simulator.get_random_gene(patient) # randomly add a gene to a patient
            if did_add_gene == 'gene_added':
                n_distractor_genes_added += 1 
    else: # we alternatively use the gene modules to add distractors
        n_tries_to_add_any_distractor = 0 

        # Loop to add distractor genes.
        # First, the gene module is chosen. This is the function what will attempt to add the gene.
        # Then, the gene module will try to add the gene if possible for module.

        while n_distractor_genes_added < n_distractor_genes: 

            # sample gene module
            if args.equal_probs: # randomly
                sampled_name = str(np.random.choice(gene_sampler_names,1)[0])
            else: # from config file probabilities
                sampled_name = str(np.random.choice(gene_sampler_names,1, p=[config.NON_SYNDROM_PHEN_PROB, \
                config.COMMON_FP_PROB, config.TISSUE_DIST_PROB, config.PATH_PHEN_PROB, config.INSUFF_EXPLAIN_PROB, \
                config.UNIVERSAL_DIST_PROB, config.PHENO_DIST_PROB])[0])
                
            gene_module, _ = simulator.gene_samplers[sampled_name]

            # In some cases, a specific gene module isn't compatible with a patient/gene. We set MAX_ADD_GENE_ATTEMPTS to limit the # of tries.
            # We also log information about whether a gene module was successfully used to add a distractor gene to the patient.
            #   - 'gene_added'                -> successfully added distractor gene
            #   - 'gene_impossible_to_add'    -> module is unable to add this particular gene
            #   - 'gene_not_added'            -> module failed this round

            for i in range(config.MAX_ADD_GENE_ATTEMPTS):  
                
                 # initialize if the gene was added to the module
                did_add_gene = gene_module(patient, n_distractor_genes_added)

                # Tabulate Simulation Outcomes
                simulation_outcome_counts[sampled_name][did_add_gene] += 1
                
                if did_add_gene == 'gene_added': 
                    if args.verbose:
                        logging.warning(f'Took {i + 1} tries to add gene using gene module {sampled_name}')
                    n_distractor_genes_added += 1
                    break
                elif did_add_gene == 'gene_impossible_to_add': 
                    if args.verbose:
                        logging.warning(f'It is impossible for gene module {sampled_name} to add any gene ') 
                    break 
                elif did_add_gene != 'gene_not_added': 
                    raise Exception('One of the gene modules is returning an invalid category')
                
                # the gene module was not able to add any genes in MAX_ADD_GENE_ATTEMPTS attempts
                if i == config.MAX_ADD_GENE_ATTEMPTS - 1:
                    if args.verbose:
                        logging.warning(f'Reached maximum number of attempts using module {sampled_name} for patient {patient}')
                    simulation_outcome_counts[sampled_name]["max_attempts"] += 1
            
            if n_tries_to_add_any_distractor > config.MAX_ADD_ANY_DISTRACTOR_ATTEMPTS:
                if args.verbose:
                    logging.warning(f'Failed to add sufficient distractors to patient: {patient} ')
                unfinished = True
                break   
            n_tries_to_add_any_distractor += 1

    # Sample noisy phenotypes & add to the patient
    logging.info('--- Sampling Noisy Phenotypes ---')
    if not args.no_phen_noise:
        simulator.sample_noisy_phenotypes(patient)

    return patient, unfinished


def run_simulation(args, filename):
    '''
    Run Rare Disease Patient Simulation
        1. Read in Orphanet Data and cerate a disease dictionary
        2. Initialize the patient simulator
        3. Simulate Patients
    '''
    logging.basicConfig(format='%(message)s', level=logging.WARNING)

    # set random seed to ensure replicability
    random.seed(SEED)
    np.random.seed(SEED)
    os.environ['PYTHONHASHSEED']=str(SEED)

    # read in orphanet data & filter to diseases from timstamp = args.timstamp
    phenotypes, genes, metadata = read_orphanet_data(args)

    # create dict mapping orphanet_id -> Disease
    disease_dict = create_disease_dict(phenotypes, genes, metadata)
    
    # initialize simulator
    logging.info('------Initializing Patient Simulator------')
    add_genemodule_distractor_phenotypes = ~ args.no_gene_module_phen
    simulator = PatientSimulator(disease_dict, config.STRONG_PHENOTYPE_THRESH, config.WEAK_PHENOTYPE_THRESH, \
        add_distractor_phenotypes=add_genemodule_distractor_phenotypes, seed=SEED)

    
    # for each disease, simulate X patients 
    logging.info('\n\n------Simulating Patients------')
    # Determine how many patients per disease will be simulated. 
    # simulate_patient() will be called for each patient and then that patient will be added to the list of finished patients.
    # The list of finished patients will then be written to a simulation output file

    patients = []
    unfinished_patients = []
    for orphanet_id, disease in tqdm(disease_dict.items()):
        if args.random_n_patients:
            # NOTE: random_n_patients will simulate PATIENTS_PER_DISEASE on average
            n_sampled_patients = 1 + np.random.poisson(config.PATIENTS_PER_DISEASE - 1) 
        else:
            # NOTE: if we take the constant approach, exactly PATIENTS_PER_DISEASE will be simulated per disease
            n_sampled_patients = config.PATIENTS_PER_DISEASE
        
        for p in range(n_sampled_patients): 
            patient = Patient(disease)
            simulated_patient, unfinished = simulate_patient(args, simulator, patient, disease)
            if unfinished:
                unfinished_patients.append(simulated_patient)
            patients.append(simulated_patient)

    print("Simulation Outcome Rates by Gene Module:")
    print(simulation_outcome_counts)

    print(f'There were {len(unfinished_patients)} patients with incomplete distractor gene sets.')

    # calculate statistics on the newly simulated data
    get_dataset_statistics(patients)

    # randomly order patients
    shuffle(patients)

    # generate save filename based on different configurations
    logging.info('\n\n------Write Patients to File------')

    #save to jsonl file
    with open(config.SIMULATED_DATA_PATH / filename, "w") as output_file:
        for n, patient in enumerate(patients):
            patient_dict = patient.to_dict()
            patient_dict['id'] = n
            
            json.dump(patient_dict, output_file)
            output_file.write('\n')

###################################################################################################
# HIGH-LEVEL SIMULATION ORCHESTRATORS

def parse_args():
    parser = argparse.ArgumentParser(description='Simulate rare disease patients')

    parser.add_argument('--random_n_patients', action='store_true', help='Whether to specify a random # of patients per orphanet disease')
    parser.add_argument('--verbose', action='store_true', help='Additional logging')

    # parameters for phenotype ablation analyses. These remove different sources of phenotypes from the pipeline
    parser.add_argument('--no_phen_noise', action='store_true', help='Remove phenotypic noise from the pipeline.')
    parser.add_argument('--no_phen_dropout', action='store_true', help='Remove phenotypic dropout from the pipeline.')
    parser.add_argument('--no_phen_corruption', action='store_true', help='Remove phenotype corruption from the pipeline.')
    parser.add_argument('--no_gene_module_phen', action='store_true', help='Remove phenotypes added through gene modules from the pipeline.')

    # parameters for ablation analyses 
    parser.add_argument('--sim_many_genes', action='store_true', help='Simulate patients with many distractor genes each. These can be down-sampled for ablation analyses.')
    parser.add_argument('--random_genes', action='store_true', help='Randomly sample distractor genes instead of using gene modules.')
    parser.add_argument('--equal_probs', action='store_true', help='Sample gene modules with equal probabilities. This is used for the ablation analysis.')

    args = parser.parse_args()
    return args

###################################################################################################
# SIMULATION ENTRY POINT

if __name__ == "__main__":
    args = parse_args()

    # get filename to save to file
    filename = get_output_filename(args)

    # run simulation
    run_simulation(args, filename)
