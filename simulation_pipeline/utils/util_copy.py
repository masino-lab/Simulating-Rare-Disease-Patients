import pandas as pd   
import obonet
import jsonlines
import sys
import config
import logging

from simulation_pipeline.modules.disease import Disease


def get_nonspecific_phenotype_list(hpo_ontology):
        '''
        Creates a list of phenotypes in the HPO ontology that are non-specific & therefore
        shouldn't be added to the simulated patients. 
        '''
        lvl1 = list(hpo_ontology.predecessors('HP:0000001')) #root HPO code
        lvl2 = []
        for p in lvl1:
            lvl2.extend(list(hpo_ontology.predecessors(p)))
        #only include level 3 for branch under phenotypic abnormality
        lvl3 = []
        for p in list(hpo_ontology.predecessors('HP:0000118')): #phenotypic abnormality
            lvl3.extend(list(hpo_ontology.predecessors(p)))
        lvl4 = []
        for p in lvl3:
            lvl4.extend(list(hpo_ontology.predecessors(p)))

        broad_phen_list = ['HP:0000001'] + lvl1 + lvl2 
        nonspecific_phen_list = broad_phen_list + lvl3 + lvl4
        return broad_phen_list, nonspecific_phen_list


def create_disease_dict(orphanet_phenotypes, orphanet_genes, orphanet_metadata):
    '''
    Creates a dictionary that maps from orphanet_id -> Disease

    '''
    hpo_ontology = obonet.read_obo(config.HPO_FILE) 
    broad_phen_list, _ = get_nonspecific_phenotype_list(hpo_ontology)

    disease_dict = {}
    for orphanet_id in orphanet_metadata['OrphaNumber'].unique(): #only loop through diseases in metadata because these have been filtered by year
        phenotypes = orphanet_phenotypes.loc[orphanet_phenotypes['OrphaNumber'] == orphanet_id]
        genes = orphanet_genes.loc[orphanet_genes['OrphaNumber'] == orphanet_id]
        metadata = orphanet_metadata.loc[orphanet_metadata['OrphaNumber'] == orphanet_id]
        hpo_ids = phenotypes['HPO_ID'].tolist()
        hpo_freqs = phenotypes['HPO_Freq'].tolist()
        disease_name = phenotypes['Disorder_Name'].tolist()[0]
        age_options = ['Onset_Infant', 'Onset_Child', 'Onset_Adolescent', 'Onset_Adult', \
            'Onset_Elderly']

        should_include = [metadata[age].values[0] for age in age_options]
        age_of_onsets = [age for age, to_include in zip(age_options, should_include) if to_include == 'True']
        
        # only include phenotypes that are not extremely broad (according to HPO hierarchy)
        phenotype_dict = {hpo_id:freq for hpo_id, freq in zip(hpo_ids, hpo_freqs) if hpo_id not in broad_phen_list}

        gene_list = genes['Ensembl_ID'].tolist()
        disease = Disease(orphanet_id, disease_name, gene_list, phenotype_dict, age_of_onsets)
        disease_dict[orphanet_id] = disease
    return disease_dict



    ##############################################
# Read in/write patients
def read_udn_patients(filename):
    patients = []
    with jsonlines.open(filename) as reader:
        for patient in reader:
            patients.append(patient)
    return patients


def read_simulated_patients(filename):
    patients = []
    with jsonlines.open(filename) as reader:
        for patient in reader:
            if type(patient['positive_phenotypes']) == dict:
                patient['positive_phenotypes'] = list(patient['positive_phenotypes'].keys())
            if type(patient['negative_phenotypes']) == dict:
                patient['negative_phenotypes'] = list(patient['negative_phenotypes'].keys())
            if type(patient['distractor_genes']) == dict:
                patient['distractor_genes'] = [gene for gene in list(patient['distractor_genes'].keys())] 
            patients.append(patient)    
    return patients

# function created before i realized you can create more than one logger
# 
# def write_logging_to_file(message: str, patients = None, path: str = "loggingOutput.txt"):
#     with open(path, "a") as file:
#         if patients is not None:
#             file.write(f"[Patient {patients.simulation_id}] message + \n")
#         else:
#             file.write(message + "\n")

def setup_logger(name, log_file, level=logging.INFO):
    # function to create multiple loggers and to write those logs to specified files

    # create a logger instance
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # create file handler and set level and formatter
    handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    handler.setFormatter(formatter)

    # add handler to logger
    logger.addHandler(handler)

    # avoid printing to console (root logger)
    logger.propagate = False 

    return logger

important_logger = setup_logger('important_info', 'loggingOutput.txt', level=logging.INFO)
