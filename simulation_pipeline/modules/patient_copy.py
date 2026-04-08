import typing
from typing import Set, List
from collections import defaultdict, OrderedDict
from enum import Enum
import random
import logging

from simulation_pipeline.modules.disease import Disease
from simulation_pipeline.utils.util_copy import important_logger

'''
This module defines the Patient class which encapsulates the data and logic needed to model, modify, and analyze a 
single simulated patient and the presentation of that simulated patient's disease.

Functions within the Patient class will initialize a patient, modify phenotypes, incorporate distractor genes, and track these
modifications. Getter methods are provided to retrieve patient information. Additional methods are provided to format the 
patient data.
'''


class Patient():

    _patient_counter = 0
###################################################################################################
# INITIALIZATION HELPERS

    def _sample_age_range(self, age_ranges: Set[str]) -> str:
        '''
        Randomly sample an age of onset for a disease from all possible 
        ages of onset associated with the disease
        '''
        # if age_ranges is empty, then no age ranges were provided for the disease. 
        # In this case, we assume uniform distribution over all ages
        if len(age_ranges) == 0:
            age_ranges = ['Onset_Infant', 'Onset_Child', 'Onset_Adolescent', 'Onset_Adult', \
            'Onset_Elderly']
           
        return random.choice(list(age_ranges))

    def _sample_gene(self, gene_list: List[str]) -> str:
        '''
        Randomly sample a gene from list of genes associated with the disease
        '''
        assert len(gene_list) > 0, 'There are no genes in the gene list.'
        return random.choice(gene_list)
    
    def __init__(self, disease: Disease):
        # dictionaries map from hpo_id -> [source_modules]
        # true = phenotype related to the disease
        # distractor = phenotype is not related to the disease
        # positive = patient has the phenotype
        # negative = patient does not have the phenotype
        self.true_positive_hpo_dict = OrderedDict()
        self.distractor_positive_hpo_dict = OrderedDict()
        self.true_negative_hpo_dict = OrderedDict()
        self.distractor_negative_hpo_dict = OrderedDict()

        #dict mapping from gene_id -> [source_modules]
        self.distractor_gene_dict = OrderedDict()

        self.true_gene = self._sample_gene(disease.get_gene_list())
        self.orphanet_id = disease.orphanet_id
        self.age = self._sample_age_range(disease.get_age_range_set())
        self.n_distractor_genes = None

        self.dropout_phen_dict = defaultdict(list)
        self.corruption_phen_dict = defaultdict(list)
        self.patient_id = Patient._patient_counter
        Patient._patient_counter += 1


###################################################################################################
# CORE MUTATION METHODS

    def _get_phenotype_dict (self, is_positive: bool, is_distractor: bool) -> None:
        # returns the HPO dictionary associated with the corresponding is_positive and is_distractor
            # is_positive = True -> patient has the phenotype
            # is_positive = False -> patient does not have the phenotype
            # is_distractor = True -> phenotype not related to the disease
            # is_distractor = False -> phenotype is related to the disease

        if is_positive:
            if is_distractor:
                return self.distractor_positive_hpo_dict
            else:
                return self.true_positive_hpo_dict
        else:
            if is_distractor:
                return self.distractor_negative_hpo_dict
            else:
                return self.true_negative_hpo_dict


    def add_phenotype(self, hpo_id: str, source_module: str, is_positive: bool, is_distractor: bool = True, allow_switch: bool = False)-> None: #TODO change
        '''
        Adds phenotype & corresponding gene sampler source to the correct phenotype dict based on 
        whether the phenotype is present in the patient (is_positive)
        & whether the phenotype is true vs a distractor (is_distractor)

        NOTE: assumes that calling method checks whether the phenotype is not an obligate when added to neg phenotypes
        & that the phenotype is not excluded when added to pos phenotypes.
        This stipulation is satisfied within the can_add_phenotype() in patient_simulator.py which calls patient.add_phenotype
        '''
        current_dict = self._get_phenotype_dict(is_positive, is_distractor)
        opposite_distractor_dict = self._get_phenotype_dict(not is_positive, True)
        opposite_true_dict = self._get_phenotype_dict(not is_positive, False)
        
        # if in the neg/distractor dict or in the neg/true dict
        if hpo_id in opposite_distractor_dict or hpo_id in opposite_true_dict:
            if allow_switch:
                self.remove_phenotype(hpo_id, not is_positive)
            else:
                raise Exception(f"Tried to add phenotype to {'positive' if is_positive else 'negative'} list,"
                    f"that is already in {'negative' if is_positive else 'negative'} phenotype list & allow_switch = False")

        # add the hpo_id to the corresponding phenotype dict
        if hpo_id in current_dict:
            current_dict[hpo_id].append(source_module)
        else:
            current_dict[hpo_id] = [source_module]

    def remove_phenotype(self, hpo_id: str, is_positive: bool) -> None:
        '''
        Loop through all phenotype dicts & remove hpo_id 
        '''
        
        logging.info(f'Removing phenotype HPO_ID: {hpo_id}')
        important_logger.info(f'{self.patient_id}: Removing phenotype HPO_ID: {hpo_id}') 

        for is_distractor in (True, False):
            current_dict = self._get_phenotype_dict(is_positive, is_distractor)
            if hpo_id in current_dict:
                current_dict.pop(hpo_id)

    def add_distractor_gene(self, gene_id: str, source_module: str) -> None:
        '''
        Adds a distractor gene & corresponding source (i.e. which gene sampler)
        to the patient's dict of distractor genes
        
        NOTE: in future iterations of the pipeline, consider checking to see if there's an association between gene & 
            disease in other data sources beyond Orphanet
        '''
        logging.info(f'Adding distractor gene {gene_id} from module {source_module}')
        important_logger.info(f'{self.patient_id}: Adding distractor gene {gene_id} from module {source_module}')

        if gene_id in self.distractor_gene_dict:
            self.distractor_gene_dict[gene_id].append(source_module)
        else: self.distractor_gene_dict[gene_id] = [source_module]

###################################################################################################
# TRACKING FUNCTIONS

    def track_dropout_phenotype(self, phenotype: str, is_positive: bool) -> None:
        '''
        This method is mainly used for error analysis & ablations to keep track of what phenotypes were dropped
        during the drop out phase.
        '''
        if is_positive:
            self.dropout_phen_dict['positive_phenotypes'].append(phenotype)
        else:
            self.dropout_phen_dict['negative_phenotypes'].append(phenotype)

    def track_corruption_phenotype(self, phenotype: str, ancestor: str, is_positive: bool) -> None:
        '''
        This method is mainly used for error analysis & ablations to keep track of what phenotypes were dropped
        during the corruption phase.
        '''
        if is_positive:
            self.corruption_phen_dict['positive_phenotypes_corrupted'].append(phenotype)
            self.corruption_phen_dict['positive_phenotypes_replacement'].append(ancestor)
        else:
            self.corruption_phen_dict['negative_phenotypes_corrupted'].append(phenotype)
            self.corruption_phen_dict['negative_phenotypes_replacement'].append(ancestor)


###################################################################################################
# GETTER METHODS

    def get_hpo_set(self, is_positive: bool) -> Set:
    
        '''
        returns set of hpo_ids
        if is_positive, returns all hpos (both true & distractor) that the patient has
        Otherwise, returns all hpos that the patient does NOT have
        
        NOTE: if remove phenotype & list is empty, then the hpo will not be in keys()
        '''

        if is_positive:
            true_hpo_set = set(self.true_positive_hpo_dict.keys())
            distractor_hpo_set = set(self.distractor_positive_hpo_dict.keys())
            return true_hpo_set.union(distractor_hpo_set)
        else:
            true_hpo_set = set(self.true_negative_hpo_dict.keys() )
            distractor_hpo_set = set(self.distractor_negative_hpo_dict.keys())
            return true_hpo_set.union(distractor_hpo_set)
            
    def get_distractor_genes(self) -> List:
        '''
        returns list of distractor genes
        '''
        return self.distractor_gene_dict.keys()

    def get_true_phenotypes(self, is_positive: bool) -> List:
        '''
        returns dict of true phenotypes {hpo_id -> [source_modules]}
        if is_positive is true, returns phenotypes that the patient has
        Otherwise, returns phenotypes the patient does NOT have
        '''
        if is_positive:
            return self.true_positive_hpo_dict.keys()
        else:
            return self.true_negative_hpo_dict.keys()

    def get_distractor_phenotypes(self, is_positive: bool) -> List:
        '''
        returns dict of distractor phenotypes {hpo_id -> [source_modules]}
        if is_positive is true, returns phenotypes that the patient has
        Otherwise, returns phenotypes the patient does NOT have
        '''
        if is_positive:
            return self.distractor_positive_hpo_dict.keys()
        else:
            return self.distractor_negative_hpo_dict.keys()

    def get_true_gene(self) -> str:
        '''
        return patient's causal gene
        '''
        return self.true_gene
    
    def get_age(self) -> str:
        '''
        return patient's age of onset
        '''
        return self.age

    def get_orphanet_id(self) -> str:
        '''
        returns the orphanet id of the patient's disease
        '''
        return self.orphanet_id

#################################################################################################
# FORMATTING/REPRESENTING PATIENT DATA METHODS

    def to_dict(self) -> dict:
        '''
        Convert the patient's data into a dict
        '''
        pos_phenotypes = {**self.true_positive_hpo_dict, **self.distractor_positive_hpo_dict}
        neg_phenotypes = {**self.true_negative_hpo_dict, **self.distractor_negative_hpo_dict}
        d = {'disease_id': str(self.orphanet_id), 'true_genes': [str(self.true_gene)], 'age':str(self.age), 'positive_phenotypes': pos_phenotypes,
        'negative_phenotypes': neg_phenotypes, 'n_distractor_genes': self.n_distractor_genes, 'distractor_genes': self.distractor_gene_dict, \
            'dropout_phenotypes': self.dropout_phen_dict, 'corruption_phenotypes': self.corruption_phen_dict}
        return d

    def __str__(self, verbose=True) -> str:
        '''
        string representation of the patient's data
        '''

        s = '''
            Patient with disease: {}, true gene: {}, age: {}  + phenotypes: {}, - phenotypes: {}, and distractor genes: {}
            '''.format(self.orphanet_id, self.true_gene, self.age,
                self.get_hpo_set(is_positive=True), self.get_hpo_set(is_positive=False),
                self.get_distractor_genes())
        return s



