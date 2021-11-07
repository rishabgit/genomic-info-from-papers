from typing import List, Tuple

import numpy as np

from genomicinfo.entity_extraction.abstract_extractor import AbstractEntityExtractor
from genomicinfo.model_blocks.ner import NERBlock
from genomicinfo.model_blocks.bow import BOWBlock
from genomicinfo.model_blocks.regex import RegexBlock


def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))


class Pipeline(AbstractEntityExtractor):

    def __init__(self, use_ner=True, use_bow=True, use_regex=True):
        self.use_ner = use_ner
        self.use_bow = use_bow
        self.use_regex = use_regex
        
    def add_details(self, span_size=150, **kwargs):
        '''
        Attributes
        ----------
        folder_name : str
            path to folder containing the saved model weights
        allele_variations : List[str]
            list of allele variations to look for
        allele_designations : List[str]
            list of allele designations to look for
        gene_names : List[str]
            list of gene names to look for
        span_size : int (default: 150)
            size of surrounding text from the mut mention
        '''
        if self.use_ner: 
            folder_name = kwargs.get('folder_name', None)
            self.ner_mut_extract = NERBlock(folder_name)

        if self.use_regex:
            allele_variations = kwargs.get('allele_variations', None)
            allele_designations = kwargs.get('allele_designations', None)
            gene_names = kwargs.get('gene_names', None)
            span_size = kwargs.get('span_size', 150)
            self.regex_mut_extract = RegexBlock(allele_variations, allele_designations, gene_names, span_size)

        if self.use_bow:
            self.bow_mut_extract = BOWBlock()


    def extract(self, text: str) -> List[Tuple[str, str]]:
        mut_and_snippets = []
        # NER 
        if self.use_ner: mut_and_snippets = mut_and_snippets + self.ner_mut_extract.extract(text)
        # Regex
        if self.use_regex: mut_and_snippets = mut_and_snippets + self.regex_mut_extract.extract(text)
        # Bag of words
        if self.use_bow: mut_and_snippets = mut_and_snippets + self.bow_mut_extract.extract(text)

        if mut_and_snippets:
            mut_and_snippets = unique_rows(mut_and_snippets).tolist()
        return mut_and_snippets
