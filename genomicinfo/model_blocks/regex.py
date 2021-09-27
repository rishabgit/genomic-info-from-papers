from typing import List, Tuple

import numpy as np

from genomicinfo.entity_extraction.abstract_extractor import AbstractEntityExtractor
from genomicinfo.entity_extraction.bow.bowdictionary import BOWdictionary
from genomicinfo.entity_extraction.regex.mutfinder import MutationFinderRegexEntityExtractor
from genomicinfo.entity_extraction.regex.tmvar import TMVarRegexEntityExtractor
from genomicinfo.entity_extraction.regex.wb_custom import WBCustomRegexEntityExtractor


def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))


class RegexBlock(AbstractEntityExtractor):

    def __init__(self, allele_variations: List[str], allele_designations: List[str], gene_names: List[str],
                 span_size: int = 150):
        self.mf_mut_extract = MutationFinderRegexEntityExtractor()
        self.tmvar_mut_extract = TMVarRegexEntityExtractor()
        self.custom_mut_extract = WBCustomRegexEntityExtractor(allele_variations=allele_variations,
                                                               allele_designations=allele_designations,
                                                               gene_names=gene_names)
        self.bow_mut_extract = BOWdictionary()
        self.span_size = span_size

    def extract(self, text: str) -> List[Tuple[str, str]]:
        mut_and_snippets = []
        # MutationFinder
        mut_and_snippets = mut_and_snippets + self.mf_mut_extract.extract(text, span_size=self.span_size)
        # tmVar
        mut_and_snippets = mut_and_snippets + self.tmvar_mut_extract.extract(text, span_size=self.span_size)
        # Custom patterns
        mut_and_snippets = mut_and_snippets + self.custom_mut_extract.extract(text, span_size=self.span_size)
        # Bag of words
        mut_and_snippets = mut_and_snippets + self.bow_mut_extract.extract(text)

        if mut_and_snippets:
            mut_and_snippets = unique_rows(mut_and_snippets).tolist()
        return mut_and_snippets
