from typing import List, Tuple

import numpy as np

from genomicinfo.entity_extraction.abstract_extractor import AbstractEntityExtractor
from genomicinfo.entity_extraction.bow.bowdictionary import BOWdictionary


class BOWBlock(AbstractEntityExtractor):

    def __init__(self):
        self.bow_mut_extract = BOWdictionary()

    def extract(self, text: str) -> List[Tuple[str, str]]:
        mut_and_snippets = []
        # BERT NER
        mut_and_snippets = mut_and_snippets + self.bow_mut_extract.extract(text)

        return mut_and_snippets