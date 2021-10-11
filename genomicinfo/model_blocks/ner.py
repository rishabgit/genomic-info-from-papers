from typing import List, Tuple

import numpy as np

from genomicinfo.entity_extraction.abstract_extractor import AbstractEntityExtractor
from genomicinfo.entity_extraction.ner.ner_bert import BERTEntityExtractor


class RegexBlock(AbstractEntityExtractor):

    def __init__(self, folder_name : str):
        self.ner_mut_extract = BERTEntityExtractor(folder_name)

    def extract(self, text: str) -> List[Tuple[str, str]]:
        mut_and_snippets = []
        # BERT NER
        mut_and_snippets = mut_and_snippets + self.ner_mut_extract.extract(text, span_size=self.span_size)

        return mut_and_snippets