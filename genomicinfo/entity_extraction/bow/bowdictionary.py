import re
from typing import List, Tuple

from genomicinfo.entity_extraction.abstract_extractor import AbstractEntityExtractor


class BOWdictionary(AbstractEntityExtractor):

    def __init__(self):
        # words whose presence would automatically tick sentence posititve without any context
        super().__init__()
        self.list_of_words = [['substitution', 'downstream', 'deletion', 'frameshift'], ]

    @staticmethod
    def tokenize_string(string):
        sentence = string
        sentence = re.sub('([0-9])([A-Za-z])', r'\1 \2', sentence)
        # separate non-ascii characters into their own tokens
        sentence = re.sub('([^\x00-\x7F])', r' \1 ', sentence)
        sentence = re.sub('([\W\-_])', r' \1 ', sentence)
        return sentence.split()  # splits by white space

    def extract(self, text: str) -> List[Tuple[str, str]]:
        final_list = []
        for single_list in self.list_of_words:
            word_set = set(single_list)
            phrase_set = set(BOWdictionary.tokenize_string(text))
            if phrase_set >= word_set:
                final_list.append(('Invalid', text))
        return final_list
