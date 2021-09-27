import os
import re
from os import listdir
from typing import List, Pattern, Tuple

from genomicinfo.entity_extraction.abstract_extractor import AbstractEntityExtractor


OPENING_CLOSING_REGEXES = [r'(', r')']
DB_VAR_REGEX = r'({designations}|m|p|ts|gf|lf|d|sd|am|cs)([0-9]+)'


class RegexEntityExtractor(AbstractEntityExtractor):
    """Extract genomic location through regexes

    Attributes
    ----------
    regex_files_folder: str
            path to folder containing the regex files to load
    """

    def __init__(self, regex_files_folder: str = None, extract_surrounding_text: bool = True,
                 surrounding_text_placeholder: str = None, min_mention_length: int = None):
        super().__init__()
        self._regular_expressions = []
        self.get_surrounding_text = extract_surrounding_text
        self.surrounding_text_placeholder = surrounding_text_placeholder
        self.min_mention_length = min_mention_length
        if regex_files_folder:
            self._regular_expressions = self._load_regexes(regex_files_folder)

    @staticmethod
    def _load_regexes(regex_files_folder) -> List[Pattern]:
        """

        Parameters
        ----------
        regex_files_folder : str
            path to folder containing the regex files to load

        Returns
        -------
        List[Pattern]
            a list of regex patterns

        """
        return [re.compile(line.strip()) for regex_file in listdir(regex_files_folder) for line in
                open(os.path.join(regex_files_folder, regex_file)) if not line.startswith("#")]

    def extract(self, text: str, span_size: int = 150) -> List[Tuple[str, str]]:
        """
        Extract entities using regexes

        Parameters
        ----------
        text : str
           text to analyze
        span_size : int
           size in number of characters of the span of text surrounding the matched entities that is returned together
           with the entity

        Returns
        -------
        List[Tuple[str, str]]
            a list of tuples containing the matched entities and their surrounding text spans

        """
        ret_list = []
        mention_surr_set = set()
        for regex in self._regular_expressions:
            for m in regex.finditer(text):
                mention, surrounding_text = self._get_mention(text=text, mention_begin=m.start(0), mention_end=m.end(0),
                                                              span_size=span_size)
                if not self.min_mention_length or len(mention) > self.min_mention_length:
                    if mention + "_" + surrounding_text not in mention_surr_set:
                        ret_list.append((mention, surrounding_text))
                        mention_surr_set.add(mention + "_" + surrounding_text)
        return ret_list

    def _get_mention(self, text, mention_begin, mention_end, span_size):
        mention = (text[mention_begin:mention_end]).strip()
        mention = mention[1:] if not mention[0].isalnum() else mention
        mention = mention[:-1] if not mention[-1].isalnum() else mention
        mention = mention.strip()
        surrounding_text = (text[max(mention_begin - span_size, 0): min(len(text), mention_end + span_size)])
        return mention, (surrounding_text if self.get_surrounding_text else self.surrounding_text_placeholder)
