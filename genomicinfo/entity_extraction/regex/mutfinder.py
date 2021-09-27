import os
import re
from os import listdir
from typing import List, Pattern

from genomicinfo.entity_extraction.regex.regex_extractor import RegexEntityExtractor


class MutationFinderRegexEntityExtractor(RegexEntityExtractor):

    def __init__(self):
        super().__init__(os.path.join(os.path.dirname(__file__), "regex_files", "mutation_finder"))

    @staticmethod
    def _load_regexes(regex_files_folder) -> List[Pattern]:
        regexes = []
        for regex_file in listdir(regex_files_folder):
            for line in open(os.path.join(regex_files_folder, regex_file)):
                line = line.strip()
                if not line.startswith("#"):
                    if line.endswith('[CASE_SENSITIVE]'):
                        regexes.append(re.compile(line[:line.rindex('[')]))
                    else:
                        regexes.append(re.compile(line, re.IGNORECASE))
        return regexes

    def extract(self, text: str, span_size: int = 150):
        final_list = []
        for regex in self._regular_expressions:
            for m in regex.finditer(text):
                span = min(m.span('wt_res')[0],
                           m.span('pos')[0],
                           m.span('mut_res')[0]),\
                       max(m.span('wt_res')[1],
                           m.span('pos')[1],
                           m.span('mut_res')[1])
                final_list.append(self._get_mention(text=text, mention_begin=span[0],
                                                    mention_end=span[1],
                                                    span_size=span_size))
        return final_list
