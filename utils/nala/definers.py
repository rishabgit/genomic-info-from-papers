import abc
import json
import csv
import re
import os
import requests
from itertools import chain


MUT_CLASS_ID = 'e_2'


class NLDefiner:
    """
    Abstract class for determining whether an annotation in the dataset is a natural language (NL) mention.
    Subclasses that inherit this class should:
    * Be named [Name]NLDefiner
    * Implement the abstract method define
    * Set the value
    """

    @abc.abstractmethod
    def define(self, dataset):
        """
        :type dataset: nalaf.structures.data.Dataset
        """
        return


class ExclusiveNLDefiner(NLDefiner):
    """NLDefiner that uses an mixed approach of max words, regexs',
    min words and a dictionary of probable nl words."""

    VERSION = "20160801"

    def __init__(self):
        self.max_words = 4

        # Before the splitter was just based on space (' ')
        # Now hyphens or slashes sourrounded by letters or space-separated brackets also produce a word
        # note, the letters are actually captured so the words are not complete. But we only care about the length
        # self.word_tokenizer = re.compile(" +|(?:[a-zA-Z])[-/](?:[a-zA-Z])")
        self.word_tokenizer = re.compile(" ")

        conventions_file = 'data/regexs/nala_regex/regex_st.json'
        tmvarregex_file = 'data/regexs/nala_regex/RegEx.NL'
        # dictionary with common English words (regarded as NL) that appear in mutation mentions
        dict_nl_words_file = 'data/regexs/nala_regex/dict_nl_words_v2.json'

        # read in file regex_st.json into conventions array
        with open(conventions_file, 'r') as f:
            conventions = json.loads(f.read())
            self.compiled_regexps_custom = [re.compile(x) for x in conventions]

        # read RegEx.NL (only codes)
        with open(tmvarregex_file) as file:
            raw_regexps = list(csv.reader(file, delimiter='\t'))
            regexps = [x[0] for x in raw_regexps if len(x[0]) < 265]
            self.compiled_regexps = [re.compile(x) for x in regexps]

        with open(dict_nl_words_file) as f:
            dict_nl_words = json.load(f)
            self.compiled_dict_nl_words = list(set([re.compile(x, re.IGNORECASE) for x in dict_nl_words]))


    def define(self, dataset):
        for ann in chain(dataset.annotations(), dataset.predicted_annotations()):
            if ann.class_id == MUT_CLASS_ID:
                ann.subclass = self.define_string(ann.text)


    def define_string(self, query):
        """
        Checks for definer but on string.
        :param query:
        :return: subclass id, which in default is 0 (standard), 1(natural language) or 2 (semi standard)
        """
        matches_tmvar = (regex.match(query) for regex in self.compiled_regexps)
        matches_custom = (regex.match(query) for regex in self.compiled_regexps_custom)

        num_nl_words_lazy = None

        def num_nl_words():
            nonlocal num_nl_words_lazy
            if num_nl_words_lazy:
                return num_nl_words_lazy
            else:
                num_nl_words_lazy = sum((regex.search(query) is not None for regex in self.compiled_dict_nl_words))
                return num_nl_words_lazy

        words = self.word_tokenizer.split(query)

        if any(matches_custom) or any(matches_tmvar):
            return 0
        elif len(words) > self.max_words or num_nl_words() >= 2:
            return 1
        elif num_nl_words() >= 1:
            return 2
        else:
            return 0
