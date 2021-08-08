import re
import os
import numpy as np
from pathlib import Path

from wbtools.db.generic import WBGenericDBManager
from wbtools.db.gene import WBGeneDBManager
from wbtools.lib.nlp.common import EntityType


# A dictionary mapping three-letter amino acids codes onto one-letter
# amino acid codes
amino_acid_three_to_one_letter_map = \
    dict([('ALA','A'),('GLY','G'),('LEU','L'),('MET','M'),\
     ('PHE','F'),('TRP','W'),('LYS','K'),('GLN','Q'),('GLU','E'),('SER','S'),\
     ('PRO','P'),('VAL','V'),('ILE','I'),('CYS','C'),('TYR','Y'),('HIS','H'),\
     ('ARG','R'),('ASN','N'),('ASP','D'),('THR','T'),('XAA','X'),('GLX','Z'),\
     ('ASX','B'), ('TER', 'X'), ('STP', 'X')])

# A dictionary mapping amino acid names to their one-letter abbreviations
amino_acid_name_to_one_letter_map = \
    dict([('ALANINE','A'),('GLYCINE','G'),('LEUCINE','L'),\
     ('METHIONINE','M'),('PHENYLALANINE','F'),('TRYPTOPHAN','W'),\
     ('LYSINE','K'),('GLUTAMINE','Q'),('GLUTAMIC ACID','E'),\
     ('GLUTAMATE','E'),('ASPARTATE','D'),('SERINE','S'),\
     ('PROLINE','P'),('VALINE','V'),('ISOLEUCINE','I'),('CYSTEINE','C'),\
     ('TYROSINE','Y'),('HISTIDINE','H'),('ARGININE','R'),\
     ('ASPARAGINE','N'),('ASPARTIC ACID','D'),('THREONINE','T'), \
     ('OCHRE', 'X'), ('AMBER', 'X'), ('OPAL', 'X'), ('UMBER', 'X'), \
     ('STOP', 'X'), ('TERM', 'X'), ('*', '*')])
     
class MutationError(Exception):
        pass

class Mutation(object):
    """ A base class for storing information about mutations """

    def __init__(self, Position):
        """ Initalize the object

            Position: the sequence position or start position of the mutation
                (must be castable to an int)
        """
        try:
            self.__position = int(Position)
        except ValueError:
            self.__position = '<See original snippet for number>'
            # NOTE: commented below lines due to inconsistency in the MF created regex
            # raise MutationError("Position must be an integer")
        # if self.__position < 1:
        #     raise MutationError("Position must be greater than 0")

    def _get_position(self):
        return self.__position
    Position = property(_get_position)

    def __str__(self):
        raise NotImplementedError('Mutation subclasses must override str()')

    def __eq__(self, other):
        raise NotImplementedError('Mutation subclasses must override ==')

    def __ne__(self, other):
        raise NotImplementedError('Mutation subclasses must override !-')

    def __hash__(self):
        raise NotImplementedError('Mutation subclasses must override hash()')

# TODO: very cluttery. Refactor later pls
# _normalize_residue_identity() fn can be used directly
class PointMutation(Mutation):
    """ 
    A class for storing information about protein point mutations
    """

    # Define a mapping for residue identity inputs to one-letter
    # abbreviations. For simplicty of the normalization procedure, a
    # one-letter to one-letter 'mapping' is also included. This
    # eliminates the need for an independent validation step, since
    # any valid identity which is passed in will be a key in this dict,
    # and it avoids having to analyze which format the input residue
    # was passed in as.
    _abbreviation_lookup = dict(zip(list('ABCDEFGHIKLMNPQRSTVWXYZ'),
                                    list('ABCDEFGHIKLMNPQRSTVWXYZ')))
    _abbreviation_lookup.update(amino_acid_three_to_one_letter_map)
    _abbreviation_lookup.update(amino_acid_name_to_one_letter_map)

    def __init__(self, Position, WtResidue, MutResidue, originalMention):
        """ Initalize the object and call the base class init

            Position: the sequence position or start position of the mutation
                (castable to an int)
            WtResidue: the wild-type (pre-mutation) residue identity (a string)
            MutReside: the mutant (post-mutation) residue identity (a string)

            Residues identities are validated to ensure that they are within
             the canonical set of amino acid residues are normalized to their
             one-letter abbreviations.
        """
        self.__wt_residue = self._normalize_residue_identity(WtResidue)
        self.__mut_residue = self._normalize_residue_identity(MutResidue)
        self.__original_mention = originalMention
        Mutation.__init__(self, Position=Position)
        
        
    @staticmethod
    def _get_abbreviation_lookup(self):
        return self._abbreviation_lookup

    
    def _normalize_residue_identity(self, residue):
        """ Normalize three-letter and full residue names to their
             one-letter abbreviations. If a residue identity is passed in
             which does not fall into the set of canonical amino acids
             a MutationError is raised.

        """
        try:
                # convert residue to its single letter abbreviation after
                # converting it to uppercase (so lookup is case-insensitive)
                return self._abbreviation_lookup[residue.upper()]
        except AttributeError:
                # if residue cannot be converted to uppercase, it is not a
                # string, so raise an error
                raise MutationError('Residue must be a string')
        except KeyError:
                # if residue is not a key in self._abbreviation_lookup, it
                # it is not a standard amino acid residue, so raise an error
                raise MutationError(\
                 'Input residue not recognized, must be a standard residue: '\
                  + residue)

    def _get_wt_residue(self):
        return self.__wt_residue
    WtResidue = property(_get_wt_residue)

    def _get_mut_residue(self):
        return self.__mut_residue
    MutResidue = property(_get_mut_residue)

    def _get_original_mention(self):
        self.__original_mention = self.__original_mention.strip()
        raw_mut = self.__original_mention[1:] if not self.__original_mention[0].isalnum() else self.__original_mention
        raw_mut = self.__original_mention[:-1] if not self.__original_mention[-1].isalnum() else self.__original_mention
        return raw_mut.strip()
    OriginalMention = property(_get_original_mention)

    def __str__(self):
        """ Return original mutation snippet"""
        return ''.join([self.__wt_residue,str(self.Position),
          self.__mut_residue])

    def __eq__(self, other):
        """ Override ==

            Two PointMutation objects are equal if their Position, WtResidue,
             and MutResidue values are all equal.
        """
        if type(self) == type(other):
          return self.Position == other.Position and \
                 self.__wt_residue == other.__wt_residue and \
                 self.__mut_residue == other.__mut_residue
        return False

    def __ne__(self,other):
        """ Override !=

            Two PointMutation obects are not equal if either their Position,
             WtResidue, or MutResidue values differ.
        """
        return not self == other

    def __hash__(self):
        """ Override hash() """
        return hash(str(type(self)) + str(self))

#######

class MutationExtractor(object):
    """ A base class for extracting Mutations from text """

    def __init__(self,ignorecase=True):
        """ Initialize the object """
        pass

class MutationFinder(MutationExtractor):

    def __init__(self,regular_expressions):
        """ Initialize the object

            regular_expressions: an interative set of regular expressions to
                be applied for extracting mutations. These are in the
                default python syntax (i.e., perl regular expressions), with
                the single exception being that regular expressions which
                should be performed in a case sensitive manner should be
                followed by the string '[CASE_SENSITIVE]', with no spaces
                between it and the regular expression.

        """
        MutationExtractor.__init__(self)
        self._regular_expressions = []

        for regular_expression in regular_expressions:
            if regular_expression.endswith('[CASE_SENSITIVE]'):
                self._regular_expressions.append(\
                 re.compile(regular_expression[:regular_expression.rindex('[')]))
            else:
                self._regular_expressions.append(\
                 re.compile(regular_expression,re.IGNORECASE))

    def _post_process(self,mutations):
        """ Perform precision increasing post-processing steps

            Remove false positives indicated by:
              -> mutant and wild-type residues being identical (e.g. A42A)
        """
        for mutation in list(mutations):
            if type(mutation) is PointMutation:
                if mutation.WtResidue == mutation.MutResidue:
                    del mutations[mutation]

    def __call__(self,raw_text, span_size=100):
        """ Extract point mutations mentions from raw_text and return them in a dict
             raw_text: a string of text

            The result of this method is a dict mapping PointMutation objects to
             a list of spans where they were identified. Spans are presented in the
             form of character-offsets in text. If counts for each mention are
             required instead of spans, apply len() to each value to convert the
             list of spans to a count.

            Example result:
             raw_text: 'We constructed A42G and L22G, and crystalized A42G.'
             result = {PointMutation(42,'A','G'):[(15,19),(46,50)],
                       PointMutation(22,'L','G'):[(24,28)]}

             Note that the spans won't necessarily be in increasing order, due
              to the order of processing regular expressions.


        """
        result = {}
        for regular_expression in self._regular_expressions:
            for m in regular_expression.finditer(raw_text):

                span = min(m.span('wt_res')[0],\
                        m.span('pos')[0],\
                        m.span('mut_res')[0]),\
                    max(m.span('wt_res')[1],\
                        m.span('pos')[1],\
                        m.span('mut_res')[1])

                current_mutation = \
                PointMutation(m.group('pos'),m.group('wt_res'),\
                            m.group('mut_res'), raw_text[span[0]:span[1]+1])

                surrounding_text = raw_text[max(span[0]-span_size, 0):\
                                      min(len(raw_text), span[1]+span_size)]

                if current_mutation not in result.keys():
                  result[current_mutation] = surrounding_text

        self._post_process(result)
        return result


def mutation_finder_from_regex_filepath(regular_expression_filepath):
    """ 
    Constructs a MutationFinder object using regular expressions in a file
    """
    regular_expressions_file = open(regular_expression_filepath)

    regular_expressions = []
    # Read in and store regular expression, ignoring lines that begin with '#'
    for line in regular_expressions_file:
        line = line.strip()
        if not line.startswith('#'):
            regular_expressions.append(line)

    return MutationFinder(regular_expressions)


class TmVar:
    def __init__(self, regex_folder):
        """ 
        regex_folder should contain the 4 regex files
        """

        self._regular_expressions = []

        regular_expressions_file = open(os.path.join(regex_folder, 'MF.RegEx.2.txt'))
        for line in regular_expressions_file:
            reg, group = line.split('\t')
            # some regex in DNAMutation group might bring tons of FP
            # but switching these off didn't have any immediate affect on test data during dev
            # manual filtering required?
            if not reg.startswith('#'):
                if group == 'DNAMutation':
                    reg = '[^0-9A-Za-z]' + reg + ' '
                else:
                    reg = '[^0-9A-Za-z]' + reg + '[^0-9A-Za-z]'
                self._regular_expressions.append(re.compile(reg))

        regular_expressions_file = open(os.path.join(regex_folder, 'SNP.RegEx.txt'))
        for reg in regular_expressions_file:
            reg = '[^0-9A-Za-z]' + reg +'[^0-9A-Za-z]'
            self._regular_expressions.append(re.compile(reg))

        regular_expressions_file = open(os.path.join(regex_folder, 'ProteinMutation.RegEx.txt'))
        for reg in regular_expressions_file:
            reg = '[^0-9A-Za-z]' + reg + '[^0-9A-Za-z]'
            self._regular_expressions.append(re.compile(reg))

        regular_expressions_file = open(os.path.join(regex_folder, 'DNAMutation.RegEx.txt'))
        for reg in regular_expressions_file:
            reg = '[^0-9A-Za-z]' + reg + '[^0-9A-Za-z]'
            self._regular_expressions.append(re.compile(reg))

    def __call__(self, text, span_size=150):
        final_list = []
        for regex in self._regular_expressions:       
            for m in regex.finditer(text):
                span = (m.start(0), m.end(0))
                # have to post process to remove mutant and wild-type residues being identical (e.g. A42A)
                # no quick way to do it tho - naive way would be manually edit the regex with WRES and MRES like in MF regex
                # TODO: manual work time? :(                
                surrounding_text = (text[max(span[0]-span_size, 0):\
                                        min(len(text), span[1]+span_size)])
                raw_mut = (text[span[0]:span[1]]).strip()
                raw_mut = raw_mut[1:] if not raw_mut[0].isalnum() else raw_mut
                raw_mut = raw_mut[:-1] if not raw_mut[-1].isalnum() else raw_mut
                final_list.append([raw_mut.strip(), surrounding_text])
        return final_list
    
    
class BOWdictionary:
    def __init__(self):
        # words whose presence would automatically tick sentence posititve without any context
        self.list_of_words = [\
            ['substitution', 'downstream', 'deletion', 'frameshift'],\
            ]
        
    @staticmethod 
    def tokenize_string(string):
        sentence = string
        sentence = re.sub('([0-9])([A-Za-z])', r'\1 \2', sentence)
        # separate non-ascii characters into their own tokens
        sentence = re.sub('([^\x00-\x7F])', r' \1 ', sentence)
        sentence = re.sub('([\W\-_])', r' \1 ', sentence)
        return sentence.split()  # splits by white space

    def __call__(self, text):
        final_list = []
        for single_list in self.list_of_words:
            word_set = set(single_list)
            phrase_set = set(BOWdictionary.tokenize_string(text))
            if phrase_set >= word_set:
                final_list.append(['Invalid', text])
        return final_list
    
    
class CustomWBregex:
    def __init__(self, db_config, extra_regex=False,locus_only=False, folder_path='data/gsoc/wbtools'):
        
        if not set(os.listdir(folder_path)) >= set(['wb_alleles_variations.npy', 'wb_allele_designations.npy', 'all_gene_names.npy']):

            db_manager = WBGenericDBManager(
                dbname=db_config['wb_database']['db_name'], user=db_config['wb_database']['db_user'],
                password=db_config['wb_database']['db_password'], host=db_config['wb_database']['db_host'])

            alleles_variations = db_manager.get_curated_entities(entity_type=EntityType.VARIATION, exclude_id_used_as_name=False)
            allele_designations = db_manager.get_allele_designations()

            db_manager = WBGeneDBManager(dbname=db_config['wb_database']['db_name'], user=db_config['wb_database']['db_user'],\
                                            password=db_config['wb_database']['db_password'], host=db_config['wb_database']['db_host'])

            all_gene_names = db_manager.get_all_gene_names()
            genes = []
            for gene in all_gene_names.values():
                if gene:
                    genes.append(gene[0])

            np.save(folder_path +'/wb_alleles_variations.npy', alleles_variations)
            np.save(folder_path +'/wb_allele_designations.npy', allele_designations)
            np.save(folder_path +'/all_gene_names.npy', genes)

            # not sure if this line even releases any memory 
            alleles_variations = None
            allele_designations = None
            all_gene_names = None
            genes = None
            
        OPENING_CLOSING_REGEXES = [r'(', r')']
        wb_genes = np.load(folder_path +'/all_gene_names.npy')
        all_genes_list = Path(folder_path +'/genes.txt').read_text().split('\n')
        for g in wb_genes: all_genes_list.append(g)
        all_genes_list = [g for g in all_genes_list if len(g) > 1]
        all_genes_list = list(set(all_genes_list))
        all_genes = OPENING_CLOSING_REGEXES[0] + '|'.join(all_genes_list) + OPENING_CLOSING_REGEXES[1]

        # the allele regex and db idea was stolen from wbtools
        allele_designations = np.load(folder_path +'/wb_allele_designations.npy').astype('U6')
        alleles_variations = np.load(folder_path +'/wb_alleles_variations.npy').astype('U6')
        DB_VAR_REGEX = r'({designations}|m|p|ts|gf|lf|d|sd|am|cs)([0-9]+)'
        var_regex_1 = OPENING_CLOSING_REGEXES[0] + DB_VAR_REGEX.format(designations="|".join(allele_designations)) + OPENING_CLOSING_REGEXES[1]
        all_var = OPENING_CLOSING_REGEXES[0] + '|'.join(alleles_variations) + '|' + var_regex_1 + OPENING_CLOSING_REGEXES[1]

        variation_regex = [\
            all_var + r'[^A-Za-z].*[^A-Za-z](bp|base pair).*([ACTG]{8,}).*([ACTG]{8,})',\
            all_var + r'[^A-Za-z].{,50}(deletes|deletion|inserts|insertion).{,50}[^A-Za-z](bp|base pair).*(flank)',\
            all_var + r'[^A-Za-z].{,50}(deletes|deletion|inserts|insertion).{,50}(exon|intron) +[0-9]+',\
            all_var + r'[^A-Za-z].{,50}[^A-Za-z](bp|base pair).{,50}(deletes|deletion|inserts|insertion)',\
            all_var + r'[^A-Za-z].*( [CISQMNPKDTFAGHLRWVEYBZJX]) *(?: *-*> *| +(in|to|into|for|of|by|with|at)) +(either +)?((an|a) +)?( *NONSENSE +)?(TERM|STOP|AMBER|OCHRE|OPAL|UMBER)',\
            all_var + r'[^A-Za-z].*( [CISQMNPKDTFAGHLRWVEYBZJX]) *(?: *-*> *| +(in|to|into|for|of|by|with|at)) +(either +)?((an|a) +)?( *NONSENSE +)?([CISQMNPKDTFAGHLRWVEYBZJX*][^0-9A-Za-z])',\
            ]
        gene_var_combo = [\
            all_var + r'[^A-Za-z]{0,2}' + all_genes + r'[^A-Za-z]',\
            all_genes + r'[^A-Za-z]{0,2}' + all_var + r'[^A-Za-z]',
            ]

        # these regexes were written after manually looking at the curator remarks
        raw_regexs = [\
            '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX]) *(\(?[1-9][0-9]*\)?)(?: *-*> *| +(in|to|into|for|of|by|with|at))? +(either +)?((an|a) +)?( *NONSENSE +)?(TERM|STOP|AMBER|OCHRE|OPAL|UMBER)',\
            '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX])(?: *-*> *| +(in|to|into|for|of|by|with|at)) +(either +)?((an|a) +)?( *NONSENSE +)?(TERM|STOP|AMBER|OCHRE|OPAL|UMBER)',\
            '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX])(?: *-*> *| +(in|to|into|for|of|by|with|at) +(either +)?((an|a) +)?)( *NONSENSE +)?([CISQMNPKDTFAGHLRWVEYBZJX*])[^0-9A-Za-z].*(flank)',\
            '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX])(?: *-*> *| +(in|to|into|for|of|by|with|at) +(either +)?((an|a) +)?)( *NONSENSE +)?([CISQMNPKDTFAGHLRWVEYBZJX*])[^0-9A-Za-z].*([ACTG]{8,}).*([ACTG]{8,})',\
            '(?:^|[\s\(\[\'"/,;\-])(flank).*([ACTG]{8,}).*([ACTG]{8,})'
            ]
        
        extra_locus_only_regexs = [\
            '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX]) *(\(?[1-9][0-9]*\)?)(?: *-*> *| +(in|to|into|for|of|by|with|at))? +(either +)?((an|a) +)?( *NONSENSE +)?(TERM|STOP|AMBER|OCHRE|OPAL|UMBER)',\
            ]
        
        if locus_only:
            self._regular_expressions = [re.compile(r,re.IGNORECASE) for r in extra_locus_only_regexs]
        else:
            self._regular_expressions = [re.compile(r,re.IGNORECASE) for r in raw_regexs + variation_regex]
        
        if extra_regex:
            self._gene_var_regex = [re.compile(r,re.IGNORECASE) for r in gene_var_combo]
            
            # repeating same set of regex is definitely not efficient
            # self._gene_var_regex and the regex below can be easily combined
            # TODO: combine
            OPENING_CLOSING_REGEXES = [r'((?:^|[\s\(\[\'"/,;\-])', r'(?:^|[\s\(\[\'"/,;\-]))']
            all_genes = OPENING_CLOSING_REGEXES[0] + '|'.join(all_genes_list) + OPENING_CLOSING_REGEXES[1]
            self._all_genes = [re.compile(r,re.IGNORECASE) for r in [all_genes]]
            

            
    def __call__(self, text, span_size=150):
        final_list = []
        for regex in self._regular_expressions:      
            for m in regex.finditer(text):
                span = (m.start(0), m.end(0))    
                surrounding_text = (text[max(span[0]-span_size, 0):\
                                        min(len(text), span[1]+span_size)])
                raw_mut = (text[span[0]:span[1]])
                raw_mut = raw_mut[1:] if not raw_mut[0].isalnum() else raw_mut
                raw_mut = raw_mut[:-1] if not raw_mut[-1].isalnum() else raw_mut
                final_list.append([raw_mut.strip(), surrounding_text])

        return final_list

    
    def var_and_gene_close(self, text, span_size=150):
        final_list = []
        for regex in self._gene_var_regex:      
            for m in regex.finditer(text):
                span = (m.start(0), m.end(0))    
                surrounding_text = (text[max(span[0]-span_size, 0):\
                                        min(len(text), span[1]+span_size)])
                raw_mut = (text[span[0]:span[1]])
                raw_mut = raw_mut[1:] if not raw_mut[0].isalnum() else raw_mut
                raw_mut = raw_mut[:-1] if not raw_mut[-1].isalnum() else raw_mut
                final_list.append([raw_mut.strip(), 'Gene & Variant'])
        return final_list
    
    
    def get_genes(self, text):
        final_list = []
        for regex in self._all_genes:      
            for m in regex.finditer(text):
                span = (m.start(0), m.end(0))    
                raw = (text[span[0]:span[1]]).strip()
                raw = raw[1:] if not raw[0].isalnum() else raw
                raw = raw[:-1] if not raw[-1].isalnum() else raw
                if len(raw.strip()) > 1:
                    final_list.append([raw.strip(), 'Just gene'])
        return final_list
    
    
# Define a mapping for residue identity inputs to one-letter
# abbreviations. For simplicty of the normalization procedure, a
# one-letter to one-letter 'mapping' is also included. This
# eliminates the need for an independent validation step, since
# any valid identity which is passed in will be a key in this dict,
# and it avoids having to analyze which format the input residue
# was passed in as.
# little modified from the ones at top of this file
amino_acid_three_to_one_letter_map = \
    dict([('ALA','A'),('GLY','G'),('LEU','L'),('MET','M'),\
     ('PHE','F'),('TRP','W'),('LYS','K'),('GLN','Q'),('GLU','E'),('SER','S'),\
     ('PRO','P'),('VAL','V'),('ILE','I'),('CYS','C'),('TYR','Y'),('HIS','H'),\
     ('ARG','R'),('ASN','N'),('ASP','D'),('THR','T'),('XAA','X'),('GLX','Z'),\
     ('ASX','B'), ('TER', 'X'), ('STP', 'X')])

# A dictionary mapping amino acid names to their one-letter abbreviations
amino_acid_name_to_one_letter_map = \
    dict([('ALANINE','A'),('GLYCINE','G'),('LEUCINE','L'),\
     ('METHIONINE','M'),('PHENYLALANINE','F'),('TRYPTOPHAN','W'),\
     ('LYSINE','K'),('GLUTAMINE','Q'),('GLUTAMIC','E'),\
     ('GLUTAMATE','E'),('ASPARTATE','D'),('SERINE','S'),\
     ('PROLINE','P'),('VALINE','V'),('ISOLEUCINE','I'),('CYSTEINE','C'),\
     ('TYROSINE','Y'),('HISTIDINE','H'),('ARGININE','R'),\
     ('ASPARAGINE','N'),('ASPARTIC','D'),('THREONINE','T'), \
     ('OCHRE', 'X'), ('AMBER', 'X'), ('OPAL', 'X'), ('UMBER', 'X'), \
     ('STOP', 'X'), ('TERM', 'X'), ('*', '*')])

amino_dict = dict(zip(list('ABCDEFGHIKLMNPQRSTVWXYZ'),
                                list('ABCDEFGHIKLMNPQRSTVWXYZ')))
amino_dict.update(amino_acid_three_to_one_letter_map)
amino_dict.update(amino_acid_name_to_one_letter_map)

def normalize_mutations(mutation):
    norm_regex_patterns = [\
        "(?P<wt_res>[A-Za-z]+)[^A-Za-z0-9]*(?P<pos>[1-9][0-9]+)(?:|( +[A-Za-z\s]*)? +(in|to|into|for|of|by|with|at) +(either +)?((an|a) +)?|[^A-Za-z0-9]*)(?P<mut_res>[A-Za-z]+)",\
        "(?P<wt_res>[A-Za-z]+)[^A-Za-z0-9\s]*(?:| *(in|to|into|for|of|by|with|at) *(either +)?((an|a) +)?([^A-Za-z0-9]*)?|[^A-Za-z0-9]*)(?P<mut_res>[A-Za-z]+)(( +.* +)|[^A-Za-z0-9]*)(?P<pos>[1-9][0-9]+)",\
        "(?P<pos>[1-9][0-9]+)[^A-Za-z0-9]*(?P<wt_res>[A-Za-z]+)[^A-Za-z0-9]+(?P<mut_res>[A-Za-z]+)",\
        "(?P<wt_res>[A-Za-z]+)[^A-Za-z0-9]+(?P<mut_res>[A-Za-z]+)[^A-Za-z0-9]*(?P<pos>[1-9][0-9]+)",\
            ]
    normalized_mutations = []
    norm_regex = [re.compile(r,re.IGNORECASE) for r in norm_regex_patterns]
    for regex in norm_regex:  
        if normalized_mutations:
            break
        for m in regex.finditer(mutation):
            # Removing false positives
            if m.group('wt_res') != m.group('mut_res'):
                wt_res = amino_dict[m.group('wt_res').upper()]
                mut_res = amino_dict[m.group('mut_res').upper()]
                # this function would be called for a KNOWN single mutation
                # so we'd be getting a single normalizaed mutation
                normalized_mutations = [wt_res, m.group('pos'), mut_res]
                break
    normalized_mutations = ''.join(normalized_mutations)
    return normalized_mutations