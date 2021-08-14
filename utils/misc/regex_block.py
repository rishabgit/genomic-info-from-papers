import re
import os
import numpy as np
from pathlib import Path

from wbtools.db.generic import WBGenericDBManager
from wbtools.db.gene import WBGeneDBManager
from wbtools.lib.nlp.common import EntityType
   
    
class MutationFinder:
    '''
    Extract mutations using modified MutationFinder regular expressions.

    ...

    Attributes
    ----------
    regex_folder : str
        should contain the modified MutationFinder regex file - seth_modified.txt

    Methods - 
    -------
    __call__(text)
        extracts mutation from the text string  
        
    '''
    
    def __init__(self, regex_path):
        """ 
        regex_folder should contain the 4 regex files
        """
        self._regular_expressions = []

        regular_expressions_file = open(regex_path)
        for line in regular_expressions_file:
            line = line.strip()
            if not line.startswith('#'):
                if line.endswith('[CASE_SENSITIVE]'):
                    self._regular_expressions.append(\
                     re.compile(line[:line.rindex('[')]))
                else:
                    self._regular_expressions.append(\
                     re.compile(line,re.IGNORECASE))
                
    def __call__(self, text, span_size=150):
        final_list = []
        for regex in self._regular_expressions:    
            for m in regex.finditer(text):
                span = min(m.span('wt_res')[0],\
                        m.span('pos')[0],\
                        m.span('mut_res')[0]),\
                    max(m.span('wt_res')[1],\
                        m.span('pos')[1],\
                        m.span('mut_res')[1])           
                surrounding_text = (text[max(span[0]-span_size, 0):\
                                        min(len(text), span[1]+span_size)])
                raw_mut = (text[span[0]:span[1]+1]).strip()
                raw_mut = raw_mut[1:] if not raw_mut[0].isalnum() else raw_mut
                raw_mut = raw_mut[:-1] if not raw_mut[-1].isalnum() else raw_mut
                final_list.append([raw_mut.strip(), surrounding_text])
        return final_list


class TmVar:
    '''
    Extract mutations using tmVar regular expressions.

    ...

    Attributes
    ----------
    regex_folder : str
        should contain the 4 tmVar regex files - 
        MF.RegEx.2.txt, SNP.RegEx.txt, ProteinMutation.RegEx.txt, DNAMutation.RegEx.txt


    Methods - 
    -------
    __call__(text)
        extracts mutation from the text string  
        
    '''
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
    '''
    Additional set of regular expressions developed after studying the manually curated remarks. 

    ...

    Attributes
    ----------
    db_config : configparser object
        contains credentials to access wbtools and textpresso api
    extra_regex : bool
        compiles regex to required in following methods -  var_and_gene_close() and get_genes()
    locus_only : bool
        compiles only the regex which can filter protein mutations with locus and can be easily normalized 
    folder_path : str
        path to store the gene and variant names extracted from wbtools in npy file 

    Methods 
    -------
    __call__(text)
        extract mutation from the text string  
    var_and_gene_close(text)
        extract var and gene mentions whenever they're close to each other in the text string  
    get_genes(text)
        extract genes from the text string  
        
    '''
    def __init__(self, db_config, extra_regex=False, locus_only=False, folder_path='data/gsoc/wbtools'):
        
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
amino_acid_three_to_one_letter_map = \
    dict([('ALA','A'),('GLY','G'),('LEU','L'),('MET','M'),\
     ('PHE','F'),('TRP','W'),('LYS','K'),('GLN','Q'),('GLU','E'),('SER','S'),\
     ('PRO','P'),('VAL','V'),('ILE','I'),('CYS','C'),('TYR','Y'),('HIS','H'),\
     ('ARG','R'),('ASN','N'),('ASP','D'),('THR','T'),('XAA','X'),('GLX','Z'),\
     ('ASX','B'), ('TER', '*'), ('STP', '*')])

# A dictionary mapping amino acid names to their one-letter abbreviations
amino_acid_name_to_one_letter_map = \
    dict([('ALANINE','A'),('GLYCINE','G'),('LEUCINE','L'),\
     ('METHIONINE','M'),('PHENYLALANINE','F'),('TRYPTOPHAN','W'),\
     ('LYSINE','K'),('GLUTAMINE','Q'),('GLUTAMIC','E'),\
     ('GLUTAMATE','E'),('ASPARTATE','D'),('SERINE','S'),\
     ('PROLINE','P'),('VALINE','V'),('ISOLEUCINE','I'),('CYSTEINE','C'),\
     ('TYROSINE','Y'),('HISTIDINE','H'),('ARGININE','R'),\
     ('ASPARAGINE','N'),('ASPARTIC','D'),('THREONINE','T'), \
     ('OCHRE', 'OCHRE'), ('AMBER', 'AMBER'), ('OPAL', 'OPAL'), ('UMBER', 'UMBER'), \
     ('STOP', 'X'), ('TERM', 'X'), ('*', 'X')])

amino_dict = dict(zip(list('ABCDEFGHIKLMNPQRSTVWXYZ'),
                                list('ABCDEFGHIKLMNPQRSTVWXYZ')))
amino_dict.update(amino_acid_three_to_one_letter_map)
amino_dict.update(amino_acid_name_to_one_letter_map)

def normalize_mutations(mutation):
    norm_regex_patterns = [\
        "(?P<wt_res>[A-Za-z]+)[^A-Za-z0-9]*(?P<pos>[1-9][0-9]+)(?:|( +[A-Za-z\s]*)? +(in|to|into|for|of|by|with|at) +(either +)?((an|a) +)?|[^A-Za-z0-9]*)(?P<mut_res>[A-Za-z]+)",\
        "(?P<wt_res>[A-Za-z]+)[^A-Za-z0-9\s]*(?:| *(in|to|into|for|of|by|with|at) *(either +)?((an|a) +)?([^A-Za-z0-9]*)?|[^A-Za-z0-9]*)(?P<mut_res>[A-Za-z]+)(( +.* +)|[^A-Za-z0-9]*)(?P<pos>[1-9][0-9]+)",\
        "(?P<wt_res>[A-Za-z]+)[^A-Za-z0-9]+(?P<mut_res>[A-Za-z]+)[^A-Za-z0-9]*(?P<pos>[1-9][0-9]+)",\
        "(?P<pos>[1-9][0-9]+)[^A-Za-z0-9]*(?P<wt_res>[A-Za-z]+)(?:| *(in|to|into|for|of|by|with|at) *(either +)?((an|a) +)?([^A-Za-z0-9]*)?|[^A-Za-z0-9]*)(?P<mut_res>[A-Za-z]+)",\
        "(?P<pos>[1-9][0-9]+)[^A-Za-z0-9]*(?P<wt_res>[A-Za-z]+)[^A-Za-z0-9]+(?P<mut_res>[A-Za-z]+)",\
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