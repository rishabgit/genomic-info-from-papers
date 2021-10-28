import re
from typing import List, Pattern

from genomicinfo.entity_extraction.regex.regex_extractor import RegexEntityExtractor, OPENING_CLOSING_REGEXES, DB_VAR_REGEX


class WBCustomRegexEntityExtractor(RegexEntityExtractor):
    """
    Additional set of regular expressions developed after studying the manually curated remarks.

    ...

    Attributes
    ----------
    locus_only : bool
        compiles only the regex which can filter protein mutations with locus and can be easily normalized

    """
    def __init__(self, allele_variations: List[str], allele_designations: List[str], gene_names: List[str],
                 locus_only=False, extract_surrounding_text: bool = True, surrounding_text_placeholder: str = None,
                 min_mention_length: int = None):
        self.gene_names = gene_names
        self.allele_variations = allele_variations
        self.allele_designations = allele_designations
        self.locus_only = locus_only
        # the allele regex and db idea was stolen from wbtools
        var_regex_1 = OPENING_CLOSING_REGEXES[0] + DB_VAR_REGEX.format(
            designations="|".join(self.allele_designations)) + \
                      OPENING_CLOSING_REGEXES[1]
        self.all_var = OPENING_CLOSING_REGEXES[0] + '|'.join(self.allele_variations) + '|' + var_regex_1 + \
                  OPENING_CLOSING_REGEXES[1]
        self.all_genes = OPENING_CLOSING_REGEXES[0] + '|'.join(self.gene_names) + OPENING_CLOSING_REGEXES[1]

        super().__init__(extract_surrounding_text=extract_surrounding_text,
                         surrounding_text_placeholder=surrounding_text_placeholder,
                         min_mention_length=min_mention_length)

    def _load_regexes(self, regex_files_folder) -> List[Pattern]:
        extra_locus_only_regexs = [
            '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX]) *(\(?[1-9][0-9]*\)?)(?: *-*> *| +(in|to|into|for|of|by|with|at))? *(either +)?((an|a) +)?( *NONSENSE +)?(TERM|STOP|AMBER|OCHRE|OPAL|UMBER)', \
            '(?:^|[\s\(\[\'"/,;\-])(?P<wt_res>[CISQMNPKDTFAGHLRWVEYBZJX])\((?P<pos>[+-]?[1-9][0-9]+(?:\s?[+-]\s?[1-9][0-9]*)?)\)(?P<mut_res>([CISQMNPKDTFAGHLRWVEYBZJX*]|[Ss]top|[Tt]erm))(?=([.,\s)\]\'":;\-?!/]|$))', \
            '(?:^|[\s\(\[\'"/,;\-])(?P<wt_res>[CISQMNPKDTFAGHLRWVEYBZJX])\((?P<pos>[1-9][0-9]*)\)(?: *-*> *| +(in|to|into|for|of|by|with|at) +)(?P<mut_res>[CISQMNPKDTFAGHLRWVEYBZJX])(?=[([.,\s)\]\'":;\-?!/]|$])', \
            '(?:^|[\s\(\[\'"/,;\-])(?P<wt_res>(?:A(?:LA(?:NINE)?|MBER|RG(?:ININE)?|S(?:P(?:AR(?:T(?:IC ACID|ATE)|AGINE))?|N|X))|MET(?:HIONINE)?|CYS(?:TEINE)?|L(?:EU(?:CINE)?|YS(?:INE)?)|O(?:CHRE|PAL)|I(?:SOLEUCINE|LE)|UMBER|T(?:ER(?:M)?|R(?:P|YPTOPHAN)|HR(?:EONINE)?|YR(?:OSINE)?)|VAL(?:INE)?|P(?:HE(?:NYLALANINE)?|RO(?:LINE)?)|S(?:T(?:P|OP)|ER(?:INE)?)|GL(?:U(?:TAM(?:ATE|I(?:C ACID|NE)))?|N|Y(?:CINE)?|X)|HIS(?:TIDINE)?|XLE))\((?P<pos>[1-9][0-9]*)\)(-*>)(?P<mut_res>(?:A(?:LA(?:NINE)?|MBER|RG(?:ININE)?|S(?:P(?:AR(?:T(?:IC ACID|ATE)|AGINE))?|N|X))|MET(?:HIONINE)?|CYS(?:TEINE)?|L(?:EU(?:CINE)?|YS(?:INE)?)|O(?:CHRE|PAL)|I(?:SOLEUCINE|LE)|UMBER|T(?:ER(?:M)?|R(?:P|YPTOPHAN)|HR(?:EONINE)?|YR(?:OSINE)?)|VAL(?:INE)?|P(?:HE(?:NYLALANINE)?|RO(?:LINE)?)|S(?:T(?:P|OP)|ER(?:INE)?)|GL(?:U(?:TAM(?:ATE|I(?:C ACID|NE)))?|N|Y(?:CINE)?|X)|HIS(?:TIDINE)?|XLE))(?=([.,\s)\]\'":;\-?!/]|$))'
        ]
        if self.locus_only:
            return [re.compile(r, re.IGNORECASE) for r in extra_locus_only_regexs]
        else:
            variation_regex = [
                self.all_var + r'[^A-Za-z].*[^A-Za-z](bp|base pair).*([ACTG]{8,}).*([ACTG]{8,})',
                self.all_var + r'[^A-Za-z].{,50}(deletes|deletion|inserts|insertion).{,50}[^A-Za-z](bp|base pair).*(flank)',
                self.all_var + r'[^A-Za-z].{,50}(deletes|deletion|inserts|insertion).{,50}(exon|intron) +[0-9]+',
                self.all_var + r'[^A-Za-z].{,50}[^A-Za-z](bp|base pair).{,50}(deletes|deletion|inserts|insertion)',
                self.all_var + r'[^A-Za-z].*( [CISQMNPKDTFAGHLRWVEYBZJX]) *(?: *-*> *| +(in|to|into|for|of|by|with|at)) +(either +)?((an|a) +)?( *NONSENSE +)?(TERM|STOP|AMBER|OCHRE|OPAL|UMBER)',
                self.all_var + r'[^A-Za-z].*( [CISQMNPKDTFAGHLRWVEYBZJX]) *(?: *-*> *| +(in|to|into|for|of|by|with|at)) +(either +)?((an|a) +)?( *NONSENSE +)?([CISQMNPKDTFAGHLRWVEYBZJX*][^0-9A-Za-z])',
            ]
            raw_regexs = [
                '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX])(?: *-*> *| +(in|to|into|for|of|by|with|at)) +(either +)?((an|a) +)?( *NONSENSE +)?(TERM|STOP|AMBER|OCHRE|OPAL|UMBER)', \
                '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX])(?: *-*> *| +(in|to|into|for|of|by|with|at) +(either +)?((an|a) +)?)( *NONSENSE +)?([CISQMNPKDTFAGHLRWVEYBZJX*])[^0-9A-Za-z].*(flank)', \
                '(?:^|[\s\(\[\'"/,;\-])([CISQMNPKDTFAGHLRWVEYBZJX])(?: *-*> *| +(in|to|into|for|of|by|with|at) +(either +)?((an|a) +)?)( *NONSENSE +)?([CISQMNPKDTFAGHLRWVEYBZJX*])[^0-9A-Za-z].*([ACTG]{8,}).*([ACTG]{8,})', \
                '(?:^|[\s\(\[\'"/,;\-])(flank).*([ACTG]{8,}).*([ACTG]{8,})'
            ]
            return [re.compile(r, re.IGNORECASE) for r in raw_regexs + extra_locus_only_regexs +
                    variation_regex]


class WBGeneVarRegexEntityExtractor(WBCustomRegexEntityExtractor):

    def __init__(self, allele_variations: List[str], allele_designations: List[str], gene_names: List[str]):
        super().__init__(allele_variations=allele_variations, allele_designations=allele_designations,
                         gene_names=gene_names, extract_surrounding_text=False,
                         surrounding_text_placeholder="Gene & Variant")

    def _load_regexes(self, regex_files_folder) -> List[Pattern]:
        gene_var_combo = [
            self.all_var + r'[^A-Za-z]{0,2}' + self.all_genes + r'[^A-Za-z]',
            self.all_genes + r'[^A-Za-z]{0,2}' + self.all_var + r'[^A-Za-z]',
        ]
        return [re.compile(r, re.IGNORECASE) for r in gene_var_combo]


class WBGeneRegexGenomicLocExtractor(WBCustomRegexEntityExtractor):

    def __init__(self, allele_variations: List[str], allele_designations: List[str], gene_names: List[str]):
        super().__init__(allele_variations=allele_variations, allele_designations=allele_designations,
                         gene_names=gene_names, extract_surrounding_text=False,
                         surrounding_text_placeholder="Just Gene", min_mention_length=1)

    def _load_regexes(self, regex_files_folder) -> List[Pattern]:
        return [re.compile(r, re.IGNORECASE) for r in [self.all_genes]]


class GenomeVersionRegexEntityExtractor(WBCustomRegexEntityExtractor):
    def __init__(self, allele_variations: List[str], allele_designations: List[str], gene_names: List[str]):
        super().__init__(allele_variations=allele_variations, allele_designations=allele_designations,
                         gene_names=gene_names, extract_surrounding_text=False,
                         surrounding_text_placeholder="Genome Version")

    def _load_regexes(self, regex_files_folder) -> List[Pattern]:
        # table 1 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6424801/pdf/nihms-1011867.pdf
        genome_vers_regex = ['ce2', 'ce4', 'ce6', 'ce8', 'ce10', 'ce11']
        genome_vers_regex = r'(' + '|'.join(genome_vers_regex) +  r')'
        genome_vers_regex = r'(?:^|[\s\(\[\'"/,;\-])' + genome_vers_regex
        return [re.compile(r) for r in [genome_vers_regex]]


class AnnotationVersionRegexEntityExtractor(WBCustomRegexEntityExtractor):
    def __init__(self, allele_variations: List[str], allele_designations: List[str], gene_names: List[str]):
        super().__init__(allele_variations=allele_variations, allele_designations=allele_designations,
                         gene_names=gene_names, extract_surrounding_text=False,
                         surrounding_text_placeholder="Annotation Version")

    def _load_regexes(self, regex_files_folder) -> List[Pattern]:
        # table 1 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6424801/pdf/nihms-1011867.pdf
        # and https://wormbase.org//about/release_schedule#0--10
        annotation_vers_regex = ['WS'+str(x) for x in range(120, 296)]
        annotation_vers_regex =  r'(' + '|'.join(annotation_vers_regex) +  r')'
        annotation_vers_regex = r'(?:^|[\s\(\[\'"/,;\-])' + annotation_vers_regex
        return [re.compile(r) for r in [annotation_vers_regex]]


