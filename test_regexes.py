from pprint import pprint
from settings import setSettings

from regex_wrapper import regex_block, extra_info_block
# from wbtools.literature.corpus import CorpusManager


def test_regexes(settings):
    # 1 Regex
    # Has high precision and low recall.
    # Consists of 3 parts - MutationFinder, tmVar and some custom patterns from WB papers
    # 1.1 Mutation Finder [Link](https://github.com/divyashan/MutationFinder),
    # Modified regex from SETH [Link](https://github.com/rockt/SETH/blob/master/resources/mutations.txt)
    # 1.2 tmVar [Link](https://www.ncbi.nlm.nih.gov/research/bionlp/Tools/tmvar/)
    # 1.3 Extra custom regexs

    mf_mut_extract = settings['mf_mut_extract']
    tmvar_mut_extract = settings['tmvar_mut_extract']
    custom_mut_extract = settings['custom_mut_extract']
    bow_mut_extract = settings['bow_mut_extract']

    result = mf_mut_extract('wa7 ) is a null allele, because a 1-bp deletion')
    pprint(result)  # -> []

    result = tmvar_mut_extract('wa7 ) is a null allele, because a 1-bp deletion')
    pprint(result)  # -> []

    result = custom_mut_extract('n2888 bp AAAAAAAAAAAAAAAAA TTTTTTTTTTTTTTTTTT')
    pprint(result)  # -> []

    result = custom_mut_extract.get_genes('Our results suggest that in C. elegans ,the trimethylation of histone H3 lysine 36 by MET-1/Set2 promotes a transcriptional repression cascade mediated by a NuRD-like complex and by the trimethylation of histone H3K9 by a SETDB1-like HMT.')
    pprint(result)  # -> [['MET-1', 'Just gene']]

    result = custom_mut_extract.get_genome_vers('The genome version used in this paper is ce4.')
    pprint(result)  # -> [['ce4', 'Genome Version']]

    result = custom_mut_extract.get_annotation_vers('When needing to refer to a specific version of the reference genome, it was therefore sufficient (and convenient) to use the WormBase release number that the genome was taken from (e.g. WS128).')
    pprint(result)  # [['MET-1', 'Just gene']] [['ce4', 'Genome Version']] [['WS128', 'Annotation Version']]

    # 1.4 Bag of words

    result = bow_mut_extract('This mutation deletes 471bp of the promoter region, the transcriptional start and 56 amino acids of the second exon.')
    pprint(result)  # -> []

    # 1.5 MF + tmVar + Custom regex + BOW

    result = regex_block(settings, ' asdf gpa-2 ::Tc1 asdf as')
    pprint(result)  # -> []

    # 1.6 * Additional details
    # These will get extracted regardless of whether the sentence has genomic information
    result = extra_info_block(custom_mut_extract, 'in the statement ced-3(n2888)')
    pprint(result)  # -> ['ced-3(n2888', 'Gene & Variant']]


if __name__ == "__main__":
    settings = setSettings()
    test_regexes(settings)
