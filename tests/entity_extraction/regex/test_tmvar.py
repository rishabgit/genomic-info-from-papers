import unittest

from genomicinfo.entity_extraction.regex.tmvar import TMVarRegexEntityExtractor


class TestTMVarRegexEntityExtractor(unittest.TestCase):

    def setUp(self) -> None:
        self.extractor = TMVarRegexEntityExtractor()

    def test_extract(self):
        ext = self.extractor.extract("Upon reanalyzing the che-1(ot73) sequence we found the molecular lesion to be a D233G change, rather than an earlier frameshift, as previ- ously reported (C hang et al .")
        self.assertTrue(ext[0][0] == 'D233G')
        ext = self.extractor.extract("The sequences of DNAs amplified from wild-type and r433. r433 contains a G to A transition at nucleotide 2518, a position coincident with the site of RNAase cleavage.")
        self.assertTrue(ext[0][0] == 'G to A transition at nucleotide 2518')
        ext = self.extractor.extract("The r387-f424-f445 mutational cluster (Pro-126- Leu) is very near Lys-128, an invariant residue that is tri- methylated in rabbit myosin.")
        self.assertTrue(ext[0][0] == 'Pro-126- Leu)')
        ext = self.extractor.extract("Mutation of lysine-48 to arginine in the yeast RAD3 protein abolishes its ATPase and DNA helicase activities but not the ability to bind ATP EMBO J. i: 3263-3269.")
        self.assertTrue(ext[0][0] == 'lysine-48 to arginine')
        ext = self.extractor.extract("Two of the four ced-3 alleles with semidominant ef- fused to a ced-3 cDNA encoding an active-site cysteine- 358-to-alanine mutant protein (see Figure 4 legend for fects, ced-3(n2871) and ced-3(n2433), alter the arginine (R359) and glycine (G360) residues, respectively, in heat-shock protocol).")
        self.assertTrue(ext[0][0] == 'cysteine- 358-to-alanine')
