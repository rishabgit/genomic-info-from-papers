import unittest

from genomicinfo.entity_extraction.regex.mutfinder import MutationFinderRegexEntityExtractor


class TestMutFinderRegexEntityExtractor(unittest.TestCase):

    def setUp(self) -> None:
        self.extractor = MutationFinderRegexEntityExtractor()

    def test_extract(self):
        ext = self.extractor.extract('Upon reanalyzing the che-1(ot73) sequence we found the molecular lesion to be a D233G change, rather than an earlier frameshift, as previ- ously reported (C hang et al .')
        self.assertTrue(ext)
        ext = self.extractor.extract("r380 (Glu-87- Lys), r342 (Arg-107-Cys), and the r381-r424-r445 muta- tional cluster (Pro-126-Leu) alter residues within a con- served, hydrophobic region that likely contributes to ATP- ase activity.")
        self.assertTrue(ext)
        ext = self.extractor.extract("The modified lysine is thought to participate in the ATPase reaction by interacting with the phosphates of bound ATP Taken together, these data suggest that the region from Lys-82 to Lys-128 contributes to the ATPase activity of the protein.")
        self.assertTrue(ext)
        ext = self.extractor.extract("This conclusion is also supported by the adventitious finding of a point mutation, u753 (and its duplication by in vitro mutagenesis), P69L, that does not produce a mutant phenotype.")
        self.assertTrue(ext)
        ext = self.extractor.extract("1996) abolished the killing activity of CED-3 in vivo. The processing of the ( n2433) mutation, which substitutes a serine for glycine at codon 360 ( Shaham and Horvitz 1996a).")
        self.assertTrue(ext)
        ext = self.extractor.extract("The processing substitutes a ala for gly at position 560.")
        self.assertTrue(ext)
