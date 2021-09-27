import unittest

from genomicinfo.entity_extraction.regex.tmvar import TMVarRegexEntityExtractor


class TestTMVarRegexEntityExtractor(unittest.TestCase):

    def setUp(self) -> None:
        self.extractor = TMVarRegexEntityExtractor()

    def test_extract(self):
        ext = self.extractor.extract('Upon reanalyzing the che-1(ot73) sequence we found the molecular lesion to be a D233G change, rather than an earlier frameshift, as previ- ously reported (C hang et al .')
        self.assertTrue(ext)
