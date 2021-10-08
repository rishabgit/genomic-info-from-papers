import unittest

from genomicinfo.entity_extraction.bow.bowdictionary import BOWdictionary


class TestBOWdictionary(unittest.TestCase):

    def setUp(self) -> None:
        self.extractor = BOWdictionary()

    def test_extract(self):
        # Bogus sentence for testing
        ext = self.extractor.extract('Upon reanalyzing the che-1(ot73) deletion we found the downstream substitution to be a D233G change, rather than an earlier frameshift, as previ- ously reported (C hang et al .')
        self.assertTrue(ext)
        