import unittest

from genomicinfo.model_blocks.bow import BOWBlock


class TestBOWEntityExtractor(unittest.TestCase):

    def setUp(self) -> None:
        self.extractor = BOWBlock()

    def test_extract(self):
        # Bogus sentence for testing
        ext = self.extractor.extract('Upon reanalyzing the che-1(ot73) deletion we found the downstream substitution to be a D233G change, rather than an earlier frameshift, as previ- ously reported (C hang et al .')
        self.assertTrue(ext)
        