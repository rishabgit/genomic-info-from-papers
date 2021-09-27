import unittest

from genomicinfo.model_blocks.regex.mutfinder import MutationFinderRegexEntityExtractor


class TestMutFinderRegexEntityExtractor(unittest.TestCase):

    def setUp(self) -> None:
        self.extractor = MutationFinderRegexEntityExtractor()

    def test_extract(self):
        ext = self.extractor.extract('Upon reanalyzing the che-1(ot73) sequence we found the molecular lesion to be a D233G change, rather than an earlier frameshift, as previ- ously reported (C hang et al .')
        self.assertTrue(ext)
