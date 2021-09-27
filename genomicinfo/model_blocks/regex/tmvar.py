import os

from genomicinfo.model_blocks.regex.regex_extractor import RegexEntityExtractor


class TMVarRegexEntityExtractor(RegexEntityExtractor):
    """Extract mutations using tmVar regular expressions."""

    def __init__(self):
        super().__init__(os.path.join(os.path.dirname(__file__), "regex_files", "tmvar"))

