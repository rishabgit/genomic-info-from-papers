import unittest

from genomicinfo.model_blocks.ner import NERBlock


class TestBERTEntityExtractor(unittest.TestCase):

    def setUp(self) -> None:
        folder_name = "nala" 
        self.extractor = NERBlock(folder_name)

    def test_extract(self):
        ext = self.extractor.extract('The fact that muta- tions affecting either the ATP or the actin binding sites lead to similar assembly defects suggests that actin-myo- sin interactions are an important part of the assembly pro- cess and that functions of the myosin ATPase are needed for these interactions.')
        self.assertTrue(ext)
        ext = self.extractor.extract('Electron microscope studies on the structure of natural and synthetic protein filaments from striated muscle. ships in the GTP binding domain of EF-Tu: mutation of Val20, the resi- due homologous to position 12 in ~21.')
        self.assertTrue(ext)
        ext = self.extractor.extract('The amber termination codon is denoted by an asterisk. Resi- dues altered in mutants are circled; residues converted to termination codons in mutants are enclosed in squares (see Table 3).')
        self.assertTrue(ext)
        ext = self.extractor.extract('The reported amino acid sequence of ceh-15 (Kamb et al., 1969) contains the substitutions Yl69T and S2OlY, pre- sumably resulting from artifacts introduced by PCR amplification.')
        self.assertTrue(ext)
        