import unittest

from genomicinfo.pipeline import Pipeline


class TestPipelineBlock(unittest.TestCase):

    def setUp(self) -> None:
        folder_name = "nala" 
        allele_variations = ['n983', 'n1550', 'n1014', 'n983', 'md2230', 'md2243',\
            'e15', 'me4', 'me9', 'me15', 'me17', 'me18', 'me19', 'me23',\
            'me24', 'me41', 'me42', 'me44', 'me45', 'n2888',
            'ok367', 'adm-4', 'mn159', 'gk329', 'mn159']
        allele_designations = ['oa', 'mh', 'qz', 'tu', 'djr', 'tt', 'shy', 'qy', 'ww', 'dm',\
            'cac', 'ex', 'reb', 'ktn', 'pen', 'uv', 'zen', 'as', 'sky', 'yqt']
        gene_names = ['unc-93','mec-1','PPA41899', 'CJA19533', \
            'sup-10', 'sup-9', 'sup-18', 'Bm9466', 'CJA41250', \
            'CRE00014', 'Bma-rgr-1', 'CRE32070', 'ok265', 'met-1', \
            'CJA02939', 'OVOC9502', 'cr01.sctg4.wum.162.1', \
            'C26F1.1', 'ceh-39', 'acr-12', 'bar-1', 'clk-2']

        self.pipeline = Pipeline(use_ner=True, use_bow=True, use_regex=True)
        self.pipeline.add_details(folder_name=folder_name, allele_designations=allele_designations,\
            gene_names=gene_names, allele_variations=allele_variations)

    def test_extract(self):
        ext = self.pipeline.extract('Upon reanalyzing the che-1(ot73) sequence we found the molecular lesion to be a D233G change, rather than an earlier frameshift, as previ- ously reported (C hang et al .')
        self.assertTrue(ext[0][0] == 'D233G')
        ext = self.pipeline.extract("r380 (Glu-87- Lys), r342 (Arg-107-Cys), and the r381-r424-r445 muta- tional cluster (Pro-126-Leu) alter residues within a con- served, hydrophobic region that likely contributes to ATP- ase activity.")
        self.assertTrue(set([t[0] for t in ext]) >= set(['Arg-107-Cys', 'Arg-107-Cys)', 'Glu-87- Lys)', 'Pro-126-Leu', 'Pro-126-Leu)']))
        ext = self.pipeline.extract("The modified lysine is thought to participate in the ATPase reaction by interacting with the phosphates of bound ATP Taken together, these data suggest that the region from Lys-82 to Lys-128 contributes to the ATPase activity of the protein.")
        self.assertTrue(ext[0][0] == 'Lys-82 to Lys')
        ext = self.pipeline.extract("This conclusion is also supported by the adventitious finding of a point mutation, u753 (and its duplication by in vitro mutagenesis), P69L, that does not produce a mutant phenotype.")
        self.assertTrue(ext[0][0] == 'P69L')
        ext = self.pipeline.extract("1996) abolished the killing activity of CED-3 in vivo. The processing of the ( n2433) mutation, which substitutes a serine for glycine at codon 360 ( Shaham and Horvitz 1996a).")
        self.assertTrue(ext[0][0] == 'serine for glycine at codon 360')
        ext = self.pipeline.extract("The processing substitutes a ala for gly at position 560.")
        self.assertTrue(ext[0][0] == 'ala for gly at position 560')
        ext = self.pipeline.extract('Table 4. sup-10 mutations Allele Mutation Effect Isolation background Mutagen Loss-of-function mutations n1017 a TGG to TAG W19 Amber sup-10(n983) EMS n1008 a CAA to TAA Q99 Opal sup-10(n983) EMS n2297 b TGG to TGA W139 Ochre sup-9(n1550); sup-18(n1014) EMS n1626 a AAA to TAA K170 Opal sup-10(n983) Gamma n619 c GAG to AAG E88K unc-93(e1500) EMS n240 c GGG to AGG G323R unc-93(e1500) DES n1007 a agAT to aaAT Exon 3 Donor sup-10(n983) EMS n250 c TTgtto TTga Exon 5 Donor unc-93(e1500) Gamma n342 c TTgtto TTga Exon 5 Donor unc-93(e1500) NTG n3558 d 85 bp del Codon 16 Fs unc-93(e1500) UV-TMP n3564 d 99 bp del Codon 34 Fs unc-93(e1500) UV-TMP n183 c 53 bp del, Codon 128 Fs unc-93(e1500) SPO e2127 e Tc1 insertion Codon 103 Fs Wild type SPO n1468 a Tc1 insertion Codon 136 Fs sup-10(n983) SPO Gain-of-function mutations n983 a gf TGG to TAG W322 Amber Wild type EMS ThesequencesofbothDNAstrandsweredeterminedforeachmutant.Forsplicemutations,intronsequenceisrepresentedinlowercase,andexonsequenceinuppercase.EMS,Ethylmethanesulfonate;SPO,spontaneous;UV-TMP,ultraviolet and trimethylpsoralen; NTG, nitrosoguanidine; DES, diethyl sulfate.')
        self.assertTrue(set([t[0] for t in ext]) >= set(['85 bp del', '99 bp del', 'G323R', 'K170 Opal', 'E88K', 'Q99 Opal', 'W19 Amber', 'W139 Ochre', 'e1500) UV-TMP n183 c 53 bp del, Codon 128 Fs unc-93(e1500) SPO e2127 e Tc1 insertion', 'W322 Amber', '53 bp del']))
        ext = self.pipeline.extract('We found six other mutations that introduce stop codons directly in the bar-1 reading frame: de6 (Q25STOP), mu350 (Q143STOP), sy324 (Q304STOP), mu236 (Q320 STOP), de7 (Q320 STOP), and ep479 (W711STOP).')
        self.assertTrue(set([t[0] for t in ext]) >= set(['Q320 STOP', 'Q143STOP', 'Q25STOP', 'Q304STOP', 'W711STOP']))
        ext = self.pipeline.extract('The fact that muta- tions affecting either the ATP or the actin binding sites lead to similar assembly defects suggests that actin-myo- sin interactions are an important part of the assembly pro- cess and that functions of the myosin ATPase are needed for these interactions.')
        self.assertTrue(ext)
        ext = self.pipeline.extract('Electron microscope studies on the structure of natural and synthetic protein filaments from striated muscle. ships in the GTP binding domain of EF-Tu: mutation of Val20, the resi- due homologous to position 12 in ~21.')
        self.assertTrue(ext)
        ext = self.pipeline.extract('The amber termination codon is denoted by an asterisk. Resi- dues altered in mutants are circled; residues converted to termination codons in mutants are enclosed in squares (see Table 3).')
        self.assertTrue(ext)
        ext = self.pipeline.extract('The reported amino acid sequence of ceh-15 (Kamb et al., 1969) contains the substitutions Yl69T and S2OlY, pre- sumably resulting from artifacts introduced by PCR amplification.')
        self.assertTrue(ext)
        ext = self.pipeline.extract('Upon reanalyzing the che-1(ot73) deletion we found the downstream substitution to be a D233G change, rather than an earlier frameshift, as previ- ously reported (C hang et al .')
        self.assertTrue(ext)
