import unittest

from genomicinfo.entity_extraction.regex.wb_custom import \
    WBCustomRegexEntityExtractor, \
    WBGeneVarRegexEntityExtractor, \
    WBGeneRegexGenomicLocExtractor


class TestMutFinderRegexEntityExtractor(unittest.TestCase):

    def setUp(self) -> None:
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

        self.custom_extractor = WBCustomRegexEntityExtractor(allele_variations, allele_designations, gene_names)
        self.custom_locus_only_extractor = WBCustomRegexEntityExtractor(allele_variations, allele_designations, gene_names, locus_only=True)
        self.var_extractor = WBGeneVarRegexEntityExtractor(allele_variations, allele_designations, gene_names)
        self.gene_extractor = WBGeneRegexGenomicLocExtractor(allele_variations, allele_designations, gene_names)

    def test_extract(self):
        ext = self.custom_extractor.extract('Table 4. sup-10 mutations Allele Mutation Effect Isolation background Mutagen Loss-of-function mutations n1017 a TGG to TAG W19 Amber sup-10(n983) EMS n1008 a CAA to TAA Q99 Opal sup-10(n983) EMS n2297 b TGG to TGA W139 Ochre sup-9(n1550); sup-18(n1014) EMS n1626 a AAA to TAA K170 Opal sup-10(n983) Gamma n619 c GAG to AAG E88K unc-93(e1500) EMS n240 c GGG to AGG G323R unc-93(e1500) DES n1007 a agAT to aaAT Exon 3 Donor sup-10(n983) EMS n250 c TTgtto TTga Exon 5 Donor unc-93(e1500) Gamma n342 c TTgtto TTga Exon 5 Donor unc-93(e1500) NTG n3558 d 85 bp del Codon 16 Fs unc-93(e1500) UV-TMP n3564 d 99 bp del Codon 34 Fs unc-93(e1500) UV-TMP n183 c 53 bp del, Codon 128 Fs unc-93(e1500) SPO e2127 e Tc1 insertion Codon 103 Fs Wild type SPO n1468 a Tc1 insertion Codon 136 Fs sup-10(n983) SPO Gain-of-function mutations n983 a gf TGG to TAG W322 Amber Wild type EMS ThesequencesofbothDNAstrandsweredeterminedforeachmutant.Forsplicemutations,intronsequenceisrepresentedinlowercase,andexonsequenceinuppercase.EMS,Ethylmethanesulfonate;SPO,spontaneous;UV-TMP,ultraviolet and trimethylpsoralen; NTG, nitrosoguanidine; DES, diethyl sulfate.')
        self.assertTrue(len(ext) == 6)
        ext = self.custom_extractor.extract('The mutation gk329 is a 1047- bp deletion in the gene ceh-39.')
        self.assertTrue(ext[0][0] == 'gk329 is a 1047- bp deletion')
        ext = self.custom_extractor.extract('ok367 has a 1,368-bp deletion and is a null allele of acr-12.')
        self.assertTrue(ext[0][0] == 'ok367 has a 1,368-bp deletion')

        ext = self.custom_locus_only_extractor.extract('We found six other mutations that introduce stop codons directly in the bar-1 reading frame: de6 (Q25STOP), mu350 (Q143STOP), sy324 (Q304STOP), mu236 (Q320 STOP), de7 (Q320 STOP), and ep479 (W711STOP).')
        self.assertTrue(len(ext) == 6)
        ext = self.custom_locus_only_extractor.extract('Table 4. sup-10 mutations Allele Mutation Effect Isolation background Mutagen Loss-of-function mutations n1017 a TGG to TAG W19 Amber sup-10(n983) EMS n1008 a CAA to TAA Q99 Opal sup-10(n983) EMS n2297 b TGG to TGA W139 Ochre sup-9(n1550); sup-18(n1014) EMS n1626 a AAA to TAA K170 Opal sup-10(n983) Gamma n619 c GAG to AAG E88K unc-93(e1500) EMS n240 c GGG to AGG G323R unc-93(e1500) DES n1007 a agAT to aaAT Exon 3 Donor sup-10(n983) EMS n250 c TTgtto TTga Exon 5 Donor unc-93(e1500) Gamma n342 c TTgtto TTga Exon 5 Donor unc-93(e1500) NTG n3558 d 85 bp del Codon 16 Fs unc-93(e1500) UV-TMP n3564 d 99 bp del Codon 34 Fs unc-93(e1500) UV-TMP n183 c 53 bp del, Codon 128 Fs unc-93(e1500) SPO e2127 e Tc1 insertion Codon 103 Fs Wild type SPO n1468 a Tc1 insertion Codon 136 Fs sup-10(n983) SPO Gain-of-function mutations n983 a gf TGG to TAG W322 Amber Wild type EMS ThesequencesofbothDNAstrandsweredeterminedforeachmutant.Forsplicemutations,intronsequenceisrepresentedinlowercase,andexonsequenceinuppercase.EMS,Ethylmethanesulfonate;SPO,spontaneous;UV-TMP,ultraviolet and trimethylpsoralen; NTG, nitrosoguanidine; DES, diethyl sulfate.')
        self.assertTrue(len(ext) == 5)
        ext = self.custom_locus_only_extractor.extract('C 123 Umber')
        self.assertTrue(ext[0][0] == 'C 123 Umber')

        ext = self.var_extractor.extract('Table 4. sup-10 mutations Allele Mutation Effect Isolation background Mutagen Loss-of-function mutations n1017 a TGG to TAG W19 Amber sup-10(n983) EMS n1008 a CAA to TAA Q99 Opal sup-10(n983) EMS n2297 b TGG to TGA W139 Ochre sup-9(n1550); sup-18(n1014) EMS n1626 a AAA to TAA K170 Opal sup-10(n983) Gamma n619 c GAG to AAG E88K unc-93(e1500) EMS n240 c GGG to AGG G323R unc-93(e1500) DES n1007 a agAT to aaAT Exon 3 Donor sup-10(n983) EMS n250 c TTgtto TTga Exon 5 Donor unc-93(e1500) Gamma n342 c TTgtto TTga Exon 5 Donor unc-93(e1500) NTG n3558 d 85 bp del Codon 16 Fs unc-93(e1500) UV-TMP n3564 d 99 bp del Codon 34 Fs unc-93(e1500) UV-TMP n183 c 53 bp del, Codon 128 Fs unc-93(e1500) SPO e2127 e Tc1 insertion Codon 103 Fs Wild type SPO n1468 a Tc1 insertion Codon 136 Fs sup-10(n983) SPO Gain-of-function mutations n983 a gf TGG to TAG W322 Amber Wild type EMS ThesequencesofbothDNAstrandsweredeterminedforeachmutant.Forsplicemutations,intronsequenceisrepresentedinlowercase,andexonsequenceinuppercase.EMS,Ethylmethanesulfonate;SPO,spontaneous;UV-TMP,ultraviolet and trimethylpsoralen; NTG, nitrosoguanidine; DES, diethyl sulfate.')
        self.assertTrue(len(ext) == 4)
        ext = self.var_extractor.extract('adm-4(ok265) is a 847 bp deletion with a extra C insertion.')
        self.assertTrue(ext[0][0] == "adm-4(ok265)")
        ext = self.var_extractor.extract('Molecular analysis of the break- points of these unc-58 deletions revealed 7- to 47-bp insertions at breakpoints for 3/6 deletions from the clk-2(mn159) background and 1- to 199-bp insertions at breakpoints for 4/11 deletions from the mrt-2 back- ground.')
        self.assertTrue(ext[0][0] == 'clk-2(mn159)')

        ext = self.gene_extractor.extract('Table 4. sup-10 mutations Allele Mutation Effect Isolation background Mutagen Loss-of-function mutations n1017 a TGG to TAG W19 Amber sup-10(n983) EMS n1008 a CAA to TAA Q99 Opal sup-10(n983) EMS n2297 b TGG to TGA W139 Ochre sup-9(n1550); sup-18(n1014) EMS n1626 a AAA to TAA K170 Opal sup-10(n983) Gamma n619 c GAG to AAG E88K unc-93(e1500) EMS n240 c GGG to AGG G323R unc-93(e1500) DES n1007 a agAT to aaAT Exon 3 Donor sup-10(n983) EMS n250 c TTgtto TTga Exon 5 Donor unc-93(e1500) Gamma n342 c TTgtto TTga Exon 5 Donor unc-93(e1500) NTG n3558 d 85 bp del Codon 16 Fs unc-93(e1500) UV-TMP n3564 d 99 bp del Codon 34 Fs unc-93(e1500) UV-TMP n183 c 53 bp del, Codon 128 Fs unc-93(e1500) SPO e2127 e Tc1 insertion Codon 103 Fs Wild type SPO n1468 a Tc1 insertion Codon 136 Fs sup-10(n983) SPO Gain-of-function mutations n983 a gf TGG to TAG W322 Amber Wild type EMS ThesequencesofbothDNAstrandsweredeterminedforeachmutant.Forsplicemutations,intronsequenceisrepresentedinlowercase,andexonsequenceinuppercase.EMS,Ethylmethanesulfonate;SPO,spontaneous;UV-TMP,ultraviolet and trimethylpsoralen; NTG, nitrosoguanidine; DES, diethyl sulfate.')
        self.assertTrue(('sup-10', 'Just Gene') in ext)
        ext = self.gene_extractor.extract('The mutation gk329 is a 1047- bp deletion in the gene ceh-39.')
        self.assertTrue(ext[0][0] == 'ceh-39')
        ext = self.gene_extractor.extract('We found six other mutations that introduce stop codons directly in the bar-1 reading frame: de6 (Q25STOP), mu350 (Q143STOP), sy324 (Q304STOP), mu236 (Q320 STOP), de7 (Q320 STOP), and ep479 (W711STOP).')
        self.assertTrue(ext[0][0] == 'bar-1')