# Setup:  
- Installation  
```  
pip3 install wheel   
pip3 install -r requirements.txt   
```  
- You'll need to put your [wbtools](https://github.com/WormBase/wbtools) credentials and TextPresso token inside utils\all_config.cfg  
  
   
# Steps:     
Follow the notebooks serially.   
1. Preparing data.ipynb: Download and prepare data for NER.   
2. Hybrid extraction.ipynb: Train NER model, run hybrid mutation extraction on C. Elegans papers by listing their WBPaper ID.   
3. Validation.ipynb: Normalize mutations and genes, get matches using the [sequence file](ftp://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS281/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS281.protein.fa.gz) and cross-check final output if paper was already manually curated.  
  
  
# References:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6371299/  
https://academic.oup.com/bioinformatics/article/33/12/1852/2991428  
https://academic.oup.com/bioinformatics/article/23/14/1862/188647  
https://www.ncbi.nlm.nih.gov/research/bionlp/Tools/tmvar/  
https://arxiv.org/pdf/2101.07450.pdf  
https://academic.oup.com/bioinformatics/article/36/4/1234/5566506  