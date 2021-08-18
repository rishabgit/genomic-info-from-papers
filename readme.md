# Extract genomic data from papers. Made during GSoC 21 for [WormBase](https://wormbase.org/)
  
![Pipeline figure](https://github.com/rishabgit/genomic-info-from-papers/blob/main/utils/misc/fig.jpg?raw=true)
  
  
# Setup:  
1. Installation  
```  
pip3 install wheel   
pip3 install -r requirements.txt   
```  
2. You'll need to put your [wbtools](https://github.com/WormBase/wbtools) credentials and TextPresso token inside utils\all_config.cfg  
  
   
# Steps:     
Follow the notebooks serially.   
1. Preparing data.ipynb: Download and prepare data for NER.   
2. Hybrid extraction.ipynb: Train NER model, run hybrid mutation extraction on C. Elegans papers by listing their WBPaper ID.   
3. Validation.ipynb: Normalize mutations and genes, get matches using the [sequence file](ftp://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS281/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS281.protein.fa.gz) and cross-check final output if paper was already manually curated.  


# Results:
In 100 papers tested (93 were in the manually curated ground truth file), gene-mutation matches were found in 53 papers.   
Total 2433 matches were present in those 53 papers. And 977 matches were found using this developed pipeline.    
Cross-checking the matches with ground truth file -  
TP: 807, FP: 170  
Precision: 82.59%  
Not all *FP* are FP. After manual verification of the final output, some were noticed to be true positive which were originally missed during the manual curation.  



# Where Next  
These ideas, while interesting, were not possible during the two-month coding period. If worked on, they might improve recall by a huge margin.   
  
Majority of the extracted mutations from notebook 2 are being ignored (almost 73% of outputs during development). This is partly due to limitation in the mutation normalization block in notebook 3 and the subpar predictions from NER due to being trained on limited mutation data in natural language form.  
1. Active learning to improve the NER model: Annotating from scratch could be very expensive, so using the existing data and the regex outputs to kick start the process, and manually verifing a few subsets of annotations with help of curators would be useful.  
2. Better mutation normalization: Sounds very vague because it is. Current approach normalizes the easiest protein mutations. This will require a bit of research, however, improved NER outputs would certainly aid a lot.  
  
More additional data which will help curators in final verification.  
1. More data: This will help in easier and quicker verification, and can also be used to fill additional details in the ace file.  
2. Smarter matching: A more complex gene-mutation matching will help in filtering out even more false positives. Especially useful when curators come across patterns that were previously unnoticed during the verification.   
   
Faster and leaner.  
1. Faster extraction: Mutation extraction in notebook 2 takes a tremendously long time - almost >1 hour for a single paper. This is mainly due to the big number of regular expressions used and while they can't be switched off, there might be a smarter way of doing it.
2. Leaner storage: Some of the data structure choices made in this project are a little questionable. While they work fine, it can be better. 



  
  
  
# References:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6371299/  
https://academic.oup.com/bioinformatics/article/33/12/1852/2991428  
https://academic.oup.com/bioinformatics/article/23/14/1862/188647  
https://www.ncbi.nlm.nih.gov/research/bionlp/Tools/tmvar/  
https://arxiv.org/pdf/2101.07450.pdf  
https://academic.oup.com/bioinformatics/article/36/4/1234/5566506  