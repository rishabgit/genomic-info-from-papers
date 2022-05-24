#!/usr/bin/env python
# coding: utf-8

# # 1 Compiling notebook 2 outputs

# In[1]:

import configparser
import glob
import json
import math
import numpy as np
import pandas as pd
import re

from utils.misc.regex_block import MutationFinder, TmVar, CustomWBregex, normalize_mutations



with open("data/model_output/processed/temp_paper_mut_count.json", "w") as outfile: 
    json.dump(paper_mut_count, outfile)

print('All', ner_count, 'NER data rows were ignored. Only', regex_count, 'regex data rows were used.')


# saving things
data = pd.DataFrame(data[:], columns=['WBPaper ID', 'Method', 'Genes', '*Gene-Variant combo ', 'Mutations', 'Normalized Mutations', 'Sentence'])
data.to_csv("data/model_output/processed/snippets_2.csv", index=False, encoding='utf-8')


# # 3 Normalizing common gene name to its WormBase ID
# And getting the gene and mutation frequency in a paper.
# In[14]:

data = pd.read_csv("data/model_output/processed/snippets_2.csv")
data = data.to_numpy() # 'WBPaper ID', 'Method', 'Genes', '*Gene-Variant combo ', 'Mutations', 'Normalized Mutations', 'Sentence'


with open("data/model_output/processed/temp_paper_wbgene_count.json", "w") as outfile: 
    json.dump(paper_wbgene_count, outfile)



# Checking if any detected gene was NOT in the WB gene dictionary

# In[18]:

data = np.array(data)
data[len(data[:,2]) != len(data[:,3])] 


# above cell takes a while to complete, so saving the data temporarily
data = pd.DataFrame(data[:], columns=['WBPaper ID', 'Method', 'Genes', 'WBGenes', '*Gene-Variant combo ', 'Mutations', 'Normalized Mutations', 'Sentence'])
data.to_csv("data/model_output/processed/snippets_3.csv", index=False, encoding='utf-8')
data = None


# # 5 Validation
# Finding the gene and mutation matches using the transcripts in c_elegans.PRJNA13758.WS281.protein.fa   
# Get the file here - ftp://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS281/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS281.protein.fa.gz


data = pd.read_csv("data/model_output/processed/snippets_3.csv")
data = data.to_numpy() # 'WBPaper ID', 'Method', 'Genes', 'WBGenes', '*Gene-Variant combo ', 'Mutations', 'Normalized Mutations', 'Sentence'


proteinfa = Path('data/gsoc/proteinfa/c_elegans.PRJNA13758.WS281.protein.fa').read_text().split('>')[1:]
wb_gene_and_prot = dict()  # {wbgene: [transcript, protein]}

for row in proteinfa:
    wbgene = re.findall("WBGene[0-9]+", row)[0]
    protein = "".join(re.findall("\n.*", row)).replace('\n', '')
    transcript = row.split(' ')[0]
    if wbgene not in wb_gene_and_prot.keys():
        wb_gene_and_prot[wbgene] = []
    wb_gene_and_prot[wbgene].append([transcript, protein])

len(wb_gene_and_prot)


# #### Create a pair of gene and mutation only when BOTH are present in same sentence.

# In[24]:


paper_raw_info_compiled = []
# 'WBPaper ID', 'Method', 'Genes', 'WBGenes', '*Gene-Variant combo ', 'Mutations', 'Normalized Mutations', 'Sentence'
for row in data:
    ppr_id = row[0]
    norm_muts = row[-2]
    wbgenes = row[3]
    sentence = row[-1]
    gene_var = row[4]
        
    # filtering out nan values
    if type(norm_muts) != float and type(wbgenes) != float:    
        norm_muts = norm_muts[1:-1].split("', '")
        wbgenes = wbgenes[1:-1].split("', '")
        for m in norm_muts:
            for g in wbgenes:
                if len(m) and len(g):
                    paper_raw_info_compiled.append([ppr_id, g, m, sentence, gene_var])


# In[25]:


matches = [] 
final_sheet = [] # ppr_id, gene, transcript

for info_from_ppr in paper_raw_info_compiled:
    ppr_id = info_from_ppr[0]
    gene = info_from_ppr[1]
    mut = info_from_ppr[2]
    sent = info_from_ppr[3]
    gene_var = info_from_ppr[4]
    if not len(mut):
        continue
    if gene not in wb_gene_and_prot.keys():
        continue
    for row in wb_gene_and_prot[gene]:
        transcript, protein_string = row
        wt_res = mut[0]
        pos = int(''.join(n for n in mut if n.isdigit()))
        mut_res = mut[-1]
        try:
            if protein_string[pos-1] == wt_res:
                matches.append([ppr_id, gene, mut, gene_var, transcript, sent])
        except IndexError:
            pass
    
for r in matches:
    p = r[0]
    p, wbg, mut, gene_var, transcript, sent = r
    # Adding gene common names column, again
    # Current code doesn't keep any link between the WB gene name and the common name
    g_common_name = all_wb_genes[wbg]
    g_common_name = ', '.join(g_common_name)
    final_sheet.append([p, wbg, g_common_name, mut, gene_var, transcript, sent])


# In[26]:


len(final_sheet)


# #### Getting metadata on genes and mutations, and adding warnings column 

# In[27]:


with open("data/model_output/processed/temp_paper_wbgene_count.json", "r") as f: 
    paper_wbgene_count = json.loads(f.read())
with open("data/model_output/processed/temp_paper_mut_count.json", "r") as f: 
    paper_mut_count = json.loads(f.read())


# In[28]:


final_sheet = np.array(final_sheet)
updated_sheet = []
for i, row in enumerate(final_sheet):
    warnings = []
    paper_id = row[0]
    wbgene = row[1]
    mut = row[3]
    sentence = row[-1]
    for ppr_mut, count in paper_mut_count[paper_id].items():
        if mut == ppr_mut and count == 1:
            warnings.append(f'{mut} mentioned only once in entire paper')
            break
            
    rows_with_same_mut = final_sheet[np.logical_and(final_sheet[:, 0] == paper_id, final_sheet[:,3] == mut)]
    same_mut_all_genes = list(set(rows_with_same_mut[:, 1]))
    # If the same variant is found in two different genes in the same paper - WARN! 
    # It is more likely to belong to the gene it is most frequently encountered
    if len(same_mut_all_genes) > 1:
        temp_warn_store = f'{mut} was paired with other genes too:'
        for ppr_gene, count in paper_wbgene_count[paper_id].items():
            if ppr_gene in same_mut_all_genes:
                temp_warn_store += (f' {ppr_gene} (seen {count} times),')
        warnings.append(temp_warn_store)
    
    cut_mut = re.sub("([A-Z])([0-9]+)([A-Za-z]+)", r'\1\2', mut)
    remaining_mut = mut.replace(cut_mut, "")
    same_cut_muts = [i for i,m in enumerate(final_sheet[:,3]) if (m[:len(cut_mut)] == cut_mut and m[len(cut_mut):] != remaining_mut)]
    if same_cut_muts:
        temp_warn_store = f'{mut} similar to:'
        for temp_i in same_cut_muts:
            temp_warn_store += (f' {final_sheet[:,3][temp_i]} (line {temp_i}),')
        warnings.append(temp_warn_store)
        
    all_muts_in_sentence = data[np.logical_and(data[:, 0] == paper_id, data[:,-1] == sentence)][:,-2]
    all_muts_in_sentence = all_muts_in_sentence[0][1:-1].split("', '")
    
    all_matched_muts_in_sentence = final_sheet[np.logical_and(final_sheet[:, 0] == paper_id, final_sheet[:,-1] == sentence)][:,3]
    all_matched_muts_in_sentence = list(set(all_matched_muts_in_sentence))
    
    unmatched_muts_in_sentence = [m for m in all_muts_in_sentence if m not in all_matched_muts_in_sentence]
    if len(unmatched_muts_in_sentence) >= 2:
        temp_warn_store = f'Sentence has multiple mutations which did not match:'
        for m in unmatched_muts_in_sentence:
            temp_warn_store += (f' {m},')
        warnings.append(temp_warn_store)
    
    all_genes_with_this_mut = final_sheet[np.logical_and(final_sheet[:, 0] == paper_id, final_sheet[:, 3] == mut)][:, 1]
    all_genes_with_this_mut = list(set(all_genes_with_this_mut))
    if len(all_genes_with_this_mut) > 3:
        temp_warn_store = f'{mut} was matched with {len(all_genes_with_this_mut)} genes:'
        for g in all_genes_with_this_mut:
            temp_warn_store += (f' {g},')
        warnings.append(temp_warn_store)
            
    if warnings:
        warnings = " || ".join(warnings)
    else:
        warnings = ""
    updated_sheet.append(np.insert(row, -1, warnings).tolist())


# In[29]:


# saving things
updated_sheet = pd.DataFrame(updated_sheet[:], columns=['WBPaper ID', 'WBGene', 'Gene', 'Mutation', 'Gene-Var combo', 'Transcript', 'Warnings', 'Sentence'])
updated_sheet.to_csv("data/model_output/processed/snippets_4.csv", index=False, encoding='utf-8')
updated_sheet = None


# # 6 Additional details

# ### 6.1 Strains

# In[30]:


data = pd.read_csv("data/model_output/processed/snippets_4.csv").to_numpy() 


# In[31]:


strains = Path('data/gsoc/Strains.txt').read_text().split('\n')
strains = [r.split('\t') for r in strains][:-1]
all_wb_strains = dict()

for row in strains:
    if row[0] not in all_wb_strains.keys():
        all_wb_strains[row[0]] = []
    for strain in row[1:]: 
        if len(strain) and strain.lower() not in all_wb_strains[row[0]]: 
            all_wb_strains[row[0]].append(strain.lower())
            
strains = [s for row in strains for s in row[1:] if len(s) and not s.isdigit()]


# In[32]:


OPENING_CLOSING_REGEXES = [r'(?:^|[^0-9A-Za-z])(', r')(?:^|[^0-9A-Za-z])']
all_strain = OPENING_CLOSING_REGEXES[0] + '|'.join(strains) + OPENING_CLOSING_REGEXES[1]
all_strain = [re.compile(r,re.IGNORECASE) for r in [all_strain]]

# 'WBPaper ID', 'WBGene', 'Gene', 'Mutation', 'Gene-Var combo', 'Transcript', 'Warnings', 'Sentence'
updated_data = []
total = len(data)
print('Total sentences: {}, processed count: '.format(total), end=' ')
for i, sent in enumerate(data[:, -1]):
    if (i+1) % 100 == 0: print(f"{i+1}", end = " ")
    paper_strains = []
    for regex in all_strain:      
        for m in regex.finditer(sent):
            span = (m.start(0), m.end(0))    
            raw = (sent[span[0]:span[1]]).strip()
            raw = raw[1:] if not raw[0].isalnum() else raw
            raw = raw[:-1] if not raw[-1].isalnum() else raw
            if len(raw.strip()) > 1 and not raw.strip().isdigit(): paper_strains.append(raw.strip())
    if paper_strains:
        paper_strains  = list(set(paper_strains))
        col_wbid = []
        for strain in paper_strains:
            for key, value in all_wb_strains.items():
                if strain.lower() in value:
                    col_wbid.append(key)
                    break
        paper_strains = "'" + "', '".join(paper_strains) + "'"
        if col_wbid:
            col_wbid = list(set(col_wbid))
            col_wbid = ", ".join(col_wbid)
        else: 
            col_wbid = ''
            # lazy way to deal with bad snippets due to special characters in the Strains.txt file
            # which are caught in regex
            paper_strains = '' 
    else:
        paper_strains = ''
        col_wbid = ''
    updated_data.append([data[i,0], data[i,1], data[i,2], col_wbid, paper_strains, data[i,3], data[i,-4], data[i,-3], data[i,-2], data[i,-1]])
data = np.array(updated_data) # 'WBPaper ID', 'WBGene', 'Gene', 'WBStrain', 'Strains', 'Mutation', 'Gene-Var combo', 'Transcript', 'Warnings', 'Sentence'
updated_data = None


# ### 6.2 Variants

# In[33]:


OPENING_CLOSING_REGEXES = [r'(?:^|[^0-9A-Za-z])(', r')(?:^|[^0-9A-Za-z])']

# the allele regex and db idea was stolen from wbtools
allele_designations = np.load('data/gsoc/wbtools/wb_allele_designations.npy').astype('U6')
alleles_variations = np.load('data/gsoc/wbtools/wb_alleles_variations.npy').astype('U6')
DB_VAR_REGEX = r'({designations}|m|p|ts|gf|lf|d|sd|am|cs)([0-9]+)'
var_regex_1 = OPENING_CLOSING_REGEXES[0] + DB_VAR_REGEX.format(designations="|".join(allele_designations)) + OPENING_CLOSING_REGEXES[1]
all_var = OPENING_CLOSING_REGEXES[0] + '|'.join(alleles_variations) + '|' + var_regex_1 + OPENING_CLOSING_REGEXES[1]
all_var = [re.compile(r,re.IGNORECASE) for r in [all_var]]

# 'WBPaper ID', 'WBGene', 'Gene', 'WBStrain', 'Strains', 'Mutation', 'Transcript', 'Warnings', 'Sentence'
updated_data = []
total = len(data)
print('Total sentences: {}, processed count: '.format(total), end=' ')
for i, sent in enumerate(data[:, -1]):
    if (i+1) % 100 == 0: print(f"{i+1}", end = " ")
    variants = []
    for regex in all_var:      
        for m in regex.finditer(sent):
            span = (m.start(0), m.end(0))    
            raw = (sent[span[0]:span[1]]).strip()
            raw = raw[1:] if not raw[0].isalnum() else raw
            raw = raw[:-1] if not raw[-1].isalnum() else raw
            if len(raw.strip()) > 1: variants.append(raw.strip())
    if variants:
        variants  = list(set(variants))
        variants = "'" + "', '".join(variants) + "'"
    else:
        variants = ''
    updated_data.append([data[i,0], data[i,1], data[i,2], data[i,3], data[i,4], variants, data[i,-5], data[i,-4], data[i,-3], data[i,-2], data[i,-1]])
    
data = np.array(updated_data) # 'WBPaper ID', 'WBGene', 'Gene', 'WBStrain', 'Strains', 'Variants', 'Mutation', 'Gene-Var combo', 'Transcript', 'Warnings', 'Sentence'
updated_data = None


# ### 6.3 Variation type
# Extraction rate would be very low as most snippets from notebook 2 are discarded due to limitation in the mutation normalization block above.

# In[34]:


Variation_type = pd.read_csv("data/gsoc/Variation_type.csv").to_numpy() 
Variation_type = [t.replace("_", " ") for t in Variation_type[:,2] if type(t)!=float]


# In[35]:


updated_sheet = []
# 'WBPaper ID', 'WBGene', 'Gene', 'WBStrain', 'Strains', 'Variants', 'Mutation', 'Gene-Var combo', 'Transcript', 'Warnings', 'Sentence'
for i, row in enumerate(data):
    sent = row[-1]
    col_var_type = []
    for sub in Variation_type:
        if re.search(sub, sent, re.IGNORECASE):
            col_var_type.append(sub)
    if col_var_type:
        col_var_type = list(set(col_var_type))
        col_var_type = ", ".join(col_var_type)
    else:
        col_var_type = ''
    updated_sheet.append(np.insert(row, -3, col_var_type).tolist())


# In[36]:


data = np.array(updated_sheet)
updated_sheet = None


# ### 6.3 Functional effect & Generation method 
# These type of data were in a few subset of papers tested during dev - expect these columns to be mostly empty.

# In[37]:


functional_effect = ['function uncertain', 'transcript function', 'translational product function',                  'decreased transcript level', 'increased transcript level', 'decreased transcript stability',                  'gain of function', 'dominant negative', 'dominant negativ', 'antimorphic',                  'hypermorphic', 'neomorphic', 'conditional activity', 'hypomorphic', 'amorphic',                  'repressible', 'misexpressed']


# In[38]:


common_gen_methods = ['CRISPR', 'ENU', 'EMS']


# In[39]:


updated_sheet = []
# 'WBPaper ID', 'WBGene', 'Gene', 'WBStrain', 'Strains', 'Variants', 'Mutation', 'Gene-Var combo', 'Variation type', 'Transcript', 'Warnings', 'Sentence'
for i, row in enumerate(data):
    sent = row[-1]
    col_functional_effect = []
    col_gen_method = [] 
    for sub in functional_effect:
        if re.search(sub, sent, re.IGNORECASE):
            col_functional_effect.append(sub)
            
    for sub in common_gen_methods:
        if re.search(sub, sent):
            col_gen_method.append(sub)
            
    if col_functional_effect:
        col_functional_effect = list(set(col_functional_effect))
        col_functional_effect = ", ".join(col_functional_effect)
    else:
        col_functional_effect = ''
        
    if col_gen_method:
        col_gen_method = list(set(col_gen_method))
        col_gen_method = ", ".join(col_gen_method)
    else:
        col_gen_method = ''
    row = np.insert(row, -3, col_functional_effect)
    row = np.insert(row, -3, col_gen_method)
    updated_sheet.append(row.tolist())
data = np.array(updated_sheet)
updated_sheet = None


# In[40]:


# saving things
updated_sheet = pd.DataFrame(data[:], columns=['WBPaper ID', 'WBGene', 'Gene', 'WBStrain', 'Strains',                                                'Variants', 'Mutation', 'Gene-Var combo', 'Variation type', 'Functional effect',                                                'Generation method', 'Transcript', 'Warnings', 'Sentence'])
updated_sheet.to_csv("data/model_output/processed/final.csv", index=False, encoding='utf-8')
updated_sheet = None


# # 7 Verification
# Finding precision by cross-checking with the manually curated data.

# In[41]:


data = pd.read_csv("data/model_output/processed/final.csv")
data = data.to_numpy() 
paper_ids_processed = np.unique(data[:,0])
paper_ids_processed = np.sort(paper_ids_processed)

temp = pd.read_csv("data/model_output/processed/snippets_1.csv")
temp = temp.to_numpy()
total_paper_ids_processed = np.unique(temp[:,0])
temp = None


# In[42]:


print('Total count of papers processed:', len(total_paper_ids_processed))
print('Count of papers:', len(paper_ids_processed))


# ### 7.1 Original ground truth

# In[43]:


ground_truth = Path('data/gsoc/variantsDB.txt').read_text().split('\n')
ground_truth = [r.split('\t') for r in ground_truth][:-1]
ground_truth = np.array(ground_truth, dtype=object)


# In[44]:


# Checking if any processed paper is not in the ground truth file
for id in total_paper_ids_processed:
    if id not in ground_truth[:,0]:
        print(id, end = ' ')
        if id in paper_ids_processed:
            print(' false positive')


# In[45]:


tp_col = []
for row in data:
    paper_id = row[0]
    gene = row[1]
    mutation = row[6]
    mutation = mutation.upper()
    transcript = row[-3]
    bool_found = False
    for label in ground_truth[ground_truth[:,0] == paper_id]:
        label[-2] = label[-2].upper()
        if transcript == label[-1] and mutation == label[-2]:
            bool_found = True
            # continue bc we're storing all the labels from a paper
            continue
    if bool_found:
        tp_col.append('True Positive')
    else:
        tp_col.append('False Positive')


# In[46]:


tp_col.count('True Positive'), tp_col.count('False Positive')


# In[47]:


print('Precision ',tp_col.count('True Positive')*100/(tp_col.count('True Positive') + tp_col.count('False Positive')), '%')


# ### 7.2 Curated ground truth  
# Updated after manually looking at the "false positives" which aren't false positives and adding them in the ground truth file.   

# In[48]:


ground_truth = Path('data/gsoc/variantsDB_curated.txt').read_text().split('\n')
ground_truth = [r.split('\t') for r in ground_truth][:-1]
ground_truth = np.array(ground_truth, dtype=object)


# In[49]:


# Checking if any processed paper is not in the ground truth file
for id in total_paper_ids_processed:
    if id not in ground_truth[:,0]:
        print(id, end = ' ')
        if id in paper_ids_processed:
            print(' false positive')


# In[50]:


tp_col = []
for row in data:
    paper_id = row[0]
    gene = row[1]
    mutation = row[6]
    mutation = mutation.upper()
    transcript = row[-3]
    bool_found = False
    for label in ground_truth[ground_truth[:,0] == paper_id]:
        label[-2] = label[-2].upper()
        if transcript == label[-1] and mutation == label[-2]:
            bool_found = True
            # continue bc we're storing all the labels from a paper
            continue
    if bool_found:
        tp_col.append('True Positive')
    else:
        tp_col.append('False Positive')


# In[51]:


tp_col.count('True Positive'), tp_col.count('False Positive')


# In[52]:


print('Precision ',tp_col.count('True Positive')*100/(tp_col.count('True Positive') + tp_col.count('False Positive')), '%')


# In[53]:


tp_col = np.array(tp_col).T.reshape(-1, 1)
final_sheet = np.hstack((data,tp_col))
# saving things
final_sheet = pd.DataFrame(final_sheet[:], columns=['WBPaper ID', 'WBGene', 'Gene', 'WBStrain', 'Strains',                                                'Variants', 'Mutation', 'Gene-Var combo', 'Variation type', 'Functional effect',                                                'Generation method', 'Transcript', 'Warnings', 'Sentence', 'Result'])
final_sheet.to_csv("data/model_output/processed/final_verified.csv", index=False, encoding='utf-8')


# ### Checking how many matches are present in the ground truth for the processed papers

# In[54]:


all_from_truth = []
for ppr in paper_ids_processed:
    for label in ground_truth[ground_truth[:,0] == ppr]:
        label[-2] = label[-2].upper()
        all_from_truth.append(label)
len(all_from_truth)


# # Results:

# Out of 100 papers, we could find gene-mutation matches on 53 papers.   
# Those 53 papers had total of 2433 matches in ground truth (which were manually curated).  
# We were able to find 977 matches.  
# TP: 472, FP: 505  
# Precision: 48.3%  
# After manually checking the false positives and updating the ground truth file -  
# TP: 807, FP: 170  
# Precision: 82.59%  
#   
# -> Not all "FP" are FP. After manual verification of the final output, some were noticied to be true positive which were originally missed during the manual curation. 

# In[ ]:




