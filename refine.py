from pprint import pprint
import numpy as np
import pandas as pd
from pathlib import Path

from regex_wrapper import point_mut_block
from settings import setSettings
from utils.misc.regex_block import normalize_mutations


def read_WB_genes():
    wb_genes_1 = Path('data/gsoc/Gene_alias.1.txt').read_text().split('\n')
    wb_genes_2 = Path('data/gsoc/Gene_alias.2.txt').read_text().split('\n')
    wb_genes_3 = Path('data/gsoc/Gene_alias.3.txt').read_text().split('\n')

    wb_genes_1 = [r.split('\t') for r in wb_genes_1]
    wb_genes_2 = [r.split(' ') for r in wb_genes_2]
    wb_genes_3 = [r.split(' ') for r in wb_genes_3]

    all_wb_genes = dict()

    for row in wb_genes_1+wb_genes_2+wb_genes_3:
        if row[0] not in all_wb_genes.keys():
            all_wb_genes[row[0]] = []
        for gene in row[1:]:
            if len(gene) and gene.lower() not in all_wb_genes[row[0]]:
                all_wb_genes[row[0]].append(gene.lower())
    return all_wb_genes


def add_normalized_mutations_column(settings, dataFrame):
    data = dataFrame.to_numpy()
    resultArray = []
    # total_count = len(data)
    ner_count = 0
    regex_count = 0
    paper_mut_count = {}

    for i, row in enumerate(data):
        if row[1] != 'Regex':
            if row[1] == 'NER':
                ner_count += 1
            resultArray.append(np.insert(data[i], -1, '').tolist())
        else:
            paper_id = row[0]
            if paper_id not in paper_mut_count.keys():
                paper_mut_count[paper_id] = {}
            regex_count += 1
            norm_mutations = []
            mutations = data[i, -2][1:-1].split("', '")
            for raw_mut in mutations:
                mut = point_mut_block(settings, raw_mut)
                if mut:
                    # helps filtering obvious ones
                    for m in mut:
                        m = m.replace(",", "")
                        if m.find(')') != -1:
                            if m.find('(') == -1:
                                continue
                        try:
                            # takes care of filtering out bad mutations where
                            # wild residues and mutants are same e.g G123G
                            norm_mut = normalize_mutations(mut[0])
                            if norm_mut:
                                if norm_mut not in paper_mut_count[paper_id].keys():
                                    paper_mut_count[paper_id][norm_mut] = 0
                                paper_mut_count[paper_id][norm_mut] += 1
                                norm_mutations.append(norm_mut)
                        except KeyError:
                            print(m)
            if norm_mutations:
                norm_mutations = list(set(norm_mutations))
                norm_mutations = "'" + "', '".join(norm_mutations) + "'"
            else:
                norm_mutations = ''
            resultArray.append(np.insert(data[i], -1, norm_mutations).tolist())

    return pd.DataFrame(resultArray[:], columns=['WBPaper ID', 'Method', 'Genes', '*Gene-Variant combo ', 'Mutations', 'Normalized Mutations', 'Sentence'])


def add_WBgenes_column(all_wb_genes, data):
    data = data.to_numpy()
    updated_data = []
    paper_wbgene_count = {}

    for i, row in enumerate(data):
        if (i+1) % 100 == 0:
            print(f"{i+1}", end=" ")
        paper_id = row[0]
        genes = row[2]
        # sentence = row[-1]
        # checking if nan
        if type(genes) == float:
            col_genes = ''
        else:
            if paper_id not in paper_wbgene_count.keys():
                paper_wbgene_count[paper_id] = {}
            genes = genes[1:-1].split("', '")
            col_genes = []

            for gene in genes:
                for key, value in all_wb_genes.items():
                    if gene.lower() in value:
                        if key not in paper_wbgene_count[paper_id]:
                            paper_wbgene_count[paper_id][key] = 0
                        paper_wbgene_count[paper_id][key] += 1
                        col_genes.append(key)
                        break
            if col_genes:
                col_genes = list(set(col_genes))
                col_genes = "'" + "', '".join(col_genes) + "'"
            else:
                col_genes = ''
        updated_data.append([data[i,0], data[i,1], data[i,2], col_genes, data[i,3], data[i,4], data[i,5], data[i,6]])

    return pd.DataFrame(updated_data[:], columns=['WBPaper ID', 'Method', 'Genes', 'WBGenes', '*Gene-Variant combo ', 'Mutations', 'Normalized Mutations', 'Sentence'])


if __name__ == "__main__":
    settings = setSettings()

    data = pd.read_csv("data/model_output/processed/snippets_1.csv")
    norm = add_normalized_mutations_column(settings, data)
    wb_genes = read_WB_genes()
    res = add_WBgenes_column(wb_genes, norm)
    res.to_csv("out.csv", index=False, encoding='utf-8')
