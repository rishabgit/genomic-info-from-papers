from pprint import pprint
import numpy as np
import pandas as pd
from pathlib import Path
import re

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


def read_WB_gene_and_prot():
    proteinfa = Path('data/gsoc/proteinfa/c_elegans.PRJNA13758.WS281.protein.fa').read_text().split('>')[1:]
    wb_gene_and_prot = dict() # {wbgene: [transcript, protein]}

    for row in proteinfa:
        wbgene = re.findall("WBGene[0-9]+", row)[0]
        protein = "".join(re.findall("\n.*", row)).replace('\n','')
        transcript = row.split(' ')[0]
        if wbgene not in wb_gene_and_prot.keys():
            wb_gene_and_prot[wbgene] = []
        wb_gene_and_prot[wbgene].append([transcript, protein])
    return wb_gene_and_prot


def add_normalized_mutations_column(settings, dataFrame):
    data = dataFrame.to_numpy()
    global paper_mut_count
    resultArray = []
    # total_count = len(data)
    ner_count = 0
    regex_count = 0

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


def create_pair_gene_mutation(dataframe):
    data = dataframe.to_numpy()
    paper_raw_info_compiled = []
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
                        paper_raw_info_compiled.append(
                            [ppr_id, g, m, sentence, gene_var])
    return paper_raw_info_compiled


def add_transcript(paper_raw_info_compiled, all_wb_genes, wb_gene_and_prot):
    matches = []
    final_sheet = []  # ppr_id, gene, transcript

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
                    matches.append(
                        [ppr_id, gene, mut, gene_var, transcript, sent])
            except IndexError:
                pass

    for r in matches:
        p = r[0]
        p, wbg, mut, gene_var, transcript, sent = r
        # Adding gene common names column, again
        # Current code doesn't keep any link between the WB gene name and the common name
        g_common_name = all_wb_genes[wbg]
        g_common_name = ', '.join(g_common_name)
        final_sheet.append(
            [p, wbg, g_common_name, mut, gene_var, transcript, sent])

    return final_sheet


def add_warnings(data):
    global paper_mut_count
    final_sheet = np.array(data)
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
    return pd.DataFrame(
        updated_sheet[:],
        columns=['WBPaper ID', 'WBGene', 'Gene', 'Mutation', 'Gene-Var combo', 'Transcript', 'Warnings', 'Sentence'])



paper_mut_count = {}
if __name__ == "__main__":
    settings = setSettings()

    data = pd.read_csv("data/model_output/processed/snippets_1.csv")
    norm = add_normalized_mutations_column(settings, data)
    wb_genes = read_WB_genes()
    withWBgenes = add_WBgenes_column(wb_genes, norm)
    wb_gene_and_prot = read_WB_gene_and_prot()
    gene_and_mutation = create_pair_gene_mutation(withWBgenes)
    with_transcript = add_transcript(
        gene_and_mutation, wb_genes, wb_gene_and_prot)
    pprint(with_transcript)
    #with_warnings = add_warnings(with_transcript)
