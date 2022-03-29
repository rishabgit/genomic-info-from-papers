#!/usr/bin/env python
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


def read_strains():
    strains = Path('data/gsoc/Strains.txt').read_text().split('\n')
    strains = [r.split('\t') for r in strains][:-1]
    all_wb_strains = dict()
    for row in strains:
        if row[0] not in all_wb_strains.keys():
            all_wb_strains[row[0]] = []
        for strain in row[1:]:
            if len(strain) and strain.lower() not in all_wb_strains[row[0]]:
                all_wb_strains[row[0]].append(strain.lower())
    strainsList = [s for row in strains for s in row[1:] if len(s) and not s.isdigit()]
    return (strainsList, all_wb_strains)


def read_variation_types():
    variation_types = pd.read_csv("data/gsoc/Variation_type.csv").to_numpy()
    variation_types = [t.replace("_", " ") for t in variation_types[:,2] if type(t)!=float]
    return variation_types


def add_normalized_mutations_column(settings, paper_mut_count, dataFrame):
    data = dataFrame.to_numpy()
    resultArray = []
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


def add_gene_wbid_column(all_wb_genes, paper_wbgene_count, data):
    data = data.to_numpy()
    updated_data = []

    for i, row in enumerate(data):
        #if (i+1) % 100 == 0:
            # print(f"{i+1}", end=" ")
        paper_id = row[0]
        genes = row[2]
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

    return pd.DataFrame(np.array(updated_data), columns=['WBPaper ID', 'Method', 'Genes', 'WBGenes', '*Gene-Variant combo ', 'Mutations', 'Normalized Mutations', 'Sentence'])


def add_pair_gene_mutation_columns(data):
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


def add_warnings(final_sheet_input, paper_mut_count, paper_wbgene_count, data):
    final_sheet = np.array(final_sheet_input)
    updated_sheet = []

    for i, row in enumerate(final_sheet):
        warnings = []
        paper_id = row[0]
        # wbgene = row[1]
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
            temp_warn_store = 'Sentence has multiple mutations which did not match:'
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


def add_strains(dataframe, strainsList, strainsDict):
    data = dataframe.to_numpy()
    OPENING_CLOSING_REGEXES = [r'(?:^|[^0-9A-Za-z])(', r')(?:^|[^0-9A-Za-z])']
    all_strain = OPENING_CLOSING_REGEXES[0] + '|'.join(strainsList) + OPENING_CLOSING_REGEXES[1]
    all_strain = [re.compile(r,re.IGNORECASE) for r in [all_strain]]
    updated_data = []

    if data:
        for i, sent in enumerate(data[:, -1]):
            # if (i+1) % 100 == 0:  # print(f"{i+1}", end=" ")
            paper_strains = []
            for regex in all_strain:
                for m in regex.finditer(sent):
                    span = (m.start(0), m.end(0))
                    raw = (sent[span[0]:span[1]]).strip()
                    raw = raw[1:] if not raw[0].isalnum() else raw
                    raw = raw[:-1] if not raw[-1].isalnum() else raw
                    if len(raw.strip()) > 1 and not raw.strip().isdigit():
                        paper_strains.append(raw.strip())
            if paper_strains:
                paper_strains = list(set(paper_strains))
                col_wbid = []
                for strain in paper_strains:
                    for key, value in strainsDict.items():
                        if strain.lower() in value:
                            col_wbid.append(key)
                            break
                paper_strains = "'" + "', '".join(paper_strains) + "'"
                if col_wbid:
                    col_wbid = list(set(col_wbid))
                    col_wbid = ", ".join(col_wbid)
                else:
                    col_wbid = ''
                    # lazy way to deal with bad snippets due to special characters
                    # in the Strains.txt file which are caught in regex
                    paper_strains = ''
            else:
                paper_strains = ''
                col_wbid = ''
            updated_data.append([data[i,0], data[i,1], data[i,2], col_wbid, paper_strains, data[i,3], data[i,-4], data[i,-3], data[i,-2], data[i,-1]])
    return np.array(updated_data)


def add_variation_type_column(data, variation_types):
    updated_sheet = []
    for i, row in enumerate(data):
        sent = row[-1]
        col_var_type = []
        for sub in variation_types:
            if re.search(sub, sent, re.IGNORECASE):
                col_var_type.append(sub)
        if col_var_type:
            col_var_type = list(set(col_var_type))
            col_var_type = ", ".join(col_var_type)
        else:
            col_var_type = ''
        updated_sheet.append(np.insert(row, -3, col_var_type).tolist())
    return np.array(updated_sheet)


def add_functional_effect_column(data):
    functional_effect = ['function uncertain', 'transcript function', 'translational product function', \
        'decreased transcript level', 'increased transcript level', 'decreased transcript stability', \
        'gain of function', 'dominant negative', 'dominant negativ', 'antimorphic', \
        'hypermorphic', 'neomorphic', 'conditional activity', 'hypomorphic', 'amorphic', \
        'repressible', 'misexpressed']
    common_gen_methods = ['CRISPR', 'ENU', 'EMS']
    updated_sheet = []

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
    return np.array(updated_sheet)


def add_variants(data):
    OPENING_CLOSING_REGEXES = [r'(?:^|[^0-9A-Za-z])(', r')(?:^|[^0-9A-Za-z])']

    # the allele regex and db idea was stolen from wbtools
    allele_designations = np.load('data/gsoc/wbtools/wb_allele_designations.npy').astype('U6')
    alleles_variations = np.load('data/gsoc/wbtools/wb_alleles_variations.npy').astype('U6')
    DB_VAR_REGEX = r'({designations}|m|p|ts|gf|lf|d|sd|am|cs)([0-9]+)'
    var_regex_1 = OPENING_CLOSING_REGEXES[0] + DB_VAR_REGEX.format(designations="|".join(allele_designations)) + OPENING_CLOSING_REGEXES[1]
    all_var = OPENING_CLOSING_REGEXES[0] + '|'.join(alleles_variations) + '|' + var_regex_1 + OPENING_CLOSING_REGEXES[1]
    all_var = [re.compile(r,re.IGNORECASE) for r in [all_var]]
    updated_data = []

    if data:
        for i, sent in enumerate(data[:, -1]):
            variants = []
            for regex in all_var:
                for m in regex.finditer(sent):
                    span = (m.start(0), m.end(0))
                    raw = (sent[span[0]:span[1]]).strip()
                    raw = raw[1:] if not raw[0].isalnum() else raw
                    raw = raw[:-1] if not raw[-1].isalnum() else raw
                    if len(raw.strip()) > 1:
                        variants.append(raw.strip())
            if variants:
                variants = list(set(variants))
                variants = "'" + "', '".join(variants) + "'"
            else:
                variants = ''
            updated_data.append([data[i,0], data[i,1], data[i,2], data[i,3], data[i,4], variants, data[i,-5], data[i,-4], data[i,-3], data[i,-2], data[i,-1]])

    return np.array(updated_data)


def refine(dataframe):
    paper_mut_count = {}
    paper_wbgene_count = {}
    settings = setSettings()
    wb_genes = read_WB_genes()
    wb_gene_and_prot = read_WB_gene_and_prot()
    strainsList, strainsDict = read_strains()
    variation_types = read_variation_types()

    norm = add_normalized_mutations_column(settings, paper_mut_count, dataframe)
    withWBgenes = add_gene_wbid_column(wb_genes, paper_wbgene_count, norm)

    data = withWBgenes.to_numpy()
    gene_and_mutation = add_pair_gene_mutation_columns(data)
    with_transcript = add_transcript(gene_and_mutation, wb_genes, wb_gene_and_prot)
    with_warnings = add_warnings(with_transcript, paper_mut_count, paper_wbgene_count, data)
    with_strains = add_strains(with_warnings, strainsList, strainsDict)
    with_variants = add_variants(with_strains)
    with_variation_type = add_variation_type_column(
        with_variants, variation_types)
    with_functional_effect = add_functional_effect_column(with_variation_type)
    out_list = []
    if with_functional_effect:
        out_list = with_functional_effect[:]
    return pd.DataFrame(out_list, columns=['WBPaper ID', 'WBGene', 'Gene', 'WBStrain', 'Strains', 'Variants',
                                           'Mutation', 'Gene-Var combo', 'Variation type', 'Functional effect',
                                           'Generation method', 'Transcript', 'Warnings', 'Sentence'])


if __name__ == "__main__":
    data = pd.read_csv("data/model_output/processed/snippets_1.csv")
    out = refine(data)
    out.to_csv("out.csv", index=False, encoding='utf-8')
