#!/usr/bin/env python
# coding: utf-8

import json
from pathlib import Path
import re
import os
import os.path
import platform

import nltk.data
import numpy as np
import pandas as pd
from regex_wrapper import regex_block
from settings import setSettings
from textpresso import textpresso_paper_text, wbtools_paper_text


def ner_mutations(tokenClassificationPipeline, stop_words, sentence):
    mutations = []
    try:
        ner_output = tokenClassificationPipeline(sentence)
        for i, grp in enumerate(ner_output):
            if grp['entity_group'] == 'LABEL_0':
                mut = grp['word']
                for j in range(i+1, len(ner_output)):
                    if ner_output[j]['entity_group'] == 'LABEL_1':
                        mut = mut + ' ' + ner_output[j]['word']
                    else:
                        # NER would be handling only data in NL form
                        if len(mut.split()) > 3 and any(word in mut.split() for word in stop_words):
                            mutations.append([mut, sentence])
                        break
    except(Exception):
        pass
    return mutations


def get_paper_sentences(wbpids, settings, store_ppr_path):
    '''
    Takes WB Paper IDs, returns a list of sentences after filtering
    Arg:
    wbpids - List of wb papers ids
        e.g. ['WBPaper00002379']
    config_path - Config file path
    store_ppr_path - Folder path to store the paper flatfiles
    retrieved from TextPresso for future use
    Returns:
    paperid_sentence_list: List of paper ID and sentence
        e.g. [['WBPaper00002379', 'First sentence'],
        ['WBPaper00002379', 'Second sentence'], ....]
    '''
    token = settings['db_config']['textpresso']['token']

    stop_words = set(nltk.corpus.stopwords.words('english'))
    stop_words = [w for w in stop_words if len(w) > 1]

    all_special_chars = []
    with open('data/nala/train_dev.json') as f:
        for jsonObj in f:
            nala_json = json.loads(jsonObj)['tokens']
            for word in nala_json:
                if not word.isalnum():
                    all_special_chars.append(word)
    # list of special characters to keep during inference
    # helps with clearing out the bad characters from old papers
    all_special_chars = list(set(all_special_chars))

    temp_paperid_sentence = np.array([])
    if os.path.isfile('data/id_and_sentence.csv'):
        temp_paperid_sentence = pd.read_csv("data/id_and_sentence.csv", lineterminator='\n', dtype=str).to_numpy()  # WBPaper ID, Sentence
    paperid_sentence_list = []

    for curr_ppr_i, id in enumerate(wbpids):
        # print(f"{curr_ppr_i+1}", end = " ")
        # textpresso_paper_text() also saves the text in flatfiles for future use
        try:
            paper_path = textpresso_paper_text(id, store_ppr_path, token)
            txt = Path(paper_path).read_text().split('\n')
            # deals with empty text files with only "0"
            if len(txt) == 2:
                if temp_paperid_sentence.size != 0:
                    txt = temp_paperid_sentence[temp_paperid_sentence[:, 0] == id[7:]][:, 1]
                    # incase the loaded numpy file didn't have the required paper
                    if len(txt) == 0 and platform.system() != 'Windows':
                        txt = wbtools_paper_text(settings, id)
                elif platform.system() != 'Windows':
                    txt = wbtools_paper_text(settings, id)

            for row in txt:
                if row.find('fifi') != -1:
                    if temp_paperid_sentence.size != 0:
                        txt = temp_paperid_sentence[temp_paperid_sentence[:, 0] == id[7:]][:, 1]
                        # incase the loaded numpy file didn't have the required paper
                        if len(txt) == 0 and platform.system() != 'Windows':
                            txt = wbtools_paper_text(settings, id)
                    elif platform.system() != 'Windows':
                        txt = wbtools_paper_text(settings, id)
                    break
        except:
            print(f'Could not download paper {id}')
            continue

        count_total_rows = len(txt)
        for current_i, row in enumerate(txt):
            if (row.lower().find("we thank") == 0 or
                row.lower().find("this work was supported") == 0 or
                row.lower().find("references") == 0 or
                row.lower().find("we also thank") == 0 or
                row.lower().find("this research was supported") == 0 or
                row.lower().find("we acknowledge") == 0 or
                row.lower().find("acknowledgments") == 0 or
                row.lower().find('literature cited') != -1):
                if current_i > count_total_rows/3:
                    break

            # usually is bad sentence
            if len(row) < 40 or not any(word in row.lower().split() for word in stop_words):
                continue
            # remove sentences with links and email ids
            if re.search('\S+@\S+\.', row) or re.search('www\.\S+\.', row) or re.search('http.?://', row):
                continue
            # filters one word sentences
            if len(row.split()) == 1:
                continue
            # sentences comprised of only single characters
            # ^ seems to be issue with wbtools extraction pipeline
            if all(len(word) < 5 for word in row.split()):
                continue
            row = re.sub("\( *cid *: *\d+ *\)", " ", row)
            # TODO: replace this block with a regex sub
            temp_row = row
            for c in temp_row:
                if (not c.isalnum() and not c == ' ') and c not in all_special_chars:
                    row = row.replace(c, "")

            # fixes bad space between each character of flanking sequence from old papers
            # Switching this off as it increases the processing time
            # also affects very small subset of old papers so not worth the extra time
            flanking_regex = re.compile('([ACTG]( +)){4,}')
            for m in flanking_regex.finditer(row):
                span = (m.start(0), m.end(0))
                span = row[span[0]:span[1]-1]
                correct_flank = re.sub('([ACTG])( +)', r'\1', row)
                row = row.replace(span, correct_flank)
            row = 'Line ' + str(current_i) + ': ' + row.strip()
            paperid_sentence_list.append((id, row))
    return paperid_sentence_list[1:]


def findVariants(settings, ids_to_extract):

    # Get the paper texts from textpresso API and wbtools
    # How it works:
    # - First paper ID is searched through textpresso API.
    # - If the recieved output is blank, then wbtools is used.

    custom_mut_extract = settings['custom_mut_extract']
    tokenClassificationPipeline = settings['nala_ner']
    stop_words = settings['stop_words']

    paperid_sentence_list = get_paper_sentences(ids_to_extract, settings, store_ppr_path='data/wbpapers')

    # remove duplicates keeping the order
    seen = set()
    paperid_sentence_list = np.array([x for x in paperid_sentence_list if x not in seen and not seen.add(x)])

    final = [
        ['temporary', 'temporary', 'temporary', 'temporary', 'temporary', 'temporary'],\
        ['WBPaper ID', 'Method', '*Genes', '*Gene-Variant combo', 'Mutation', 'Sentence']]
    total_sentences = len(paperid_sentence_list)
    print_snips = False

    for ppr_sen_idx, row in enumerate(paperid_sentence_list):
        if (ppr_sen_idx+1) % 50 == 0:
            print(f"{ppr_sen_idx+1}>{len(final)-1}", end=" ")
        paper_id = row[0]
        sentence = str()
        # traverse upto 2 sentences at a time (if they aren't super long)
        # this should be removed if table content thing from wbtools is ever fixed
        # or else, it might group lines of table which will result in higher false positive count
        limit = min(ppr_sen_idx+2, total_sentences)
        # some sentences - mostly table content are all in a single sentence
        # temp fix, need to have a nice sentence splitter to minimize manual verification time
        not_single_sentence = False
        for i in range(ppr_sen_idx, limit):

            sentence = sentence + paperid_sentence_list[i][1] + ' '
            raw_sentence = sentence

            if (len(sentence) > 250 and not_single_sentence):
                break
            if paper_id != paperid_sentence_list[i][0]:
                break

            var_plus_genes = ''
            # Look for the special data e.g. gene-variant combo (e.g 'ced-3(n2888)') only on single sentences
            if not not_single_sentence:
                var_plus_genes = []
                all_genes = []
                genome_versions = []
                annotation_versions = []

                for data_and_cat in custom_mut_extract.var_and_gene_close(sentence.strip()):
                    var_plus_genes.append(data_and_cat[0])
                if var_plus_genes:
                    var_plus_genes = list(set(var_plus_genes))
                    var_plus_genes = "'" + "', '".join(var_plus_genes) + "'"
                else:
                    var_plus_genes = ''

                for data_and_cat in custom_mut_extract.get_genes(sentence.strip()):
                    all_genes.append(data_and_cat[0])
                if all_genes:
                    all_genes = list(set(all_genes))
                    # Removing gene mentions from sentence
                    # e.g. "lysine 36 by MET-1/Set2" regex will classify this as protein mutation
                    # but MET-1 is a gene name so the mutation isn't valid
                    for gene in all_genes:
                        sentence = sentence.replace(gene, "")
                    all_genes = "'" + "', '".join(all_genes) + "'"
                else:
                    all_genes = ''

            output = regex_block(settings, sentence.strip())
            if output:
                mutations = []
                for mut_and_snip in output:
                    # temp fix to deal with same mutation getting detected due to stiching multiple sentences
                    if ((mut_and_snip[0] not in final[-1][4][1:-1].split(", ") and
                         mut_and_snip[0] not in final[-2][4][1:-1].split(", ")) and
                         mut_and_snip[0] not in mutations):
                        mutations.append(mut_and_snip[0])
                if mutations:
                    mutations = "'" + "', '".join(mutations) + "'"
                    if print_snips:
                        print(1, mutations)
                    final.append([paper_id, 'Regex', all_genes, var_plus_genes, mutations, raw_sentence.strip()])
                break

            output = ner_mutations(tokenClassificationPipeline, stop_words, sentence.strip())
            if output:
                mutations = []
                for mut_and_snip in output:
                    # temp fix to deal with same mutation getting detected due to stiching multiple sentences
                    if ((mut_and_snip[0] not in final[-1][4][1:-1].split(", ") and
                         mut_and_snip[0] not in final[-2][4][1:-1].split(", ")) and
                         not all(len(word) < 4 for word in mut_and_snip[0].split()) and
                         mut_and_snip[0] not in mutations):
                        mutations.append(mut_and_snip[0])
                if mutations:
                    mutations = "'" + "', '".join(mutations) + "'"
                    if print_snips:
                        print(2, mutations)
                    final.append([paper_id, 'NER', all_genes, var_plus_genes, mutations, raw_sentence.strip()])
                break

            # these data, if found, are going to be important if no mutations are in that sentence
            if var_plus_genes or all_genes:
                final.append([paper_id, '', all_genes, var_plus_genes, '', raw_sentence.strip()])

            not_single_sentence = True

    temp = final[2:]  # removing the temporary first row and header

    # this sheet will contain high number of duplicates
    # which will get filtered in the refinement process
    # columns with asterisk contain data which are useful
    # regardless of whether the sentence has  mutation info
    return pd.DataFrame(
        temp[:],
        columns=['WBPaper ID', 'Method', '* Genes', '* Gene-Variant combo',
                 'Mutation', 'Sentence'])


if __name__ == "__main__":
    settings = setSettings()

    # test: papers mentioned the remarks ace file in data/gsoc
    ids_to_extract = np.load('data/top100.npy').tolist()[98:]
    findVariants(settings, ids_to_extract)
