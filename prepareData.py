#!/usr/bin/env python
# coding: utf-8

import json
import math
import os
import pandas as pd
from pprint import pprint

from utils.nala.readers import HTMLReader
from utils.nala.annotation_readers import AnnJsonAnnotationReader, AnnJsonMergerAnnotationReader
from utils.nala.definers import ExclusiveNLDefiner
from utils.nala.tokenizers import TmVarTokenizer
from utils.nala.spliters import NLTK_SPLITTER
from utils.nala.labelers import BIOLabeler


def prepareData():

    MUT_CLASS_ID = 'e_2'
    dir_path = os.path.dirname(os.path.realpath(__file__))
    base_folder = os.path.join(dir_path, 'data/nala/tagtog_IDP4+_anndoc/tagtog_IDP4')
    pprint(base_folder)
    html_folder = os.path.join(base_folder, 'html')
    annjson_folder = os.path.join(base_folder, 'annjson')
    dataset = HTMLReader(html_folder).read()

    AnnJsonMergerAnnotationReader(
        os.path.join(annjson_folder, 'members'),
        read_only_class_id=MUT_CLASS_ID,
        strategy='union',
        entity_strategy='priority',
        priority=['Ectelion', 'abojchevski', 'sanjeevkrn', 'Shpendi'],
        delete_incomplete_docs=True).annotate(dataset)

    base_folder = os.path.join(dir_path, 'data/nala/tagtog_IDP4+_anndoc/tagtog_nala_anndoc')
    html_folder = os.path.join(base_folder, 'nala_plain_html', 'pool')
    annjson_folder = os.path.join(base_folder, 'nala_members_json')
    nala_anndoc_dataset = HTMLReader(html_folder).read()
    AnnJsonMergerAnnotationReader(
        os.path.join(annjson_folder, 'pool'),
        read_only_class_id=MUT_CLASS_ID,
        strategy='union',
        entity_strategy='priority',
        priority=['abojchevski', 'cuhlig', 'jmcejuela'],
        delete_incomplete_docs=True).annotate(nala_anndoc_dataset)
    dataset.extend_dataset(nala_anndoc_dataset)
    nala_anndoc_dataset = None

    base_folder = os.path.join(dir_path, 'data/nala/tagtog_IDP4+_anndoc/tagtog_nala_discoveries')
    html_folder = os.path.join(base_folder, 'html')
    annjson_folder = os.path.join(base_folder, 'annjson')
    nala_dis_dataset = HTMLReader(html_folder).read()
    AnnJsonAnnotationReader(
        annjson_folder,
        read_only_class_id=MUT_CLASS_ID,
        delete_incomplete_docs=True).annotate(nala_dis_dataset)
    dataset.extend_dataset(nala_dis_dataset)
    nala_dis_dataset = None

    definer = ExclusiveNLDefiner()
    definer.define(dataset)
    NLTK_SPLITTER.split(dataset)
    tokenizer = TmVarTokenizer()
    tokenizer.tokenize(dataset=dataset)

    # 0 (standard), 1(natural language) or 2 (semi standard)
    remove_subclasses = [0]
    dataset.delete_subclass_annotations(subclasses=remove_subclasses)
    labeler = BIOLabeler()
    labeler.label(dataset)

    # Run this cell to create file with sentence and label of true if it contains mutation, false if it doesn't.
    # Used in notebook 2 during testing.
    # NOTE: Do not run delete_subclass_annotations() above to store mutations of all three kinds - ST, SST and NL

    final_binary = []
    for i, part in enumerate(dataset.parts()):
        for tokenized_sent, raw_sent in zip(part.sentences, part.sentences_):
            mark_postitive = 0
            # sanity check, checking if first character same, very dumb way
            assert str(tokenized_sent[0])[0] == raw_sent[0], f'{str(tokenized_sent[0])[0], raw_sent[0]}'
            for token in tokenized_sent:
                if mark_postitive:
                    break
                for ann in part.annotations:
                    start = ann.offset
                    end = ann.offset + len(ann.text)
                    if start <= token.start < end:
                        mark_postitive = 1
                        break
            final_binary.append([raw_sent, mark_postitive])

    print('Sentences count : ', len(final_binary))

    data = pd.DataFrame(final_binary[:], columns=["Sentence", "Contains mutation?"])
    data.to_csv(os.path.join(dir_path, "data/nala/nala_binary.csv"), index=False, encoding='utf-8')

    dataset.prune_sentences(0.1)  # tried with 0.8, no change in NER acc
    final = []
    for doc_id, doc in dataset.documents.items():
        for part_id, part in doc.parts.items():
            for sentence in part.sentences:
                for token in sentence:
                    final.append([token.word, token.original_labels[0].value.split('-')[0]])
                final.append([])

    chunk_size = 10000
    total = len(final)
    devel_thres = math.ceil(total*0.8)

    for i in range(math.ceil(len(final)/chunk_size)):
        temp = final[i*chunk_size:(i+1)*chunk_size]

        if i > math.ceil(len(final)//10000)*0.8:
            file = open(os.path.join(dir_path, "data/nala/devel.txt"), "a", encoding="utf-8")

        file = open(os.path.join(dir_path, "data/nala/train_dev.txt"), "a", encoding="utf-8")

        ok_to_switch = False
        for index in range(len(temp)):
            if i*chunk_size + index > devel_thres and ok_to_switch:
                file.close()
                file = open(os.path.join(dir_path, "data/nala/devel.txt"), "a", encoding="utf-8")
            if temp[index]:
                file.write(str(temp[index][0]) + " " + str(temp[index][1]) + "\n")
                ok_to_switch = False
            else:
                file.write("\n")
                ok_to_switch = True
        file.close()

    # Convert BIO2 text files to JSON
    # Deals with [Issue #8698](https://github.com/huggingface/transformers/issues/8698)

    data = []
    with open(os.path.join(dir_path, 'data/nala/devel.txt'), 'r', encoding="utf-8") as f_in:
        for line in f_in:
            line = line.split()
            data.append(line)

    token = []
    label = []
    for row in data:
        if row:
            token.append(row[0])
            label.append(row[1])
        else:
            assert len(token) == len(label)
            # for l in label:
            #     if l not in ['B', 'O', 'I']:
            #         print('Error')
            #         break
            dictionary = {
                "tokens": token,
                "tags": label,
            }
            with open(os.path.join(dir_path, "data/nala/devel.json"), "a", encoding="utf-8") as outfile:
                json.dump(dictionary, outfile)
                outfile.write('\n')
            token = []
            label = []

    data = []
    with open(os.path.join(dir_path, 'data/nala/train_dev.txt'), 'r', encoding="utf-8") as f_in:
        for line in f_in:
            line = line.split()
            data.append(line)

    token = []
    label = []
    for row in data:
        if row:
            token.append(row[0])
            label.append(row[1])
        else:
            assert len(token) == len(label)
            # for l in label:
            #     if l not in ['B', 'O', 'I']:
            #         print('Error')
            #         break
            dictionary = {
                "tokens": token,
                "tags": label,
            }
            with open(os.path.join(dir_path, "data/nala/train_dev.json"), "a", encoding="utf-8") as outfile:
                json.dump(dictionary, outfile)
                outfile.write('\n')
            token = []
            label = []


if __name__ == "__main__":
    prepareData()
