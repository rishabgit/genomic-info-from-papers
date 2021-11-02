from argparse import ArgumentParser
import re
import os
from typing import List, Tuple
import nltk
import math
import pandas as pd
import json 
import requests
from tqdm import tqdm
import zipfile


from genomicinfo.entity_extraction.ner.nala.readers import HTMLReader
from genomicinfo.entity_extraction.ner.nala.annotation_readers import AnnJsonAnnotationReader, AnnJsonMergerAnnotationReader
from genomicinfo.entity_extraction.ner.nala.definers import ExclusiveNLDefiner
from genomicinfo.entity_extraction.ner.nala.tokenizers import TmVarTokenizer
from genomicinfo.entity_extraction.ner.nala.spliters import NLTK_SPLITTER
from genomicinfo.entity_extraction.ner.nala.labelers import BIOLabeler


parser = ArgumentParser()
parser.add_argument('--folder', type=str, default='IDP4+',\
    help='Path to folder containing the data to use for training NER.\nDefault: IDP4+\nNote: Currently, only data stored in the IDP4 format is supported.')

args = parser.parse_args()


class NERDataPrep():
    def __init__(self, folder_name: str):
        data_directory = os.path.join(os.path.dirname(__file__), "training_data")
        os.makedirs(data_directory, exist_ok=True)
        if folder_name == 'IDP4+':
            idp4_directory = os.path.join(data_directory, "tagtog_IDP4+_anndoc")
            if not os.path.isdir(idp4_directory):
                print('Downloading and unzipping main data zip file...')
                zip_path = os.path.join(data_directory, 'tagtog_IDP4+_anndoc.zip')
                self.download(
                    url='https://s3.amazonaws.com/net.tagtog.public/resources/corpora/tagtog_IDP4%2B_anndoc.zip', 
                    fname=zip_path
                    )
                self.unzip(zip_path, data_directory)
                os.remove(zip_path)

            self.data_setup_using_nala(idp4_directory)
            self.convert_bio2_text_to_json(idp4_directory)
            print('Data convertion complete.')
        
    def data_setup_using_nala(self, path: str):
        MUT_CLASS_ID = 'e_2'
        base_folder = os.path.join(path, "tagtog_IDP4")
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

        base_folder = os.path.join(path, "tagtog_nala_anndoc")
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

        base_folder = os.path.join(path, "tagtog_nala_discoveries")
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

        dataset.prune_sentences(0.1) # tried with 0.8, no change in NER acc 

        self.create_bio2_text_files(dataset, path)

    def create_bio2_text_files(self, dataset, path):
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
            print(len(temp))

            if i > math.ceil(len(final)//10000)*0.8:
                file = open(os.path.join(path, "devel.txt"), "a", encoding="utf-8")
            file = open(os.path.join(path, "train_dev.txt"), "a", encoding="utf-8")

            for index in range(len(temp)):
                if i*chunk_size + index > devel_thres and ok_to_switch:
                    file.close()
                    file = open(os.path.join(path, "devel.txt"), "a", encoding="utf-8")
                if temp[index]:
                    file.write(str(temp[index][0]) + " " + str(temp[index][1]) + "\n")
                    ok_to_switch = False
                else:
                    file.write("\n")
                    ok_to_switch = True
            file.close()

    @staticmethod
    def convert_bio2_text_to_json(path: str):

        data = []
        with open(os.path.join(path, "devel.txt"), 'r', encoding="utf-8") as f_in:
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
                    "tokens" : token, 
                    "tags" : label, 
                }  
                with open(os.path.join(path, "devel.json"), "a", encoding="utf-8") as outfile: 
                    json.dump(dictionary, outfile)
                    outfile.write('\n')
                token = []
                label = []

        data = []
        with open(os.path.join(path, "train_dev.txt"), 'r', encoding="utf-8") as f_in:
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
                    "tokens" : token, 
                    "tags" : label, 
                }  
                with open(os.path.join(path, "train_dev.json"), "a", encoding="utf-8") as outfile: 
                    json.dump(dictionary, outfile)
                    outfile.write('\n')
                token = []
                label = []

    @staticmethod
    def download(url: str, fname: str):
        resp = requests.get(url, stream=True)
        total = int(resp.headers.get('content-length', 0))
        with open(fname, 'wb') as file, tqdm(
            desc='Progress',
            total=total,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for data in resp.iter_content(chunk_size=1024):
                size = file.write(data)
                bar.update(size)

    @staticmethod
    def unzip(path: str, unzip_path: str):
        with zipfile.ZipFile(path,"r") as zip_ref:
            zip_ref.extractall(unzip_path)

if __name__ == '__main__':
    data_prep = NERDataPrep(args.folder)