import configparser
import os
import sys
import os.path
import argparse
import subprocess
import requests
import csv
import re
import json
import nltk.data
import os
from xml.dom import minidom
from bs4 import BeautifulSoup
from wbtools.literature.corpus import CorpusManager
import platform
import numpy as np
from pathlib import Path
from nltk.corpus import stopwords


def textpresso_paper_text(wbpid, path, token):
    """This sub takes a wbpid eg WBPaper00056731 and returns the fulltext paper in sentences"""
    
    ft=[0];
    # Check that wbpid is a valid WBPaper
    if not re.match( 'WBPaper', wbpid):
        print (wbpid, "is not a valid WBPaper ID")
        return ft
    # Download paper if it doesn't exist
    fn = path + '/temp/' + wbpid + '.json'

    if os.path.exists(fn) and os.path.getsize(fn) > 16:
        pass
    else:
        com1 = '-o '+fn +'\n-k '+ '\n'+'-d "{\\"token\\":\\"'+ token + '\\", \\"query\\": {\\"accession\\": \\"' + wbpid +'\\", \\"type\\": \\"document\\", \\"corpora\\": [\\"C. elegans\\"]}, \\"include_fulltext\\": true}"'
        configf= path +'/temp/' + wbpid + '.tmp.config'
        curlf = open(configf,'w')
        print (com1, file=curlf)
        curlf.close()
        command = 'curl -o '+ fn +' -K '+ configf+' https://textpressocentral.org:18080/v1/textpresso/api/search_documents' 
        comlist = command.split()
        os.system(command)

    # Read the paper, and split into sentences
    if os.path.exists(fn) and os.path.getsize(fn) > 20:
        # Open our JSON file and load it into python
        input_file = open (fn)
        json_array = json.load(input_file)
        for item in json_array:
            abs = item["abstract"]
            fullt =  item["fulltext"]
            tokenizer = nltk.data.load('tokenizers/punkt/english.pickle')
            ft = tokenizer.tokenize(abs)
            ftt=tokenizer.tokenize(fullt)
            ft = ft +ftt
    else:
        # some paper texts are blank for some reason
        # pipeline uses wbtools to get sentences in such case
        pass

    outfilen = os.path.join(path, 'text_flatfiles', wbpid+'.txt')
    outf = open(outfilen, 'w')
    for sen in ft:
        sen =str(sen)
        print(sen, file=outf)
    outf.close()
    return outfilen


'''
Sentences out of wbtools are sometimes weird, especially for the old papers
This increases the false neg rate so this is used only when the textpresso api provides bad text
Most of the issues have band aid fixes in the text preprocessing cell below but there are some 
with no easy fix or need to be worked on - 
1. Table content is extracted column wise, not row wise 
2. Some sentences have white space between every character and are also somehow inverted (???)
    e.g. Check out the sentences of WBPaper00002018 (1994)
    Line 154 -  '0 7 1 u ( 7 - c e m    ) F (   .'
    inverted and without the white space is (mec-7)u170 which is extremely useful and will get missed 
    by the pipeline unless processed correctly.
TODO: Not sure how to solve point 1 but point 2 is easy to solve and also helps a LOT.
    Rishab, work on this after you complete the project. 
    Not high priority as this might be only for the >10 year old papers (which are already manually curated)
'''
cm = CorpusManager()
def wbtools_paper_text(wbpid, db_name, db_user, db_password, db_host, ssh_host,\
    ssh_user, ssh_passwd):
    # sectioning might not be always correct, text processing is done separately in the pipeline
    # remove_sections = [PaperSections.ACKNOWLEDGEMENTS, PaperSections.REFERENCES, PaperSections.RELATED_WORK, PaperSections.INTRODUCTION]
    remove_sections = []
    paper_id = wbpid[7:]
    cm.load_from_wb_database(db_name=db_name, db_user=db_user, db_password=db_password,
        db_host=db_host, paper_ids=[paper_id],
        ssh_host=ssh_host, ssh_user=ssh_user, ssh_passwd=ssh_passwd,
        load_bib_info=False, load_afp_info=False, load_curation_info=False)
    sentences = cm.get_paper(paper_id).get_text_docs(remove_sections=remove_sections,split_sentences=True)
    return sentences


def get_paper_sentences(wbpids, config_path, store_ppr_path):
    '''
    Takes WB Paper IDs and returns a list of sentences from those papers after filtering
    Arg:
    wbpids - List of wb papers ids 
        e.g. ['WBPaper00002379']
    config_path - Config file path
    store_ppr_path - Folder path to store the paper flatfiles retrieved from TextPresso for future use
    Returns:
    paperid_sentence_list: List of paper ID and sentence
        e.g. [['WBPaper00002379', 'First sentence'], ['WBPaper00002379', 'Second sentence'], ....]
    '''
    stop_words = set(stopwords.words('english'))
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

    config = configparser.ConfigParser()
    config.read(config_path)
    token = config['textpresso']['token']
    db_name=config['wb_database']['db_name']
    db_user=config['wb_database']['db_user']
    db_password=config['wb_database']['db_password']
    db_host=config['wb_database']['db_host']
    ssh_host=config['wb_database']['ssh_host']
    ssh_user=config['wb_database']['ssh_user']
    ssh_passwd=config['wb_database']['ssh_passwd']

    temp_paperid_sentence = np.array([])
    if os.path.isfile('data/id_and_sentence.npy'):
        temp_paperid_sentence = np.load('data/id_and_sentence.npy')
    paperid_sentence_list = np.array([['WBPaperID', 'Sentence']])

    for id in wbpids:
        # textpresso_paper_text() also saves the text in flatfiles for future use 
        paper_path = textpresso_paper_text(id, store_ppr_path, token)
        txt = Path(paper_path).read_text().split('\n')
        # deals with empty text files with only "0"
        if len(txt) == 2:
            if platform.system() != 'Windows':
                txt = wbtools_paper_text(id[7:], db_name, db_user, db_password, db_host, ssh_host,\
                    ssh_user, ssh_passwd)
            elif temp_paperid_sentence.size != 0:
                txt = temp_paperid_sentence[temp_paperid_sentence[:, 0] == id[7:]][:, 1]
            
        for row in txt: 
            if row.find('fifi') == -1:
                if platform.system() != 'Windows':
                    txt = wbtools_paper_text(id[7:], db_name, db_user, db_password, db_host, ssh_host,\
                        ssh_user, ssh_passwd)
                elif temp_paperid_sentence.size != 0:
                    txt = temp_paperid_sentence[temp_paperid_sentence[:, 0] == id[7:]][:, 1]
                break
            
        count_total_rows = len(txt)
        for current_i, row in enumerate(txt):
            if row.lower().find("we thank") == 0 or row.lower().find("this work was supported") == 0 \
                or row.lower().find("references") == 0 or row.lower().find("we also thank") == 0 \
                or row.lower().find("this research was supported") == 0 or row.lower().find("we acknowledge") == 0 \
                or row.lower().find("acknowledgments") == 0 or row.lower().find('literature cited') != -1:
                if current_i > count_total_rows/2:
                    break

            # usually is bad sentence
            if len(row) < 40 or not any(word in row.lower().split() for word in stop_words):
                continue
            # remove sentences with links and email ids
            if re.search('\S+@\S+\.', row) or re.search('www.\S+\.', row):
                continue
            # filters one word sentences
            if len(row.split()) == 1:
                continue
            # sentences comprised of only single characters 
            # ^ seems to be issue with wbtools extraction pipeline 
            if all(len(word) < 5 for word in row.split()):
                continue
            row = re.sub("\( *cid *: *\d+ *\)", " ", row)
            temp_row = row
            for c in temp_row:
                if (not c.isalnum() and not c == ' ') and c not in all_special_chars:
                        row = row.replace(c, "")
            # fixes bad space between each character of flanking sequence from old papers
            flanking_regex = re.compile('([ACTG]( +)){4,}', re.IGNORECASE)
            for m in flanking_regex.finditer(row):
                span = (m.start(0), m.end(0))   
                span = row[span[0]:span[1]-1]
                correct_flank = re.sub('([ACTG])( +)', r'\1', row, flags=re.I)
                row = row.replace(span, correct_flank)

            # filters out repeated lines, e.g. check out WBPaper00028727.txt in flatfiles folder
            if row not in paperid_sentence_list[paperid_sentence_list[:,0]==id][:,1]:
                paperid_sentence_list = np.vstack((paperid_sentence_list, [id, row]))
    return paperid_sentence_list[1:]
