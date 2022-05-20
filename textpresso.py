import json
import logging
import requests
from typing import List

import nltk
from wbtools.literature.corpus import CorpusManager

logger = logging.getLogger(__name__)


def textpresso_paper_text(wbpid, token):
    ''' Takes a wbpid eg WBPaper00056731 and returns the fulltext paper
        in sentences
    '''
    tokenizer = nltk.data.load('tokenizers/punkt/english.pickle')
    url = 'https://textpressocentral.org:18080/v1/textpresso/api/search_documents'
    headers = {'Accept': 'application/json',
               'Content-Type': 'application/json'}
    body = json.dumps({'token': token,
                       'query': {'accession': wbpid,
                                 'type': 'document',
                                 'corpora': ['C. elegans']},
                       'include_fulltext': True})
    response = requests.post(url, data=body, headers=headers, verify=False)
    if response.status_code == 200:
        if response.json() is None:
            return []
        else:
            paper = response.json()[0]
        abstract = tokenizer.tokenize(paper['abstract'])
        fulltext = tokenizer.tokenize(paper['fulltext'])
        return abstract + fulltext
    else:
        return []


def wbtools_paper_text(settings, wbpid):
    db_name = settings['wb_database']['db_name']
    db_user = settings['wb_database']['db_user']
    db_password = settings['wb_database']['db_password']
    db_host = settings['wb_database']['db_host']
    file_server_host = settings['wb_database']['file_server_host']
    file_server_user = settings['wb_database']['file_server_user']
    file_server_passwd = settings['wb_database']['file_server_passwd']

    cm = CorpusManager()
    # sectioning might not be always correct, text processing is done separately in the pipeline
    # remove_sections = [PaperSections.ACKNOWLEDGEMENTS, PaperSections.REFERENCES, PaperSections.RELATED_WORK, PaperSections.INTRODUCTION]
    remove_sections = []
    paper_id = wbpid[7:]
    cm.load_from_wb_database(db_name=db_name, db_user=db_user, db_password=db_password,
                             db_host=db_host, paper_ids=[paper_id],
                             file_server_host=file_server_host, file_server_user=file_server_user,
                             file_server_passwd=file_server_passwd,
                             load_bib_info=False, load_afp_info=False, load_curation_info=False)
    sentences = cm.get_paper(paper_id).get_text_docs(remove_sections=remove_sections, split_sentences=True)
    return sentences


def wbtools_get_papers_last_month(settings, day=None, max_num_papers: int = None):
    ''' List of paper Ids since the last day of previous month'''
    if day is None:
        day = datetime.now()

    logger.info("getting papers from " + day.strftime("%Y-%m-%d"))
    if day.month == 1:
        previous_month = 12
        year = day.year - 1
    else:
        previous_month = day.month-1
        year = day.year
    first_day, last_day = calendar.monthrange(
        year, previous_month)
    query_date = datetime(
        year, previous_month, last_day)

    db_name = settings['wb_database']['db_name']
    db_user = settings['wb_database']['db_user']
    db_password = settings['wb_database']['db_password']
    db_host = settings['wb_database']['db_host']
    file_server_host = settings['wb_database']['file_server_host']
    file_server_user = settings['wb_database']['file_server_user']
    file_server_passwd = settings['wb_database']['file_server_passwd']
    cm = CorpusManager()
    cm.load_from_wb_database(
        db_name=db_name, db_user=db_user, db_password=db_password,
        db_host=db_host, from_date=query_date.strftime("%Y-%m-%d"),
        file_server_host=file_server_host, file_server_user=file_server_user, file_server_passwd=file_server_passwd,
        max_num_papers=max_num_papers)

    return [paper.paper_id for paper in cm.get_all_papers()]


def wbtools_get_papers(settings, paper_ids: List[str]):
    db_name = settings['wb_database']['db_name']
    db_user = settings['wb_database']['db_user']
    db_password = settings['wb_database']['db_password']
    db_host = settings['wb_database']['db_host']
    file_server_host = settings['wb_database']['file_server_host']
    file_server_user = settings['wb_database']['file_server_user']
    file_server_passwd = settings['wb_database']['file_server_passwd']

    cm = CorpusManager()
    cm.load_from_wb_database(
        db_name=db_name, db_user=db_user, db_password=db_password,
        db_host=db_host, paper_ids=paper_ids,
        file_server_host=file_server_host, file_server_user=file_server_user, file_server_passwd=file_server_passwd,
        load_bib_info=False, load_afp_info=False, load_curation_info=False)
    return cm


if __name__ == "__main__":
    from settings import setSettings
    settings = setSettings()
    print(textpresso_paper_text('WBPaper00002627', settings['db_config']['textpresso']['token']))
    # print(wbtools_get_papers_last_month(settings['db_config']))
    # feb_day = datetime.strptime('2022-02-15', '%Y-%m-%d')
    # print(wbtools_get_papers_last_month(settings['db_config'], day=feb_day))

