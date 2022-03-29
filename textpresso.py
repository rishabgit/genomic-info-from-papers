import calendar
from datetime import datetime
import json
import os
import re

import nltk.data
from wbtools.literature.corpus import CorpusManager


def textpresso_paper_text(wbpid, path, token):
    '''
    This sub takes a wbpid eg WBPaper00056731 and
    returns the fulltext paper in sentences
    '''
    ft = [0]
    # Check that wbpid is a valid WBPaper
    if not re.match('WBPaper', wbpid):
        print(wbpid, "is not a valid WBPaper ID")
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
        os.system(command)

    # Read the paper, and split into sentences
    if os.path.exists(fn) and os.path.getsize(fn) > 20:
        # Open our JSON file and load it into python
        input_file = open(fn)
        json_array = json.load(input_file)
        for item in json_array:
            abs = item["abstract"]
            fullt = item["fulltext"]
            tokenizer = nltk.data.load('tokenizers/punkt/english.pickle')
            ft = tokenizer.tokenize(abs)
            ftt = tokenizer.tokenize(fullt)
            ft = ft + ftt
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


def wbtools_paper_text(settings, wbpid):
    db_name = settings['wb_database']['db_name']
    db_user = settings['wb_database']['db_user']
    db_password = settings['wb_database']['db_password']
    db_host = settings['wb_database']['db_host']
    ssh_host = settings['wb_database']['ssh_host']
    ssh_user = settings['wb_database']['ssh_user']
    ssh_passwd = settings['wb_database']['ssh_passwd']

    cm = CorpusManager()
    # sectioning might not be always correct, text processing is done separately in the pipeline
    # remove_sections = [PaperSections.ACKNOWLEDGEMENTS, PaperSections.REFERENCES, PaperSections.RELATED_WORK, PaperSections.INTRODUCTION]
    remove_sections = []
    paper_id = wbpid[7:]
    cm.load_from_wb_database(db_name=db_name, db_user=db_user, db_password=db_password,
                             db_host=db_host, paper_ids=[paper_id],
                             ssh_host=ssh_host, ssh_user=ssh_user, ssh_passwd=ssh_passwd,
                             load_bib_info=False, load_afp_info=False, load_curation_info=False)
    sentences = cm.get_paper(paper_id).get_text_docs(remove_sections=remove_sections, split_sentences=True)
    return sentences


def wbtools_get_papers_last_month(settings, day=None):
    ''' List of paper Ids since the last day of previous month'''
    if day is None:
        day = datetime.now()
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
    ssh_host = settings['wb_database']['ssh_host']
    ssh_user = settings['wb_database']['ssh_user']
    ssh_passwd = settings['wb_database']['ssh_passwd']
    cm = CorpusManager()
    cm.load_from_wb_database(
        db_name=db_name, db_user=db_user, db_password=db_password,
        db_host=db_host, from_date=query_date,
        ssh_host=ssh_host, ssh_user=ssh_user, ssh_passwd=ssh_passwd)

    return [paper.paper_id for paper in cm.get_all_papers()]


if __name__ == "__main__":
    from settings import setSettings
    settings = setSettings()
    print(wbtools_get_papers_last_month(settings['db_config']))
    feb_day = datetime.strptime('2022-02-15', '%Y-%m-%d')
    print(wbtools_get_papers_last_month(settings['db_config'], day=feb_day))
