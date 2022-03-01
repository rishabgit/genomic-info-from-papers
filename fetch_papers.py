import numpy as np
import pandas as pd

from settings import setSettings
from wbtools.literature.corpus import CorpusManager


def fetchPapers(db_config):
    '''
    Create a temporary cvs file to store the paper sentences
    to avoid spending time on pulling sentences while in development
    Format - [paper_id, sentence]
    '''

    final_data = []

    remove_sections = []
    # random 100 papers mentioned the remarks ace file in data/gsoc
    paper_ids = np.load('data/top100.npy').tolist()
    cm = CorpusManager()
    for i, paper_id in enumerate(paper_ids):
        paper_id = paper_id[7:]
        cm.load_from_wb_database(
            db_name=db_config['wb_database']['db_name'],
            db_user=db_config['wb_database']['db_user'],
            db_password=db_config['wb_database']['db_password'],
            db_host=db_config['wb_database']['db_host'],
            paper_ids=[paper_id],
            ssh_host=db_config['wb_database']['ssh_host'],
            ssh_user=db_config['wb_database']['ssh_user'],
            ssh_passwd=db_config['wb_database']['ssh_passwd'],
            load_bib_info=False,
            load_afp_info=False,
            load_curation_info=False)
        sentences = cm.get_paper(paper_id).get_text_docs(
            remove_sections=remove_sections, split_sentences=True)
        for sent in sentences:
            final_data.append([paper_id, sent])
        print(i, end=" ")
    final_data = pd.DataFrame(final_data[:], columns=['WBPaper ID', 'Sentence'])
    final_data.to_csv('data/id_and_sentence.csv',
                      index=False,
                      encoding='utf-8')


if __name__ == "__main__":
    settings = setSettings()
    fetchPapers(settings['db_config'])

