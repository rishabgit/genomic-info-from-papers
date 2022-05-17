
import configparser

import nltk
from transformers import AutoTokenizer, AutoConfig, AutoModelForTokenClassification
from transformers import TokenClassificationPipeline
from utils.misc.regex_block import MutationFinder, TmVar, BOWdictionary, CustomWBregex


def setSettings():
    '''
    Returns a dictionary with the DBconfig, Regex objects and TC pipeline
    '''
    settings = {}

    db_config = configparser.ConfigParser()
    db_config.read('utils/all_config.cfg')

    #settings['textpresso'] = db_config['textpresso']
    settings['db_config'] = db_config

    mf_mut_extract = MutationFinder('data/regexs/mutationfinder_regex/seth_modified.txt')
    custom_mut_extract = CustomWBregex(db_config, extra_regex=True)
    bow_mut_extract = BOWdictionary()
    tmvar_mut_extract = TmVar('data/regexs/tmvar_regex/final_regex_path')
    model_name_or_path = 'models/nala'
    settings['mf_mut_extract'] = mf_mut_extract
    settings['custom_mut_extract'] = custom_mut_extract
    settings['tmvar_mut_extract'] = tmvar_mut_extract
    settings['bow_mut_extract'] = bow_mut_extract
    settings['model_name_or_path'] = model_name_or_path

    config = AutoConfig.from_pretrained(model_name_or_path)
    tokenizer = AutoTokenizer.from_pretrained(model_name_or_path)
    model = AutoModelForTokenClassification.from_pretrained(
        model_name_or_path,
        from_tf=bool(".ckpt" in model_name_or_path),
        config=config,
    )

    # LABEL_0 - B-mut, LABEL_1 - I-mut, LABEL_2 - O
    nala_ner = TokenClassificationPipeline(model=model, tokenizer=tokenizer, task='ner', aggregation_strategy='first')
    settings['nala_ner'] = nala_ner
    stop_words = set(nltk.corpus.stopwords.words('english'))
    stop_words = [w for w in stop_words if len(w) > 1]
    settings['stop_words'] = stop_words
    return settings


if __name__ == "__main__":
    settings = setSettings()
    print(settings['db_config']['textpresso']['token'])
