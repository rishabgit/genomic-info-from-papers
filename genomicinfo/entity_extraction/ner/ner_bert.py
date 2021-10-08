import re
import os
from typing import List, Tuple
import nltk
import torch
from transformers import AutoTokenizer, AutoConfig, AutoModelForTokenClassification
from transformers import TokenClassificationPipeline

from genomicinfo.entity_extraction.abstract_extractor import AbstractEntityExtractor

nltk.download('stopwords')
nltk.download('punkt')

class BERTEntityExtractor(AbstractEntityExtractor):

    def __init__(self):
        super().__init__()
        # Load the model
        self.bert_ner  = self._load_model(os.path.join(os.path.dirname(__file__), "saved_models", "nala"))
        # Loading stopwords from nltk 
        self.stop_words = set(nltk.corpus.stopwords.words('english'))
        self.stop_words = [w for w in self.stop_words if len(w) > 1]

    @staticmethod
    def _load_model(model_name_or_path):
        config = AutoConfig.from_pretrained(model_name_or_path)
        tokenizer = AutoTokenizer.from_pretrained(model_name_or_path)
        model = AutoModelForTokenClassification.from_pretrained(
            model_name_or_path,
            from_tf=bool(".ckpt" in model_name_or_path),
            config=config,
        )
        return TokenClassificationPipeline(model=model, tokenizer=tokenizer, task='ner', aggregation_strategy='first')

    def extract(self, sentence: str()) -> List[Tuple[str, str]]:
        final_list = []
        try:
            ner_output = self.bert_ner(sentence)
            for i, grp in enumerate(ner_output):
                if grp['entity_group'] == 'LABEL_0':
                    mut = grp['word']
                    for j in range(i+1, len(ner_output)):
                        if ner_output[j]['entity_group'] == 'LABEL_1':
                            mut  = mut + ' ' + ner_output[j]['word']
                        else:
                            # NER would be handling only data in NL form
                            if len(mut.split()) > 3 and any(word in mut.split() for word in self.stop_words):
                                final_list.append([mut, sentence])
                            break
        except:
            pass
        return final_list
