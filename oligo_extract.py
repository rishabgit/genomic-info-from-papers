import json
import re
import time

import nltk
nltk.download("stopwords")
nltk.download("punkt")

import numpy as np
import pandas as pd
from settings import setSettings

from regex_wrapper import regex_block
from textpresso import textpresso_paper_text


def get_paper_sentences_with_TE(wbpids, settings):
    """
    Takes WB Paper IDs, returns a list of sentences after filtering

    Arg:
    wbpids - List of wb papers ids e.g. ['WBPaper00002379']
    settings - Dictionary with db_config properties and texpresso token

    Returns:
    paperid_sentence_list: List of paper ID and sentence e.g. [['WBPaper00002379', 'First sentence'],
        ['WBPaper00002379', 'Second sentence'], ....]
    """
    # Creates a list of commonly found stop words
    stop_words = set(nltk.corpus.stopwords.words("english"))
    stop_words = [w for w in stop_words if len(w) > 1]


    # Creates a list "all_special_chars" having special charaters
    # like -, +, = 
    # The list got it by checking which words of "train_dev.json" 
    # file has neither alphabets nor numeric characters, i.e it's
    # a special chracter like -, + or = 
    # Also avoids again appending of characters in the list
    all_special_chars = []
    with open("data/nala/train_dev.json") as f:
        for jsonObj in f:
            nala_json = json.loads(jsonObj)["tokens"]
            for word in nala_json:
                if not word.isalnum() and word not in all_special_chars:
                    all_special_chars.append(word)
    all_special_chars = list(set(all_special_chars))


    # Gets text from wb paper(s) using "textpresso_paper_text" and 
    # stores it in "txt" variable
    textpresso_token = settings["db_config"]["textpresso"]["token"]
    paperid_sentence_list = []
    for id in wbpids:
        txt = textpresso_paper_text(id, textpresso_token)
        #with open("texttt.txt", "w") as g:
          #g.write(str(txt))
        count_total_rows = len(txt)

        for current_i, row in enumerate(txt):
            # Limiting content inside "txt" variable to only "abstract"
            # and "main content" of paper not "reference" or 
            # "acknowledgments" sections, by checking if current row of
            # "txt" variable dont has any of the following things which 
            # is usually found in "reference" or "acknowledgments" sections
            # If it finds any of these words to be in the current "row",
            # it take it to reached the "reference" or "acknowledgments"
            # section and thus stops any further iteration through "rows"
            # of "txt" variable by using "break".
            # To avoid breaking/stopping process if these words comes in 
            # main content, it uses "current_i > count_total_rows / 3"
            if any(
                    word in row.lower().split() for word in ("literature", "cited",
                    "supported", "references", "thank", "acknowledge", 
                    "acknowledgments", "thank", "contributions")  
            ):
                if current_i > count_total_rows / 3:
                    break

            # Remove sentence if number of characters (not words) are less 
            # than 40 or there is only 1 word in "row"
            if len(row) < 40 or len(row.split()) == 1:
                continue

            # Cleans stop words from current sentence
            clean_row = ""
            word_bank = row.split()
            for word in word_bank:
                if word.lower() not in stop_words:
                    clean_row+=word+" "

            # Remove sentence with links and email ids
            if (re.search("\S+@\S+\.", clean_row) 
                or re.search("www\.\S+\.", clean_row) 
                or re.search("http.?://", clean_row)
            ):
                continue

            # Remove special characters like "@", except those present in "all_special_chars" list
            for c in clean_row:
                if (not c.isalnum() and not c == " ") and c not in all_special_chars:
                    clean_row = clean_row.replace(c, "")
                    clean_row = re.sub(' +', ' ', clean_row)

            paperid_sentence_list.append((id, clean_row))
    return paperid_sentence_list 


if __name__ == "__main__":
    start_time = time.time()
    settings = setSettings()
    paper_ids = ["WBPaper00050743"]
    paperid_sentence_list = get_paper_sentences_with_TE(paper_ids, settings)
    with open("text_final.txt", "w") as f:
        f.write(str(paperid_sentence_list))
    end_time = time.time()
    print(end_time-start_time)

