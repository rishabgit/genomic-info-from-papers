from sklearn.metrics import precision_recall_fscore_support


def precision(mf_mut_extract, tmvar_mut_extract, custom_mut_extract, bow_mut_extract):
    ## # 3 Testing
    ## Using the IDP4+ dataset downloaded and setup in notebook 1.

    nala_db = pd.read_csv('data/nala/nala_binary.csv').to_numpy()


    ## ### 3.1 RegEx


    y_true = []
    y_pred = []
    print('Total sentences to process: ', len(nala_db))
    for i, row in enumerate(nala_db):
        if (i+1) % 500 == 0: print(f"{i+1}", end = " ")
        sentence = row[0]
        true = row[1]
        if regex_block(settings, sentence):
            pred = 1
        else:
            pred = 0
        y_true.append(true)
        y_pred.append(pred)

    precision_recall_fscore_support(y_true, y_pred, average='binary')
    # 3.3 RegEx (custom) + NER

    #y_true = []
    #y_pred = []
    #start_time = time.time()
    #print('Total sentences to process: ', len(nala_db))
    #for i, row in enumerate(nala_db):
        #if (i+1) % 500 == 0: print(f"{i+1}", end = " ")
        #if (i+1) % 1000 == 0:
            #print('Time for 1000 lines:', int((time.time() - start_time)/60))
            #start_time = time.time()
        #sentence = row[0]
        #true = row[1]
        #if regex_block(sentence) or ner_mutations(sentence):
            #pred = 1
        #else:
            #pred = 0
        #y_true.append(true)
        #y_pred.append(pred)

    #precision_recall_fscore_support(y_true, y_pred, average='binary')


    # 3.2 BioBERT NER
    y_true = []
    y_pred = []
    print('Total sentences to process: ', len(nala_db))
    for i, row in enumerate(nala_db):
        if (i+1) % 500 == 0: print(f"{i+1}", end = " ")
        sentence = row[0]
        true = row[1]
        if ner_mutations(nala_ner, sentence):
            pred = 1
        else:
            pred = 0
        y_true.append(true)
        y_pred.append(pred)

    precision_recall_fscore_support(y_true, y_pred, average='binary')
    # taking too long to execute - not tested
