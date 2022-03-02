import numpy as np


def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))


def extra_info_block(custom_mut_extract, sentence, span_size=150):
    info_and_snippets = []

    # look for gene and variant combo
    info_and_snippets += custom_mut_extract.var_and_gene_close(sentence, span_size=span_size)

    if info_and_snippets:
        info_and_snippets = unique_rows(info_and_snippets).tolist()
    return info_and_snippets


def regex_block(settings, sentence, span_size=150):
    mf_mut_extract = settings['mf_mut_extract']
    tmvar_mut_extract = settings['tmvar_mut_extract']
    custom_mut_extract = settings['custom_mut_extract']
    bow_mut_extract = settings['bow_mut_extract']
    mut_and_snippets = []

    # MutationFinder
    mut_and_snippets += mf_mut_extract(sentence, span_size=span_size)
    # tmVar
    mut_and_snippets += tmvar_mut_extract(sentence, span_size=span_size)
    # Custom patterns
    mut_and_snippets += custom_mut_extract(sentence, span_size=span_size)
    # Bag of words
    mut_and_snippets += bow_mut_extract(sentence)

    if mut_and_snippets:
        mut_and_snippets = unique_rows(mut_and_snippets).tolist()
    return mut_and_snippets


def point_mut_block(settings, sentence, span_size=150):
    mf_mut_extract = settings['mf_mut_extract']
    tmvar_mut_extract = settings['tmvar_mut_extract']
    custom_mut_extract = settings['custom_mut_extract']
    mut_and_snippets = []
    # MutationFinder
    mut_and_snippets += mf_mut_extract(sentence, span_size=span_size)
    # tmVar
    mut_and_snippets += tmvar_mut_extract(sentence, span_size=span_size)
    # Custom patterns
    mut_and_snippets += custom_mut_extract(sentence, span_size=span_size)

    if mut_and_snippets:
        mut_and_snippets = np.array(mut_and_snippets)
        mut_and_snippets = mut_and_snippets[:, 0].tolist()
        mut_and_snippets = list(set(mut_and_snippets))
    return mut_and_snippets
