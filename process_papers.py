import argparse
import logging

from hybrid_extraction import findVariants
from refine import refine
from settings import setSettings


logger = logging.getLogger(__name__)


def load_ids(paper_ids):
    with open(paper_ids) as f:
        return [line.strip() for line in f.read().splitlines() if line != '']
    return []


def main():
    parser = argparse.ArgumentParser(description='Extract genomic location info from papers')
    parser.add_argument('-l', '--log-file', metavar='log_file', dest='log_file', type=str, default=None,
                        help='path to the log file to generate. Default ./afp_pipeline.log')
    parser.add_argument('-L', '--log-level', dest='log_level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR',
                                                                        'CRITICAL'], default='INFO',
                        help='set the logging level')
    parser.add_argument('-i', '--paper-ids', metavar='paper_ids', dest='paper_ids', type=str,
                        help='text file containing paper ids, one per line')
    parser.add_argument('-o', '--output_file', metavar='output_file', dest='output_file', type=str,
                        help='output file path')
    parser.add_argument('-m', '--method', metavar='method', dest='method',  choices=['textpresso', 'wbtools'],
                        help='textpresso or wbtools')
    args = parser.parse_args()
    logging.basicConfig(filename=args.log_file, level=args.log_level,
                        format='%(asctime)s - %(name)s - %(levelname)s:%(message)s')

    logger.info('init settings')
    settings = setSettings()
    logger.info('finished init settings')
    ids = load_ids(args.paper_ids)
    variants = findVariants(settings, ids, args.method)
    df = refine(variants)
    df.to_csv(args.output_file, index=False, encoding='utf-8')


if __name__ == '__main__':
    settings = setSettings()
    #import numpy as np
    #ids_to_extract = np.load('data/top100.npy').tolist()[90:]
    ids_to_extract = load_ids('data/exampleInput.txt')
    variants = findVariants(settings, ids_to_extract, 'textpresso')
    df = refine(variants)
    df.to_csv('output.csv', index=False, encoding='utf-8')

