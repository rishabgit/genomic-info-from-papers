import argparse
import logging

from hybrid_extraction import findVariants
from refine import refine
from settings import setSettings
from textpresso import wbtools_get_papers_last_month
from dateutil import parser as date_parser


logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Extract genomic location info from papers")
    parser.add_argument("-l", "--log-file", metavar="log_file", dest="log_file", type=str, default=None,
                        help="path to the log file to generate. Default ./afp_pipeline.log")
    parser.add_argument("-L", "--log-level", dest="log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR',
                                                                        'CRITICAL'], default="INFO",
                        help="set the logging level")
    parser.add_argument("-m", "--max-papers", metavar="max_papers", dest="max_papers", type=int, default=10,
                        help="maximum number of papers to process")
    parser.add_argument("-i", "--paper-ids", metavar="paper_ids", dest="paper_ids", type=str, nargs="+",
                        help="process the provided list of papers instead of reading them from db")
    parser.add_argument("-f", "--from-date", metavar="from_date", dest="from_date", default=None,
                        help="process papers added to WormBase or modified on or after the provided date")
    args = parser.parse_args()
    logging.basicConfig(filename=args.log_file, level=args.log_level,
                        format='%(asctime)s - %(name)s - %(levelname)s:%(message)s')

    logger.info("init settings")
    settings = setSettings()
    logger.info("finished init settings")
    if args.paper_ids:
        paper_ids = args.paper_ids
    else:
        logger.info("getting paper ids from DB")
        paper_ids = wbtools_get_papers_last_month(settings['db_config'], day=date_parser.parse(args.from_date),
                                                  max_num_papers=args.max_papers)
    wb_paper_ids = [f'WBPaper{id}' for id in paper_ids]
    variants = findVariants(settings, wb_paper_ids)
    df = refine(variants)
    df.to_csv("output.csv", index=False, encoding='utf-8')


if __name__ == "__main__":
    main()
