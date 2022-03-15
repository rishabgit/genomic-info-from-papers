
from datetime import datetime

from hybrid_extraction import findVariants
from refine import refine
from settings import setSettings
from textpresso import wbtools_get_papers_last_month


def main():
    settings = setSettings()
    paper_ids = wbtools_get_papers_last_month(settings['db_config'], day=datetime.now())
    # feb_day = datetime.strptime('2022-02-15', '%Y-%m-%d')
    # paper_ids = wbtools_get_papers_last_month(settings['db_config'], day=feb_day)
    wb_paper_ids = [f'WBPaper{id}' for id in paper_ids]
    variants = findVariants(settings, wb_paper_ids)
    df = refine(variants)
    df.to_csv("output.csv", index=False, encoding='utf-8')


if __name__ == "__main__":
    main()
