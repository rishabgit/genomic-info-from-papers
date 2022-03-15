
from datetime import datetime

from hybrid_extraction import findVariants
from refine import refine
from settings import setSettings
from textpresso import wbtools_get_papers_last_month


def main():
    settings = setSettings()
    paperIds = wbtools_get_papers_last_month(settings, day=datetime.now())
    variants = findVariants(settings, paperIds)
    df = refine(variants)
    df.to_csv("output.csv", index=False, encoding='utf-8')


if __name__ == "__main__":
    main()
