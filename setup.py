import setuptools

with open("readme.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="genomicinfo",
    version="0.0.1",
    author="Rishab Mallick",
    author_email="",
    description="Extract genomic data from scientific articles",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rishabgit/genomic-info-from-papers",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'torch==1.9.0',
        'transformers==4.9.1',
        'nervaluate==0.1.8',
        'accelerate==0.3.0',
        'seqeval==1.2.2',
        'datasets==1.11.0',
        'bs4',
        'wbtools==1.2.0'
    ]
)