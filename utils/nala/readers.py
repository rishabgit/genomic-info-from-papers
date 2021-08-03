import abc
from bs4 import BeautifulSoup
import re
import glob
import csv
import os
import xml.etree.ElementTree as ET
import warnings

from utils.nala.data import Dataset, Document, Part, Entity


class Reader:
    """
    Abstract class for reading in a dataset in some format.
    Subclasses that inherit this class should:
    * Be named [Name]Reader
    * Implement the abstract method read that returns an object of type Dataset
    """

    @abc.abstractmethod
    def read(self):
        """
        :returns: nalaf.structures.data.Dataset
        """
        return


class HTMLReader(Reader):
    """
    Reader for a local file system of tagtog plain.html files.

    It reads either a single file or a directory (the contained .html's)
    """

    def __init__(self, path):
        self.path = path
        """an html file or a directory containing .html files"""
        self.whole_basename_as_docid = None

    def __read_directory_localfs(self):
        dataset = Dataset()

        filenames = glob.glob(str(self.path + "/**/*.html"), recursive=True) + glob.glob(str(self.path + "/**/*.xml"), recursive=True)
        for filename in filenames:
            dataset = self.__read_file_path_localfs(filename, dataset)

        return dataset

    def __read_file_path_localfs(self, filename, dataset=None):
        if dataset is None:
            dataset = Dataset()

        with open(filename, 'rb') as a_file:
            HTMLReader.read_file(a_file, filename, dataset, self.whole_basename_as_docid)

        return dataset


    @staticmethod
    def read_file(a_file, filename, dataset=None, whole_basename_as_docid=False):
        if dataset is None:
            dataset = Dataset()

        soup = BeautifulSoup(a_file, "html.parser")
        document = Document()

        for part in soup.find_all(id=re.compile('^s')):
            if re.match(r'^s[3-9]', part['id']):
                is_abstract = False
            else:
                is_abstract = True
            document.parts[part['id']] = Part(str(part.string), is_abstract=is_abstract)

        doc_id = os.path.basename(filename).replace('.plain.html', '').replace('.html', '').replace('.xml', '')
        if not whole_basename_as_docid and '-' in doc_id:
            doc_id = doc_id.split('-')[-1]

        dataset.documents[doc_id] = document

        return dataset

    def read(self):
        if os.path.isdir(self.path):
            return self.__read_directory_localfs()