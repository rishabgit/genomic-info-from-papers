import abc
import logging
from typing import List, Tuple

logger = logging.getLogger(__name__)


class AbstractEntityExtractor(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def __init__(self):
        pass

    @abc.abstractmethod
    def extract(self, text: str) -> List[Tuple[str, str]]:
        pass
