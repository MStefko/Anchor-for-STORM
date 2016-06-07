import numpy as np
import pandas as pd
import logging, sys


class STORM:
    """Encapsulator for STORM data. Not much functionality for now. """
    def __init__(self, filename = None):
        if filename:
            self._loadLocalizationData(filename)
        else:
            self.data = None

    def _loadLocalizationData(self, filename):
        self.data = pd.read_csv(filename)
        return

