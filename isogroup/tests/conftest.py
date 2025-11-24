"""
Variable configuration for testing
"""

import pytest
from decimal import Decimal as D
import pandas as pd

@pytest.fixture
def dataset_df():
    return pd.DataFrame(
        {'id': ['F1', 'F2', 'F3',  'F4',  'F5', 'F6', 'F7', 'F8', 'F9'], 
         'mz': [119.025753, 120.0291332, 191.0191775654, 119.0232843, 137.0275004, 136.024129, 135.0208168, 134.0174803, 133.0140851], 
         'rt': [667.779067, 667.9255408, 679.9930235, 678.1606593, 676.4604364, 676.5620229, 676.6045604, 676.8898827, 676.8952154], 
         'Sample_1': [1571414706.0, 1059554882.0, 31398195.78, 0.0, 529223407.9, 2090662547.0, 3105587268.0, 2077278842.0, 543216118.8], 
         'Sample_2': [266171108.6, 129533534.2, 5324316.124, 0.0, 28994270.58, 97127965.25, 154077393.8, 218743897.0, 155940888.7]}
    )

@pytest.fixture
def database_df():
    return pd.DataFrame(
        {'metabolite': ['Fumarate', 'Succinate', 'Citrate', 'Isocitrate', 'Malate', 'a-KG', 'G6P', 'ADP'], 
         'rt': [979, 668, 680, 680, 676, 883, 890, 2050], 
         'formula': ['C4H4O4', 'C4H6O4', 'C6H8O7', 'C6H8O7', 'C4H6O5', 'C5H6O5', 'C6H13O9P', 'C10H15N5O10P2'], 
         'charge': [-1, -1, -1, -1, -1, -1, -1, -1]})