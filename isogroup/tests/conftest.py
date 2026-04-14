"""
Variable configuration for testing
"""

import pytest
from isogroup.base.feature import Feature
import pandas as pd

@pytest.fixture
def dataset_df():
    """
    Dataset for testing IsoGroup (Targeted and Untargeted experiment)
    """
    return pd.DataFrame(
        {'id': ['F1', 'F2', 'F3',  'F4',  'F5', 'F6', 'F7', 'F8', 'F9'], 
         'mz': [119.025753, 120.0291332, 191.0191775654, 119.0232843, 137.0275004, 136.024129, 135.0208168, 134.0174803, 133.0140851], 
         'rt': [667.779067, 667.9255408, 679.9930235, 678.1606593, 676.4604364, 676.5620229, 676.6045604, 676.8898827, 676.8952154], 
         'Sample_1': [1571414706.0, 1059554882.0, 31398195.78, 0.0, 529223407.9, 2090662547.0, 3105587268.0, 2077278842.0, 543216118.8], 
         'Sample_2': [266171108.6, 129533534.2, 5324316.124, 0.0, 28994270.58, 97127965.25, 154077393.8, 218743897.0, 155940888.7]}
    )

@pytest.fixture
def database_df():
    """
    Database for testing IsoGroup (Targeted experiment)
    """
    return pd.DataFrame(
        {'metabolite': ['Fumarate', 'Succinate', 'Citrate', 'Isocitrate', 'Malate', 'a-KG', 'G6P', 'ADP'], 
         'rt': [979, 668, 680, 680, 676, 883, 890, 2050], 
         'formula': ['C4H4O4', 'C4H6O4', 'C6H8O7', 'C6H8O7', 'C4H6O5', 'C5H6O5', 'C6H13O9P', 'C10H15N5O10P2'], 
         'charge': [-1, -1, -1, -1, -1, -1, -1, -1]})

@pytest.fixture
def features_dict():
    """
    """
    return {'Sample_1': {'F1': Feature(feature_id = "F1", rt=667.779067, tracer="13C", tracer_element="C", mz=119.025753, intensity=1571414706.0, sample="Sample_1"), 
                        'F2': Feature(feature_id = "F2", rt=667.9255408, tracer="13C", tracer_element="C", mz=120.0291332, intensity=1059554882.0, sample="Sample_1"), 
                        'F3': Feature(feature_id = "F3", rt=679.9930235, tracer="13C", tracer_element="C", mz=191.0191775654, intensity=31398195.78, sample="Sample_1"), 
                        'F4': Feature(feature_id = "F4", rt=678.1606593, tracer="13C", tracer_element="C", mz=119.0232843, intensity=0.0, sample="Sample_1"), 
                        'F5': Feature(feature_id = "F5", rt=676.4604364, tracer="13C", tracer_element="C", mz=137.0275004, intensity=529223407.9, sample="Sample_1"), 
                        'F6': Feature(feature_id = "F6", rt=676.5620229, tracer="13C", tracer_element="C", mz=136.024129, intensity=2090662547.0, sample="Sample_1"), 
                        'F7': Feature(feature_id = "F7", rt=676.6045604, tracer="13C", tracer_element="C", mz=135.0208168, intensity=3105587268.0, sample="Sample_1"), 
                        'F8': Feature(feature_id = "F8", rt=676.8898827, tracer="13C", tracer_element="C", mz=134.0174803, intensity=2077278842.0, sample="Sample_1"), 
                        'F9': Feature(feature_id = "F9", rt=676.8952154, tracer="13C", tracer_element="C", mz=133.0140851, intensity=543216118.8, sample="Sample_1")}, 
            'Sample_2': {'F1': Feature(feature_id = "F1", rt=667.779067, tracer="13C", tracer_element="C", mz=119.025753, intensity=266171108.6, sample="Sample_2"), 
                        'F2': Feature(feature_id = "F2", rt=667.9255408, tracer="13C", tracer_element="C", mz=120.0291332, intensity=129533534.2, sample="Sample_2"), 
                        'F3': Feature(feature_id = "F3", rt=679.9930235, tracer="13C", tracer_element="C", mz=191.0191775654, intensity=5324316.124, sample="Sample_2"), 
                        'F4': Feature(feature_id = "F4", rt=678.1606593, tracer="13C", tracer_element="C", mz=119.0232843, intensity=0.0, sample="Sample_2"), 
                        'F5': Feature(feature_id = "F5", rt=676.4604364, tracer="13C", tracer_element="C", mz=137.0275004, intensity=28994270.58, sample="Sample_2"), 
                        'F6': Feature(feature_id = "F6", rt=676.5620229, tracer="13C", tracer_element="C", mz=136.024129, intensity=97127965.25, sample="Sample_2"), 
                        'F7': Feature(feature_id = "F7", rt=676.6045604, tracer="13C", tracer_element="C", mz=135.0208168, intensity=154077393.8, sample="Sample_2"), 
                        'F8': Feature(feature_id = "F8", rt=676.8898827, tracer="13C", tracer_element="C", mz=134.0174803, intensity=218743897.0, sample="Sample_2"), 
                        'F9': Feature(feature_id = "F9", rt=676.8952154, tracer="13C", tracer_element="C", mz=133.0140851, intensity=155940888.7, sample="Sample_2")}}

@pytest.fixture
def dataset_df_duplicates():
    """
    Dataset with duplicates for testing IsoGroup (Untargeted experiment)
    """

    return pd.DataFrame(
        {'id': ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11'], 
         'mz': [119.025753, 120.0291332, 191.0191775654, 119.0232843, 137.0275004, 136.024129, 137.0275104, 136.024229, 135.0208168, 134.0174803, 133.0140851], 
         'rt': [667.779067, 667.9255408, 679.9930235, 678.1606593, 676.4604364, 676.5620229, 676.4604364, 676.5620229, 676.6045604, 676.8898827, 676.8952154], 
         'Sample_1': [1571414706.0, 1059554882.0, 31398195.78, 0.0, 529223407.9, 2090662547.0, 529223407.9, 2090662547.0, 3105587268.0, 2077278842.0, 543216118.8], 
         'Sample_2': [266171108.6, 129533534.2, 5324316.124, 0.0, 28994270.58, 97127965.25, 28994270.58, 97127965.25, 154077393.8, 218743897.0, 155940888.7]}
    )
