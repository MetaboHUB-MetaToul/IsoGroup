import pandas as pd
from isogroup.base.feature import Feature


class Sample:
    def __init__(self, dataset: pd.DataFrame, sample_name : str, sample_type = None):
        self.data = dataset 
        self.sample_name = sample_name # Name of the sample
        self.sample_type = sample_type # Type of the sample : control, TP, fully labelled, etc.
        self.features = [] # List of features in the sample

    def initialize_features(self, annotated_data):
        """
        Initialize the features for this sample from the annotated data.
        For each feature, create a new instance with the intensity of the sample.
        """
        pass

    def __repr__(self):
        """
        Return a string representation of the sample.
        :return: str
        """
        return f"Sample({self.sample_name}, {self.sample_type}, {self.features})"
    