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
        for feature in annotated_data:
            # Extract the intensity of the sample
            feature_intensity = {self.sample_name: feature.intensity[self.sample_name]}

            # Create a new feature for this sample (copy the feature and update the intensity)
            new_feature = Feature(
                rt=feature.rt,
                mz=feature.mz,
                feature_id=feature.feature_id,
                intensity=feature_intensity,
                formula=feature.formula,
                metabolite=feature.metabolite,
                isotopologue=feature.isotopologue,
                mz_error=feature.mz_error,
                rt_error=feature.rt_error
            )

            # Add the new feature to the sample
            self.features.append(new_feature)

    def __repr__(self):
        """
        Return a string representation of the sample.
        :return: str
        """
        return f"Sample({self.sample_name}, {self.sample_type}, {self.features})"


