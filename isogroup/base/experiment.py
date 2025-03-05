import pandas as pd
from isogroup.base.database import Database
from isogroup.base.feature import Feature
from isogroup.base.cluster import Cluster

class Experiment:

    def __init__(self, dataset: pd.DataFrame, database: 'Database' = None):
        self.dataset = dataset
        self.database = database
        self.samples: dict = {} # Dictionary to store the samples
        self.mz_tol: None | float = None
        self.rt_tol: None | float = None
        self.tracer: None | str = None
        self.tracer_element: None | str = None
        self.experimental_features : list = [] # List of experimental features
        self.annotated_features: list = [] # List of experimental features after annotation # Inutile ? Modification de l'objet feature directement dans la liste experimental_features
        self.annotated_clusters: list = []   # List of annotated clusters 


    def initialize_experimental_features(self):
        """
        Initialize the experimental features from the dataset
        """
        for idx, _ in self.dataset.iterrows():
            mz = idx[0]
            rt = idx[1]
            identity = idx[2]

            # Extract the intensity for each sample in the dataset
            for sample in self.dataset.columns:
                if sample not in ["mz", "rt", "identity"]:
                    intensity = self.dataset.loc[idx, sample]

                    # Initialize the experimental features for each sample
                    feature = Feature(
                        rt=rt, mz=mz, 
                        feature_id=identity, 
                        intensity=intensity,
                        metabolite=[],
                        name=[],
                        isotopologue=[],
                        mz_error=[],
                        rt_error=[],
                        sample=sample
                        ) 
                    
                    # Add the feature in the list corresponding to the sample
                    if sample not in self.samples:
                        self.samples[sample] = []
                    self.samples[sample].append(feature)
                    
                    # Store all experimental features
                    self.experimental_features.append(feature)


    def annotate_features(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        Calculate the mz error and the rt error
        """
        for feature in self.experimental_features:
            mz = feature.mz
            rt = feature.rt

            for th_feature in self.database.features:
                mz_db = float(th_feature.mz)
                rt_db = float(th_feature.rt)

                # Calculate the exact mz and rt errors
                mz_error = (mz_db - mz)
                rt_error = (rt_db - rt)

                # Check if the experimental feature is within tolerance
                if abs(mz_error) <= mz_tol and abs(rt_error) <= rt_tol:
                    feature.metabolite.append(th_feature.metabolite)
                    feature.isotopologue.append(th_feature.isotopologue)
                    feature.name.append(th_feature.metabolite.label)
                    feature.mz_error.append(mz_error)
                    feature.rt_error.append(rt_error)

            # Store all experimental features after annotation
            self.annotated_features.append(feature)  # Inutile ? Modification de l'objet feature directement dans la liste experimental_features, sample mis à jour

        self.mz_tol = mz_tol
        self.rt_tol = rt_tol


    def annotate_experiment(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        MultiIndex DataFrame
        """
        # Initialize the experimental features from the dataset
        self.initialize_experimental_features()

        # Annotate the experimental features
        self.annotate_features(mz_tol, rt_tol)





