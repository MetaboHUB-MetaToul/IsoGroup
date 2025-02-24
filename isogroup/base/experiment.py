import pandas as pd
from isogroup.base.database import Database
from isogroup.base.sample import Sample
from isogroup.base.feature import Feature
from isogroup.base.cluster import Cluster

class Experiment:

    def __init__(self, dataset: pd.DataFrame, database: 'Database' = None):
        self.experiment = dataset
        self.database = database
        self.samples: dict = {} # Dictionary to store the samples
        self.annotated_df: None | pd.DataFrame = None
        self.clusters_summary: None | pd.DataFrame = None # Summary of the clusters of annotated features in the experiment
        self.mz_tol: None | float = None
        self.rt_tol: None | float = None

        self.annotated_data: list = [] # List of experimental features
        self.annotated_clusters: list = []   # List of annotated clusters


    def annotate_experiment(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        And calculate the mz error and the rt error.
        MultiIndex DataFrame
        """

        for idx, _ in self.experiment.iterrows():
            mz = idx[0]
            rt = idx[1]
            identity = idx[2]

            # Initialize the experiment features from the dataset
            exp_feature = Feature(
                rt=rt, mz=mz, 
                Fid=identity, 
                intensity=None, 
                metabolite=[], 
                isotopologue=[],
                mz_error=[], 
                rt_error=[]
                ) 
            
            for th_feature in self.database.features:
                mz_db = float(th_feature.mz)
                rt_db = float(th_feature.rt)

                # Calculate the exact mz and rt errors
                mz_error = (mz_db - exp_feature.mz)
                rt_error = (rt_db - exp_feature.rt)

                # Check if the feature is within tolerance
                if abs(mz_error) <= mz_tol and abs(rt_error) <= rt_tol:
                
                    exp_feature.metabolite.append(th_feature.metabolite)
                    exp_feature.isotopologue.append(th_feature.isotopologue)
                    exp_feature.mz_error.append(mz_error)
                    exp_feature.rt_error.append(rt_error)

            # Store all experimental features
            self.annotated_data.append(exp_feature)

        self.mz_tol = mz_tol
        self.rt_tol = rt_tol

        # Create a DataFrame to summarize the annotated data
        self.to_dataframe()


    def initialize_samples(self):
        """
        Initialize the samples from the dataset.
        :return:
        """
        data = self.annotated_experiment if (self.annotated_experiment is not None) else self.experiment
        for sample in data.columns:
            self.samples[sample] = Sample(dataset=data[[sample]], sample_type="test")

 

    def get_metabolite_clusters(self):
        """
        Create clusters of annotated features based on the metabolite.
        """
        # Check if the experiment has been annotated
        if not self.annotated_data:
            raise ValueError("The experiment has not been annotated yet.")

        metabolite_clusters = {}

        # Iterate over annotated_data to retrieve all features with their annotations
        for feature in self.annotated_data:
            if feature.metabolite:  # Only consider annotated features
                key = tuple(sorted(feature.metabolite))  # Sort and convert to tuple to use as a dictionary key
                if key not in metabolite_clusters:
                    metabolite_clusters[key] = []
                metabolite_clusters[key].append(feature)

        # Create Cluster objects
        self.annotated_clusters = [
            Cluster(features=features, cluster_id=idx) 
            for idx, features in enumerate(metabolite_clusters.values())
        ]

        # Create a dataframe to summarize the annotated clusters
        cluster_data = []

        # Iterate over the annotated clusters
        for cluster in self.annotated_clusters:
            for feature in cluster.features:
                # Append the feature data to the cluster_data list
                cluster_data.append(
                        [
                        cluster.cluster_id, 
                        cluster.metabolite, 
                        feature.isotopologue, 
                        feature.feature_id, 
                        feature.mz, 
                        feature.rt, 
                        feature.mz_error, 
                        feature.rt_error
                        ]
                    )
            
        # Create a DataFrame to summarize the annotated clusters

        self.clusters_summary = pd.DataFrame(
            cluster_data,
            columns=["cluster_id", "metabolite", "isotopologues", "identity", "mz", "rt", "mz_error", "rt_error"]
        ).set_index(["mz", "rt", "identity"])

        # Export the cluster summary to a tsv file
        self.clusters_summary.to_csv("cluster_summary.tsv", sep="\t", index=False)

# class AnnotationError(Exception):
#     pass


    def to_dataframe(self):
        """
        Create a DataFrame to summarize the annotated data
        """
        # Check if the experiment has been annotated
        if not self.annotated_data:
            raise ValueError("The experiment has not been annotated yet.")

        data = []

        # Iterate over the annotated data
        for feature in self.annotated_data:
            data.append(
                [
                    feature.feature_id, 
                    feature.metabolite, 
                    feature.isotopologue, 
                    feature.mz, 
                    feature.rt, 
                    feature.mz_error, 
                    feature.rt_error
                ]
            )

        # Create a DataFrame to summarize the annotated data
        self.annotated_df = pd.DataFrame(
            data,
            columns=["identity", "metabolite", "isotopologue", "mz", "rt", "mz_error", "rt_error"]
        ).set_index(["mz", "rt", "identity"])

        # Export the annotated data to a tsv file
        self.annotated_df.to_csv("annotated_data.tsv", sep="\t", index=True)




