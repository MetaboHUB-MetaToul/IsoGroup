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
        self.annotated_experiment: None | pd.DataFrame = None
        self.clusters_summary: None | pd.DataFrame = None # Summary of the clusters of annotated features in the experiment
        self.mz_tol: None | float = None
        self.rt_tol: None | float = None

        self.feature_exp: list = [] # List of experiment features
        self.annotated_clusters: list = []   # List of annotated clusters


    def annotate_experiment(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        And calculate the mz error and the rt error.
        MultiIndex DataFrame
        """
        annotated_data = []  # List to store annotated rows

        for idx, _ in self.experiment.iterrows():
            mz = idx[0]
            rt = idx[1]
            identity = idx[2]

            # Initialize the experiment features from the dataset
            feature_data = Feature(rt=rt, mz=mz, Fid=identity, intensity=None, metabolite=None, isotopologue=None, mz_error=None, rt_error=None) 

            # Initialize lists to store the metabolites, isotopologues, mz errors and rt errors
            metabolites, isotopologues, mz_errors, rt_errors = [], [], [], []
            
            for feature_db in self.database.features:
                mz_db = float(feature_db.mz)
                rt_db = float(feature_db.rt)

                # Calculate the exact mz and rt errors
                mz_error = (mz_db - feature_data.mz)
                rt_error = (rt_db - feature_data.rt)

                # Check if the feature is within tolerance
                if abs(mz_error) <= mz_tol and abs(rt_error) <= rt_tol:
                    
                    metabolites.append(feature_db.metabolite)
                    isotopologues.append(feature_db.isotopologue)
                    mz_errors.append(mz_error)
                    rt_errors.append(rt_error)
 
            annotated_data.append([mz, rt, identity, metabolites, isotopologues, mz_errors, rt_errors])
            
            # Store all annotations in feature_data
            feature_data.metabolite = metabolites
            feature_data.isotopologue = isotopologues
            feature_data.mz_error = mz_errors
            feature_data.rt_error = rt_errors
            self.feature_exp.append(feature_data)
                
        self.mz_tol = mz_tol
        self.rt_tol = rt_tol
        self.annotated_experiment = pd.DataFrame(
            annotated_data, 
            columns=["mz", "rt", "id", "metabolite", "isotopologues", "mz_error", "rt_error"]
        ).set_index(["mz", "rt", "id"])

        # Export the annotated experiment to a tsv file
        self.annotated_experiment.to_csv("annotated_experiment.tsv", sep="\t")


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
        if not self.feature_exp:
            raise ValueError("The experiment has not been annotated yet.")

        metabolite_clusters = {}

        # Iterate over feature_exp to retrieve all features with their annotations
        for feature in self.feature_exp:
            if feature.metabolite:  # Only consider annotated features
                key = tuple(sorted(feature.metabolite))  # Sort and convert to tuple to use as a dictionary key
                if key not in metabolite_clusters:
                    metabolite_clusters[key] = []
                metabolite_clusters[key].append(feature)

        # Create Cluster objects
        self.annotated_clusters = [Cluster(features) for features in metabolite_clusters.values()]

        # Assign a unique identifier Cid to each cluster
        for idx, cluster in enumerate(self.annotated_clusters):
            cluster.Cid = idx

        # Create a dataframe to summarize the annotated clusters
        cluster_data = []

        # Iterate over the annotated clusters
        for cluster in self.annotated_clusters:
            for feature in cluster.features:
                # Get the mz_error, rt_error and identity of the feature
                isotopologue = feature.isotopologue
                mz = feature.mz
                rt = feature.rt
                mz_error = feature.mz_error
                rt_error = feature.rt_error
                identity = feature.Fid

                # Append the feature data to the cluster_data list
                cluster_data.append([cluster.Cid, cluster.metabolite, isotopologue, identity, mz, rt, mz_error, rt_error])
            
        # Create a DataFrame to summarize the annotated clusters

        self.clusters_summary = pd.DataFrame(
            cluster_data,
            columns=["Cid", "metabolite", "isotopologues", "identity", "mz", "rt", "mz_error", "rt_error"]
        ).set_index(["mz", "rt", "identity"])

        # Export the cluster summary to a tsv file
        self.clusters_summary.to_csv("cluster_summary.tsv", sep="\t", index=False)
