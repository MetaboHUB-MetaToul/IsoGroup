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

        self.annotated_experiment: None | pd.DataFrame = None # Rename


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
                        annotation=[],
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
                    feature.annotation.append(th_feature.metabolite.label)
                    feature.mz_error.append(mz_error)
                    feature.rt_error.append(rt_error)

            # Store all experimental features after annotation
            self.annotated_features.append(feature)  # Inutile ? Modification de l'objet feature directement dans la liste experimental_features, sample mis à jour

        self.mz_tol = mz_tol
        self.rt_tol = rt_tol


    def annotate_experiment(self, mz_tol, rt_tol, tracer, tracer_element):
        """
        Annotate the experiment features with the database within a given tolerance
        MultiIndex DataFrame
        """
        # Initialize the experimental features from the dataset
        self.initialize_experimental_features()

        # Annotate the experimental features
        self.annotate_features(mz_tol, rt_tol)

        self.tracer = tracer
        self.tracer_element = tracer_element


    def features_summary(self, filename = None, sample_name = None):
        """
        Create a DataFrame to summarize the annotated data
        Optionnal: Export the DataFrame to a tsv file if a filename is provided
        Optionnal: Export the Dataframe of only one sample if a sample name is provided
        """

        # Create a DataFrame to summarize the experimental features
        df = pd.DataFrame([vars(f) for f in self.experimental_features])
        df = df.drop(columns=["is_adduct", "in_cluster"], errors="ignore")

        self.annotated_experiment = df

        # Export the DataFrame to a csv file if a filename is provided
        if filename:
            df = df.drop(columns=["metabolite"], errors="ignore")
            df.to_csv(filename, sep="\t", index=False)

        # Export the Dataframe of only one sample if a sample name is provided
        if filename and sample_name:
            # Check if the sample name is in the DataFrame
            if sample_name not in df["sample"].unique():
                raise ValueError(f"The sample {sample_name} is not in the DataFrame")
            
            df = df[df["sample"] == sample_name]
            df = df.drop(columns=["metabolite"], errors="ignore")
            df.to_csv(filename, sep="\t", index=False)

        return df

    
    def get_annotated_clusters(self):
        """
        Create clusters by metabolite from the annotated features
        Duplicate the features if they have multiple annotations to create a cluster for each annotation
        """
        # Check if the experiment has been annotated
        if not self.annotated_features:
            raise ValueError("The experiment has not been annotated yet")
        
        # Create a dictionary to store the clusters
        clusters = {}
        for feature in self.experimental_features:
            for metabolite_name in feature.annotation: # A feature can have multiple annotations
                if metabolite_name not in clusters:
                    clusters[metabolite_name] = []
                clusters[metabolite_name].append(feature)
                
        # Create a Cluster object for each cluster and add a cluster id
        self.annotated_clusters = [
            Cluster(features, cluster_id=i, cluster_annotation=metabolite_name, tracer=self.tracer, tracer_element=self.tracer_element)
            for i, (metabolite_name, features) in enumerate(clusters.items(), start=0)
            ]

        # self.annotated_clusters = [
        #     Cluster(features, cluster_id=i)
        #     for i, (metabolite_name, features) in enumerate(clusters.items(), start=0)
        #     ]


    def cluster_summary(self, filename = None, sample_name = None):
        """
        Create a DataFrame to summarize the annotated clusters
        Optionnal: Export the DataFrame to a tsv file if a filename is provided
        Optionnal: Export the Dataframe of only one sample if a sample name is provided
        """
        # Check if the experiment has been annotated
        if not self.annotated_clusters:
            raise ValueError("No annotated clusters found. Run get_annotated_clusters() first")
        
        # Check if the sample name is in the DataFrame
        all_samples = {f.sample for c in self.annotated_clusters for f in c.features}
        if sample_name and sample_name not in all_samples:
            raise ValueError(f"Sample {sample_name} not found in annotated clusters. Available samples: {', '.join(all_samples)}")
        
        cluster_data = [] # List to store the cluster data

        for cluster in self.annotated_clusters:
            for feature in cluster.features:
            
                # Apply filter on the sample name if provided
                if sample_name and feature.sample != sample_name:
                    continue

                cluster_data.append({
                    "cluster_id": cluster.cluster_id,
                    "annotation": cluster.cluster_annotation,
                    "feature_id": feature.feature_id,
                    "mz": feature.mz,
                    "rt": feature.rt,
                    "feature_possible_annotation": feature.annotation,
                    "isolopologue": feature.isotopologue,
                    "mz_error": feature.mz_error,
                    "rt_error": feature.rt_error,
                    "sample": feature.sample,
                    "intensity": feature.intensity,
                    "is_complete": cluster.is_complete,
                    "missing_isotopologue": cluster.missing_isotopologue
                })

        # Create a DataFrame to summarize the annotated clusters
        df = pd.DataFrame(cluster_data)

        # Check if the DataFrame is empty
        # if df.empty:
        #     raise ValueError("No data found for the provided sample") 

        # Export the DataFrame to a tsv file if a filename is provided
        if filename:
            df.to_csv(filename, sep="\t", index=False)

        return df



    def feature_by_id(self):
        """
        Return all the features with the same feature_id
        """
        pass

    def feature_by_sample(self):
        """
        Return all the features for a specific sample
        """
        pass