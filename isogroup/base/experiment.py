import pandas as pd
from isogroup.base.database import Database
from isogroup.base.feature import Feature
from isogroup.base.cluster import Cluster
from isogroup.base.misc import Misc


class Experiment:

    def __init__(self, dataset: pd.DataFrame, database: 'Database' = None, tracer=None):
        self.dataset = dataset
        self.database = database
        self.samples: dict = {} 
        self._mz_tol: None | float = None
        self._rt_tol: None | float = None

        self._tracer: None | str = tracer
        self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer)
        self.clusters = None
        

    @property
    def rt_tol(self):
        return self._rt_tol

    @property
    def tracer(self):
        return self._tracer

    @property
    def tracer_element(self):
        return self._tracer_element

    @property
    def mz_tol(self):
        return self._mz_tol
    

    def initialize_experimental_features(self):
        """
        Initialize the experimental features from the dataset
        """
        for idx, _ in self.dataset.iterrows():
            mz = idx[0]
            rt = idx[1]
            id = idx[2]

            # Extract the intensity for each sample in the dataset
            for sample in self.dataset.columns:
                if sample not in ["mz", "rt", "id"]:
                    intensity = self.dataset.loc[idx, sample]

                    # Initialize the experimental features for each sample
                    feature = Feature(
                        rt=rt, mz=mz, tracer=self.tracer,
                        feature_id=id, 
                        intensity=intensity,
                        sample=sample
                        )
                    
                    # Add the feature in the list corresponding to the sample
                    if sample not in self.samples:
                        self.samples[sample] = {}
                    self.samples[sample][id] = feature


    def annotate_features(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        Calculate the mz error and the rt error
        """

        for sample in self.samples.values():
            for feature in sample.values():
                    
                for db_feature in self.database.features:

                    # Calculate the exact mz and rt errors
                    mz_error = (db_feature.mz - feature.mz)
                    rt_error = (db_feature.rt - feature.rt)

                    # Covert mz_error to ppm 
                    mz_error = (mz_error / feature.mz) * 1e6

                    # Check if the experimental feature is within tolerance
                    if abs(mz_error) <= mz_tol and abs(rt_error) <= rt_tol:
                        feature.chemical.append(db_feature.chemical[0])
                        feature.isotopologue.append(db_feature.isotopologue[0])
                        feature.metabolite.append(db_feature.chemical[0].label)
                        feature.formula.append(db_feature.chemical[0].formula)
                        feature.mz_error.append(mz_error)
                        feature.rt_error.append(rt_error)

        self._mz_tol = mz_tol
        self._rt_tol = rt_tol


    def annotate_experiment(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        MultiIndex DataFrame
        """
        # Initialize the experimental features from the dataset
        self.initialize_experimental_features()

        # Annotate the experimental features
        self.annotate_features(mz_tol, rt_tol)



    def export_features(self, filename = None, sample_name = None):
        """
        Create a DataFrame to summarize the annotated data
        Optionnal: Export the DataFrame to a tsv file if a filename is provided with samples in column
        Optionnal: Export the Dataframe of only one sample if a sample name is provided
        """

        # Create a DataFrame to summarize the experimental features
        feature_data = []
        for sample in self.samples.values():
            for feature in sample.values():
                feature_data.append({
                    "feature_id": feature.feature_id,
                    "mz": feature.mz,
                    "rt": feature.rt,
                    "metabolite": feature.metabolite,
                    "isotopologue": feature.isotopologue,
                    "mz_error": feature.mz_error,
                    "rt_error": feature.rt_error,
                    "sample": feature.sample,
                    "intensity": feature.intensity
                })

        # Create a DataFrame to summarize the annotated data
        df = pd.DataFrame(feature_data)

        # Export the DataFrame to a tsv file if a filename is provided
        if filename:
            df.to_csv(filename, sep="\t", index=False)

            # Export the Dataframe of only one sample if a sample name is provided
            if sample_name:
                df = df[df["sample"] == sample_name] # Filter the DataFrame by sample name
                df.to_csv(filename, sep="\t", index=False)

        return df

    def clusterize(self):
        """
        Create unique clusters from annotated features based on their names.
        """
        cluster_names = []
        
        # Group features by metabolite
        for sample in self.samples.values():
            for feature in sample.values():
                cluster_names += feature.metabolite

        cluster_names = set(cluster_names)

        # Create unique clusters
        self.clusters = {}

        for sample in self.samples.keys():
            self.clusters[sample] = {}
            for i,c in enumerate(cluster_names):
                features = self.get_features_from_name(c, sample)

                # Sort features by isotopologues
                features.sort(key=lambda f: f.isotopologue)

                # Assign the cluster_id to the features in the cluster
                for feature in features:
                    if not hasattr(feature, "in_cluster") or feature.in_cluster is None:
                        feature.in_cluster = [] 
                    feature.in_cluster.append(f"C{i}")  

                self.clusters[sample][c] = Cluster(features=features, cluster_id=f"C{i}", name=c)


    def get_features_from_name(self, name, sample_name:str):
        """
        Get a feature from the experiment by its name, in a given sample if provided
        """
        features = []
        for feature in self.samples[sample_name].values():
            if name in feature.metabolite:
                features.append(feature)
        return features

    
    def export_clusters(self, filename = None, sample_name = None):
        """
        Create a DataFrame to summarize the annotated clusters
        Optionnal: Export the DataFrame to a tsv file if a filename is provided
        Optionnal: Export the Dataframe of only one sample if a sample name is provided
        """
        
        # Check if the sample name is in the DataFrame
        all_samples = list(self.samples.keys())
        if sample_name is not None:
            if sample_name not in all_samples:
                raise ValueError(f"Sample {sample_name} not found in annotated clusters. Available samples: {', '.join(all_samples)}")
        
        cluster_data = []
        for sample, clusters in self.clusters.items():
            if sample_name is None or sample_name == sample: # Filter the DataFrame by sample name if provided
                for cname, cluster in clusters.items():
                    for feature in cluster.features:
                        idx = [i for i,j in enumerate(feature.metabolite) if j == cname][0]
                        # Get the cluster_id of the features in another cluster
                        other_clusters = [c.cluster_id for cluster_name, c in clusters.items() if feature in c.features and c.cluster_id != cluster.cluster_id]
                        cluster_data.append({
                            "cluster_id": cluster.cluster_id,
                            "metabolite": cluster.name,
                            "feature_id": feature.feature_id,
                            "mz": feature.mz,
                            "rt": feature.rt,
                            "feature_potential_metabolite": feature.metabolite,
                            "isotopologue": feature.isotopologue[idx],
                            "mz_error": feature.mz_error[idx],
                            "rt_error": feature.rt_error[idx],
                            "sample": feature.sample,
                            "intensity": feature.intensity,
                            "status": cluster.status,
                            "missing_isotopologue": cluster.missing_isotopologues,
                            "duplicated_isotopologue": cluster.duplicated_isotopologues,
                            # "in_cluster": feature.in_cluster,
                            "in_another_cluster": other_clusters
                        })

        # Create a DataFrame to summarize the annotated clusters
        df = pd.DataFrame(cluster_data)

        # Export the DataFrame to a tsv file if a filename is provided
        if filename:
            df.to_csv(filename, sep="\t", index=False)

        return df


    def get_clusters_from_name(self, name, sample_name:str):
        """
        Get a cluster from the experiment by its name, in a given sample if provided
        """
        for cluster in self.clusters[sample_name].values():
            if cluster.name == name:
                return cluster
        return None
    

    def clusters_summary(self, filename = None):
        """
        Export a tsv file with a summary of the clusters
        """
        # List to store the cluster summary data
        cluster_summary = []
        cluster_id_unique = set() # To store unique cluster_id

        for sample, clusters in self.clusters.items():
            for cluster in clusters.values():

                # Check if the cluster_id is unique
                if cluster.cluster_id not in cluster_id_unique:
                    cluster_id_unique.add(cluster.cluster_id)

                    summary = cluster.cluster_summary

                    # Retrieve the samples in which the cluster is present
                    samples_in_cluster = {sample for sample, clusters in self.clusters.items() if cluster.cluster_id in [c.cluster_summary["cluster_id"] for c in clusters.values()]}
                    summary["samples"] = len(samples_in_cluster)

                    cluster_summary.append(summary)

        # Create a DataFrame with the collected information
        df = pd.DataFrame(cluster_summary)

        # Export the DataFrame to a tsv file if a filename is provided
        if filename:
            df.to_csv(filename, sep="\t", index=False)

        return df
