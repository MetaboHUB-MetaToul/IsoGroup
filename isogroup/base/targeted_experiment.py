from __future__ import annotations
import pandas as pd
from isogroup.base.experiment import Experiment
from isogroup.base.cluster import Cluster
from isogroup.base.database import Database
import logging
import time

logger = logging.getLogger(f"IsoGroup")

class TargetedExperiment(Experiment):
    """
    Represents a targeted mass spectrometry experiment.
    Used to group and annotate detected features from an experimental dataset using a reference database with isotopic tracer information.
    """

    def __init__(self, dataset:pd.DataFrame, tracer:str, ppm_tol:float, rt_tol:float, database:pd.DataFrame):
        """
        :param dataset: DataFrame containing experimental data with columns for m/z, retention time (RT), feature ID and sample intensities.
        :param tracer: Tracer code used in the experiment (e.g. "13C").
        :param ppm_tol: m/z tolerance (in ppm).
        :param rt_tol: Retention time tolerance (in seconds).
        :param database: DataFrame containing theoretical features with columns retention time (RT), metabolite names, and formulas.
        """
        super().__init__(dataset = dataset, tracer=tracer, ppm_tol=ppm_tol, rt_tol=rt_tol, database=database)
        self.database = Database(dataset=database, 
                                 tracer=self._tracer,
                                 tracer_element=self.tracer_element)
        # self.ppm_tol = ppm_tol
        # self.rt_tol = rt_tol

        # self.tracer = tracer
        # self.cluster = cluster
        # self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer)

    def run_targeted_pipeline(self):
        """
        Run the full targeted annotation pipeline for the experiment.
        
        This includes:
        - Initializing Feature objects from the dataset.
        - Matching experimental features to the database within specified tolerances.
        - Clustering features by metabolite names.
        """
        start_time = time.time()
        
        self.initialize_experimental_features()
        self.annotate_features()
        self.clusterize()

        total_time = time.time() - start_time

        logger.info(f"Targeted grouping completed in {total_time:.2f} seconds.")

    def annotate_features(self):
        """
        Annotate experimental features by matching them with the database 
        features within specified m/z and retention time tolerances.

        """
        logger.info("Find matches between experimental features and database features...")
        
        nb_features_annotated = 0 

        for features_id in self.features.values():
            for feature in features_id.values():    
                for db_feature in self.database.theoretical_features:
                    # Calculate the exact mz and rt errors
                    mz_error = (db_feature.mz - feature.mz)
                    rt_error = (db_feature.rt - feature.rt)
                    # Covert mz_error to ppm 
                    mz_error = (mz_error / feature.mz) * 1e6

                    # Check if the experimental feature is within tolerance
                    if abs(mz_error) <= self.ppm_tol and abs(rt_error) <= self.rt_tol:
                        feature.chemical.append(db_feature.chemical[0])
                        # feature.isotopologue.append(db_feature.isotopologue[0])
                        feature.cluster_isotopologue[db_feature.chemical[0].label] = db_feature.cluster_isotopologue[db_feature.chemical[0].label]
                        feature.metabolite.append(db_feature.chemical[0].label)
                        feature.formula.append(db_feature.chemical[0].formula)
                        feature.mz_error.append(mz_error)
                        feature.rt_error.append(rt_error)
                        nb_features_annotated += 1
                        logger.debug(f"Feature {feature.feature_id} in sample {feature.sample} annotated with {db_feature.chemical[0].label} (isotopologue: {db_feature.cluster_isotopologue[db_feature.chemical[0].label]})")
                        logger.debug(f" - mz error (ppm): {mz_error}, rt error (sec): {rt_error}")
        
        logger.info(f"    => {nb_features_annotated} experimental features matched with database features.\n")
        

    def clusterize(self):
        """
        Group features by metabolite names within each sample and assign a unique cluster ID to each group.
        Populates `self.clusters` as a dictionary of the form:
        {sample_name: {cluster_id: Cluster object}}
        """
        # cluster_names = []
        
        # # Group features by metabolite
        # for sample in self.features.values():
        #     for feature in sample.values():
        #         cluster_names += feature.metabolite
        
        # # cluster_names = set(cluster_names)

         # Create unique clusters
        # # # self.clusters = {}
        logger.info("Grouping features by metabolite names...")
        
        cluster_names = []

        for _, features in self.features.items():
            for feature in features.values():
                cluster_names += [metabolite_name for metabolite_name in feature.metabolite 
                                  if metabolite_name not in cluster_names]

        for sample in self.features.keys():
            self.clusters[sample] = {}
            for index, clusters in enumerate(cluster_names):
                features = self.get_features_from_name(clusters, sample)
                
                # Sort features by isotopologues
                # features.sort(key=lambda f: f.isotopologue)
                features.sort(key=lambda f: f.cluster_isotopologue[clusters])
                # Assign the cluster_id to the features in the cluster
                for feature in features:
                    if not hasattr(feature, "in_cluster") or feature.in_cluster is None:
                        feature.in_cluster = [] 
                    feature.in_cluster.append(f"C{index}")  

                self.clusters[sample][clusters] = Cluster(features=features, cluster_id=f"C{index}", name=clusters)
                logger.debug(f"Cluster C{index} ({clusters}) identified with {len(features)} features in sample {sample}.")
                logger.debug(f"    {[features.feature_id for features in features]} ")
        
        logger.info(f"    => {len(cluster_names)} clusters identified.\n")
    
    def get_features_from_name(self, name:str, sample_name:str):
        """
        Retrieve all features in a given sample that are annotated with a specific metabolite name.

        :param name: Name of the metabolite for which to retrieve features
        :param sample_name: Name of the sample from which to retrieve features

        :return: List of Feature objects that match the metabolite name in the specified sample
        """
        features = []
        for feature in self.features[sample_name].values():
            if name in feature.metabolite:
                features.append(feature)
        return features

    def get_clusters_from_name(self, name, sample_name:str):
        """
        Get a cluster from the experiment by its name, in a given sample if provided

        :param name: Name of the cluster to retrieve
        :param sample_name: Name of the sample to retrieve the cluster from

        :return: Cluster object if found, None otherwise
        """
        for cluster in self.clusters[sample_name].values():
            if cluster.name == name:
                return cluster
        return None
    
# if __name__ == "__main__":
#     from isogroup.base.io import IoHandler
#     from isogroup.base.database import Database
#     from pathlib import Path
#     io= IoHandler()
#     data= io.read_dataset(Path(r"..\..\data\dataset_test_XCMS.txt"))
#     database = io.read_database(Path(r"..\..\data\database.csv"))
#     experiment = TargetedExperiment(data, tracer="13C", mz_tol=5, rt_tol=15, database=database)
    
#     experiment.run_targeted_pipeline()
#     for sample, clusters in experiment.clusters.items():
#         for cluster in clusters.values():
#             print(cluster.expected_isotopologues_in_cluster)
    
###############################################################################
        # @property
    # def rt_tol(self):
    #     """
    #     Returns the retention time tolerance used for feature annotation.
    #     :return: float        
    #     """
    #     return self._rt_tol

    # @property
    # def tracer(self):
    #     """
    #     Returns the tracer used for the experiment.
    #     :return: str | None
    #     """
    #     return self._tracer

    # @property
    # def tracer_element(self):
    #     """
    #     Returns the tracer element used in the experiment.
    #     :return: str | None
    #     """
    #     return self._tracer_element

    # @property
    # def mz_tol(self):
    #     """
    #     Returns the m/z tolerance used for feature annotation.
    #     :return: float | None
    #     """
    #     return self._mz_tol
    

    # def initialize_experimental_features(self):
    #     """
    #     Initialize Feature objects from the dataset and organize them by sample.
    #     Each feature is created with its retention time, m/z, tracer, intensity, and sample name.
    #     Populates `self.samples` as a dictionary of the form:
    #     {sample_name: {feature_id: Feature object}}
    #     """
    #     for idx, _ in self.dataset.iterrows():
    #         mz = idx[0]
    #         rt = idx[1]
    #         id = idx[2]

    #         # Extract the intensity for each sample in the dataset
    #         for sample in self.dataset.columns:
    #             if sample not in ["mz", "rt", "id"]:
    #                 intensity = self.dataset.loc[idx, sample]

    #                 # Initialize the experimental features for each sample
    #                 feature = Feature(
    #                     rt=rt, mz=mz, tracer=self.tracer,
    #                     feature_id=id, 
    #                     intensity=intensity,
    #                     sample=sample
    #                     )
                    
    #                 # Add the feature in the list corresponding to the sample
    #                 if sample not in self.samples:
    #                     self.samples[sample] = {}
    #                 self.samples[sample][id] = feature

    # def annotate_experiment(self, mz_tol, rt_tol):
    #     """
    #     Run the full annotation process for the experiment.
    #     This includes:
    #     - Initializing Feature objects from the dataset.
    #     - Matching experimental features to the database within specified tolerances.
    #     :param mz_tol: m/z tolerance in ppm
    #     :param rt_tol: retention time tolerance in seconds
    #     """
    #     # Initialize the experimental features from the dataset
    #     self.initialize_experimental_features()

    #     # Annotate the experimental features
    #     self.annotate_features(mz_tol, rt_tol)



    # def export_features(self, filename = None, sample_name = None):
    #     """
    #     Summarize annotated features into a DataFrame and optionally export it to a tsv file.
    #     :param filename: Name of the file to export the summary to
    #     :param sample_name: Name of the sample to filter the DataFrame by, if provided
    #     :return: pd.DataFrame with the summary of the annotated features
    #     """

    #     # Create a DataFrame to summarize the experimental features
    #     feature_data = []
    #     for sample in self.samples.values():
    #         for feature in sample.values():
    #             feature_data.append({
    #                 "feature_id": feature.feature_id,
    #                 "mz": feature.mz,
    #                 "rt": feature.rt,
    #                 "metabolite": feature.metabolite,
    #                 "isotopologue": feature.isotopologue,
    #                 "mz_error": feature.mz_error,
    #                 "rt_error": feature.rt_error,
    #                 "sample": feature.sample,
    #                 "intensity": feature.intensity
    #             })

    #     # Create a DataFrame to summarize the annotated data
    #     df = pd.DataFrame(feature_data)

    #     # Export the DataFrame to a tsv file if a filename is provided
    #     if filename:
    #         df.to_csv(filename, sep="\t", index=False)

    #         # Export the Dataframe of only one sample if a sample name is provided
    #         if sample_name:
    #             df = df[df["sample"] == sample_name] # Filter the DataFrame by sample name
    #             df.to_csv(filename, sep="\t", index=False)

    #     return df


    
    
    # def export_clusters(self, filename = None, sample_name = None):
    #     """
    #     Summarize annotated clusters into a DataFrame and optionally export it to a tsv file.
    #     :param filename: Name of the file to export the summary to
    #     :param sample_name: Name of the sample to filter the DataFrame by, if provided
    #     :return: pd.DataFrame with the summary of the annotated clusters
    #     """
        
    #     # Check if the sample name is in the DataFrame
    #     all_samples = list(self.samples.keys())
    #     if sample_name is not None:
    #         if sample_name not in all_samples:
    #             raise ValueError(f"Sample {sample_name} not found in annotated clusters. Available samples: {', '.join(all_samples)}")
        
    #     cluster_data = []
    #     for sample, clusters in self.clusters.items():
    #         if sample_name is None or sample_name == sample: # Filter the DataFrame by sample name if provided
    #             for cname, cluster in clusters.items():
    #                 for feature in cluster.features:
    #                     idx = [i for i,j in enumerate(feature.metabolite) if j == cname][0]
    #                     # Get the cluster_id of the features in another cluster
    #                     other_clusters = [c.cluster_id for cluster_name, c in clusters.items() if feature in c.features and c.cluster_id != cluster.cluster_id]
    #                     cluster_data.append({
    #                         "cluster_id": cluster.cluster_id,
    #                         "metabolite": cluster.name,
    #                         "feature_id": feature.feature_id,
    #                         "mz": feature.mz,
    #                         "rt": feature.rt,
    #                         "feature_potential_metabolite": feature.metabolite,
    #                         "isotopologue": feature.isotopologue[idx],
    #                         "mz_error": feature.mz_error[idx],
    #                         "rt_error": feature.rt_error[idx],
    #                         "sample": feature.sample,
    #                         "intensity": feature.intensity,
    #                         "status": cluster.status,
    #                         "missing_isotopologue": cluster.missing_isotopologues,
    #                         "duplicated_isotopologue": cluster.duplicated_isotopologues,
    #                         # "in_cluster": feature.in_cluster,
    #                         "in_another_cluster": other_clusters
    #                     })

    #     # Create a DataFrame to summarize the annotated clusters
    #     df = pd.DataFrame(cluster_data)

    #     # Export the DataFrame to a tsv file if a filename is provided
    #     if filename:
    #         df.to_csv(filename, sep="\t", index=False)

    #     return df


    # def clusters_summary(self, filename = None):
    #     """
    #     Export a tsv file with a summary of the clusters
    #     :param filename: Name of the file to export the summary to
    #     :return: pd.DataFrame with the summary of the clusters
    #     """
    #     # List to store the cluster summary data
    #     cluster_summary = []
    #     cluster_id_unique = set() # To store unique cluster_id

    #     for sample, clusters in self.clusters.items():
    #         for cluster in clusters.values():

    #             # Check if the cluster_id is unique
    #             if cluster.cluster_id not in cluster_id_unique:
    #                 cluster_id_unique.add(cluster.cluster_id)

    #                 summary = cluster.cluster_summary

    #                 # Retrieve the samples in which the cluster is present
    #                 samples_in_cluster = {sample for sample, clusters in self.clusters.items() if cluster.cluster_id in [c.cluster_summary["cluster_id"] for c in clusters.values()]}
    #                 summary["samples"] = len(samples_in_cluster)

    #                 cluster_summary.append(summary)

    #     # Create a DataFrame with the collected information
    #     df = pd.DataFrame(cluster_summary)

    #     # Export the DataFrame to a tsv file if a filename is provided
    #     if filename:
    #         df.to_csv(filename, sep="\t", index=False)

    #     return df
