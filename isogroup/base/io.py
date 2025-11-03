from isogroup.base.database import Database
from isogroup.base.feature import Feature
from isogroup.base.misc import Misc

import pandas as pd
from pathlib import Path
import logging



class IoHandler:
    """
    Handles input and output operations for both targeted and untargeted experiments.

    Args:
        dataset (str): Path to the dataset file containing raw mass spectrometry data.
        tracer (str): Tracer code used in the experiment (e.g. "13C").
        mz_tol (float|None): m/z tolerance (in ppm).
        rt_tol (float|None): Retention time tolerance (in sec).
        database (str|None): Path to the database file.
        rt_window (float|None): rt tolerance for clustering (e.g. "10") (Untargeted case).
        ppm_tol (float|None): m/z tolerance in ppm (Untargeted case).
        max_atoms (int|None): Maximum number of tracer atoms to consider in a molecule (e.g. "20").
        outputs_path (str|None): Path to the output directory.    
    """
    def __init__(self, dataset, tracer, mz_tol, rt_tol, database=None, max_atoms=None, outputs_path=None): 
        self.dataset = dataset
        self.dataset_name = None
        self.database = database
        self.tracer = tracer
        self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer)
        self.mz_tol = mz_tol
        self.rt_tol = rt_tol
        self.max_atoms = max_atoms
        self.outputs_path = outputs_path
        
        self.features = {} # {sample_name: {feature_id: Feature object}}

    
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

    # # @property
    # # def tracer_element(self):
    # #     """
    # #     Returns the tracer element used in the experiment.
    # #     :return: str | None
    # #     """
    # #     return self._tracer_element

    # @property
    # def mz_tol(self):
    #     """
    #     Returns the m/z tolerance used for feature annotation.
    #     :return: float | None
    #     """
    #     return self._mz_tol
    
    
    def read_dataset(self):
        """
        Reads the dataset from the specified file path and loads it into a pandas DataFrame.
        """

        inputdata = Path(self.dataset)

        if not inputdata.exists():
            raise FileNotFoundError(f"File {self.dataset} does not exist")

        self.dataset_name = inputdata.stem
        self.dataset = pd.read_csv(inputdata, sep="\t").set_index(["mz", "rt", "id"])

        # logging.info(f"Dataset loaded from {inputdata} with shape {data.shape}")    
    
    def read_database(self):
        """
        Reads the database from the specified file path and initializes a Database object.
        """
        if not isinstance(self.database, str):
            raise FileNotFoundError(f"File {self.database} does not exist")
        
        # Load database
        db_data = pd.read_csv(self.database, sep=";")
        self.database = Database(dataset=db_data, tracer=self.tracer)


    def initialize_experimental_features(self):
        """
        Initialize Feature objects from the dataset and organize them by sample.
        Each feature is created with its retention time, m/z, tracer, intensity, and sample name.
        Populates `self.features` as a dictionary of the form:
        {sample_name: {feature_id: Feature object}}

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
                    if sample not in self.features:
                        self.features[sample] = {}
                    self.features[sample][id] = feature

    
    def create_output_directory(self):
        """
        Create an output directory for saving results.
        """
        # if self.outputs_path is None:
        #     raise ValueError("Output path is not set. Please read a dataset first.")
        
        res_dir = Path(f"{self.outputs_path}/{self.dataset_name}_res")
        res_dir.mkdir(parents=True, exist_ok=True)
        self.outputs_path = res_dir
        
        # logging.info(f"Results will be saved to: {self.outputs_path}")

   
    def export_annotated_features(self, sample_name = None):
        """
        Summarize and export annotated features into a DataFrame and export it to a tsv file.
        
        :param sample_name: Name of the sample to filter the DataFrame by, if provided
        """

        # Create a DataFrame to summarize the experimental features
        feature_data = []
        for sample in self.features.values():
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
        df.to_csv(f"{self.outputs_path}/{self.dataset_name}_features.tsv", sep="\t", index=False)

        # Export the Dataframe of only one sample if a sample name is provided
        if sample_name:
            df = df[df["sample"] == sample_name] # Filter the DataFrame by sample name
            df.to_csv(f"{self.outputs_path}/{self.dataset_name}_features.tsv", sep="\t", index=False)
        

        # return df

    
    def export_unannotated_features(self):
        """
        Export all features to a TSV file (Untargeted case).
        
        """
        records = []
        for _, features in self.features.items():
            for f in features.values():
                # If not in any cluster, mark accordingly
                cluster_ids = f.in_cluster if f.in_cluster else ["None"]
                iso_labels = [f.cluster_isotopologue.get(cid, "N/A") for cid in cluster_ids]

                records.append({
                    "FeatureID": f.feature_id,
                    "RT": f.rt,
                    "m/z": f.mz,
                    "sample": f.sample,
                    "Intensity": f.intensity,
                    "InClusters": cluster_ids,
                    "Isotopologues": iso_labels
                })

        df = pd.DataFrame.from_records(records)
        df.to_csv(f"{self.outputs_path}/{self.dataset_name}_unannotated_features.tsv", sep="\t", index=False)

    
    def export_clusters(self, clusters_to_export, sample_name = None):
        """
        Summarize annotated clusters into a DataFrame and export it to a tsv file.
        :param clusters_to_export: dict containing clusters to export
        :param sample_name: Name of the sample to filter the DataFrame by, if provided
        """
        
        # Check if the sample name is in the DataFrame
        all_samples = list(self.features.keys())
        if sample_name is not None:
            if sample_name not in all_samples:
                raise ValueError(f"Sample {sample_name} not found in annotated clusters. Available samples: {', '.join(all_samples)}")
        
        cluster_data = []
        for sample, clusters in clusters_to_export.items():
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
        # if filename:
        df.to_csv(f"{self.outputs_path}/{self.dataset_name}_annotated_clusters.tsv", sep="\t", index=False)

        # return df
    
    def clusters_summary(self, clusters_to_summarize):
        """
        Export a tsv file with a summary of the clusters
        :param clusters_to_summarize: dict containing clusters to summarize
        :return: pd.DataFrame with the summary of the clusters
        """
        # List to store the cluster summary data
        cluster_summary = []
        cluster_id_unique = set() # To store unique cluster_id

        for _, clusters in clusters_to_summarize.items():
            for cluster in clusters.values():

                # Check if the cluster_id is unique
                if cluster.cluster_id not in cluster_id_unique:
                    cluster_id_unique.add(cluster.cluster_id)

                    summary = cluster.cluster_summary

                    # Retrieve the samples in which the cluster is present
                    samples_in_cluster = {sample for sample, clusters in clusters_to_summarize.items() if cluster.cluster_id in [c.cluster_summary["cluster_id"] for c in clusters.values()]}
                    summary["samples"] = len(samples_in_cluster)

                    cluster_summary.append(summary)

        # Create a DataFrame with the collected information
        df = pd.DataFrame(cluster_summary)

        # Export the DataFrame to a tsv file if a filename is provided
        # if filename:
        df.to_csv(f"{self.outputs_path}/{self.dataset_name}_summary.tsv", sep="\t", index=False)

        # return df
    
    def clusters_to_dataframe(self, cluster_to_export):
        """
        Convert the clusters into a pandas DataFrame for easier analysis and export (Untargeted case).

        :param cluster_to_export: dict containing clusters to export
        """
        records = []
        for _, clusters in cluster_to_export.items():
            for cluster in clusters.values():
                sorted_features = sorted(cluster.features, key=lambda f: f.mz)

                for idx, f in enumerate(sorted_features):
                    iso_label = f.cluster_isotopologue.get(cluster.cluster_id, "Mx")
                    records.append({
                        "ClusterID": cluster.cluster_id,
                        "FeatureID": f.feature_id,
                        "RT": f.rt,
                        "m/z": f.mz,
                        "sample": f.sample,
                        "Intensity": f.intensity,
                        "Isotopologue": iso_label,
                        "InClusters": f.in_cluster,
                        "AlsoIn": f.also_in
                    })

        df = pd.DataFrame.from_records(records)
        df.to_csv(f"{self.outputs_path}/{self.dataset_name}_unannotated_clusters.tsv", sep="\t", index=False)
        # return pd.DataFrame.from_records(records)


    # def export_clusters_to_tsv(self, filepath: str):
    #     """
    #     Export the clusters to a CSV file.
    #     :param filepath: str
    #     """
    #     df = self.clusters_to_dataframe()
    #     df.to_csv(filepath, sep="\t", index=False)


# if __name__ == "__main__":
#     test = IoHandler(dataset=r"C:\Users\kouakou\Documents\IsoGroup_test\data\dataset_test_XCMS.txt",
#                      tracer="13C",
#                     )
#     test.read_dataset()
#     # test.create_output_directory(r"C:\Users\kouakou\Documents\IsoGroup_test")
#     test.initialize_experimental_features()
#     print(test.outputs_path)
#     print(test.tracer)
#     print(test._tracer_element)
#     # test.export_annotated_features()
#     # print(test.samples)