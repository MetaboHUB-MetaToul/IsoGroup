# TODO : Refactor all export methods

from isogroup.base.database import Database
import pandas as pd
from pathlib import Path



class IoHandler:
    """
    Handles input and output operations.
 
    """

    def __init__(self):
        self.dataset_path = None
        self.dataset_name = None
        self.database_path = None
        self.outputs_path = None

    def read_dataset(self, dataset):
        """
        Reads the dataset from the specified file path and loads it into a pandas DataFrame.
        """
        self.dataset_path = dataset

        if not self.dataset_path.exists():
            raise FileNotFoundError(f"File {self.dataset_path} does not exist.")
        
        self.dataset_name = self.dataset_path.stem
        
        return pd.read_csv(self.dataset_path, sep="\t")
         
        # logging.info(f"Dataset loaded from {inputdata} with shape {data.shape}")    
    
    def read_database(self, database):
        """
        Reads the database from the specified file path and loads it into a pandas DataFrame.
        """
        self.database_path = database

        if not self.database_path.exists():
            raise FileNotFoundError(f"File {self.database_path} does not exist.")
        
        return pd.read_csv(self.database_path, sep=";")
    

    def export_theoretical_database(self, database: Database):
        """
        Summarize theoretical features into a DataFrame and optionally export it to a tsv file.
        """

        # Create a DataFrame to summarize the theoretical features
        feature_data = {
            "mz": [],
            "rt": [],
            "metabolite": [],
            "isotopologue": [],
            "formula": []
        }
        for feature in database.theoretical_features:
            feature_data["mz"].append(feature.mz)
            feature_data["rt"].append(feature.rt)
            feature_data["metabolite"].append(', '.join(feature.metabolite))
            feature_data["isotopologue"].append(', '.join(map(str, feature.isotopologue)))
            feature_data["formula"].append(feature.formula)
       
        pd.DataFrame.from_dict(feature_data).to_csv(Path(f"{self.outputs_path}/{self.dataset_name}_theoretical_db.tsv"), 
                                          sep="\t", 
                                          index=False)

    def create_output_directory(self, outputs_path):
        """
        Create an output directory for saving results.
        """
        res_dir = Path(f"{outputs_path}/{self.dataset_name}_res")
        res_dir.mkdir(parents=True, exist_ok=True)
        self.outputs_path = res_dir
        
        # logging.info(f"Results will be saved to: {self.outputs_path}")

   
    def targ_export_features(self, features_to_export, sample_name = None):
        """
        Summarize and export annotated features into a DataFrame and export it to a tsv file.
        
        :param sample_name: Name of the sample to filter the DataFrame by, if provided
        """

        # Create a DataFrame to summarize the experimental features
        feature_data = []
        for sample in features_to_export.values():
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

    def targ_export_clusters(self, features, clusters_to_export, sample_name = None):
        """
        Summarize annotated clusters into a DataFrame and export it to a tsv file.
        :param clusters_to_export: dict containing clusters to export
        :param sample_name: Name of the sample to filter the DataFrame by, if provided
        """
        
        # Check if the sample name is in the DataFrame
        all_samples = list(features.keys())
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
        df.to_csv(f"{self.outputs_path}/{self.dataset_name}_clusters.tsv", sep="\t", index=False)

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

    def untarg_export_features(self, features_to_export):
        """
        Export all features to a TSV file (Untargeted case).
        
        """
        records = []
        for _, features in features_to_export.items():
            for f in features.values():
                # If not in any cluster, mark accordingly
                cluster_ids = f.in_cluster if f.in_cluster else ["None"]
                # iso_labels = [f.cluster_isotopologue.get(cid, "N/A") for cid in cluster_ids]

                records.append({
                    "FeatureID": f.feature_id,
                    "RT": f.rt,
                    "m/z": f.mz,
                    "sample": f.sample,
                    "Intensity": f.intensity,
                    "InClusters": cluster_ids,
                    "Isotopologues": [f.cluster_isotopologue.get(cid, "N/A") for cid in cluster_ids]
                })

        df = pd.DataFrame(records)
        df.to_csv(f"{self.outputs_path}/{self.dataset_name}_features.tsv", sep="\t", index=False)

    def untarg_export_clusters(self, cluster_to_export):
        """
        Convert the clusters into a pandas DataFrame for easier analysis and export (Untargeted case).

        :param cluster_to_export: dict containing clusters to export
        """
        records = []
        for _, clusters in cluster_to_export.items():
            for cluster in clusters.values():
                sorted_features = sorted(cluster.features, key=lambda f: f.mz)

                for _, f in enumerate(sorted_features):
                    # iso_label = f.cluster_isotopologue.get(cluster.cluster_id, "Mx")
                    records.append({
                        "ClusterID": cluster.cluster_id,
                        "FeatureID": f.feature_id,
                        "RT": f.rt,
                        "m/z": f.mz,
                        "sample": f.sample,
                        "Intensity": f.intensity,
                        "Isotopologue": f.cluster_isotopologue[cluster.cluster_id],
                        # "InClusters": f.in_cluster,
                        "AlsoIn": f.also_in[cluster.cluster_id]
                    })

        df = pd.DataFrame(records)
        df.to_csv(f"{self.outputs_path}/{self.dataset_name}_clusters.tsv", sep="\t", index=False)
        # return pd.DataFrame.from_records(records)


    # def export_clusters_to_tsv(self, filepath: str):
    #     """
    #     Export the clusters to a CSV file.
    #     :param filepath: str
    #     """
    #     df = self.clusters_to_dataframe()
    #     df.to_csv(filepath, sep="\t", index=False)


# if __name__ == "__main__":
#     test = IoHandler(
#                     )
    # print(test.read_dataset(r"C:\Users\kouakou\Documents\IsoGroup_test\data\dataset_test_XCMS.txt"))

    # print(test.outputs_path)
    # print(test.tracer)
    # print(test._tracer_element)
    # test.export_annotated_features()
    # print(test.samples)
