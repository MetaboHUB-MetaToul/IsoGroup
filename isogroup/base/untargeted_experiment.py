import pandas as pd
import numpy as np
import bisect
from isogroup.base.feature import Feature
from isogroup.base.cluster import Cluster
from isogroup.base.misc import Misc


class UntargetedExperiment:
    """
    Represents an untargeted mass spectrometry experiment.
    An untargeted experiment contains a collection of features and clusters, along with associated metadata.

    Args:
        dataset (pd.DataFrame|None): DataFrame containing raw mass spectrometry data.
        tracer (str|None): Tracer code used in the experiment (e.g. "13C").
    """

    def __init__(self, dataset: pd.DataFrame|None = None, tracer: str|None = None):
        """
        Initialize the untargeted experiment with raw data and tracer information.
        """
        self.dataset = dataset

        self._tracer = tracer
        self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer) if tracer is not None else (None, None)
        self._RTwindow: float|None = None
        self._ppm_tolerance: float|None = None
        self.mzshift_tracer = float(Misc.calculate_mzshift(self._tracer)) if tracer is not None else None

        self.features: dict = {}  # {sample_name: {feature_id: Feature object}}
        self.clusters: dict = {}  # {sample_name: [Cluster objects]}

    @property
    def RTwindow(self) -> float|None:
        """
        Returns the retention time tolerance (RT window) used to group features into clusters.
        :return: float|None
        """
        return self._RTwindow
    
    @property
    def ppm_tolerance(self) -> float|None:
        """
        Returns the mass-to-charge ratio (m/z) tolerance in parts per million (ppm) used to group features into clusters.
        :return: float|None
        """
        return self._ppm_tolerance
    
    @property
    def tracer(self) -> str|None:
        """
        Returns the tracer code used in the experiment.
        :return: str|None
        """
        return self._tracer
    
    @property
    def tracer_element(self) -> str|None:
        """
        Returns the tracer element extracted from the tracer code.
        :return: str|None
        """
        return self._tracer_element
    
    def build_final_clusters(self, RTwindow: float, ppm_tolerance: float):
        """
        Build and deduplicate clusters of features based on retention time and m/z tolerances.
        This is a convenience method that combines the initialization of features,
        building of clusters, and deduplication into a single step.
        :param RTwindow: float
        :param ppm_tolerance: float
        :return: None
        """
        self.initialize_experimental_features()
        self.build_clusters(RTwindow=RTwindow, ppm_tolerance=ppm_tolerance)
        self.deduplicate_clusters()


    def initialize_experimental_features(self):
        """
        Initialize Feature objects from the dataset and organize them by sample.
        Each feature is created with its retention time, m/z, tracer, intensity, and sample name.
        Populates `self.samples` as a dictionary of the form:
        {sample_name: {feature_id: Feature object}}
        """
        if self.dataset is None:
            raise ValueError("Dataset is not provided.")
        
        # required_columns = {"mz", "rt", "id"}
        # if not required_columns.issubset(self.dataset.columns):
        #     raise ValueError(f"Dataset must contain the following columns: {required_columns}")
        
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
                    if sample not in self.features:   ## Interrogation : features ou samples pour le dictionnaire ? en ciblé : sample
                        self.features[sample] = {}
                    self.features[sample][id] = feature


    def build_clusters(self, RTwindow: float, ppm_tolerance: float):
        """
        Build clusters of features based on retention time and m/z tolerances.
        Features within the specified RT window and ppm tolerance are grouped into clusters.
        Populates `self.clusters` with Cluster objects.
        :param RTwindow: float
        :param ppm_tolerance: float
        """
        self._RTwindow = RTwindow
        self._ppm_tolerance = ppm_tolerance

        if not self.features:
            raise ValueError("Features are not initialized. Please run 'initialize_experimental_features()' first.")

        self.clusters = {}

        for sample_name, features in self.features.items():
            all_features = sorted(features.values(), key=lambda f: f.rt)
            rts = [f.rt for f in all_features]

            clusters = {}
            cluster_id = 0

            # For each feature, find potential isotopologues within the RT window
            for i, base_feature in enumerate(all_features):
                left_bound = bisect.bisect_left(rts, base_feature.rt - RTwindow)
                right_bound = bisect.bisect_right(rts, base_feature.rt + RTwindow)

                candidates = all_features[left_bound:right_bound]
                group = {base_feature}

                for candidate in candidates:
                    if candidate == base_feature:
                        continue

                    iso_number = round((candidate.mz - base_feature.mz) / self.mzshift_tracer)
                    max_iso = Misc.get_max_isotopologues_for_mz(base_feature.mz, self.tracer_element)
                    
                    if abs(iso_number) > max_iso:
                        continue

                    expected_mz = base_feature.mz + iso_number * self.mzshift_tracer
                    if abs(expected_mz - candidate.mz) <= (expected_mz * ppm_tolerance / 1e6):
                        group.add(candidate)

                if len(group) > 1:
                    cluster_id = f"C{cluster_id}"
                    group_sorted = sorted(list(group), key=lambda f: f.mz)
                    min_mz = group_sorted[0].mz

                    for f in group_sorted:
                        f.cluster_isotopologue[cluster_id] = round((f.mz - min_mz) / self.mzshift_tracer) # Specific to clusters

                        if cluster_id not in f.in_cluster:
                            f.in_cluster.append(cluster_id)
                    clusters[cluster_id] = Cluster(cluster_id=cluster_id, features=group_sorted)
                    cluster_id = int(cluster_id[1:]) + 1

            self.clusters[sample_name] = clusters


    def deduplicate_clusters(self):
        """
        Deduplicate clusters to remove redundancies after initial clustering.
        Identical clusters (same features) are merged into a single cluster.
        Updates `self.clusters` to contain only unique clusters with new cluster IDs and updates features' cluster memberships.
        """
        unique_clusters = {}
        new_cluster_id = 0
        signature_to_cid = {} # Map cluster signatures to new cluster IDs

        for sample_name, clusters in self.clusters.items():
            unique_clusters[sample_name] = {}
            for cluster in clusters.values():

                # Create a signature for the cluster based on its feature IDs
                signature = frozenset(f.feature_id for f in cluster.features)

                if signature in signature_to_cid:
                    cid = signature_to_cid[signature]
                else:
                    cid = f"C{new_cluster_id}"
                    signature_to_cid[signature] = cid
                    new_cluster_id += 1

                unique_clusters[sample_name][cid] = Cluster(cluster_id=cid, features=list(cluster.features))

                # Construct mapping of features to their new cluster IDs
        features_to_clusters = {}
        for sample_clusters in unique_clusters.values():
            for cluster in sample_clusters.values():
                for f in cluster.features:
                    features_to_clusters.setdefault(f.feature_id, set()).add(cluster.cluster_id)

                # Update each feature's cluster membership
        for sample_clusters in unique_clusters.values():
            for cluster in sample_clusters.values():
                cluster.features.sort(key=lambda f: f.mz)
                min_mz = cluster.features[0].mz

                for f in cluster.features:
                    f.in_cluster = list(features_to_clusters.get(f.feature_id, []))
                    f.also_in = [cid for cid in f.in_cluster if cid != cluster.cluster_id]
                    f.cluster_isotopologue[cluster.cluster_id] = round((f.mz - min_mz) / self.mzshift_tracer) # Specific to clusters

        self.clusters = unique_clusters


    def clusters_to_dataframe(self) -> pd.DataFrame:
        """
        Convert the clusters into a pandas DataFrame for easier analysis and export.
        Each row in the DataFrame represents a feature within a cluster, along with relevant metadata.
        Missing isotopologues are represented as empty rows (Nan)
        :return: pd.DataFrame
        """
        records = []
        for sample_name, clusters in self.clusters.items():
            for cluster in clusters.values():
                min_mz = cluster.lowest_mz
                max_mz = cluster.highest_mz
                max_iso = round((max_mz - min_mz) / self.mzshift_tracer)

                iso_map = {f.cluster_isotopologue.get(cluster.cluster_id, None): f for f in cluster.features}
                for iso in range(max_iso + 1):
                    if iso in iso_map and iso_map[iso] is not None:
                        feature = iso_map[iso]
                        record = {
                            "ClusterID": cluster.cluster_id,
                            "FeatureID": feature.feature_id,
                            "RT": feature.rt,
                            "m/z": feature.mz,
                            "sample": feature.sample,
                            "Intensity": feature.intensity,
                            #"Isotopologue": feature.isotopologue,
                            "Isotopologue": iso,
                            "InClusters": feature.in_cluster,
                            "AlsoIn": feature.also_in
                        }
                    else:
                        record = {
                            "ClusterID": cluster.cluster_id,
                            "FeatureID": None,
                            "RT": None,
                            "m/z": None,
                            "sample": sample_name,
                            "Intensity": None,
                            "Isotopologue": iso,
                            "InClusters": [],
                            "AlsoIn": []
                        }

                    records.append(record)
        
        df = pd.DataFrame.from_records(records)
        df.fillna(np.nan, inplace=True)
        return df
    

    def export_clusters_to_tsv(self, filepath: str):
        """
        Export the clusters to a CSV file.
        :param filepath: str
        """
        df = self.clusters_to_dataframe()
        df.to_csv(filepath, sep="\t", index=False)
                    


