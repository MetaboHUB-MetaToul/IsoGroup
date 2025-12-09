from __future__ import annotations
from isogroup.base.experiment import Experiment
import bisect
from collections import defaultdict
from isogroup.base.cluster import Cluster
from isogroup.base.misc import Misc
import logging
import time
from datetime import datetime

class UntargetedExperiment(Experiment):
    """
    Represents an untargeted mass spectrometry experiment.
    An untargeted experiment contains a collection of features and clusters, along with associated metadata.

    """

    def __init__(self, dataset, tracer:str, mz_tol:float, rt_tol:float, max_atoms:int = None, #  keep_best_candidate: bool = False, #  keep_richest: bool = False,
                 log_file: str = "untargeted_experiment_log.txt"):
        """
        :param tracer: Tracer code used in the experiment (e.g. "13C").
        :param mz_tol: m/z tolerance in ppm.
        :param rt_tol: Retention time tolerance in seconds.
        :param max_atoms: Maximum number of tracer atoms to consider for isotopologues. If None, it will be determined based on m/z.
        """

        super().__init__(dataset= dataset, tracer=tracer, mz_tol=mz_tol, rt_tol=rt_tol, max_atoms=max_atoms)

        # self.dataset = dataset
        # self.features = features
        self.log_file = log_file

        # self.tracer = tracer
        # self._tracer_element, self._tracer_idx = tracer_element, tracer_idx
        # self.RTwindow = rt_window
        # self.ppm_tolerance = ppm_tolerance
        # self.max_atoms = max_atoms
        self.mzshift_tracer = float(Misc.calculate_mzshift(self.tracer)) 

        # self.keep_best_candidate = keep_best_candidate
        # self.keep_richest = keep_richest

        self.unclustered_features: dict = {}  # {sample_name: [Feature objects]}

        # --- Set up logging ---
        self.log_file = log_file
        logging.basicConfig(filename=self.log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger("IsoGroup.UntargetedExperiment")
        self.logger.info(f"Tracer: {self.tracer}, Tracer element: {self.tracer_element}, m/z shift: {self.mzshift_tracer}")


    def build_final_clusters(self, keep_best_candidate: bool = False, keep_richest: bool = True, verbose: bool = False):
        """
        Complete pipeline to build and deduplicate clusters from the dataset with logging and timing.
        Parameters:
            RTwindow (float): Retention time window for clustering.
            ppm_tolerance (float): m/z tolerance in parts per million for clustering.
            max_atoms (int|None): Maximum number of tracer atoms to consider for isotopologues. If None, it will be determined based on m/z.
            keep_best_candidate (bool): If True, retain only the feature with the highest intensity for each isotopologue within a cluster.
            keep_richest (bool): If True, retain only the largest cluster when multiple clusters share features.
        """
        start_time = time.time()
        start_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.logger.info(f"Starting untargeted clustering pipeline")


        # --- Initialization of features ---
        print(" Initializing features...", end=" ", flush=True)
        # t0 = time.time()
        # self.initialize_experimental_features()
        features_count = len(next(iter(self.features.values())))
        nb_samples = len(self.features)
        print(f" done ({features_count} features per sample)")
        self.logger.info(f"Initialized {features_count} features for {nb_samples} samples")


        # --- Construction of clusters ---
        print(" Building clusters without filtration...", end=" ", flush=True)
        # t0 = time.time()
        self.build_clusters(self.rt_tol, self.mz_tol, self.max_atoms)
        clusters_count = len(next(iter(self.clusters.values())))  
        print(f" done ({clusters_count} clusters per sample)")

        # --- Deduplication and cleaning of clusters ---
        print(" Cleaning clusters...", end=" ", flush=True)
        # t0 = time.time()
        merged, subset_removed, final, unclustered = self.deduplicate_clusters(keep_best_candidate=keep_best_candidate, keep_richest=keep_richest)
        print(f"→ {merged} merged, {subset_removed} subsets removed, {final} final clusters remained/sample")
        self.logger.info(
            f"Deduplication completed: merged clusters={merged}, removed subsets={subset_removed}, final cleaned clusters={final}, unclustered features={unclustered}"
        )

        total_time = time.time() - start_time
        print(f"[IsoGroup] Untargeted clustering completed in {total_time:.2f} seconds.")
        self.logger.info(f"Pipeline completed in {total_time:.2f} seconds.")

        # --- Verbose logging to file ---
        if verbose:
            summary = [
                ("Start Time", start_dt),
                ("Tracer", self.tracer),
                ("Number of samples", nb_samples),
                ("Features/sample", features_count),
                ("RT window (s)", self.RTwindow),
                ("m/z tolerance (ppm)", self.ppm_tolerance),
                ("Clusters before cleaning", clusters_count),
                ("Clusters merged", merged),
                ("Subset clusters removed", subset_removed),
                ("Final isotopic clusters/sample", final),
                ("Unclustered features", unclustered),
                ("Total time (s)", f"{total_time:.2f}")
            ]
            with open(self.log_file, "a") as f:
                f.write("\n" + "=" * 80 + "\nUntargeted Isotopic Clustering Summary\n" + "=" * 80 + "\n")
                for key, value in summary:
                    f.write(f"{key}: {value}\n")

        # TODO: Add Dataset name in summary
        # TODO: Erase previous log file content if any

    def build_clusters(self, RTwindow: float, ppm_tolerance: float, max_atoms: int = None):
        """
        Group features into potential isotopologue clusters based on retention time proximity and m/z differences.
        Parameters:
            RTwindow (float): Retention time window for clustering.
            ppm_tolerance (float): m/z tolerance in parts per million for clustering.
            max_atoms (int|None): Maximum number of tracer atoms to consider for isotopologues. If None, it will be determined based on m/z.
        """
        # self._RTwindow = RTwindow
        # self._ppm_tolerance = ppm_tolerance

        if not self.features:
            raise ValueError("Features not initialized.")

        # self.clusters = {}

        for sample_name, features in self.features.items():
            all_features = sorted(features.values(), key=lambda f: f.rt)
            
            rts = [f.rt for f in all_features]
    
            clusters = {}
            cluster_id_local = 0
        
            # For each feature, find potential isotopologues within the RT window
            for base_feature in all_features:
                # --- Find candidates within the RT window ---
                left_bound = bisect.bisect_left(rts, base_feature.rt - RTwindow)
                right_bound = bisect.bisect_right(rts, base_feature.rt + RTwindow)

                candidates = all_features[left_bound:right_bound]
                potential_group = {base_feature}
                # --- Identification of candidates for isotopologues ---
                for candidate in candidates:
                    if candidate == base_feature:
                        continue
                    
                    # iso_index = round((candidate.mz - base_feature.mz) / self.mzshift_tracer)
                    iso_index = Misc.calculate_isotopologue_index(candidate.mz, base_feature.mz, self.mzshift_tracer)
                    # Define a maximum number of tracer atoms if specified
                    max_iso = Misc.get_max_isotopologues_for_mz(base_feature.mz, self.tracer_element) if max_atoms is None else max_atoms

                    if abs(iso_index) > max_iso:
                        continue
                    
                    expected_mz = base_feature.mz + iso_index * self.mzshift_tracer
                    delta_ppm = abs(expected_mz - candidate.mz) / expected_mz * 1e6

                    if delta_ppm <= ppm_tolerance:
                        potential_group.add(candidate)                
                
                # --- If a group of isotopologues is found, create a cluster ---
                if len(potential_group) > 1:
                    cluster_id = f"C{cluster_id_local}"
                    group_sorted = sorted(list(potential_group), key=lambda f: f.mz)
                
                    for f in group_sorted:
                        # iso_index = round((f.mz - group_sorted[0].mz) / self.mzshift_tracer) 
                        iso_index = Misc.calculate_isotopologue_index(f.mz, group_sorted[0].mz, self.mzshift_tracer) # Theoretical isotopologue index
                        iso_label_tmp = "Mx" if iso_index == 0 else f"M+{iso_index}"
                
                        f.cluster_isotopologue[cluster_id] = iso_label_tmp # Specific to clusters
                        # print(f.cluster_isotopologue)
                        if cluster_id not in f.in_cluster:
                            f.in_cluster.append(cluster_id)
                            

                    clusters[cluster_id] = Cluster(cluster_id=cluster_id, features=group_sorted)
                    cluster_id_local += 1

            self.clusters[sample_name] = clusters          
            


    def _keep_longest_cluster(self, cluster):
        """
        Retain only the largest cluster.
        """
        subsets_to_remove = []
        signatures = {cid: set(f.feature_id for f in c.features) for cid, c in cluster.items()}
        # print("signatures:", signatures)
        sorted_clusters = sorted(signatures.items(), key=lambda x: len(x[1]), reverse=True)
        # print("sorted clusters:", sorted_clusters)
        kept = []
        # Compare from largest to smallest cluster to identify subsets
        # If a smaller cluster is a subset of any kept larger cluster, mark it for removal
        for cid, sig1 in sorted_clusters:
            is_subset = False
            for c_id, sig2 in kept:
                if sig1 < sig2:
                    is_subset = True
                    subsets_to_remove.append(f"{sig1} is subset of {sig2}")
                    del cluster[cid]
                    break
            if not is_subset:
                kept.append((cid, sig1))

    def _keep_closest_mz_candidate(self, cluster):
        """
        Keep only the feature closest to the expected m/z for each isotopologue in the cluster.
        """
        for cluster in cluster.values():
            iso_to_candidate  = defaultdict(list)
            base_mz = cluster.lowest_mz
            
            for feature in cluster.features:
                iso_index = Misc.calculate_isotopologue_index(feature.mz, base_mz, self.mzshift_tracer)
                iso_to_candidate[iso_index].append(feature)
                cluster.features = [min(candidates, key=lambda f: abs(f.mz - (base_mz + iso * self.mzshift_tracer))) for iso, candidates in iso_to_candidate.items()]

    def deduplicate_clusters(self, keep: str=None):
        """
        Clean up and deduplicate clusters by :
        - Merging clusters with identical feature compositions.
        - Removing clusters that are subsets of larger clusters (if keep is "longest").
        - Keeping only the best candidate feature for each isotopologue (if keep is "best_candidate").
        - Updating each feature's cluster memberships, isotopologue numbers, and also_in lists.

        Parameters:
            keep (str): Strategy for deduplication. Options are "longest" to keep the largest cluster,
                        "best_candidate" to retain only the feature with the highest intensity for each isotopologue within a cluster,
                        or "both" to apply both strategies.
        """
        final_clusters = {}
        merged = 0
        for sample, clusters in self.clusters.items():
            # --- Merge identical clusters ---
            final_clusters[sample] = {}
            seen_signatures = {}

            for cluster in clusters.values():
                signature = frozenset(f.feature_id for f in cluster.features)
                if signature not in seen_signatures:
                    seen_signatures[signature] = cluster.cluster_id
                    final_clusters[sample][cluster.cluster_id] = cluster
                else:
                    merged += 1
            # --- Remove subset clusters if keep_richest is True ---
            if keep == "longest":
                self._keep_longest_cluster(final_clusters[sample])
            elif keep =="closest_mz":
                self._keep_closest_mz_candidate(final_clusters[sample])
            elif keep == "both":
                self._keep_longest_cluster(final_clusters[sample])
                self._keep_closest_mz_candidate(final_clusters[sample])
            
            # --- Assign final cluster_id, isotopologues label, in_cluster and also_in to features ---
            new_clusters = {}
            features_to_clusters = defaultdict(set)
            for new_index, cluster in enumerate(final_clusters[sample].values()):
                cluster.cluster_id = f"C{new_index}"
                cluster.features.sort(key=lambda f: f.mz)
                min_mz=cluster.lowest_mz
                new_clusters[cluster.cluster_id] = cluster
                for f in cluster.features:
                    features_to_clusters[f.feature_id].add(cluster.cluster_id)
                    iso_index = round((f.mz - min_mz) / self.mzshift_tracer)
                    iso_label = "Mx" if iso_index == 0 else f"Mx+{iso_index}"
                    f.cluster_isotopologue[cluster.cluster_id] = iso_label
                    f.in_cluster = list(features_to_clusters[f.feature_id])
                    # f.also_in = [c for c in f.in_cluster if c != cluster.cluster_id]
                
            final_clusters[sample] = new_clusters
            self.clusters = final_clusters

            self.unclustered_features = {}
            for sample, features in self.features.items():
                self.unclustered_features[sample] = [f for f in features.values() if not f.in_cluster]

            final = len(next(iter(self.clusters.values()))) if self.clusters else 0
            unclustered = sum(1 for f in next(iter(self.features.values())).values() if not f.in_cluster) if self.features else 0
            # return merged, subset_removed, final, unclustered



# if __name__ == "__main__":
#     from isogroup.base.io import IoHandler
#     import pandas as pd
#     # from isogroup.base.database import Database
#     io = IoHandler()
#     # data= io.read_dataset(r"..\..\data\dataset_test_XCMS.txt")
#     untargeted = UntargetedExperiment(dataset=data, tracer="13C", mz_tol=5, rt_tol=15)
#     untargeted.initialize_experimental_features()
#     untargeted.build_clusters(RTwindow=15, ppm_tolerance=5)
#     # print(untargeted.clusters)
#     untargeted.deduplicate_clusters()


##################################################

    # def deduplicate_clusters(self, keep_best_candidate: bool = False, keep_richest: bool = True):
    #     """
    #     Clean up and deduplicate clusters by :
    #     - Merging clusters with identical feature compositions.
    #     - Removing clusters that are subsets of larger clusters (if keep_richest is True).
    #     - Keeping only the best candidate feature for each isotopologue (if keep_best_candidate is True).
    #     - Updating each feature's cluster memberships, isotopologue numbers, and also_in lists.

    #     Parameters:
    #         keep_best_candidate (bool): If True, retain only the feature with the highest intensity for each isotopologue within a cluster.
    #         keep_richest (bool): If True, retain only the largest cluster when multiple clusters share features.
    #     """
    #     merged = 0
    #     subset_removed = 0
    #     final_clusters = {}
    
    #     for sample, clusters in self.clusters.items():

    #         # --- Merge identical clusters ---
    #         final_clusters[sample] = {}
            
    #         seen_signatures = {}
    #         next_cluster_id = 0

    #         for cluster in clusters.values():
    #             signature = frozenset(f.feature_id for f in cluster.features)
    #             if signature not in seen_signatures:
    #                 cluster.cluster_id = f"C{next_cluster_id}"
    #                 seen_signatures[signature] = cluster.cluster_id
    #                 final_clusters[sample][cluster.cluster_id] = cluster
    #                 next_cluster_id += 1
    #             else:
    #                 merged += 1

    #         # --- Remove subset clusters if keep_richest is True ---
    #         if keep_richest:
    #             signatures = {cid: set(f.feature_id for f in c.features) for cid, c in final_clusters[sample].items()}
    #             sorted_clusters = sorted(signatures.items(), key=lambda x: len(x[1]), reverse=True)
    #             to_remove = set()
    #             kept = []
    #             # Compare from largest to smallest cluster to identify subsets
    #             # If a smaller cluster is a subset of any kept larger cluster, mark it for removal
    #             for cid, sig1 in sorted_clusters:
    #                 if any(sig1 < sig2 for _, sig2 in kept):
    #                     to_remove.add(cid)
    #                     print(to_remove)
    #                 else:
    #                     kept.append((cid, sig1))
    #             subset_removed += len(to_remove)
    #             final_clusters[sample] = {cid: c for cid, c in final_clusters[sample].items() if cid not in to_remove}

    #         # --- Keep only the best candidate for each isotopologue (based on the the closest m/z to expected) if keep_best_candidate is True ---
    #         if keep_best_candidate:
    #             for cluster in final_clusters[sample].values():
    #                 iso_to_candidate  = defaultdict(list)
    #                 base_mz = cluster.lowest_mz
    #                 for f in cluster.features:
    #                     iso_index = round((f.mz - base_mz) / self.mzshift_tracer)
    #                     iso_to_candidate[iso_index].append(f)
    #                     cluster.features = [min(candidates, key=lambda f: abs(f.mz - (base_mz + iso * self.mzshift_tracer))) for iso, candidates in iso_to_candidate.items()]
    #         # --- Assign final cluster_id, isotopologues label, in_cluster and also_in to features ---
    #         features_to_clusters = defaultdict(set)
    #         for cluster in final_clusters[sample].values():
    #             for f in cluster.features:
    #                 features_to_clusters[f.feature_id].add(cluster.cluster_id)
    #         for cluster in final_clusters[sample].values():
    #             cluster.features.sort(key=lambda f: f.mz)
    #             min_mz=cluster.lowest_mz
    #             for f in cluster.features:
    #                 iso_index = round((f.mz - min_mz) / self.mzshift_tracer)
    #                 iso_label = "Mx" if iso_index == 0 else f"Mx+{iso_index}"
    #                 f.cluster_isotopologue[cluster.cluster_id] = iso_label
    #                 f.in_cluster = list(features_to_clusters[f.feature_id])
    #                 f.also_in = [c for c in f.in_cluster if c != cluster.cluster_id]
        
    #     self.clusters = final_clusters
        
    #     # Keep unclustered features for reference
    #     self.unclustered_features = {}
    #     for sample, features in self.features.items():
    #         self.unclustered_features[sample] = [f for f in features.values() if not f.in_cluster]

    #     final = len(next(iter(self.clusters.values()))) if self.clusters else 0
    #     unclustered = sum(1 for f in next(iter(self.features.values())).values() if not f.in_cluster) if self.features else 0
    #     return merged, subset_removed, final, unclustered

    # def clusters_to_dataframe(self) -> pd.DataFrame:
    #     """
    #     Convert the clusters into a pandas DataFrame for easier analysis and export.
    #     :return: pd.DataFrame
    #     """
    #     records = []
    #     for sample_name, clusters in self.clusters.items():
    #         for cluster in clusters.values():
    #             sorted_features = sorted(cluster.features, key=lambda f: f.mz)

    #             for idx, f in enumerate(sorted_features):
    #                 iso_label = f.cluster_isotopologue.get(cluster.cluster_id, "Mx")
    #                 records.append({
    #                     "ClusterID": cluster.cluster_id,
    #                     "FeatureID": f.feature_id,
    #                     "RT": f.rt,
    #                     "m/z": f.mz,
    #                     "sample": f.sample,
    #                     "Intensity": f.intensity,
    #                     "Isotopologue": iso_label,
    #                     "InClusters": f.in_cluster,
    #                     "AlsoIn": f.also_in
    #                 })

    #     return pd.DataFrame.from_records(records)


    # def export_clusters_to_tsv(self, filepath: str):
    #     """
    #     Export the clusters to a CSV file.
    #     :param filepath: str
    #     """
    #     df = self.clusters_to_dataframe()
    #     df.to_csv(filepath, sep="\t", index=False)


    # def export_features(self, filename: str):
    #     """
    #     Export all features to a TSV file.
    #     :param filename: str
    #     """
    #     records = []
    #     for sample_name, features in self.features.items():
    #         for f in features.values():
    #             # If not in any cluster, mark accordingly
    #             cluster_ids = f.in_cluster if f.in_cluster else ["None"]
    #             iso_labels = [f.cluster_isotopologue.get(cid, "N/A") for cid in cluster_ids]

    #             records.append({
    #                 "FeatureID": f.feature_id,
    #                 "RT": f.rt,
    #                 "m/z": f.mz,
    #                 "sample": f.sample,
    #                 "Intensity": f.intensity,
    #                 "InClusters": cluster_ids,
    #                 "Isotopologues": iso_labels
    #             })

    #     df = pd.DataFrame.from_records(records)
    #     df.to_csv(filename, sep="\t", index=False)
