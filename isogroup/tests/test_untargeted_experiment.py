from isogroup.base.untargeted_experiment import UntargetedExperiment
import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch

@pytest.mark.parametrize("cluster_id, nb_features, features_id",
                         [("C0", 2, ["F1", "F2"]),
                          ("C1", 2,["F1", "F2"]),
                          ("C2", 5, ['F9', 'F8', 'F7', 'F6', 'F5']),
                          ("C3", 5, ['F9', 'F8', 'F7', 'F6', 'F5']),
                          ("C4", 5, ['F9', 'F8', 'F7', 'F6', 'F5']),
                          ("C5", 5, ['F9', 'F8', 'F7', 'F6', 'F5']),
                          ("C6", 5, ['F9', 'F8', 'F7', 'F6', 'F5'])])

def test_build_cluster(features_dict, cluster_id, nb_features, features_id):
    """
    Test the build_clusters method of the UntargetedExperiment class.
    It checks if the clusters are built correctly based on the provided dataset.

    :param dataset_df: DataFrame containing the dataset for the experiment.
    :param cluster_id: Expected cluster ID to check.
    :param nb_features: Expected number of features in the cluster.
    :param features_id: List of expected feature IDs in the cluster.
    """
    with patch.object(UntargetedExperiment, 'initialize_experimental_features', return_value=None):
        untargeted_experiment = UntargetedExperiment(dataset=pd.DataFrame(),
                                                tracer="13C",
                                                ppm_tol=5,
                                                rt_tol=15,
                                                max_atoms=None)
        untargeted_experiment.features = features_dict
   
    untargeted_experiment.build_clusters(rt_tol=15, ppm_tol=5, max_atoms=None)
    assert len(untargeted_experiment.clusters["Sample_1"]) == 7
    assert len(untargeted_experiment.clusters["Sample_2"]) == 7
    assert cluster_id in untargeted_experiment.clusters["Sample_1"]
    assert cluster_id in untargeted_experiment.clusters["Sample_2"]
    assert len(untargeted_experiment.clusters["Sample_1"][cluster_id].features) == nb_features
    assert len(untargeted_experiment.clusters["Sample_2"][cluster_id].features) == nb_features
    for f_id in features_id:
        assert f_id in [f.feature_id for f in untargeted_experiment.clusters["Sample_1"][cluster_id].features]
        assert f_id in [f.feature_id for f in untargeted_experiment.clusters["Sample_2"][cluster_id].features]

def test_no_cluster(features_dict):
    """
    Test the build_clusters method of the UntargetedExperiment class with strict tolerances.
    It checks that no clusters are formed when the tolerances are too strict.

    :param dataset_df: DataFrame containing the dataset for the experiment.
    """
    # Using very strict tolerances to ensure no clusters
    untargeted_experiment = UntargetedExperiment(dataset=pd.DataFrame(),
                                                tracer="13C",
                                                ppm_tol=1,
                                                rt_tol=1,
                                                max_atoms=None)
    
    with patch.object(UntargetedExperiment, 'initialize_experimental_features', return_value=None):
        untargeted_experiment.features = features_dict

    untargeted_experiment.build_clusters(rt_tol=0.01, ppm_tol=0.01)
    assert len(untargeted_experiment.clusters["Sample_1"]) == 0
    assert len(untargeted_experiment.clusters["Sample_2"]) == 0

@pytest.mark.parametrize("deduplication_method, cluster_id, nb_features, features_id",
                         [(None,"C0", 2, ["F1", "F2"]),
                          (None, "C1", 7, ['F11', 'F10', 'F9', 'F6', 'F5', 'F7', 'F8']),
                          ("closest_mz", "C0", 2, ["F1", "F2"]),
                          ("closest_mz", "C1", 5, ['F11', 'F10', 'F9', 'F6', 'F5'])])


def test_deduplicate_clusters(dataset_df_duplicates, deduplication_method, cluster_id, nb_features, features_id):
    """
    Test the deduplicate_clusters method of the UntargetedExperiment class.

    :param dataset_df: DataFrame containing the dataset for the experiment.
    :param cluster_id: Expected cluster ID to check after deduplication.
    :param nb_features: Expected number of features in the cluster after deduplication.
    :param features_id: List of expected feature IDs in the cluster after deduplication.
    """
    untargeted_experiment = UntargetedExperiment(dataset=dataset_df_duplicates,
                                                tracer="13C",
                                                ppm_tol=5,
                                                rt_tol=15,
                                                max_atoms=None)
    untargeted_experiment.initialize_experimental_features()
    untargeted_experiment.build_clusters(rt_tol=15, ppm_tol=5, max_atoms=None)
    untargeted_experiment.deduplicate_clusters(deduplication_method)
    assert len(untargeted_experiment.clusters["Sample_1"]) == 2
    assert len(untargeted_experiment.clusters["Sample_2"]) == 2
    assert cluster_id in untargeted_experiment.clusters["Sample_1"]
    assert cluster_id in untargeted_experiment.clusters["Sample_2"]
    assert len(untargeted_experiment.clusters["Sample_1"][cluster_id].features) == nb_features
    assert len(untargeted_experiment.clusters["Sample_2"][cluster_id].features) == nb_features
    for f_id in features_id:
        assert f_id in [f.feature_id for f in untargeted_experiment.clusters["Sample_1"][cluster_id].features]
        assert f_id in [f.feature_id for f in untargeted_experiment.clusters["Sample_2"][cluster_id].features]

def test_unlabeled_enhancer():
    """
    Test the unlabeled_enhancer method of the UntargetedExperiment class with the "unlabeled" enhancer.
    It checks that the enhancer correctly identifies features found in the unlabeled sample and calculates 
    the Mx+1/Mx ratio for the clusters.
    """
    untargeted_experiment = UntargetedExperiment(dataset=pd.DataFrame(),
                                                tracer="13C",
                                                ppm_tol=5,
                                                rt_tol=15,
                                                max_atoms=None)
    
    cluster_df = pd.DataFrame({'ClusterID': ['C0', 'C0', 'C1', 'C1', 'C1', 'C1', 'C1', 'C0', 'C0', 'C1', 'C1', 'C1', 'C1', 'C1'], 
                                'FeatureID': ['F1', 'F2', 'F9', 'F8', 'F7', 'F6', 'F5', 'F1', 'F2', 'F9', 'F8', 'F7', 'F6', 'F5'], 
                                'RT': [667.779067, 667.9255408, 676.8952154, 676.8898827, 676.6045604, 676.5620229, 676.4604364, 667.779067, 667.9255408, 676.8952154, 676.8898827, 676.6045604, 676.5620229, 676.4604364], 
                                'm/z': [119.025753, 120.0291332, 133.0140851, 134.0174803, 135.0208168, 136.024129, 137.0275004, 119.025753, 120.0291332, 133.0140851, 134.0174803, 135.0208168, 136.024129, 137.0275004], 
                                'sample': ['Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2'], 
                                'Intensity': [1571414706.0, 0.0, 543216118.8, 2077278842.0, 3105587268.0, 2090662547.0, 529223407.9, 266171108.6, 129533534.2, 155940888.7, 218743897.0, 154077393.8, 97127965.25, 28994270.58], 
                                'Isotopologue': ['Mx', 'Mx+1', 'Mx', 'Mx+1', 'Mx+2', 'Mx+3', 'Mx+4', 'Mx', 'Mx+1', 'Mx', 'Mx+1', 'Mx+2', 'Mx+3', 'Mx+4'], 
                                'AlsoIn': ['[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]']})
    
    untargeted_experiment.unlabeled_enhancer(cluster_df, sample_name="Sample_1")
    assert "Found in unlabeled sample (Sample_1)" in untargeted_experiment.all_clusters_df.columns
    assert (untargeted_experiment.all_clusters_df.loc[untargeted_experiment.all_clusters_df["FeatureID"] == "F2",  "Found in unlabeled sample (Sample_1)"] == "No").all()
    
    assert "Mx+1/Mx ratio" in untargeted_experiment.all_clusters_df.columns
    assert (untargeted_experiment.all_clusters_df.loc[(untargeted_experiment.all_clusters_df["ClusterID"] == "C0") 
                                                      & (untargeted_experiment.all_clusters_df["sample"] == "Sample_1"), "Mx+1/Mx ratio"] == "ND").all()
    assert np.allclose(untargeted_experiment.all_clusters_df.loc[(untargeted_experiment.all_clusters_df["ClusterID"] == "C1") 
                                                                 & (untargeted_experiment.all_clusters_df["sample"] == "Sample_1"), "Mx+1/Mx ratio"], 3.8240, rtol=0, atol=1e-4)
    
    assert (untargeted_experiment.all_clusters_df.loc[(untargeted_experiment.all_clusters_df["ClusterID"] == "C0") 
                                                      & (untargeted_experiment.all_clusters_df["sample"] == "Sample_2"), "Mx+1/Mx ratio"]).isna().all()
    assert (untargeted_experiment.all_clusters_df.loc[(untargeted_experiment.all_clusters_df["ClusterID"] == "C1") & 
                                                      (untargeted_experiment.all_clusters_df["sample"] == "Sample_2"), "Mx+1/Mx ratio"]).isna().all()

     
def test_fully_labeled_enhancer():
    """
    Test the fully_labeled_enhancer method of the UntargetedExperiment class with the "fully_labeled" enhancer.
    It checks that the enhancer correctly identifies features found in the fully labeled sample.
    """
    untargeted_experiment = UntargetedExperiment(dataset=pd.DataFrame(),
                                                tracer="13C",
                                                ppm_tol=5,
                                                rt_tol=15,
                                                max_atoms=None)
    
    cluster_df = pd.DataFrame({'ClusterID': ['C0', 'C0', 'C1', 'C1', 'C1', 'C1', 'C1', 'C0', 'C0', 'C1', 'C1', 'C1', 'C1', 'C1'], 
                                'FeatureID': ['F1', 'F2', 'F9', 'F8', 'F7', 'F6', 'F5', 'F1', 'F2', 'F9', 'F8', 'F7', 'F6', 'F5'], 
                                'RT': [667.779067, 667.9255408, 676.8952154, 676.8898827, 676.6045604, 676.5620229, 676.4604364, 667.779067, 667.9255408, 676.8952154, 676.8898827, 676.6045604, 676.5620229, 676.4604364], 
                                'm/z': [119.025753, 120.0291332, 133.0140851, 134.0174803, 135.0208168, 136.024129, 137.0275004, 119.025753, 120.0291332, 133.0140851, 134.0174803, 135.0208168, 136.024129, 137.0275004], 
                                'sample': ['Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2'], 
                                'Intensity': [1571414706.0, 0.0, 543216118.8, 2077278842.0, 3105587268.0, 2090662547.0, 529223407.9, 266171108.6, 129533534.2, 155940888.7, 218743897.0, 154077393.8, 97127965.25, 28994270.58], 
                                'Isotopologue': ['Mx', 'Mx+1', 'Mx', 'Mx+1', 'Mx+2', 'Mx+3', 'Mx+4', 'Mx', 'Mx+1', 'Mx', 'Mx+1', 'Mx+2', 'Mx+3', 'Mx+4'], 
                                'AlsoIn': ['[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]']})
    
    untargeted_experiment.fully_labeled_enhancer(cluster_df, sample_name="Sample_1")
    assert "Found in fully labeled sample (Sample_1)" in untargeted_experiment.all_clusters_df.columns
    assert (untargeted_experiment.all_clusters_df.loc[untargeted_experiment.all_clusters_df["FeatureID"] == "F2",  "Found in fully labeled sample (Sample_1)"] == "No").all()
   
    


