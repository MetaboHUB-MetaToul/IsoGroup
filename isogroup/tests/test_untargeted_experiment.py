from isogroup.base.untargeted_experiment import UntargetedExperiment
import math
import pytest

@pytest.mark.parametrize("cluster_id, nb_features, features_id",
                         [("C0", 2, ["F1", "F2"]),
                          ("C1", 2,["F1", "F2"]),
                          ("C2", 5, ['F9', 'F8', 'F7', 'F6', 'F5']),
                          ("C3", 5, ['F9', 'F8', 'F7', 'F6', 'F5']),
                          ("C4", 5, ['F9', 'F8', 'F7', 'F6', 'F5']),
                          ("C5", 5, ['F9', 'F8', 'F7', 'F6', 'F5']),
                          ("C6", 5, ['F9', 'F8', 'F7', 'F6', 'F5'])])


def test_build_cluster(dataset_df, cluster_id, nb_features, features_id):
    untargeted_experiment = UntargetedExperiment(dataset=dataset_df,
                                                tracer="13C",
                                                mz_tol=5,
                                                rt_tol=15,
                                                max_atoms=None)
   
    untargeted_experiment.initialize_experimental_features()
    untargeted_experiment.build_clusters(RTwindow=15, ppm_tolerance=5, max_atoms=None)
    assert len(untargeted_experiment.clusters["Sample_1"]) == 7
    assert len(untargeted_experiment.clusters["Sample_2"]) == 7
    assert cluster_id in untargeted_experiment.clusters["Sample_1"]
    assert cluster_id in untargeted_experiment.clusters["Sample_2"]
    assert len(untargeted_experiment.clusters["Sample_1"][cluster_id].features) == nb_features
    assert len(untargeted_experiment.clusters["Sample_2"][cluster_id].features) == nb_features
    for f_id in features_id:
        assert f_id in [f.feature_id for f in untargeted_experiment.clusters["Sample_1"][cluster_id].features]
        assert f_id in [f.feature_id for f in untargeted_experiment.clusters["Sample_2"][cluster_id].features]

def test_no_cluster(dataset_df):
    # Using very strict tolerances to ensure no clusters
    untargeted_experiment = UntargetedExperiment(dataset=dataset_df,
                                                tracer="13C",
                                                mz_tol=1,
                                                rt_tol=1,
                                                max_atoms=None)
    untargeted_experiment.initialize_experimental_features()
    untargeted_experiment.build_clusters(RTwindow=0.01, ppm_tolerance=0.01)
    assert len(untargeted_experiment.clusters["Sample_1"]) == 0
    assert len(untargeted_experiment.clusters["Sample_2"]) == 0