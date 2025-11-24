from isogroup.base.targeted_experiment import TargetedExperiment
import math
import pytest


@pytest.mark.parametrize("metabolite, feature_id",
                         [(["Succinate"], "F1"),
                          (["Succinate"], "F2"),
                          (["Citrate", "Isocitrate"], "F3"),
                          ([], "F4"),
                          (["Malate"], "F5"),
                          (["Malate"], "F6"),
                          (["Malate"], "F7"),
                          (["Malate"], "F8"),
                          (["Malate"], "F9")])

def test_annotate_features(dataset_df, database_df, metabolite, feature_id):
    targeted_experiment = TargetedExperiment(dataset=dataset_df,
                                            tracer="13C",
                                            mz_tol=5,
                                            rt_tol=15,
                                            database=database_df)
    targeted_experiment.initialize_experimental_features()
    targeted_experiment.annotate_features()
    assert targeted_experiment.features["Sample_1"][feature_id].metabolite == metabolite
    assert targeted_experiment.features["Sample_2"][feature_id].metabolite == metabolite
    assert math.isclose(targeted_experiment.features["Sample_1"]["F2"].rt_error[0], 0.07445919999997841)
    assert math.isclose(targeted_experiment.features["Sample_2"]["F3"].rt_error[1], 0.006976499999950647)
    assert math.isclose(targeted_experiment.features["Sample_1"]["F9"].mz_error[0], -2.9082559166254343)

def test_no_annotation(dataset_df, database_df):
    # Using very strict tolerances to ensure no matches
    targeted_experiment = TargetedExperiment(dataset=dataset_df,
                                            tracer="13C",
                                            mz_tol=1,
                                            rt_tol=1,
                                            database=database_df)
    targeted_experiment.initialize_experimental_features()
    targeted_experiment.annotate_features()
    assert targeted_experiment.features["Sample_1"]["F1"].metabolite == []
    assert targeted_experiment.features["Sample_2"]["F5"].metabolite == []

@pytest.mark.parametrize("Cluster, features_nb",
                         [("Succinate", 2),
                          ("Malate", 5),
                          ("Citrate", 1),
                          ("Isocitrate", 1)])

def test_clusterize(dataset_df, database_df, Cluster, features_nb):
    targeted_experiment = TargetedExperiment(dataset=dataset_df,
                                            tracer="13C",
                                            mz_tol=5,
                                            rt_tol=15,
                                            database=database_df)
    targeted_experiment.initialize_experimental_features()
    targeted_experiment.annotate_features()
    targeted_experiment.clusterize()
    assert len(targeted_experiment.clusters["Sample_1"]) == 4
    assert len(targeted_experiment.clusters["Sample_2"]) == 4
    assert Cluster in targeted_experiment.clusters["Sample_1"]
    assert Cluster in targeted_experiment.clusters["Sample_2"]
    assert len(targeted_experiment.clusters["Sample_1"][Cluster]) == features_nb
    assert len(targeted_experiment.clusters["Sample_2"][Cluster]) == features_nb