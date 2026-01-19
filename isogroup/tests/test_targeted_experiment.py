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
    """
    Test the annotation of features in a TargetedExperiment.

    :param dataset_df: DataFrame containing the dataset features.
    :param database_df: DataFrame containing the database of known metabolites.
    :param metabolite: Expected metabolite annotation for the feature.
    :param feature_id: ID of the feature to be tested.
    """
    targeted_experiment = TargetedExperiment(dataset=dataset_df,
                                            tracer="13C",
                                            ppm_tol=5,
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
    """
    Test that strict tolerances result in no feature annotations.

    :param dataset_df: DataFrame containing the dataset features.
    :param database_df: DataFrame containing the database of known metabolites.
    """
    # Using very strict tolerances to ensure no matches
    targeted_experiment = TargetedExperiment(dataset=dataset_df,
                                            tracer="13C",
                                            ppm_tol=1,
                                            rt_tol=1,
                                            database=database_df)
    targeted_experiment.initialize_experimental_features()
    targeted_experiment.annotate_features()
    assert targeted_experiment.features["Sample_1"]["F1"].metabolite == []
    assert targeted_experiment.features["Sample_2"]["F5"].metabolite == []

@pytest.mark.parametrize("cluster, features_nb",
                         [("Succinate", 2),
                          ("Malate", 5),
                          ("Citrate", 1),
                          ("Isocitrate", 1)])

def test_clusterize(dataset_df, database_df, cluster, features_nb):
    """
    Test the clustering of features in a TargetedExperiment.

    :param dataset_df: DataFrame containing the dataset features.
    :param database_df: DataFrame containing the database of known metabolites.
    :param cluster: Expected cluster name.
    :param features_nb: Expected number of features in the cluster.
    """
    targeted_experiment = TargetedExperiment(dataset=dataset_df,
                                            tracer="13C",
                                            ppm_tol=5,
                                            rt_tol=15,
                                            database=database_df)
    targeted_experiment.initialize_experimental_features()
    targeted_experiment.annotate_features()
    targeted_experiment.clusterize()
    assert len(targeted_experiment.clusters["Sample_1"]) == 4
    assert len(targeted_experiment.clusters["Sample_2"]) == 4
    assert cluster in targeted_experiment.clusters["Sample_1"]
    assert cluster in targeted_experiment.clusters["Sample_2"]
    assert len(targeted_experiment.clusters["Sample_1"][cluster]) == features_nb
    assert len(targeted_experiment.clusters["Sample_2"][cluster]) == features_nb