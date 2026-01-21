from isogroup.base.experiment import Experiment
import math
import pytest

def test_initialize_experimental_features(dataset_df):
    """
    Test the initialization of experimental features in the Experiment class.

    :param dataset_df: DataFrame containing the dataset for testing.
    """
    experiment = Experiment(dataset_df, tracer="13C", ppm_tol=5, rt_tol=15)
    experiment.initialize_experimental_features()
    assert len(experiment.features) == 2
    assert "Sample_1" in experiment.features
    assert "Sample_2" in experiment.features
    assert len(experiment.features["Sample_1"]) == 9
    assert len(experiment.features["Sample_2"]) == 9
    assert math.isclose(experiment.features["Sample_1"]["F1"].mz, 119.025753)
    assert math.isclose(experiment.features["Sample_2"]["F3"].rt, 679.9930235)
    assert math.isclose(experiment.features["Sample_1"]["F4"].intensity, 0.0)


@pytest.mark.parametrize("tracer_code, expected_element, expected_index", [
    ("13C", "C", 1),
    ("15N", "N", 1),
    ("12C", "C", 0)])

def test_tracer_properties(dataset_df, tracer_code, expected_element, expected_index):
    """
    Test the tracer properties in the Experiment class.

    :param dataset_df: DataFrame containing the dataset for testing.
    :param tracer_code: Tracer code to test.
    :param expected_element: Expected tracer element.
    :param expected_index: Expected tracer index.
    """
    experiment = Experiment(dataset_df,tracer=tracer_code, ppm_tol=5, rt_tol=15)
    assert experiment.tracer_element == expected_element
    assert experiment.tracer_idx == expected_index

@pytest.mark.parametrize("tracer_code", [
    "XyZ",
    "123",
    " "])

def test_wrong_tracer(dataset_df, tracer_code):
    """
    Test the handling of wrong tracer codes in the Experiment class.

    :param dataset_df: DataFrame containing the dataset for testing.
    :param tracer_code: Invalid tracer code to test.
    """
    with pytest.raises(ValueError):
        Experiment(dataset_df, tracer=tracer_code, ppm_tol=5, rt_tol=15)   


@pytest.mark.parametrize("rt_tol", [8.5, 11.0, 15.5])
@pytest.mark.parametrize("ppm_tol", [2.0, 5.0, 10.0])

def test_feature_tol_setters(dataset_df, rt_tol, ppm_tol):
    """
    Test the feature tolerance setters in the Experiment class.

    :param dataset_df: DataFrame containing the dataset for testing.
    :param rt_tol: Retention time tolerance to set.
    :param ppm_tol: ppm tolerance to set.
    """
    experiment = Experiment(dataset_df, tracer="13C", ppm_tol=ppm_tol, rt_tol=rt_tol)
    assert experiment.rt_tol == rt_tol
    assert experiment.ppm_tol == ppm_tol

@pytest.mark.parametrize("rt_tol", ["a", None, [10, 9]])
@pytest.mark.parametrize("ppm_tol", ["b", None, [5, 6]])

def test_wrong_tol_setters(dataset_df, rt_tol, ppm_tol):
    """
    Test the handling of wrong tolerance values in the Experiment class.

    :param dataset_df: DataFrame containing the dataset for testing.
    :param rt_tol: Invalid retention time tolerance to test.
    :param ppm_tol: Invalid ppm tolerance to test.
    """
    experiment = Experiment(dataset_df, tracer="13C", ppm_tol=ppm_tol, rt_tol=rt_tol)
    with pytest.raises(ValueError):
        experiment.rt_tol = rt_tol
    with pytest.raises(ValueError):
        experiment.ppm_tol = ppm_tol

# def test_database_initialization(database_df):
#     experiment = Experiment(tracer="13C", mz_tol=5, rt_tol=10, database=database_df)
#     assert len(experiment.database.theoretical_features) == 39


    