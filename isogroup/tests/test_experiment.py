from isogroup.base.experiment import Experiment
import pandas as pd
import math
import pytest

# {'Sample_1': {'F1': Feature(ID = F1, RT=667.779067, Metabolite=[], Isotopologue={}, mz=119.025753, intensity=1571414706.0), 
#               'F2': Feature(ID = F2, RT=667.9255408, Metabolite=[], Isotopologue={}, mz=120.0291332, intensity=1059554882.0), 
#               'F3': Feature(ID = F3, RT=679.9930235, Metabolite=[], Isotopologue={}, mz=191.0191775654, intensity=31398195.78), 
#               'F4': Feature(ID = F4, RT=678.1606593, Metabolite=[], Isotopologue={}, mz=119.0232843, intensity=0.0), 
#               'F5': Feature(ID = F5, RT=676.4604364, Metabolite=[], Isotopologue={}, mz=137.0275004, intensity=529223407.9), 
#               'F6': Feature(ID = F6, RT=676.5620229, Metabolite=[], Isotopologue={}, mz=136.024129, intensity=2090662547.0), 
#               'F7': Feature(ID = F7, RT=676.6045604, Metabolite=[], Isotopologue={}, mz=135.0208168, intensity=3105587268.0), 
#               'F8': Feature(ID = F8, RT=676.8898827, Metabolite=[], Isotopologue={}, mz=134.0174803, intensity=2077278842.0), 
#               'F9': Feature(ID = F9, RT=676.8952154, Metabolite=[], Isotopologue={}, mz=133.0140851, intensity=543216118.8)}, 
# 'Sample_2': {'F1': Feature(ID = F1, RT=667.779067, Metabolite=[], Isotopologue={}, mz=119.025753, intensity=266171108.6), 
#               'F2': Feature(ID = F2, RT=667.9255408, Metabolite=[], Isotopologue={}, mz=120.0291332, intensity=129533534.2), 
#               'F3': Feature(ID = F3, RT=679.9930235, Metabolite=[], Isotopologue={}, mz=191.0191775654, intensity=5324316.124), 
#               'F4': Feature(ID = F4, RT=678.1606593, Metabolite=[], Isotopologue={}, mz=119.0232843, intensity=0.0), 
#               'F5': Feature(ID = F5, RT=676.4604364, Metabolite=[], Isotopologue={}, mz=137.0275004, intensity=28994270.58), 
#               'F6': Feature(ID = F6, RT=676.5620229, Metabolite=[], Isotopologue={}, mz=136.024129, intensity=97127965.25), 
#               'F7': Feature(ID = F7, RT=676.6045604, Metabolite=[], Isotopologue={}, mz=135.0208168, intensity=154077393.8), 
#               'F8': Feature(ID = F8, RT=676.8898827, Metabolite=[], Isotopologue={}, mz=134.0174803, intensity=218743897.0), 
#               'F9': Feature(ID = F9, RT=676.8952154, Metabolite=[], Isotopologue={}, mz=133.0140851, intensity=155940888.7)}}

@pytest.mark.parametrize("f_id, mz, rt, intensity",
        [("F1", 119.0257, 667.7790, 1571414706.0),
         ("F2", 120.0291, 667.9255, 1059554882.0),
         ("F3", 191.0191, 679.9930, 31398195.78),
         ("F4", 119.0232, 678.1606, 0.0),
         ("F5", 137.0275, 676.4604, 529223407.9),
         ("F6", 136.0241, 676.5620, 2090662547.0),
         ("F7", 135.0208, 676.6045, 3105587268.0),
         ("F8", 134.0174, 676.8898, 2077278842.0),
         ("F9", 133.0140, 676.8952, 543216118.8)])

def test_initialize_experimental_features(dataset_df, f_id, mz, rt, intensity):
    """
    Test the initialization of experimental features in the Experiment class.

    :param dataset_df: DataFrame containing the dataset for testing.
    :param f_id: Feature ID to test.
    :param mz: m/z value to test.
    :param rt: Retention time value to test.
    :param intensity: Intensity value to test.
    """
    experiment = Experiment(dataset_df, tracer="13C", ppm_tol=5, rt_tol=15)
    experiment.initialize_experimental_features()
    assert len(experiment.features) == 2
    assert "Sample_1" in experiment.features
    assert "Sample_2" in experiment.features
    assert len(experiment.features["Sample_1"]) == 9
    assert len(experiment.features["Sample_2"]) == 9
    assert f_id in experiment.features["Sample_1"]
    assert f_id in experiment.features["Sample_2"]
    assert math.isclose(experiment.features["Sample_1"][f_id].mz, mz, rel_tol=0, abs_tol=1e-4 )
    assert math.isclose(experiment.features["Sample_1"][f_id].rt, rt, rel_tol=0, abs_tol=1e-4)
    assert math.isclose(experiment.features["Sample_1"][f_id].intensity, intensity, rel_tol=0, abs_tol=1e-4)


@pytest.mark.parametrize("wrong_input_dataframe", [
    pd.DataFrame({"mz": [119.0257, 120.0291], "rt": [667.7790, 667.9255], "id": ["F1", "F2"]}),
    pd.DataFrame({"m/z": [119.0257, 120.0291], "rt": [667.7790, 667.9255], "id": ["F1", "F2"], "Sample_1": [1571414706.0, 1059554882.0], "Sample_2": [266171108.6, 129533534.2]}),
    pd.DataFrame({"mz": [119.0257, 120.0291], "rt": [667.7790, 667.9255], "Sample_1": [1571414706.0, 1059554882.0], "Sample_2": [266171108.6, 129533534.2]})])


def test_wrong_column_entry(wrong_input_dataframe):
    """
    Test the handling of wrong column entries in the Experiment class.
    """
    experiment = Experiment(wrong_input_dataframe, tracer="13C", ppm_tol=5, rt_tol=15)
    with pytest.raises(ValueError):
        experiment.initialize_experimental_features()

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
    
# @pytest.mark.parametrize("rt_tol", [8.5, 11.0, 15.5])
# @pytest.mark.parametrize("ppm_tol", [2.0, 5.0, 10.0])

# def test_feature_tol_setters(dataset_df, rt_tol, ppm_tol):
#     """
#     Test the feature tolerance setters in the Experiment class.

#     :param dataset_df: DataFrame containing the dataset for testing.
#     :param rt_tol: Retention time tolerance to set.
#     :param ppm_tol: ppm tolerance to set.
#     """
#     experiment = Experiment(dataset_df, tracer="13C", ppm_tol=ppm_tol, rt_tol=rt_tol)
#     assert experiment.rt_tol == rt_tol
#     assert experiment.ppm_tol == ppm_tol

# @pytest.mark.parametrize("rt_tol", ["a", None, [10, 9]])
# @pytest.mark.parametrize("ppm_tol", ["b", None, [5, 6]])

# def test_wrong_tol_setters(dataset_df, rt_tol, ppm_tol):
#     """
#     Test the handling of wrong tolerance values in the Experiment class.

#     :param dataset_df: DataFrame containing the dataset for testing.
#     :param rt_tol: Invalid retention time tolerance to test.
#     :param ppm_tol: Invalid ppm tolerance to test.
#     """
#     experiment = Experiment(dataset_df, tracer="13C", ppm_tol=ppm_tol, rt_tol=rt_tol)
#     with pytest.raises(ValueError):
#         experiment.rt_tol = rt_tol
#     with pytest.raises(ValueError):
#         experiment.ppm_tol = ppm_tol

# def test_database_initialization(database_df):
#     experiment = Experiment(tracer="13C", mz_tol=5, rt_tol=10, database=database_df)
#     assert len(experiment.database.theoretical_features) == 39


    