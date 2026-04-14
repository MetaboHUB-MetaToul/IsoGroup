import pytest
import pandas as pd
import isogroup.enhancer.labeled_enhancer as labeled_enhancer


def test_annotate_feature_found():
    """
    Test the annotate_feature_found function of the labeled_enhancer module.

    :param dataset_df: DataFrame containing the dataset for the experiment.
    """
    clusters_df = pd.DataFrame({'ClusterID': ['C0', 'C0', 'C1', 'C1', 'C1', 'C1', 'C1', 'C0', 'C0', 'C1', 'C1', 'C1', 'C1', 'C1'], 
                                'FeatureID': ['F1', 'F2', 'F9', 'F8', 'F7', 'F6', 'F5', 'F1', 'F2', 'F9', 'F8', 'F7', 'F6', 'F5'], 
                                'RT': [667.779067, 667.9255408, 676.8952154, 676.8898827, 676.6045604, 676.5620229, 676.4604364, 667.779067, 667.9255408, 676.8952154, 676.8898827, 676.6045604, 676.5620229, 676.4604364], 
                                'm/z': [119.025753, 120.0291332, 133.0140851, 134.0174803, 135.0208168, 136.024129, 137.0275004, 119.025753, 120.0291332, 133.0140851, 134.0174803, 135.0208168, 136.024129, 137.0275004], 
                                'sample': ['Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2', 'Sample_2'], 
                                'Intensity': [1571414706.0, 0.0, 543216118.8, 2077278842.0, 3105587268.0, 2090662547.0, 529223407.9, 266171108.6, 129533534.2, 155940888.7, 218743897.0, 154077393.8, 97127965.25, 28994270.58], 
                                'Isotopologue': ['Mx', 'Mx+1', 'Mx', 'Mx+1', 'Mx+2', 'Mx+3', 'Mx+4', 'Mx', 'Mx+1', 'Mx', 'Mx+1', 'Mx+2', 'Mx+3', 'Mx+4'], 
                                'AlsoIn': ['[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]', '[]']})
    
    result_df = labeled_enhancer.annotate_feature_found(clusters_df, sample_name="Sample_1")
    assert "Found in fully labeled sample (Sample_1)" in result_df.columns
    assert result_df["Found in fully labeled sample (Sample_1)"].tolist() == ["Yes", "No", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "No", "Yes", "Yes", "Yes", "Yes", "Yes"]

def test_sample_not_found():
    """
    Test the annotate_feature_found function of the labeled_enhancer module when the specified sample is not found.
    """
    clusters_df = pd.DataFrame({'ClusterID': ['C0', 'C0', 'C1', 'C1', 'C1', 'C1', 'C1'], 
                                'FeatureID': ['F1', 'F2', 'F9', 'F8', 'F7', 'F6', 'F5'], 
                                'RT': [667.779067, 667.9255408, 676.8952154, 676.8898827, 676.6045604, 676.5620229, 676.4604364], 
                                'm/z': [119.025753, 120.0291332, 133.0140851, 134.0174803, 135.0208168, 136.024129, 137.0275004], 
                                'sample': ['Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1', 'Sample_1'], 
                                'Intensity': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                                'Isotopologue': ['Mx', 'Mx+1', 'Mx', 'Mx+1', 'Mx+2', 'Mx+3', 'Mx+4'], 
                                'AlsoIn': ['[]', '[]', '[]', '[]', '[]', '[]', '[]']})
    
    with pytest.raises(ValueError):
        labeled_enhancer.annotate_feature_found(clusters_df, sample_name="Sample_2")
