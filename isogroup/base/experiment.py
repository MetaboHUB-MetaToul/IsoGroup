from isogroup.base.feature import Feature
from isogroup.base.misc import Misc
import pandas as pd
import logging

logger = logging.getLogger(f"IsoGroup")

class Experiment:
    """
    Represents a mass spectrometry experiment with experimental features.

    Args:
        dataset (pd.DataFrame): DataFrame containing experimental data with columns for m/z, retention time (RT), feature ID, and sample intensities.
        tracer (str): Tracer code used in the experiment (e.g. "13C").
        mz_tol (float) : m/z tolerance (in ppm).
        rt_tol(float) : Retention time tolerance (in sec).
        max_atoms (int|None):  Maximum number of tracer atoms to consider for isotopologues.
        database (None): DataFrame containing theoretical features with columns retention time (RT), metabolite names, and formulas.
    """
    def __init__(self, dataset : pd.DataFrame, tracer, mz_tol, rt_tol, max_atoms=None, database=None): 
        self.dataset = dataset 
        self._tracer = tracer
        self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer)
        self._mz_tol = mz_tol
        self._rt_tol = rt_tol
        self.max_atoms = max_atoms
        self.database = database
        self.features = {} # {sample_name: {feature_id: Feature object}}
        self.clusters = {} # {sample_name: {cluster_id: Cluster object}}
        
    @property
    def rt_tol(self):
        """
        Returns the retention time tolerance used for feature annotation.
        :return: float        
        """
        return self._rt_tol
    
    @rt_tol.setter
    def rt_tol(self, value):
        if not isinstance(value, (float)):
            raise ValueError("RT tolerance must be a number.")
        if self._rt_tol is None:
            raise ValueError("RT tolerance must be provided.") 
        self._rt_tol = value

    @property
    def tracer(self):
        """
        Returns the tracer used for the experiment.
        :return: str 
        """
        return self._tracer

    @property
    def mz_tol(self):
        """
        Returns the m/z tolerance used for feature annotation.
        :return: float 
        """
        return self._mz_tol
    
    @mz_tol.setter
    def mz_tol(self, value):
        if not isinstance(value, (float)):
            raise ValueError("mz tolerance must be a number.")
        if self._mz_tol is None:
            raise ValueError("mz tolerance must be provided.") 
        self._mz_tol = value
        

    @property
    def tracer_element(self):
        """
        Returns the tracer element used in the experiment.
        :return: str 
        """
        return self._tracer_element
    
    @property
    def tracer_idx(self):
        """
        Returns the tracer index used in the experiment.
        :return: int 
        """
        return self._tracer_idx

    def initialize_experimental_features(self):
        """
        Initialize Feature objects from the dataset and organize them by sample.
        Each feature is created with its retention time, m/z, tracer, intensity, and sample name.

        :param dataset: DataFrame containing experimental data with columns for m/z, retention time (RT), 
                        feature ID, and sample intensities.
        """
        dataset = self.dataset.set_index(["mz", "rt", "id"])
        
        for idx, _ in dataset.iterrows():
            mz = idx[0]
            rt = idx[1]
            id = idx[2]
            
            for sample in dataset.columns:
                # Extract the intensity for each sample in the dataset
                intensity = dataset.loc[idx, sample]

                # Initialize the experimental features for each sample
                feature = Feature(
                    rt=rt, mz=mz, tracer=self.tracer,
                    feature_id=id, 
                    intensity=intensity,
                    sample=sample,
                    tracer_element=self.tracer_element,
                    )
                
                # Add the feature in the list corresponding to the sample
                if sample not in self.features:
                    self.features[sample] = {}
                self.features[sample][id] = feature
        
        features_count = len(next(iter(self.features.values())))
        logger.info(f"Initialized {features_count} features for {len(self.features)} samples")

# if __name__ == "__main__":
#     # from isogroup.base.io import IoHandler
#     from isogroup.base.targeted_experiment import TargetedExperiment
#     # io= IoHandler()
#     # data= io.read_dataset(r"..\..\data\dataset_test_XCMS.txt")
    
#     # database = io.read_database(r"..\..\data\database.csv")
    
#     test = TargetedExperiment(data, tracer="13C", mz_tol=5, rt_tol=10, database=database)
#     test.initialize_experimental_features()
#     print(test.database.theoretical_features)
#     io.export_theoretical_database(theoretical_db)
    # test.initialize_experimental_features()
    # print(test.database.theoretical_features)
    # # print(test.features["C13_WT_2"])