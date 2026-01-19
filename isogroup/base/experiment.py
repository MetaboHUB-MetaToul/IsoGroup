from isogroup.base.feature import Feature
from isogroup.base.misc import Misc
import pandas as pd
import logging

logger = logging.getLogger(f"IsoGroup")

class Experiment:
    """
    Represents a mass spectrometry experiment with experimental features.
        
    """
    def __init__(self, dataset : pd.DataFrame, tracer:str, ppm_tol:float, rt_tol:float, max_atoms:int=None, database:pd.DataFrame=None): 
        """
        :param dataset: DataFrame containing experimental data with columns for m/z, retention time (RT), feature ID, and sample intensities.
        :param tracer: Tracer code used in the experiment (e.g. "13C").
        :param ppm_tol: m/z tolerance (in ppm).
        :param rt_tol: Retention time tolerance (in sec).
        :param max_atoms: Maximum number of tracer atoms to consider for isotopologues. If None, IsoGroup automatically estimates the maximum number of isotopologues based on the feature m/z and tracer element. 
        :param database: DataFrame containing theoretical features with columns retention time (RT), metabolite names, and formulas.
        """
        self.dataset = dataset 
        self._tracer = tracer
        self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer)
        self._ppm_tol = ppm_tol
        self._rt_tol = rt_tol
        self.max_atoms = max_atoms
        self.database = database
        self.features = {} # {sample_name: {feature_id: Feature object}}
        self.clusters = {} # {sample_name: {cluster_id: Cluster object}}
        
    @property
    def rt_tol(self) -> float:
        """
        Returns the retention time tolerance used for feature annotation.      
        """
        return self._rt_tol
    
    @rt_tol.setter
    def rt_tol(self, value) -> float:
        """
        Sets the retention time tolerance used for feature annotation.
        """
        if not isinstance(value, (float)):
            raise ValueError("RT tolerance must be a number.")
        if self._rt_tol is None:
            raise ValueError("RT tolerance must be provided.") 
        self._rt_tol = value

    @property
    def tracer(self) -> str:
        """
        Returns the tracer used for the experiment.
        """
        return self._tracer

    @property
    def ppm_tol(self) -> float:
        """
        Returns the m/z tolerance (in ppm) used for feature annotation.
        """
        return self._ppm_tol
    
    @ppm_tol.setter
    def ppm_tol(self, value):
        """
        Sets the m/z tolerance (in ppm) used for feature annotation.
        """
        if not isinstance(value, (float)):
            raise ValueError("mz tolerance must be a number.")
        if self._ppm_tol is None:
            raise ValueError("mz tolerance must be provided.") 
        self._ppm_tol = value
        

    @property
    def tracer_element(self) -> str:
        """
        Returns the tracer element used in the experiment.
        """
        return self._tracer_element
    
    @property
    def tracer_idx(self) -> int:
        """
        Returns the tracer index used in the experiment.
        """
        return self._tracer_idx

    def initialize_experimental_features(self):
        """
        Initialize Feature objects from the dataset and organize them by sample.
        Each feature is created with its retention time, m/z, tracer, intensity, and sample name.
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
        logger.info(f"{features_count} features loaded per sample ({len(self.features)} sample(s)).\n")

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