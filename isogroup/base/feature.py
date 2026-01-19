from __future__ import annotations
from isogroup.base.misc import Misc

class Feature:
    """
    Represents a mass spectrometry feature in the dataset.
    A feature is characterized by its retention time (RT), mass-to-charge ratio (m/z), intensity.
    It can also have associated chemical information, isotopologues, and other metadata.
    """

    def __init__(self, rt:float, mz:float, tracer:str, intensity:float, feature_id:str=None, tracer_element = None, formula:list=None, sample:str=None,
                 chemical:list=None, metabolite:list=None, mz_error:list=None, rt_error: list|None=None, **extra_dims:dict):
        """
        :param rt: Retention time tolerance (in sec).
        :param mz: Mass-to-charge ratio of the feature.
        :param tracer: Tracer code (e.g. "13C").
        :param intensity: Intensity of the feature.
        :param feature_id: Unique identifier for the feature.
        :param formula: Formula of the feature.
        :param sample: Name of the sample the feature belongs to.
        :param chemical: List of chemical objects (LabelledChemical) associated with the feature.
        :param metabolite: List of metabolite names associated with the feature.
        :param isotopologue: List of isotopologues for the annotated feature.
        :param mz_error: List of m/z errors for the annotated feature.
        :param rt_error: List of retention time errors for the annotated feature.
        :param extra_dims: Additional dimensions to be added to the feature.
        """
        self.rt = float(rt)
        self.mz = float(mz)
        self.tracer = tracer
        self._tracer_element = tracer_element 
        # self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer) 
        self.intensity = intensity
        self.feature_id = feature_id
        self.chemical = chemical if chemical is not None else []
        self.counter_formula = [i.formula for i in self.chemical] if self.chemical is not None else formula # Counter formula of the feature
        self.formula = formula if formula is not None else []
        self.sample = sample
        self.mz_error = mz_error if mz_error is not None else []
        self.rt_error = rt_error if rt_error is not None else []
        self.metabolite = [i.label for i in self.chemical] if self.chemical is not None else metabolite #metabolite ou [] ?
        # self.isotopologue = isotopologue if isotopologue is not None else [] # Targeted version
        # self.cluster_isotopologue = isotopologue if isotopologue is not None else {} # Test 10/09 for Untargeted, dict, iso associated to the cluster -> check impact on Targeted version
        self.cluster_isotopologue = {} # Store the isotopologue number per cluster {cluster_name: isotopologue_number}
        self.__dict__.update(extra_dims)
        self.is_adduct: tuple[bool, str] = (False, "")
        self.in_cluster = []
        self.also_in = {}


    def __repr__(self) -> str:
        """
        Return a string representation of the feature.
        :return: str
        """
        return (f"Feature(ID = {self.feature_id}, RT={self.rt}, Metabolite={self.metabolite}, Isotopologue={self.cluster_isotopologue}, "
                f"mz={self.mz}, "
                f"intensity={self.intensity})")
    
    # @property
    # def in_cluster(self):
    #     """
    #     Check if the feature is in another cluster
    #     Return the cluster_id if the feature is in it
    #     """
    #     pass


