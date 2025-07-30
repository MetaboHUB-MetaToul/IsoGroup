from __future__ import annotations
from isogroup.base.misc import Misc

class Feature:
    """
    Represents a mass spectrometry feature in the dataset.
    A feature is characterized by its retention time (RT), mass-to-charge ratio (m/z), intensity.
    It can also have associated chemical information, isotopologues, and other metadata.
    
    Args:
        rt (float): Retention time of the feature.
        mz (float): Mass-to-charge ratio of the feature.
        tracer (str): Tracer code (e.g. "13C").
        intensity (float|None): Intensity of the feature.
        feature_id (int|None): Unique identifier for the feature.
        counter_formula (list|None): Counter formula of the feature.
        formula (list|None): Formula of the feature.
        sample (str|None): Name of the sample where the feature was detected.
        chemical (list|None): List of chemical objects associated with the feature.
        metabolite (list|None): List of metabolite names associated with the feature.
        isotopologue (list|None): List of isotopologues for the annotated feature.
        mz_error (list|None): List of m/z errors for the annotated feature.
        rt_error (list|None): List of retention time errors for the annotatedfeature.
        extra_dims (dict): Additional dimensions to be added to the feature.
        
        """

    def __init__(self, rt: float, mz: float, tracer: str, intensity:float|None, feature_id = None, counter_formula: list|None=None, formula: list|None=None, sample: str|None=None,
                 chemical: list|None=None, metabolite: list|None=None, isotopologue: list|None=None, mz_error: list|None=None, rt_error: list|None=None, **extra_dims):
        """
        Initialize a Feature instance with mass spectrometry data and annotated information.
        """
        self.rt = float(rt)
        self.mz = float(mz)
        self.tracer = tracer
        self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer) 
        self.intensity = intensity
        self.feature_id = feature_id
        self.chemical = chemical if chemical is not None else []
        self.counter_formula = [i.formula for i in self.chemical] if self.chemical is not None else formula #formula ou [] ?
        self.formula = formula if formula is not None else []
        self.sample = sample
        self.mz_error = mz_error if mz_error is not None else []
        self.rt_error = rt_error if rt_error is not None else []
        self.metabolite = [i.label for i in self.chemical] if self.chemical is not None else metabolite #metabolite ou [] ?
        self.isotopologue = isotopologue if isotopologue is not None else []
        self.__dict__.update(extra_dims)
        self.is_adduct: tuple[bool, str] = (False, "")
        self.in_cluster = []


    def __repr__(self) -> str:
        """
        Return a string representation of the feature.
        :return: str
        """
        return (f"Feature(ID = {self.feature_id}, RT={self.rt}, Metabolite={self.metabolite}, Isotopologue={self.isotopologue}, "
                f"mz={self.mz}, "
                f"intensity={self.intensity})")
    
    # @property
    # def in_cluster(self):
    #     """
    #     Check if the feature is in another cluster
    #     Return the cluster_id if the feature is in it
    #     """
    #     pass


