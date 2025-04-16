from isogroup.base.misc import Misc

class Feature:

    def __init__(self, rt: float, mz: float, tracer: str, intensity:float|None, feature_id = None, counter_formula: list|None=None, formula: list|None=None, sample: str|None=None,
                 chemical: list|None=None, metabolite: list|None=None, isotopologue: list|None=None, mz_error: list|None=None, rt_error: list|None=None, **extra_dims):
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


