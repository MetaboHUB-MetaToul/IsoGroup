class Feature:

    def __init__(self, rt: float, mz: float, intensity: dict[float]|None, feature_id = None, formula: list|None=None, sample: str|None=None,
                 metabolite: list|None=None, annotation: list|None=None, isotopologue: list|None=None, mz_error: list|None=None, rt_error: list|None=None, **extra_dims):
        self.rt = rt
        self.mz = mz
        self.intensity = intensity
        self.feature_id = feature_id
        self.formula = formula
        self.sample = sample
        self.mz_error = mz_error
        self.rt_error = rt_error
        self.metabolite = metabolite
        self.annotation = self.metabolite.label if self.metabolite else annotation
        self.isotopologue = isotopologue
        self.__dict__.update(extra_dims)
        self.is_adduct: tuple[bool, str] = (False, "")
        # self.in_cluster: bool = False

    def __repr__(self) -> str:
        """
        Return a string representation of the feature.
        :return: str
        """
        return (f", Feature(Identity = {self.feature_id}, RT={self.rt}, Annotation={self.annotation}, Isotopologue={self.isotopologue}, "
                f"mz={self.mz}, "
                f"intensity={self.intensity})")
    
    @property
    def in_cluster(self, clusters):
        """
        Check if the feature is in another cluster
        Return the cluster_id if the feature is in it
        """
    pass


# class AnnotatedFeature(Feature):

#     def __init__(self, mz_error: float, rt_error: float, **kwargs):
#         super().__init__(**kwargs)
#         self.mz_error = mz_error
#         self.rt_error = rt_error

