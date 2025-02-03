class Feature:

    def __init__(self, rt: float, mz: float, intensity: float|None,
                 metabolite=None, isotopologue=None, **extra_dims):

        self.rt = rt
        self.mz = mz
        self.intensity = intensity
        self.metabolite = metabolite
        self.isotopologue = isotopologue
        self.__dict__.update(extra_dims)
        self.is_adduct: tuple[bool, str] = (False, "")
        self.in_cluster: bool = False

    def __repr__(self) -> str:
        """
        Return a string representation of the feature.
        :return: str
        """
        return (f", Feature(RT={self.rt}, Metabolite={self.metabolite} "
                f"mz={self.mz}, "
                f"intensity={self.intensity})")



