class Feature:

    def __init__(self, rt: float, mz: float, intensity: float|None, Fid = None,
                 mz_error=None, rt_error=None,
                 metabolite=None, isotopologue=None, **extra_dims):

        self.rt = rt
        self.mz = mz
        self.intensity = intensity
        self.Fid = Fid
        self.mz_error = mz_error
        self.rt_error = rt_error
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
        return (f", Feature(Fid = {self.Fid}, RT={self.rt}, Metabolite={self.metabolite}, Isotopologue={self.isotopologue}, "
                f"mz={self.mz}, "
                f"intensity={self.intensity})")



