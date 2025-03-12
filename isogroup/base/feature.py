import re
from isocor.base import LabelledChemical

class Feature:

    def __init__(self, rt: float, mz: float, tracer: str, intensity:float|None, feature_id = None, formula: list|None=None, sample: str|None=None,
                 metabolite: list|None=None, annotation: list|None=None, isotopologue: list|None=None, mz_error: list|None=None, rt_error: list|None=None, **extra_dims):
        self.rt = float(rt)
        self.mz = float(mz)
        self.tracer = tracer
        self._tracer_element, self._tracer_idx = self._parse_strtracer(tracer)
        self.intensity = intensity
        self.feature_id = feature_id
        self.formula = formula
        self.sample = sample
        self.mz_error = mz_error
        self.rt_error = rt_error
        self.metabolite = metabolite
        self.annotation = [i.label for i in self.metabolite] if self.metabolite is not None else annotation
        self.isotopologue = isotopologue
        self.__dict__.update(extra_dims)
        self.is_adduct: tuple[bool, str] = (False, "")
        # self.in_cluster: bool = False

    @staticmethod
    def _parse_strtracer(str_tracer):
        """Parse the tracer code.

        Args:
            str_tracer (str): tracer code (e.g. "13C")

        Returns:
            tuple
                - (str) tracer element (e.g. "C")
                - (int) tracer index in :py:attr:`~data_isotopes`
        """
        try:
            tracer = re.search(r'(\d*)([A-Z][a-z]*)', str_tracer)
            count = int(tracer.group(1))
            tracer_el = tracer.group(2)
        except (ValueError, AttributeError):
            raise ValueError("Invalid tracer code: '{}'."
                             " Please check your inputs.".format(str_tracer))
        best_diff = float("inf")
        idx_tracer = None
        unexpected_msg = "Unexpected tracer code. Are you sure this isotope is "\
                         "in data_isotopes? '{}'".format(str_tracer)
        assert tracer_el in LabelledChemical.DEFAULT_ISODATA, unexpected_msg
        for i, mass in enumerate(LabelledChemical.DEFAULT_ISODATA[tracer_el]["mass"]):
            test_diff = abs(mass - count)
            if test_diff < best_diff:
                best_diff = test_diff
                idx_tracer = i
        assert best_diff < 0.5, unexpected_msg
        return (tracer_el, idx_tracer)

    def __repr__(self) -> str:
        """
        Return a string representation of the feature.
        :return: str
        """
        return (f"Feature(ID = {self.feature_id}, RT={self.rt}, Metabolite={self.annotation}, Isotopologue={self.isotopologue}, "
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

