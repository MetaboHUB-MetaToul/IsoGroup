from isogroup.base.feature import Feature
from isocor.base import LabelledChemical
from isogroup.base.misc import Misc
import pandas as pd


class Database:

    def __init__(self, dataset: pd.DataFrame, tracer="13C"):
        self.dataset = dataset
        self.features: list = []
        self.tracer: str = tracer
        self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer)
        self.clusters: list = []

        _isodata: dict = LabelledChemical.DEFAULT_ISODATA
        self._delta_mz_tracer: float = _isodata[self._tracer_element]["mass"][1] - _isodata[
            self._tracer_element]["mass"][0]
        self._delta_mz_hydrogen: float = _isodata["H"]["mass"][0]

        self.initialize_theoretical_features()

    def __len__(self) -> int:
        return len(self.dataset)

    def initialize_theoretical_features(self):

        """
        Creates chemical labelled from isocor functions
        then initializes the theoretical features from a database file
        """
        for _, line in self.dataset.iterrows():
            chemical = LabelledChemical(
                formula=line["formula"],
                tracer=self.tracer,
                derivative_formula="",
                tracer_purity=[1.0, 0.0],
                correct_NA_tracer=False,
                data_isotopes=None,
                charge=line["charge"],
                label=line["metabolite"],
            )
            for i in range(chemical.formula[self._tracer_element] + 1):
                mz = (chemical.molecular_weight + i * self._delta_mz_tracer
                      + line["charge"] * self._delta_mz_hydrogen)
                feature = Feature(
                    rt=line["rt"],
                    mz=mz,
                    tracer=self.tracer,
                    intensity=None,
                    chemical=[chemical],
                    isotopologue=[i],
                    metabolite=[chemical.label]
                )
                self.features.append(feature)
