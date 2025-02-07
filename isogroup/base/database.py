from isogroup.base.feature import Feature
from isocor.base import LabelledChemical
import pandas as pd


class Database:

    def __init__(self, dataset: pd.DataFrame, tracer="13C", tracer_element="C"):
        self.dataset = dataset
        self.features: list = []
        self.tracer: str = tracer
        self.tracer_element: str = tracer_element

        _isodata: dict = LabelledChemical.DEFAULT_ISODATA
        self._delta_mz_tracer: float = _isodata["C"]["mass"][1] - _isodata[
            "C"]["mass"][0]
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
                label=line["metabolite"]
            )
            for i in range(chemical.formula[self.tracer_element] + 1):
                mz = (chemical.molecular_weight + i * self._delta_mz_tracer
                      + line["charge"] * self._delta_mz_hydrogen)
                feature = Feature(
                    rt=line["rt"],
                    mz=mz,
                    intensity=None,
                    metabolite=chemical.label,
                    isotopologue=i
                )
                self.features.append(feature)
