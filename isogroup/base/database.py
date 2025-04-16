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
        self.export_database(filename="Database_isotopic.tsv")

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
                    metabolite=[chemical.label],
                    formula = line["formula"],
                )
                self.features.append(feature)


    def export_database(self, filename = None):
        """
        Create a DataFrame to summarize the database
        Optionnal: Export the DataFrame to a tsv file if a filename is provided with samples in column
        """

        # Create a DataFrame to summarize the theoretical features
        feature_data = []
        for feature in self.features:
            feature_data.append({
                "mz": feature.mz,
                "rt": feature.rt,
                "metabolite": ', '.join(feature.metabolite),
                "isotopologue": ', '.join(map(str, feature.isotopologue)),
                "formula": feature.formula,
                })

        df = pd.DataFrame(feature_data)

        # Export the DataFrame to a tsv file if a filename is provided
        if filename:
            df.to_csv(filename, sep="\t", index=False)

            return df