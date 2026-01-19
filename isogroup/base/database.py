from __future__ import annotations
from isogroup.base.feature import Feature
from isocor.base import LabelledChemical
# from isogroup.base.misc import Misc
import pandas as pd


class Database:
    """
    Represents a database of theoretical features for a specific tracer.

    """

    def __init__(self, dataset: pd.DataFrame, tracer: str, tracer_element: str):
        """
        :param dataset: DataFrame containing theoretical features with columns retention time (RT), metabolite names, and formulas.
        :param tracer: Tracer code (e.g. "13C") used to initialize the database.
        :param tracer_element:  Tracer element (e.g. "C") used.
        """
        self.dataset = dataset
        self.theoretical_features = []
        self.tracer = tracer
        self._tracer_element = tracer_element
        # self._tracer_element, self._tracer_idx = Misc._parse_strtracer(tracer)
        self.clusters = []

        _isodata: dict = LabelledChemical.DEFAULT_ISODATA
        self._delta_mz_tracer: float = _isodata[self._tracer_element]["mass"][1] - _isodata[
            self._tracer_element]["mass"][0]
        self._delta_mz_hydrogen: float = _isodata["H"]["mass"][0]

        self.initialize_theoretical_features()
        # self.export_database(filename="isotopic_db_export.tsv")

    def __len__(self) -> int:
        return len(self.dataset)

    def initialize_theoretical_features(self):
        """
        Creates chemical labelled objects from the dataset and initializes theoretical features.
        For each chemical, it generates features with isotopologues based on the tracer.
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
            for isotopologue in range(chemical.formula[self._tracer_element] + 1):
                mz = (chemical.molecular_weight + isotopologue * self._delta_mz_tracer
                      + line["charge"] * self._delta_mz_hydrogen)
                feature = Feature(
                    rt=line["rt"],
                    mz=mz,
                    tracer=self.tracer,
                    intensity=None,
                    chemical=[chemical],
                    # isotopologue=[isotopologue],
                    cluster_isotopologue={chemical.label: isotopologue},
                    metabolite=[chemical.label],
                    formula = line["formula"],
                )
                self.theoretical_features.append(feature)


    # def export_database(self, filename = None):
    #     """
    #     Summarize theoretical features into a DataFrame and optionally export it to a tsv file.
    #     :param filename: Name of the file to export the summary to
    #     :return: pd.DataFrame with the summary of the theoretical features
    #     """

    #     # Create a DataFrame to summarize the theoretical features
    #     feature_data = []
    #     for feature in self.features:
    #         feature_data.append({
    #             "mz": feature.mz,
    #             "rt": feature.rt,
    #             "metabolite": ', '.join(feature.metabolite),
    #             "isotopologue": ', '.join(map(str, feature.isotopologue)),
    #             "formula": feature.formula,
    #             })

    #     df = pd.DataFrame(feature_data)

    #     # Export the DataFrame to a tsv file if a filename is provided
    #     if filename:
    #         df.to_csv(filename, sep="\t", index=False)

    #         return df

# if __name__ == "__main__":
#     from isogroup.base.io import IoHandler
#     from pathlib import Path
#     io= IoHandler()
#     database_df= io.read_database(Path(r"..\..\data\database.csv"))
#     test_db = Database(dataset=database_df, tracer="13C", tracer_element="C")
#     test_db.initialize_theoretical_features()
#     for feature in test_db.theoretical_features:
#         print(feature.metabolite)
#         print(feature.cluster_isotopologue)
        
    