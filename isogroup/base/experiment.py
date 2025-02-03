import pandas as pd
from isogroup.base.database import Database
from isogroup.base.sample import Sample
from isogroup.base.feature import Feature
class Experiment:

    def __init__(self, dataset: pd.DataFrame, database: Database=None):
        self.experiment = dataset
        self.database = database
        self.samples = {}
        self.annotated_experiment = None
        self.initialize_samples()
        self.mz_tol = None
        self.rt_tol = None

    def annotate_experiment(self, mz_tol, rt_tol):

        annotated_index = []
        for idx, _ in self.experiment.iterrows():
            mz = idx[0]
            rt = idx[1]
            id = idx[2]
            dummy_feature = Feature(rt=rt, mz=mz, intensity=None)
            annotated_row = [mz, rt, id, [], []]
            for feature in self.database.features:
                if ((abs(
                        float(feature.mz) - dummy_feature.mz) <=
                    mz_tol and abs(float(feature.rt) - dummy_feature.rt) <=
                        rt_tol)):
                    annotated_row[3].append(feature.metabolite)
                    annotated_row[4].append(feature.isotopologue)
            annotated_row[3] = str(annotated_row[3])
            annotated_row[4] = str(annotated_row[4])
            annotated_index.append(tuple(annotated_row))
        print(annotated_index)
        self.mz_tol = mz_tol
        self.rt_tol = rt_tol
        self.annotated_experiment = self.experiment.copy(deep=True)
        self.annotated_experiment.index = pd.MultiIndex.from_tuples(
            tuple(annotated_index),
            names=["mz", "rt", "id", "metabolite", "isotopologue"])

    def initialize_samples(self):
        """
        Initialize the samples from the dataset.
        :return:
        """
        for sample in self.experiment.columns:
            self.samples[sample] = Sample(dataset=self.experiment[[sample]],
                                          sample_type="test")
        print(self.samples)

if __name__ == "__main__":

    data = pd.read_csv(
        r"C:\Users\legregam\PycharmProjects\IsoGroup\data\dataset_test.txt",
        sep="\t"
    )
    data = data.set_index(["mz", "rt", "id"])
    db_data = pd.read_csv(r"C:\Users\legregam\PycharmProjects\IsoGroup\data"
                       r"\database_test.csv", sep=";")
    database = Database(dataset=db_data)
    experiment = Experiment(dataset=data, database=database)
    experiment.annotate_experiment(mz_tol=0.01, rt_tol=10)
    print(experiment.annotated_experiment)


