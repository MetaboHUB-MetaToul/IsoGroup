import pandas as pd
from isogroup.base.database import Database
from isogroup.base.sample import Sample
from isogroup.base.feature import Feature

class Experiment:

    def __init__(self, dataset: pd.DataFrame, database: Database = None):
        self.experiment = dataset
        self.database = database
        self.samples: dict = {}
        self.annotated_experiment: None | pd.DataFrame = None
        self.mz_tol: None | float = None
        self.rt_tol: None | float = None

    def annotate_experiment(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        And calculate the mz error and the rt error
        """
        annotated_index = []  # Store index

        for idx, _ in self.experiment.iterrows():
            mz = idx[0]
            rt = idx[1]
            identity = idx[2]
            dummy_feature = Feature(rt=rt, mz=mz, intensity=None)
            annotated_row = [mz, rt, identity, [], [], [], []]

            
            for feature in self.database.features:
                mz_db = float(feature.mz)
                rt_db = float(feature.rt)

                # Calculate the exact mz and rt error
                mz_error = abs(mz_db - dummy_feature.mz)
                rt_error = abs(rt_db - dummy_feature.rt)

                if (mz_error <= mz_tol and
                     rt_error <= rt_tol):

                    annotated_row[3].append(feature.metabolite)
                    annotated_row[4].append(feature.isotopologue)
                    
                    # Store the exact mz and rt error
                    annotated_row[5].append(mz_error)
                    annotated_row[6].append(rt_error)
            
            annotated_row[3] = str(annotated_row[3])
            annotated_row[4] = str(annotated_row[4])
            annotated_row[5] = str(annotated_row[5])
            annotated_row[6] = str(annotated_row[6])
            annotated_index.append(tuple(annotated_row))
        
        self.mz_tol = mz_tol
        self.rt_tol = rt_tol
        self.annotated_experiment = self.experiment.copy(deep=True)
        self.annotated_experiment.index = pd.MultiIndex.from_tuples(
            tuple(annotated_index),
            names=["mz", "rt", "id", "metabolite", "isotopologue", "mz_error", "rt_error"])

    def initialize_samples(self):
        """
        Initialize the samples from the dataset.
        :return:
        """
        data = self.annotated_experiment if (self.annotated_experiment is not None) else self.experiment
        for sample in data.columns:
            self.samples[sample] = Sample(dataset=data[[sample]], sample_type="test")
