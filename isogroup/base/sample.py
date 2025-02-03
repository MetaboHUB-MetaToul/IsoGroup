import pandas as pd


class Sample:
    def __init__(self, dataset: pd.DataFrame, sample_type: str):
        self.dataset = dataset
        self.sample_type = sample_type
        self.features = []
        self.initialize_features()

    def initialize_features(self):
        pass

#    @classmethod
        #    def from_file(cls, path: str, sample_type: str):
        #dataset = pd.read_csv(path, sep=";")
        #return cls(dataset, sample_type)
