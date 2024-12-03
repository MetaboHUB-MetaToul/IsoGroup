import pandas as pd


class Sample:
    def __init__(self, dataset: pd.DataFrame, type: str):
        self.dataset = dataset
        self.type = type
        self.clusters = []
        self.nb_clusters = len(self.clusters)

    def initialize_clusters(self):
        pass
