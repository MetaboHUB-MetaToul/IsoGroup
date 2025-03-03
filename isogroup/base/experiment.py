import pandas as pd
from isogroup.base.database import Database
from isogroup.base.sample import Sample
from isogroup.base.feature import Feature
from isogroup.base.cluster import Cluster

class Experiment:

    def __init__(self, dataset: pd.DataFrame, database: 'Database' = None):
        self.experiment = dataset
        self.database = database
        self.samples: dict = {} # Dictionary to store the samples
        self.annotated_df: None | pd.DataFrame = None
        self.clusters_summary: None | pd.DataFrame = None # Summary of the clusters of annotated features in the experiment
        self.mz_tol: None | float = None
        self.rt_tol: None | float = None

        self.annotated_data: list = [] # List of experimental features
        self.annotated_clusters: list = []   # List of annotated clusters


    def annotate_experiment(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        And calculate the mz error and the rt error.
        MultiIndex DataFrame
        """

        for idx, _ in self.experiment.iterrows():
            mz = idx[0]
            rt = idx[1]
            identity = idx[2]

            # 28/02 NBU : Recover the intensities of the samples pour les colonnes qui ne sont pas mz, rt et identity
            intensities = {sample: self.experiment.loc[idx, sample] for sample in self.experiment.columns if sample not in ["mz", "rt", "identity"]}

            # Initialize the experiment features from the dataset
            exp_feature = Feature(
                rt=rt, mz=mz, 
                feature_id=identity, 
                intensity=intensities, # NBU 28/02 : Add the intensities of the samples as a dictionary
                formula=[], 
                metabolite=[], 
                isotopologue=[],
                mz_error=[], 
                rt_error=[]
                ) 
            
            for th_feature in self.database.features:
                mz_db = float(th_feature.mz)
                rt_db = float(th_feature.rt)

                # Calculate the exact mz and rt errors
                mz_error = (mz_db - exp_feature.mz)
                rt_error = (rt_db - exp_feature.rt)

                # Check if the feature is within tolerance
                if abs(mz_error) <= mz_tol and abs(rt_error) <= rt_tol:
                
                    exp_feature.metabolite.append(th_feature.metabolite)
                    exp_feature.isotopologue.append(th_feature.isotopologue)
                    exp_feature.formula = th_feature.formula
                    exp_feature.mz_error.append(mz_error)
                    exp_feature.rt_error.append(rt_error)

            # Store all experimental features
            self.annotated_data.append(exp_feature)

        self.mz_tol = mz_tol
        self.rt_tol = rt_tol

        # Initialize the samples from the dataset
        self.initialize_samples()

        # Create a DataFrame to summarize the annotated data
        self.to_df()


    def initialize_samples(self):
        """
        Initialize Sample objects for each sample in the experiment.
        For each sample, create a Sample object and initialize its features
        """
        for sample_name in self.experiment.columns:
            if sample_name not in ["mz", "rt", "identity"]:
                sample = Sample(dataset=self.experiment[[sample_name]], sample_name=sample_name)
                sample.initialize_features(self.annotated_data)
                self.samples[sample_name] = sample


    def to_df(self):
        """
        Return a DataFrame of the annotated data
        """
        if not self.annotated_data:
            raise ValueError("The experiment has not been annotated yet.")
        
        # Create a DataFrame
        df = pd.DataFrame([vars(feature) for feature in self.annotated_data])

        # Supprimer certaines colonnes
        df = df.drop(columns=["is_adduct", "in_cluster"], errors="ignore")  # colonnes à exclure pour l'instant

        # Export the DataFrame to a tsv file
        df.to_csv("annotated_data.tsv", sep="\t", index=False)

        return df