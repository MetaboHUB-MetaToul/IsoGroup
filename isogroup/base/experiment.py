import pandas as pd
from isogroup.base.database import Database
from isogroup.base.feature import Feature
from isogroup.base.cluster import Cluster
from isocor.base import LabelledChemical
import re

class Experiment:

    def __init__(self, dataset: pd.DataFrame, database: 'Database' = None, tracer=None):
        self.dataset = dataset
        self.database = database
        self.samples: dict = {} # Dictionary to store the samples
        self._mz_tol: None | float = None
        self._rt_tol: None | float = None

        self._tracer: None | str = tracer
        self._tracer_element, self._tracer_idx = self._parse_strtracer(tracer)
        #self.experimental_features : list = [] # List of experimental features
        #self.annotated_features: list = [] # List of experimental features after annotation # Inutile ? Modification de l'objet feature directement dans la liste experimental_features
        #self.annotated_clusters: list = []   # List of annotated clusters 
        #self.annotated_sample_clusters: list = [] # List of annotated custers for each sample
        self.clusters = None
        #self.annotated_experiment: None | pd.DataFrame = None # Rename
        

    @property
    def rt_tol(self):
        return self._rt_tol

    @property
    def tracer(self):
        return self._tracer

    @property
    def tracer_element(self):
        return self._tracer_element

    @property
    def mz_tol(self):
        return self._mz_tol
    
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

    def initialize_experimental_features(self):
        """
        Initialize the experimental features from the dataset
        """
        for idx, _ in self.dataset.iterrows():
            mz = idx[0]
            rt = idx[1]
            id = idx[2]

            # Extract the intensity for each sample in the dataset
            for sample in self.dataset.columns:
                if sample not in ["mz", "rt", "id"]:
                    intensity = self.dataset.loc[idx, sample]

                    # Initialize the experimental features for each sample
                    feature = Feature(
                        rt=rt, mz=mz, tracer=self.tracer,
                        feature_id=id, 
                        intensity=intensity,
                        metabolite=[],
                        annotation=[],
                        isotopologue=[],
                        mz_error=[],
                        rt_error=[],
                        formula=[],
                        sample=sample
                        ) 
                    
                    # Add the feature in the list corresponding to the sample
                    if sample not in self.samples:
                        self.samples[sample] = {}
                    self.samples[sample][id] = feature
                    
                    # Store all experimental features
                    #self.experimental_features.append(feature)


    def annotate_features(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        Calculate the mz error and the rt error
        """

        for sample in self.samples.values():
            for feature in sample.values():
                
                for db_feature in self.database.features:

                    # Calculate the exact mz and rt errors
                    mz_error = (db_feature.mz - feature.mz)
                    rt_error = (db_feature.rt - feature.rt)

                    # Check if the experimental feature is within tolerance
                    if abs(mz_error) <= mz_tol and abs(rt_error) <= rt_tol:
                        feature.metabolite.append(db_feature.metabolite[0])
                        feature.isotopologue.append(db_feature.isotopologue[0])
                        feature.annotation.append(db_feature.metabolite[0].label)
                        feature.formula.append(db_feature.metabolite[0].formula)
                        feature.mz_error.append(mz_error)
                        feature.rt_error.append(rt_error)

        self._mz_tol = mz_tol
        self._rt_tol = rt_tol


    def annotate_experiment(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        MultiIndex DataFrame
        """
        # Initialize the experimental features from the dataset
        self.initialize_experimental_features()

        # Annotate the experimental features
        self.annotate_features(mz_tol, rt_tol)



    def features_export(self, filename = None, sample_name = None):
        """
        Create a DataFrame to summarize the annotated data
        Optionnal: Export the DataFrame to a tsv file if a filename is provided with samples in column
        Optionnal: Export the Dataframe of only one sample if a sample name is provided
        """

        # Create a DataFrame to summarize the experimental features
        df = pd.DataFrame([vars(f) for f in self.experimental_features])
        df = df.drop(columns=["is_adduct", "in_cluster", "metabolite", "formula"], errors="ignore")

        columns_to_fix = ["annotation", "mz_error", "rt_error", "isotopologue"]
        for col in columns_to_fix:
            df[col] = df[col].apply(lambda x: ", ".join(map(str, x)) if isinstance(x, list) else str(x))
        df = df.pivot_table(index=["mz", "rt", "feature_id", "annotation", "mz_error", "rt_error", "isotopologue"],
                                  columns="sample", values="intensity",
                                  aggfunc="first").reset_index()
        
        self.annotated_experiment = df

        # Export the DataFrame to a csv file if a filename is provided
        if filename:
            df.to_csv(filename, sep="\t", index=False)

        # Export the Dataframe of only one sample if a sample name is provided
        if filename and sample_name:
            # Check if the sample name is in the DataFrame
            if sample_name not in df.columns:
                raise ValueError(f"The sample {sample_name} is not in the DataFrame")
            
            df = df[["mz", "rt", "feature_id", "annotation", "mz_error", "rt_error", "isotopologue", sample_name]]
            df.to_csv(filename, sep="\t", index=False)

        return df
    
    def clusterize(self):
        """
        Create unique clusters from annotated features based on their names.
        """
        cluster_names = []
        
        # Group features by metabolite
        for sample in self.samples.values():
            for feature in sample.values():
                cluster_names += feature.annotation

        cluster_names = set(cluster_names)
        
        # Create unique clusters
        self.clusters = {}
        for sample in self.samples.keys():
            self.clusters[sample] = {}
            for i,c in enumerate(cluster_names):
                features = self.get_features_from_name(c, sample)
                self.clusters[sample][c] = Cluster(features=features, cluster_id=f"C{i}", name=c)


    def get_features_from_name(self, name, sample_name:str):
        """
        Get a feature from the experiment by its name, in a given sample if provided
        """
        features = []
        for feature in self.samples[sample_name].values():
            if name in feature.annotation:
                features.append(feature)
        return features

    
    def export_clusters(self, filename = None, sample_name = None):
        """
        Create a DataFrame to summarize the annotated clusters
        Optionnal: Export the DataFrame to a tsv file if a filename is provided
        Optionnal: Export the Dataframe of only one sample if a sample name is provided
        """
        
        # Check if the sample name is in the DataFrame
        all_samples = list(self.samples.keys())
        if sample_name is not None:
            if sample_name not in all_samples:
                raise ValueError(f"Sample {sample_name} not found in annotated clusters. Available samples: {', '.join(all_samples)}")
        
        cluster_data = []
        for sample, clusters in self.clusters.items():
            for cname, cluster in clusters.items():
                for feature in cluster.features:
                    idx = [i for i,j in enumerate(feature.annotation) if j == cname][0]
                    cluster_data.append({
                        "cluster_id": cluster.cluster_id,
                        "annotation": cluster.name,
                        "feature_id": feature.feature_id,
                        "mz": feature.mz,
                        "rt": feature.rt,
                        "feature_possible_annotation": feature.annotation,
                        "isotopologue": feature.isotopologue[idx],
                        "mz_error": feature.mz_error[idx],
                        "rt_error": feature.rt_error[idx],
                        "sample": feature.sample,
                        "intensity": feature.intensity,
                        "is_complete": "TODO",#cluster.is_complete,
                        "missing_isotopologue": "TODO"# cluster.missing_isotopologue
                    })

        # Create a DataFrame to summarize the annotated clusters
        df = pd.DataFrame(cluster_data)

        # Export the DataFrame to a tsv file if a filename is provided
        if filename:
            df.to_csv(filename, sep="\t", index=False)

        return df


