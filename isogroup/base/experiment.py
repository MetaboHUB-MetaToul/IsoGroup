# import pandas as pd
# from isogroup.base.database import Database
# from isogroup.base.sample import Sample
# from isogroup.base.feature import Feature

# class Experiment:

#     def __init__(self, dataset: pd.DataFrame, database: Database = None):
#         self.experiment = dataset
#         self.database = database
#         self.samples: dict = {}
#         self.annotated_experiment: None | pd.DataFrame = None
#         self.mz_tol: None | float = None
#         self.rt_tol: None | float = None
#         # self.annotated_clusters : dict = {}

#     def annotate_experiment(self, mz_tol, rt_tol):
#         """
#         Annotate the experiment features with the database within a given tolerance
#         And calculate the mz error and the rt error
#         """
#         annotated_index = []  # Store index

#         for idx, _ in self.experiment.iterrows():
#             mz = idx[0]
#             rt = idx[1]
#             identity = idx[2]
#             dummy_feature = Feature(rt=rt, mz=mz, intensity=None)
#             annotated_row = [mz, rt, identity, [], [], [], []]

            
#             for feature in self.database.features:
#                 mz_db = float(feature.mz)
#                 rt_db = float(feature.rt)

#                 # Calculate the exact mz and rt error
#                 mz_error = abs(mz_db - dummy_feature.mz)
#                 rt_error = abs(rt_db - dummy_feature.rt)

#                 if (mz_error <= mz_tol and
#                      rt_error <= rt_tol):

#                     annotated_row[3].append(feature.metabolite)
#                     annotated_row[4].append(feature.isotopologue)
                    
#                     # Store the exact mz and rt error
#                     annotated_row[5].append(mz_error)
#                     annotated_row[6].append(rt_error)
            
#             annotated_row[3] = str(annotated_row[3])
#             annotated_row[4] = str(annotated_row[4])
#             annotated_row[5] = str(annotated_row[5])
#             annotated_row[6] = str(annotated_row[6])
#             annotated_index.append(tuple(annotated_row))
        
#         self.mz_tol = mz_tol
#         self.rt_tol = rt_tol
#         self.annotated_experiment = self.experiment.copy(deep=True)
#         self.annotated_experiment.index = pd.MultiIndex.from_tuples(
#             tuple(annotated_index),
#             names=["mz", "rt", "id", "metabolite", "isotopologue", "mz_error", "rt_error"])

#     def initialize_samples(self):
#         """
#         Initialize the samples from the dataset.
#         :return:
#         """
#         data = self.annotated_experiment if (self.annotated_experiment is not None) else self.experiment
#         for sample in data.columns:
#             self.samples[sample] = Sample(dataset=data[[sample]], sample_type="test")

# # fonction qui prend un dataframe contenant des features avec des listes de métabolites et les sépare en plusieurs lignes pour chaque métabolite
#     def separate_metabolites(self):
#         if self.annotated_experiment is None:
#             raise ValueError("The experiment is not annotated. Please annotate the experiment first.")
#         data = self.annotated_experiment
#         new_data = pd.DataFrame(columns=data.columns)
#         for idx, row in data.iterrows():
#             metabolites = eval(row[3])
#             isotopologues = eval(row[4])
#             mz_errors = eval(row[5])
#             rt_errors = eval(row[6])
#             for i in range(len(metabolites)):
#                 new_row = [row[0], row[1], row[2], metabolites[i], isotopologues[i], mz_errors[i], rt_errors[i]]
#                 new_data.loc[len(new_data)] = new_row
#         self.annotated_experiment = new_data
#         print(self.annotated_experiment)


# # Créer un dictionnaire de cluster par métabolite unique

#     def create_cluster(self):

#         # Check if the experiment is annotated
#         if self.annotated_experiment is None:
#             raise ValueError("The experiment is not annotated. Please annotate the experiment first.")
        
#         dict_cluster = {}
#         for idx, _ in self.annotated_experiment.iterrows():
#             mz = idx[0]
#             rt = idx[1]
#             identity = idx[2]
#             metabolite = idx[3]
#             isotopologue = idx[4]
#             mz_error = idx[5]
#             rt_error = idx[6]
#             if metabolite not in dict_cluster:
#                 dict_cluster[metabolite] = []
#             dict_cluster[metabolite].append([mz, rt, identity, isotopologue, mz_error, rt_error])

#         print(dict_cluster)

# # Pour chaque feature, créer un cluster de features en fonction du métabolite
# # Pour chaque cluster, ajouter les données correspondantes dans le dictionnaire
# # Pour chaque métabolite, créer un cluster de features
# # Pour chaque feature, ajouter le feature au cluster si le métabolite est le même
# # Sinon, créer un nouveau cluster
# # Pour chaque cluster, ajouter le cluster au dictionnaire avec le métabolite comme clé
# # Pour chaque cluster, ajouter les données correspondantes dans le dictionnaire

"""
New version test
"""

import pandas as pd
from isogroup.base.database import Database
from isogroup.base.sample import Sample
from isogroup.base.feature import Feature

class Experiment:

    def __init__(self, dataset: pd.DataFrame, database: 'Database' = None):
        self.experiment = dataset
        self.database = database
        self.samples: dict = {}
        self.annotated_experiment: None | pd.DataFrame = None
        self.mz_tol: None | float = None
        self.rt_tol: None | float = None

    def annotate_experiment(self, mz_tol, rt_tol):
        """
        Annotate the experiment features with the database within a given tolerance
        And calculate the mz error and the rt error.
        """
        annotated_data = []  # List to store annotated rows

        # Loop over the experiment features
        for idx, _ in self.experiment.iterrows():
            mz, rt, identity = idx
            dummy_feature = Feature(rt=rt, mz=mz, intensity=None)
            metabolites, isotopologues, mz_errors, rt_errors = [], [], [], []

            # Loop over the database features
            for feature in self.database.features:
                mz_db = float(feature.mz)
                rt_db = float(feature.rt)

                # Calculate the exact mz and rt error
                mz_error = abs(mz_db - dummy_feature.mz)
                rt_error = abs(rt_db - dummy_feature.rt)

                # Check if the feature is within the tolerance
                if mz_error <= mz_tol and rt_error <= rt_tol:
                    metabolites.append(feature.metabolite)
                    isotopologues.append(feature.isotopologue)
                    mz_errors.append(mz_error)
                    rt_errors.append(rt_error)

            # Append annotated data
            annotated_data.append([mz, rt, identity, metabolites, isotopologues, mz_errors, rt_errors])

        # Create DataFrame with MultiIndex
        self.mz_tol = mz_tol
        self.rt_tol = rt_tol
        self.annotated_experiment = pd.DataFrame(
            annotated_data, 
            columns=["mz", "rt", "id", "metabolite", "isotopologues", "mz_error", "rt_error"]
        ).set_index(["mz", "rt", "id"])

        self.annotated_experiment.to_excel("annotated_experiment.xlsx", index=True)

    def initialize_samples(self):
        """
        Initialize the samples from the dataset.
        :return:
        """
        data = self.annotated_experiment if (self.annotated_experiment is not None) else self.experiment
        for sample in data.columns:
            self.samples[sample] = Sample(dataset=data[[sample]], sample_type="test")        

    def export_to_excel(self, filename="annotated_experiment.xlsx"):
        """ Exporte le DataFrame annoté au format Excel. """
        if self.annotated_experiment is not None:
            df_export = self.annotated_experiment.copy()

            # Convertir les listes en chaînes de caractères pour éviter les problèmes dans Excel
            #for col in ["metabolite", "isotopologues", "mz_error", "rt_error"]:
            #    df_export[col] = df_export[col].apply(lambda x: ", ".join(map(str, x)) if isinstance(x, list) else str(x))

            # Exporter en Excel
            df_export.to_excel(filename, index=True)
            print(f"✅ Exporté avec succès sous {filename}")
        else:
            print("⚠️ Erreur : Aucune annotation disponible. Exécute d'abord annotate_experiment()")


