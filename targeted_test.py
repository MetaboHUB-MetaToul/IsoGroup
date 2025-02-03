from isogroup.base.feature import Feature
from isogroup.base.sample import Sample
import pandas as pd

# data = pd.read_csv("data/data_test.csv", sep=";")
# echantillon = data[['id', 'mz', 'rt', 'C12_WT_1']]
# # print(echantillon)
# #
# # feature_list = []
# #
# # for index, row in echantillon.iterrows():
# #     feature_list.append(
# #         Feature(
# #             rt=row['rt'],
# #             mz=row['mz'],
# #             intensity=row['C12_WT_1'],
# #             id=row['id'])
# #     )
# #
# # print(feature_list)
#
# sample = Sample(dataset=echantillon, sample_type="C12")
# print(sample.dataset)

# WORKPLAN

# But: faire une analyse ciblé d'un échantillon à partir d'une database

# 1. Importer les données
# 1.1 Importer librairie pandas
import pandas as pd
# 1.2 utiliser méthode read_csv
data = pd.read_csv("data/data_test.csv", sep=";")
data_sample = data[["id", "mz", "rt", "C12_WT_1"]]
sample = Sample(dataset=data_sample, sample_type="C12")
print(sample.dataset)
# 1.3 Vérifier l'intégrité des données


# 2. Importer la database

# 3. Créer les features de la database, donc les théoretical features
# 4. Créer les features de l'échantillon
# 5. Faire le matching

