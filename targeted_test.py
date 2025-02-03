from isogroup.base.database import Database
from isogroup.base.experiment import Experiment
import pandas as pd


# WORKPLAN

# But: faire une analyse ciblé d'un échantillon à partir d'une database

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
experiment.initialize_samples()
for key in experiment.samples.keys():
    print(experiment.samples[key].dataset)
