import argparse
from isogroup.base.database import Database
from isogroup.base.experiment import Experiment
from pathlib import Path
import pandas as pd


def process(args):

    # load data file
    inputdata = Path(args.inputdata)

    if not inputdata.exists():
        msg = f"File {inputdata} does not exist"
        raise FileNotFoundError(msg)

    # load database file
    if hasattr(args, 'D'):
        database = Path(args.D)
    else:
        msg = "No database file provided"
        raise ValueError(msg)
    if not database.exists():
        msg = f"File {database} does not exist"
        raise FileNotFoundError(msg)

    mztol = float(getattr(args, 'mztol', None))
    rttol = float(getattr(args, 'rttol', None))
    tracer = getattr(args, 'tracer', None)

    db_data = pd.read_csv(database, sep=";")

    database = Database(dataset=db_data, tracer="13C")
    data = pd.read_csv(inputdata, sep="\t")
    data = data.set_index(["mz", "rt", "id"])

    experiment = Experiment(dataset=data, database=database, tracer=tracer)
    experiment.annotate_experiment(mz_tol=mztol, rt_tol=rttol)
    experiment.clusterize()

    # Set working directory from output path
    if args.output:
        output = Path(args.output)
        experiment.export_features(filename=output.with_suffix('.features.tsv'))
        experiment.export_clusters(filename=output.with_suffix('.clusters.tsv'))
        experiment.clusters_summary(filename=output.with_suffix('.clusters_summary.tsv'))
    else:
        msg = "No output file provided"
        raise ValueError(msg)

def parseArgs():
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description='annotation of isotopic datasets')

    parser.add_argument("inputdata", help="measurements file to process")
    parser.add_argument("-D", type=str, help="path to database file (csv)") 

    parser.add_argument("-t", "--tracer", type=str, required=True,
                        help='the isotopic tracer (e.g. "13C")')
    parser.add_argument("--mztol", type=float, required=True,
                        help='mz tolerance in ppm (e.g. "5")')
    parser.add_argument("--rttol", type=float, required=True,
                        help='rt tolerance (e.g. "10")')
    parser.add_argument("-o", "--output", type=str,
                        help='output file for the clusters')
    # parser.add_argument("-v", "--verbose",
    #                     help="flag to enable verbose logs", action='store_true')
    return parser


def start_cli():
    parser = parseArgs()
    args = parser.parse_args()
    process(args)