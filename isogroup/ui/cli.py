import argparse
from isogroup.base.database import Database
from isogroup.base.experiment import Experiment
from pathlib import Path
import pandas as pd


def process(args):

    # create environment
    inputdata = Path(args.inputdata)
    database = Path(args.D)

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
        experiment.export_clusters(filename=output)
        experiment.clusters_summary(filename=output.with_suffix('.summary.tsv'))
    else:
        experiment.export_clusters()
        experiment.clusters_summary()

def parseArgs():
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description='annotation of isotopic datasets')

    parser.add_argument("inputdata", help="measurements file to process")
    parser.add_argument("-D", type=str, help="path to database")

    parser.add_argument("-t", "--tracer", type=str, required=True,
                        help='the isotopic tracer (e.g. "13C")')
    parser.add_argument("--mztol", type=float,
                        help='mz tolerance in ppm (e.g. "5")')
    parser.add_argument("--rttol", type=float,
                        help='rt tolerance in seconds (e.g. "10")')
    parser.add_argument("-o", "--output", type=str,
                        help='output file for the clusters')
    parser.add_argument("-v", "--verbose",
                        help="flag to enable verbose logs", action='store_true')
    return parser


def start_cli():
    parser = parseArgs()
    args = parser.parse_args()
    process(args)