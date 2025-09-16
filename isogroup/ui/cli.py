import argparse
from isogroup.base.database import Database
from isogroup.base.targeted_experiment import Experiment
from isogroup.base.untargeted_experiment import UntargetedExperiment
from pathlib import Path
import pandas as pd


# -------------------
# Targeted processing
# -------------------

def process_targeted(args):

    # load data file
    inputdata = Path(args.inputdata)

    if not inputdata.exists():
        raise FileNotFoundError(f"File {inputdata} does not exist")

    # load database file
    if args.D is None:
        raise ValueError("No database file provided. Use -D <file.csv>")
    database = Path(args.D)
    if not database.exists():
        msg = f"File {database} does not exist"
        raise FileNotFoundError(msg)
    
    #Check tolerances
    if args.mztol is None or args.rttol is None:
        raise ValueError("Both --mztol and --rttol must be provided for targeted mode.")
    
    db_data = pd.read_csv(database, sep=";")
    database = Database(dataset=db_data, tracer=args.tracer)

    data = pd.read_csv(inputdata, sep="\t").set_index(["mz", "rt", "id"])

    experiment = Experiment(dataset=data, database=database, tracer=args.tracer)
    experiment.annotate_experiment(mz_tol=args.mztol, rt_tol=args.rttol)
    experiment.clusterize()

    # Set working directory from output path
    if args.output:
        output = Path(args.output).resolve()

        # Create the output directory 
        res_dir = output.parent / "res"
        res_dir.mkdir(parents=True, exist_ok=True)
        output = res_dir / output.name
        print(f"Results will be saved to: {res_dir}")

        experiment.export_features(filename=output.with_suffix('.features.tsv'))
        experiment.export_clusters(filename=output.with_suffix('.annotated_clusters.tsv'))
        experiment.clusters_summary(filename=output.with_suffix('.clusters_summary.tsv'))
    else:
        raise ValueError("No output file provided")


# ---------------------
# Untargeted processing
# ---------------------

def process_untargeted(args):
    # Check if arguments are provided
    if args.ppm_tol is None or args.rt_window is None:
        raise ValueError("Both --ppm_tol and --rt_window must be provided for untargeted mode.")
    
    inputdata = Path(args.inputdata)
    if not inputdata.exists():
        raise FileNotFoundError(f"File {inputdata} does not exist")
    data = pd.read_csv(inputdata, sep="\t").set_index(["mz", "rt", "id"])
    experiment = UntargetedExperiment(dataset=data, tracer=args.tracer)
    experiment.build_final_clusters(RTwindow=args.rt_window, ppm_tolerance=args.ppm_tol)

    res_dir = Path(args.inputdata).parent / "res"
    res_dir.mkdir(parents=True, exist_ok=True)


    if args.output:
        # If user provided an output name
        output = res_dir / Path(args.output).name
    else:
        # If no output name provided, generate one
        base = Path(args.inputdata).stem
        output_name = f"{base}_clusters_RT{args.rt_window}_ppm{args.ppm_tol}.tsv"
        output = res_dir / output_name

    print(f"Results will be saved to: {res_dir}")

    experiment.export_clusters_to_tsv(filepath=output)


# -------------------
# Argument parsing
# -------------------
def parseArgs():
    parser = argparse.ArgumentParser(description='Annotation or clustering of isotopic datasets')
    parser.add_argument("inputdata", help="measurements file to process")
    parser.add_argument("-t", "--tracer", type=str, required=True,
                        help='the isotopic tracer (e.g. "13C")')
    parser.add_argument("--mode", type=str, choices=['targeted', 'untargeted'], required=True,
                        help='mode of operation: "targeted" or "untargeted"')
    parser.add_argument("-o", "--output", type=str,
                        help='output file for the clusters')
    
    # Targeted specific arguments
    parser.add_argument("-D", type=str,
                        help="path to database file (csv)")
    parser.add_argument("--mztol", type=float,
                        help='mz tolerance in ppm (e.g. "5"), required for targeted')
    parser.add_argument("--rttol", type=float,
                        help='rt tolerance (e.g. "10"), required for targeted')
    
    # Untargeted specific arguments
    parser.add_argument("--ppm_tol", type=float,
                        help='mz tolerance in ppm for clustering (e.g. "5"), required for untargeted')
    parser.add_argument("--rt_window", type=float,
                        help='rt tolerance for clustering (e.g. "10"), required for untargeted')
    
    return parser

# -------------------
# CLI entry point
# -------------------
def start_cli():
    parser = parseArgs()
    args = parser.parse_args()

    # Define the mode
    if args.mode == 'targeted':
        # Check required args for targeted
        if not all(hasattr(args, attr) for attr in ['D', 'mztol', 'rttol']):
            parser.error("For targeted mode, --D, --mztol and --rttol are required.")
        process_targeted(args)
    elif args.mode == 'untargeted':
        # Check required args for untargeted
        if not all(hasattr(args, attr) for attr in ['ppm_tol', 'rt_window']):
            parser.error("For untargeted mode, --ppm_tol and --rt_window are required.")
        process_untargeted(args)
    else:
        parser.error("Mode must be either 'targeted' or 'untargeted'.")

if __name__ == "__main__":
    start_cli()