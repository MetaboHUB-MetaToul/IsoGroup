import argparse
from isogroup.base.database import Database
from isogroup.base.targeted_experiment import TargetedExperiment
from isogroup.base.untargeted_experiment import UntargetedExperiment
from isogroup.base.io import IoHandler
from pathlib import Path
import pandas as pd


# -------------------
# Targeted processing
# -------------------

# def process_targeted(args):
#     """Processing function for targeted mode."""

#     # load data file
#     inputdata = Path(args.inputdata)

#     if not inputdata.exists():
#         raise FileNotFoundError(f"File {inputdata} does not exist")

#     # load database file
#     if args.D is None:
#         raise ValueError("No database file provided. Use -D <file.csv>")
#     database = Path(args.D)
#     if not database.exists():
#         msg = f"File {database} does not exist"
#         raise FileNotFoundError(msg)
    
#     #Check tolerances
#     if args.mztol is None or args.rttol is None:
#         raise ValueError("Both --mztol and --rttol must be provided for targeted mode.")
    
#     # Load data and database
#     db_data = pd.read_csv(database, sep=";")
#     database = Database(dataset=db_data, tracer=args.tracer)

#     data = pd.read_csv(inputdata, sep="\t").set_index(["mz", "rt", "id"])

#     experiment = TargetedExperiment(dataset=data, database=database, tracer=args.tracer)
#     experiment.annotate_experiment(mz_tol=args.mztol, rt_tol=args.rttol)
#     experiment.clusterize()

#     # Set working directory from output path)
#     if args.output:
#         output = Path(args.output).resolve()

#         # Create the output directory 
#         res_dir = output.parent / "res"
#         res_dir.mkdir(parents=True, exist_ok=True)
#         output = res_dir / output.name
#         print(f"Results will be saved to: {res_dir}")

#         experiment.export_features(filename=output.with_suffix('.features.tsv'))
#         experiment.export_clusters(filename=output.with_suffix('.annotated_clusters.tsv'))
#         experiment.clusters_summary(filename=output.with_suffix('.clusters_summary.tsv'))
#     else:
#         raise ValueError("No output file provided")

def process_targeted(args):

    # load data file
    io = IoHandler(dataset=args.inputdata, tracer=args.tracer, database=args.D, mz_tol=args.mztol, rt_tol=args.rttol, outputs_path=args.output)
    io.read_dataset()
    io.initialize_experimental_features()
    io.read_database()

    io.create_output_directory()

    targeted_experiment= TargetedExperiment(
        features=io.features,
        database=io.database,
        mz_tol=io.mz_tol,
        rt_tol=io.rt_tol)
    targeted_experiment.annotate_features()
    targeted_experiment.clusterize()
    
    io.export_annotated_features()
    io.export_clusters(targeted_experiment.clusters)
    io.clusters_summary(targeted_experiment.clusters)
    print(f"Results will be saved to: {io.outputs_path}")

# ---------------------
# Untargeted processing
# ---------------------

# def process_untargeted(args):
#     """Processing function for untargeted mode."""

#     # load data file
#     inputdata = Path(args.inputdata)
#     if not inputdata.exists():
#         raise FileNotFoundError(f"File {inputdata} does not exist")

#     # Check if arguments are provided
#     if args.ppm_tol is None or args.rt_window is None:
#         raise ValueError("Both --ppm_tol and --rt_window must be provided for untargeted mode.")
    
#     data = pd.read_csv(inputdata, sep="\t").set_index(["mz", "rt", "id"])

#     # Output directory
#     res_dir = inputdata.parent / "res"
#     res_dir.mkdir(parents=True, exist_ok=True)

#     log_path = (Path(args.output).with_suffix('.log') if args.output else res_dir / f"{inputdata.stem}_untargeted.log")
#     print(f"Log file will be saved to: {log_path}")

#     experiment = UntargetedExperiment(dataset=data, tracer=args.tracer, log_file=str(log_path))
    
#     experiment.build_final_clusters(
#         RTwindow=args.rt_window,
#         ppm_tolerance=args.ppm_tol,
#         max_atoms=args.max_atoms,
#         verbose=args.verbose,
#         keep_best_candidate=args.kbc,
#         keep_richest=args.kr,    
#         )

#     if args.output:
#         # If user provided an output name
#         output = res_dir / Path(args.output).name
#     else:
#         # If no output name provided, generate one
#         base = inputdata.stem
#         output_name = f"{base}_clusters_RT{args.rt_window}_ppm{args.ppm_tol}.tsv"
#         output = res_dir / output_name

#     experiment.export_clusters_to_tsv(filepath=output)
#     experiment.export_features(filename=output.with_suffix('.features.tsv'))
    
#     print(f"Results will be saved to: {res_dir}")

def process_untargeted(args):
    io= IoHandler(dataset=args.inputdata, tracer=args.tracer, mz_tol=args.ppm_tol, rt_tol=args.rt_window, max_atoms=args.max_atoms, outputs_path=args.output)
    io.read_dataset()
    io.initialize_experimental_features()
    io.create_output_directory()


    untargeted_experiment= UntargetedExperiment(
        features=io.features,
        tracer=io.tracer,
        tracer_element=io._tracer_element,
        tracer_idx=io._tracer_idx,
        ppm_tolerance=io.mz_tol,
        rt_window=io.rt_tol,
        max_atoms=io.max_atoms,
        log_file=None)
    

    untargeted_experiment.build_final_clusters(
        verbose=args.verbose,
        keep_best_candidate=args.kbc,
        keep_richest=args.kr,)
    
    io.export_unannotated_features()
    io.clusters_to_dataframe(untargeted_experiment.clusters)

# -------------------
# CLI setup
# -------------------
def build_parser_targeted():
    parser = argparse.ArgumentParser(
        prog='isogroup_targeted',
        description='Annotation of isotopic datasets',
    )

    parser.add_argument("inputdata", help="measurements file to process")
    parser.add_argument("-t", "--tracer", type=str, required=True,
                        help='the isotopic tracer (e.g. "13C")')
    parser.add_argument("-D", type=str, required=True,
                        help="path to database file (csv)")
    parser.add_argument("--mztol", type=float, required=True,
                        help='mz tolerance in ppm (e.g. "5")')
    parser.add_argument("--rttol", type=float, required=True,
                        help='rt tolerance (e.g. "10")')
    parser.add_argument("-o", "--output", type=str, required=True,
                        help='output file for the clusters')
    parser.set_defaults(func=process_targeted)
    return parser

def build_parser_untargeted():
    parser = argparse.ArgumentParser(
        prog='isogroup_untargeted',
        description='Clustering of isotopic datasets',
    )
    parser.add_argument("inputdata", help="measurements file to process")
    parser.add_argument("-t", "--tracer", type=str, required=True,
                        help='the isotopic tracer (e.g. "13C")')
    parser.add_argument("--ppm_tol", type=float, required=True,
                        help='mz tolerance in ppm for clustering (e.g. "5")')
    parser.add_argument("--rt_window", type=float, required=True,
                        help='rt tolerance for clustering (e.g. "10")')
    parser.add_argument("--max_atoms", type=int, default=None,
                        help='maximum number of tracer atoms in a molecule (e.g. "20"), optional')
    parser.add_argument("--kbc", type=bool, default=False,
                        help='keep only the best candidate among overlapping clusters during clustering (default: False)')
    parser.add_argument("--kr", type=bool, default=True,
                        help='keep only the richest cluster among overlapping clusters during clustering (default: True)')
    parser.add_argument("-o", "--output", type=str,
                        help='output file for the clusters')
    parser.add_argument("-v", "--verbose", action="store_true",
                        help='enable verbose logging')
    parser.set_defaults(func=process_untargeted)
    return parser

# ---------------------
# CLI entry point
# ---------------------
def main_targeted():
    parser = build_parser_targeted()
    args = parser.parse_args()
    args.func(args)

def main_untargeted():
    parser = build_parser_untargeted()
    args = parser.parse_args()
    args.func(args)


# TODO: Homogeneize the output files
