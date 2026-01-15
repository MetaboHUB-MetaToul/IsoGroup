import argparse
import isogroup
from isogroup.base.targeted_experiment import TargetedExperiment
from isogroup.base.untargeted_experiment import UntargetedExperiment
from isogroup.base.io import IoHandler
import logging
from pathlib import Path

def _build_logger(args, output_path):
    """
    Build the logger object
    
    Args:
        args: command line arguments
    """
    _logger = logging.getLogger("IsoGroup")
    log_file = logging.FileHandler(f"{output_path}/log.txt", mode='w')
    cli_handle = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - '
                                  '%(message)s')
    log_file.setFormatter(formatter)
    cli_handle.setFormatter(formatter)
    if args.verbose :
        _logger.setLevel(logging.DEBUG)
    else:
        _logger.setLevel(logging.INFO)

    _logger.addHandler(cli_handle)
    _logger.addHandler(log_file)
    return _logger


# -------------------
# Targeted processing
# -------------------

def targeted_process(args):
    
    # load data file
    io = IoHandler()
    dataset = io.read_dataset(Path(args.inputdata))
    io.create_output_directory(Path(args.output))

    _logger = _build_logger(args, io.outputs_path)
    _logger.info("=============================================")
    _logger.info(f"Using IsoGroup targeted version: {isogroup.__version__}")    
    _logger.info("=============================================\n")
    _logger.info(f"Dataset loaded from {args.inputdata}")
    database = io.read_database(Path(args.database))
    _logger.info(f"Database loaded from {args.database}\n")

    targeted_experiment= TargetedExperiment(
        dataset=dataset,
        tracer=args.tracer,
        mz_tol=args.ppm_tol,
        rt_tol=args.rt_tol,
        database=database)
    
    _logger.info("-------- Parameters --------")
    _logger.info(f"  Tracer : {args.tracer}")
    _logger.info(f"  m/z tolerance (ppm) : {args.ppm_tol}")
    _logger.info(f"  RT tolerance (sec) : {args.rt_tol}")
    _logger.info("----------------------------\n")

    io.export_theoretical_database(targeted_experiment.database)

    targeted_experiment.run_targeted_pipeline()
    
    io.targ_export_features(targeted_experiment.features)
    io.targ_export_clusters(targeted_experiment.features, targeted_experiment.clusters)
    io.clusters_summary(targeted_experiment.clusters)
    _logger.info(f"Path to results files: {io.outputs_path}")

# ---------------------
# Untargeted processing
# ---------------------

def untargeted_process(args):
    io= IoHandler()
    dataset = io.read_dataset(Path(args.inputdata))
    io.create_output_directory(Path(args.output))

    _logger=_build_logger(args, io.outputs_path)
    _logger.info("=============================================")
    _logger.info(f"Using IsoGroup untargeted version: {isogroup.__version__}")    
    _logger.info("=============================================\n")
    _logger.info(f"Dataset loaded from {args.inputdata}\n")

    _logger.info("-------- Parameters --------")
    _logger.info(f"  Tracer : {args.tracer}")
    _logger.info(f"  m/z tolerance (ppm) : {args.ppm_tol}")
    _logger.info(f"  RT tolerance (sec) : {args.rt_tol}")
    _logger.info(f"  Max atoms : {args.max_atoms}")
    _logger.info("----------------------------\n")

    untargeted_experiment= UntargetedExperiment(
        dataset=dataset,
        tracer=args.tracer,
        mz_tol=args.ppm_tol,
        rt_tol=args.rt_tol,
        max_atoms=args.max_atoms,
        keep=args.keep)
    
    # untargeted_experiment.build_final_clusters(
    #     verbose=args.verbose,
    #     keep_best_candidate=args.kbc,
    #     keep_richest=args.kr,)
    untargeted_experiment.run_untargeted_pipeline()
    
    io.untarg_export_features(untargeted_experiment.features)
    io.untarg_export_clusters(untargeted_experiment.clusters)
    _logger.info(f"Path to results files: {io.outputs_path}")

# -------------------
# CLI setup
# -------------------
def build_parser_targeted():
    parser = argparse.ArgumentParser(
        prog='isogroup_targeted',
        description='Annotation of isotopic datasets',
    )

    parser.add_argument("inputdata", help="input dataset file")
    parser.add_argument("-t", "--tracer", type=str, required=True,
                        help='the isotopic tracer (e.g. "13C")')
    parser.add_argument("-D", "--database", type=str, required=True,
                        help="path to database file (csv)")
    parser.add_argument("-ppm", "--ppm_tol", type=float, required=True,
                        help='m/z tolerance in ppm (e.g. "5")')
    parser.add_argument("-rt", "--rt_tol", type=float, required=True,
                        help='retention time tolerance in sec (e.g. "10")')
    parser.add_argument("-o", "--output", type=str, required=True,
                        help='path to generate the output files')
    parser.add_argument("-v", "--verbose",
                        help='enable verbose logging', action="store_true")
    parser.set_defaults(func=targeted_process)
    return parser

def build_parser_untargeted():
    parser = argparse.ArgumentParser(
        prog='isogroup_untargeted',
        description='Clustering of isotopic datasets',
    )
    parser.add_argument("inputdata", help="input dataset file")
    parser.add_argument("-t", "--tracer", type=str, required=True,
                        help='the isotopic tracer (e.g. "13C")')
    parser.add_argument("-ppm", "--ppm_tol", type=float, required=True,
                        help='m/z tolerance in ppm for clustering (e.g. "5")')
    parser.add_argument("-rt","--rt_tol", type=float, required=True,
                        help='rt tolerance in sec for clustering (e.g. "10")')
    parser.add_argument("--max_atoms", type=int, default=None,
                        help='maximum number of tracer atoms in a molecule (e.g. "20"). OPTIONAL')
    # parser.add_argument("--kbc", type=bool, default=False,
    #                     help='keep only the best candidate among overlapping clusters during clustering (default: False)')
    # parser.add_argument("--kr", type=bool, default=True,
    #                     help='keep only the richest cluster among overlapping clusters during clustering (default: True)')
    parser.add_argument("-k","--keep", type=str, default="All",
                        help='strategy to deduplicate overlapping clusters: "longest", "closest_mz", "both", "All". OPTIONAL')
    parser.add_argument("-o", "--output", type=str, required=True,
                        help='path to generate the output files')
    parser.add_argument("-v", "--verbose", action="store_true",
                        help='enable verbose logging')
    parser.set_defaults(func=untargeted_process)
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

# -------------------
# Old Targeted processing
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

# ---------------------
# old Untargeted processing
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