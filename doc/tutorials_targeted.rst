..  _Tutorials:

################################################################################
Tutorials : Targeted mode
################################################################################
.. .. seealso:: If you have a question that is not covered in the tutorials, have a look
             at the :ref:`faq`.


.. _Required input data files:

********************************************************************************
Required input data files
********************************************************************************

.. ..  _`Input data`:

.. Input data
.. --------------------------------------------------------------------------------

IsoGroup takes as input a list of MS features (m/z, retention time, intensity) and a database of
elemental formulas of metabolites, in order to extract and annotate :ref:`isotopic clusters <isotopic cluster>`.
The upstream untargeted extraction of MS features is not handled by IsoGroup,
but can be performed using any MS data processing software such as XCMS, MS-DIAL, or MZMine.
The input data should be provided in a specific format, as detailed below.

..  _`Measurements file`:

Measurements file
--------------------------------------------------------------------------------

**This file contains the MS features of each sample**,
i.e. the mass-to-charge ratios (m/z), retention times, and intensities of the detected peaks.

The measurement file is a TXT file with one row by feature and the following columns:

:id: The feature identifier, as it is referred in the MS data processing software; e.g. "feature_1".
:mz: The mass-to-charge ratio of the feature; e.g. "123.456".
:rt: The retention time of the feature, in seconds; e.g. "12.34".
:sample: The sample name, as it is referred in the MS data processing software; e.g. "Sample_1". The sample column **contains the intensity** (or peak area) of the detected feature in the corresponding sample.

:download:`Example file <../data/dataset_test_XCMS.txt>`.

..  _`Database file`:

Database file
--------------------------------------------------------------------------------

A database file is required to annotate the isotopic clusters. 

This file stores **elemental formulas of the metabolites**.

It is a CSV file with the following columns:

:metabolite: Metabolite name or abbreviation; e.g. "pyruvic acid" or "PYR".
:rt: Theoretical retention time of the metabolite using your analytical method, in seconds; e.g. "12.34".
:formula: Elemental formula of the metabolite moiety of the molecular entity that
          gives rise to the measured :ref:`isotopic cluster <isotopic cluster>`; e.g. "C\ :sub:`3`\ H\ :sub:`4`\ O\ :sub:`3`\ ".
:charge: Charge state of the detected ion; e.g. "-1".

:download:`Example file <../data/database.csv>`.


..  _`Annotation parameters`:

********************************************************************************
Annotation parameters
********************************************************************************
.. Annotation parameters
.. ================================================================================

IsoGroup provides flexible options to adapt to various experimental conditions, such as isotopic tracer, mass accuracy, or data processing quality.

:Measurements file: Path to the :ref:`Measurements file`.
:Database file: Path to the :ref:`Database file` containing the elemental formulas of metabolites.
:Isotopic tracer: The tracer used for your experiment.
:ppm tolerance: The mass accuracy allowed for the annotation of isotopic clusters, in ppm (parts per million).
:rt tolerance: The retention time tolerance for the annotation of isotopic clusters compared to the theoretical retention time of the metabolite, in seconds.
:Output data path: Path to the :ref:`Output data`. A log file with the same name will be created in the same directory, with a ‘.log’ extension.
:Verbose logs: If set, the console and the log-file will contain all information necessary to check intermediate results of the annotation process.


..  _`Output data`:

********************************************************************************
Output files
********************************************************************************

.. Output files
.. ================================================================================

A ``res`` directory is created, which includes the following files:

..  _`Feature file`:
Feature file (``.features.tsv``)
--------------------------------------------------------------------------------

Contains all the annotated or still unknown feature of the measurement file, with the following columns:

:feature_id: Identifier of the feature, as it was provided in the :ref:`Measurements file`.
:mz: Mass-to-charge ratio of the feature, as it was provided in the :ref:`Measurements file`.
:rt: Retention time of the feature, as it was provided in the :ref:`Measurements file`.
:metabolite: Name of the metabolite corresponding to the annotated feature, as it was provided in the :ref:`Database file`.
:isotopologue: The index of the isotopologues of the metabolite, as it was annotated. 
:mz_error: The mass error between the annotated feature and the theoretical m/z of the metabolite, in ppm (parts per million).
:rt_error: The retention time error between the annotated feature and the theoretical rt of the metabolite, in seconds.
:sample: Name of the sample, as it was provided in the :ref:`Measurements file`.
:intensity: The intensity of the feature in the sample, as it was provided in the :ref:`Measurements file`.

.. warning:: A single feature may be annotated with multiple metabolites.
             In such cases, the annotation-related columns ('metabolite', 'isotopologue', 'mz_error', 'rt_error') will contain multiple values separated by commas.
             It is up to the user to handle these cases carefully when interpreting the data, using the calculated error values as guidance.


Cluster file (``.clusters.tsv``)
--------------------------------------------------------------------------------

Contains the annotated isotopic clusters, including status information (completeness, duplications, etc):

:cluster_id: Identifier of the isotopic cluster (generated by IsoGroup).
:metabolite: Name of the metabolite corresponding to the cluster.
:feature_id: Identifier of the feature, as it was provided in the :ref:`Measurements file`.
:mz: Mass-to-charge ratio of the feature, as it was provided in the :ref:`Measurements file`.
:rt: Retention time of the feature, as it was provided in the :ref:`Measurements file`.
:feature_potential_metabolite: Potential annotation of the feature, as it was provided in the :ref:`Feature file`.
:isotopologues: The index of the isotopologues of the metabolite corresponding to the cluster.
:mz_error: The mass error between the annotated feature and the theoretical m/z of the metabolite, in ppm (parts per million).
:rt_error: The retention time error between the annotated feature and the theoretical rt of the metabolite, in seconds.
:sample: Name of the sample, as it was provided in the :ref:`Measurements file`.
:intensity: The intensity of the feature in the sample, as it was provided in the :ref:`Measurements file`.
:status: Status of the cluster, it can be "Ok", "Incomplete" and/or "Duplicated isotopologues"
:missing_isotopologues: List of missing isotopologues in the cluster, if any.
:duplicated_isotopologues: List of duplicated isotopologues in the cluster, if any.
:in_another_cluster: List of other clusters in which a feature is also present, if any.

.. warning:: A single feature can be retrieved in multiple clusters if it is annotated with multiple metabolites.
           

Cluster summary file (``.summary.tsv``)
--------------------------------------------------------------------------------

Provides an overview of each cluster, with the following columns:

:cluster_id: Identifier of the isotopic cluster (generated by IsoGroup).
:name: Name of the metabolite corresponding to the cluster.
:number_of_features: Number of features in the cluster.
:isotopologues: The isotopologues that make up the cluster, listed by their indices or identifiers.
:status: Status of the cluster, it can be "Ok", "Incomplete", "Duplicated isotopologues".
:missing_isotopologues: List of missing isotopologues in the cluster, if any
:duplicated_isotopologues: List of duplicated isotopologues in the cluster, if any.
number_of_samples: number of samples in which the cluster is present.

.. warning:: The summary file is intended to provide a quick overview of the clusters and their status.
             It does not contain all the details of the features, but rather a high-level summary of the clusters.


Theoretical database (``.theoretical_db.tsv``)
--------------------------------------------------------------------------------

Contains the theoretical isotopologues generated from the :ref:`Database file`, with the following columns:

:mz: Theoretical mass-to-charge ratio of the isotopologue.
:rt: Theoretical retention time of the isotopologue, as provided in the :ref:`Database file`.
:metabolite: Name of the metabolite corresponding to the isotopologue, as provided in the :ref:`Database file`.
:isotopologue: The theoretical index of the isotopologue.
:formula: Elemental formula of the isotopologue, as provided in the :ref:`Database file`.


Log file (``.log``)
--------------------------------------------------------------------------------

Extensive information on the annotation process can be found in the log file if ‘Verbose logs’ option has been checked.


Warning and error messages
--------------------------------------------------------------------------------

Error messages are explicit. You should examine carefully any warning/error message.
After correcting the problem, you might have to restart IsoGroup and re-run the analysis.