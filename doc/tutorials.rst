..  _Tutorials:

################################################################################
Tutorials
################################################################################

.. seealso:: If you have a question that is not covered in the tutorials, have a look
             at the :ref:`faq`.


.. _First time using IsoGroup:

********************************************************************************
First time using IsoGroup
********************************************************************************

..  _`Input data`:

Input data
================================================================================

IsoGroup takes as input a list of MS features (m/z, retention time, intensity) and a database of
elemental formulas of metabolites, in order to extract and annotate :ref:`isotopic clusters <isotopic clusters>`.
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
:mz: The mass-to-charge ratio of the feature, as a float; e.g. "123.456".
:rt: The retention time of the feature, as a float, in seconds or minutes; e.g. "12.34".
:sample: The sample name, as it is referred in the MS data processing software; e.g. "Sample_1". The Sample column contains the intensity (or peak area) of the detected feature in the corresponding sample.

:download:`Example file <../isogroup/data/dataset_test.txt>`.

..  _`Database file`:

Database file
--------------------------------------------------------------------------------

A database file is required to annotate the isotopic clusters. 

This file stores **elemental formulas of the metabolites**.

It is a CSV file with the following columns:

:metabolite: Metabolite name or abbreviation; e.g. "pyruvic acid" or "PYR".
:rt: Theoretical retention time of the metabolite using your analytical method, as a float, in seconds or minutes; e.g. "12.34".
:formula: Elemental formula of the metabolite moiety of the molecular entity that
          gives rise to the measured :ref:`isotopic cluster <isotopic cluster>`; e.g. "C\ :sub:`3`\ H\ :sub:`4`\ O\ :sub:`3`\ ". See also :ref:`Formulas`.
:charge: Charge state of the detected ion; e.g. "-1".

:download:`Example file <../isogroup/data/database_test.csv>`.



Annotation parameters
================================================================================

IsoGroup provides options to adapt to many situations that can be encountered
in terms of isotopic tracer, mass accuracy, quality of data processing, etc.

:Measurements file: Path to the :ref:`Measurements file`.
:Database file: Path to the :ref:`Database file` containing the elemental formulas of metabolites.
:Isotopic tracer: The tracer used for your experiment.
:m/z tolerance: The mass accuracy allowed for the annotation of isotopic clusters, in ppm (parts per million).
:rt tolerance: The retention time tolerance for the annotation of isotopic clusters compared to the theoretical retention time of the metabolite, in seconds or minutes.
:Output data path: Path to the :ref:`Output data`. A res folder will be created in the same directory.


..  _`Output data`:

Output files
================================================================================

Feature file
--------------------------------------------------------------------------------

A file containing all the annotated or still unknown feature of the measurement file is created in the res folder.

The feature file is a TSV file with the suffix ".features" and contains the following columns:

:feature_id: Identifier of the feature, as it was provided in the :ref:`Measurements file`.
:mz: Mass-to-charge ratio of the feature, as it was provided in the :ref:`Measurements file`.
:rt: Retention time of the feature, as it was provided in the :ref:`Measurements file`.
:metabolite: Name of the metabolite corresponding to the annotated feature, as it was provided in the :ref:`Database file`.
:isotopologues: The index of the isotopologues of the metabolite, as it was annotated. 
:mz_error: The mass error between the annotated feature and the theoretical m/z of the metabolite, in ppm (parts per million).
:rt_error: The retention time error between the annotated feature and the theoretical rt of the metabolite, in seconds or minutes.
:sample: Name of the sample, as it was provided in the :ref:`Measurements file`.
:intensity: The intensity of the feature in the sample, as it was provided in the :ref:`Measurements file`.

.. warning:: A single feature may be annotated with multiple metabolites
             In such cases, the annotation-related columns ('metabolite', 'isotopologue', 'mz_error', 'rt_error') will contain multiple values separated by commas.
             It is up to the user to handle these cases carefully when interpreting the data, using the calculated error values as guidance.


Cluster file
--------------------------------------------------------------------------------

A file containing the annotated isotopic clusters is created in the same directory as the Feature file.

The cluster file is a TSV file with the suffix ".annotated_clusters" and contains the following columns:

:cluster_id: Identifier of the isotopic cluster, as it was generated by IsoGroup.
:metabolite: Name of the metabolite corresponding to the cluster.
:feature_id: Identifier of the feature, as it was provided in the :ref:`Measurements file`.
:mz: Mass-to-charge ratio of the feature, as it was provided in the :ref:`Measurements file`.
:rt: Retention time of the feature, as it was provided in the :ref:`Measurements file`.
:feature_potential_metabolite: Potential annotation of the feature, as it was provided in the :ref:`Feature file`.
:isotopologues: The index of the isotopologues of the metabolite corresponding to the cluster.
:mz_error: 
:rt_error: 
:sample: Name of the sample, as it was provided in the :ref:`Measurements file`.
:intensity: The intensity of the feature in the sample, as it was provided in the :ref:`Measurements file`.
:status: Status of the cluster, it can be "Ok", "Incomplete", ...
:missing_isotopologues: List of missing isotopologues in the cluster, if any.
:duplicated_isotopologues: List of duplicated isotopologues in the cluster, if any.
:in_another_cluster: List of other clusters in which a feature is also present, if any.

.. warning:: A single feature can be retrieved in multiple clusters if it is annotated with multiple metabolites.
           

Cluster summary file
--------------------------------------------------------------------------------

A summary file is created in the same directory with a '_summary' extension.

This file summarizes the consistency of annotated isotopic clusters. It contains the following columns:

:cluster_id: Identifier of the isotopic cluster, as it was generated by IsoGroup.
:name: Name of the metabolite corresponding to the cluster.
:number_of_features: Number of features in the cluster.
:isotopologues: The isotopologues that make up the cluster, listed by their indices or identifiers.
:status: Status of the cluster, it can be "Ok", "Incomplete", "Duplicated", etc.
:missing_isotopologues: List of missing isotopologues in the cluster, if any
:duplicated_isotopologues: List of duplicated isotopologues in the cluster, if any.
:samples: List of samples in which the cluster is present.

.. warning:: The summary file is intended to provide a quick overview of the clusters and their status.
             It does not contain all the details of the features, but rather a high-level summary of the clusters.
