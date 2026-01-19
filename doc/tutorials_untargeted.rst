..  _Tutorials:

################################################################################
Tutorials : Untargeted mode
################################################################################

.. .. seealso:: If you have a question that is not covered in the tutorials, have a look
             at the :ref:`faq`.


.. _Required input data files:

********************************************************************************
Required input data files
********************************************************************************

.. ..  _`Input data`:

.. Input data
.. ================================================================================

IsoGroup takes as input a list of MS features (m/z, retention time, intensity) in order to extract :ref:`isotopic clusters <isotopic cluster>` without prior knowledge. 
The upstream untargeted extraction of MS features is not handled by IsoGroup,
but can be performed using any MS data processing software such as XCMS, MS-DIAL, or MZMine.
The input data should be provided in a specific format, as detailed below.

..  _`Measurements file (untargeted)`:

Measurements file
--------------------------------------------------------------------------------

**This file contains the MS features of each sample**,
i.e. the mass-to-charge ratios (m/z), retention times, and intensities of the detected peaks.

The measurement file is a TXT file with one row by feature and the following columns:

:id: The feature identifier, as it is referred in the MS data processing software; e.g. "feature_1".
:mz: The mass-to-charge ratio of the feature, as a float; e.g. "123.456".
:rt: The retention time of the feature, as a float, in seconds; e.g. "12.34".
:sample: The sample name, as it is referred in the MS data processing software; e.g. "Sample_1". The Sample column contains the intensity (or peak area) of the detected feature in the corresponding sample.

:download:`Example file <../data/dataset_test_XCMS.txt>`.

..  _`Grouping parameters`:

********************************************************************************
Grouping parameters
********************************************************************************

IsoGroup provides flexible options to adapt to various experimental conditions, such as isotopic tracer, mass accuracy, or data processing quality

:Measurements file: Path to the :ref:`Measurements file (untargeted)`.
:Isotopic tracer: The tracer used for your experiment.
:ppm tolerance: The mass accuracy allowed for the grouping of isotopologues into isotopic clusters, in ppm (parts per million).
:rt tolerance: The retention time tolerance for the grouping of isotopologues into isotopic clusters, in seconds.
:Verbose: If set, the console and the log-file will contain all information necessary to check intermediate results of the annotation process.

Additional optional parameters
----------------------------------------------------------------------------------

IsoGroup also provides additional optional parameters to refine the grouping of isotopologues into isotopic clusters:

:Max atoms: The maximum number of tracer atoms expected for any molecule in your dataset. Restricting this parameter reduces the search space and thus the computation time. 
            By default, IsoGroup automatically estimates the maximum number of isotopologues based on the feature m/z and tracer element.
:keep: Strategy to deduplicate overlapping clusters. Possible values are:
  - ``longest``: When multiple clusters share subsets of features, this option keeps only the **largest** cluster and removes its strict subsets.
  - ``closest_mz``: If, in a cluster, an isotopologue have multiple feature candidates, this option keep only the feature with the **closest m/z to the expected theoretical m/z** (minimizing Δppm) for each isotopologue in overlapping clusters.
  - ``both``: applies both previous strategies.
  **If this parameter is not set, all clusters are kept, even if they share features.**

.. :Keep best candidate: *(bool, default = False)* If set to ``True``, only the best candidate feature is retained for each isotopologue in a cluster. The best candidate is defines as the one **closest to the expected theoretical m/z** (minimizing Δppm).
.. :Keep richest: *(bool, default = True)* When multiple clusters share subsets of features, this option keeps only the **largest (richest)** cluster and removes its strict subsets. If set to ``False``, all clusters are kept, even if they share features.


..  _`Output data`:

********************************************************************************
Output files
********************************************************************************

A ``res`` directory is created, which includes the following files:

Feature file (``.features.tsv``)
--------------------------------------------------------------------------------
Contains a summary of the features included or not in isotopic clusters, with the following columns:

:FeatureID: Identifier of the feature, as it was provided in the :ref:`Measurements file (untargeted)`.
:RT: Retention time of the feature, as it was provided in the :ref:`Measurements file (untargeted)`.
:m/z: Mass-to-charge ratio of the feature, as it was provided in the :ref:`Measurements file (untargeted)`.
:sample: Name of the sample, as it was provided in the :ref:`Measurements file (untargeted)`.
:intensity: The intensity of the feature in the sample, as it was provided in the `Measurements file (untargeted)`:
:inClusters: The cluster identifiers in which the feature is included. If the feature is not included in any cluster, the value is "None".
:Isotopologues: The index of the isotopologues of the metabolite corresponding to the cluster. If the feature is not included in any cluster, the value is "None".

.. warning:: **A single feature may be included in multiple clusters.**
             In such cases, the 'inClusters' and 'isotopologues' columns will contain multiple values separated by commas.
             It is up to the user to handle these cases carefully when interpreting the data.



Cluster file (``.clusters.tsv``)
--------------------------------------------------------------------------------

Contains the isotopic clusters formed without prior knowledge based on the tracer element, with the following columns:

:ClusterID: Identifier of the isotopic cluster (generated by IsoGroup).
:FeatureID: Identifier of the feature, as it was provided in the :ref:`Measurements file (untargeted)`.
:RT: Retention time of the feature, as it was provided in the :ref:`Measurements file (untargeted)`.
:m/z: Mass-to-charge ratio of the feature, as it was provided in the :ref:`Measurements file (untargeted)`.
:sample: Name of the sample, as it was provided in the :ref:`Measurements file (untargeted)`.
:Intensity: The intensity of the feature in the sample, as it was provided in the :ref:`Measurements file (untargeted)`.
:Isotopologue: The index of the isotopologue of the metabolite corresponding to the cluster.
:AlsoIn: Identifiers of other clusters in which the feature is included. If the feature is not included in any other cluster, the list is empty.

Log file (``.log``)
--------------------------------------------------------------------------------

Extensive information on the annotation process can be found in the log file if ‘Verbose logs’ option has been checked.

Warning and error messages
--------------------------------------------------------------------------------

Error messages are explicit. You should examine carefully any warning/error message.
After correcting the problem, you might have to restart IsoGroup and re-run the analysis.