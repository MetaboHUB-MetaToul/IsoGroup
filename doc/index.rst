IsoGroup: Isotopic Grouping for mass spectrometry labeling experiments
********************************************************************************

Welcome to IsoGroup documentation!
----------------------------------------

**IsoGroup is a scientific software tool dedicated to the processing of isotopic data from untargeted mass spectrometry (MS) labeling experiments.**
IsoGroup groups mass features (m/z, retention time, intensity) identified as :ref:`isotopologues <isotopologues>` into :ref:`isotopic clusters <isotopic clusters>` based on the tracer element used in the experiment.
It also generates a database of isotopologues specific to the selected tracer, which is used to annotate the clusters.
The output of IsoGroup is a list of annotated isotopic clusters (m/z, retention time, intensity, isotopologues). 

.. It supports diverse applications, including annotation, isotopic profiling, quantification or fluxomics. # Pour version non-ciblée
.. It automatically extracts and annotation of isotopic information from a list of MS features (m/z, retention time, intensity).

The code is open-source, and available on `GitHub <https://github.com/MetaboHUB-MetaToul/IsoGroup>`_ under a :ref:`GPLv3 license <license>`.

.. This documentation is available on Read the Docs (`https://isogroup.readthedocs.io <https://isogroup.readthedocs.io/>`_) # MAJ
.. and can be downloaded as a `PDF file <https://readthedocs.org/projects/isocor/downloads/pdf/latest/>`_.


.. rubric:: Key features

* group MS features into **isotopic clusters**
* generate a database with isotopologues for a given tracer element,
* **annotate** isotopic clusters,
* calculate **exact mass and retention time errors** for annotated clusters,
* can be used with any tracer element,
* open-source, free and easy to install everywhere where Python 3 and pip run,
.. * handle multiple types of samples (unlabeled, fully labeled, Pascal Triangle, etc)
.. * work in natural abundance & enriched settings
.. * integrate the possibility of using in-house databases for targeted analysis


.. seealso:: We strongly encourage you to read the :ref:`Tutorials` before using IsoGroup.


.. toctree::
   :maxdepth: 2
   :caption: Usage

   quickstart.rst
   tutorials.rst
   .. cite.rst

.. toctree::
   :maxdepth: 1
   :caption: Miscellaneous

   definitions.rst
   license.rst
   library_doc.rst
   .. faq.rst

.. todolist::
