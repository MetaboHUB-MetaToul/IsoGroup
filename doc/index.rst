IsoGroup: Isotopic Grouping for mass spectrometry labeling experiments
********************************************************************************

Welcome to IsoGroup documentation!
----------------------------------------

**IsoGroup is a scientific software tool dedicated to the processing of isotopic data from untargeted mass spectrometry (MS) labeling experiments.**
It is suitable for diverse applications including annotation, isotopic profiling, quantification, and fluxomics experiments.

IsoGroup supports two complementary modes of data processing: 

* In **targeted** mode, IsoGroup performs **annotation** of isotopic clusters based on a user-provided database of compounds.
  Mass features (m/z, retention time, intensity) are matched to known metabolites and isotopologues according to their exact mass and retention time.
  A database of isotopologues specific to the selected tracer is generated and used for the annotation.
* In **untargeted** mode, IsoGroup performs **clustering** by grouping MS features identified as :ref:`isotopologues <isotopologues>` into :ref:`isotopic clusters <isotopic clusters>` **without prior knowledge**, based solely on the characteristic mass differences induced by the tracer element.

The output of IsoGroup is a list of isotopic clusters with their m/z, retention time, intensity, and isotopologues. 

The code is open-source, and available on `GitHub <https://github.com/MetaboHUB-MetaToul/IsoGroup>`_ under a :ref:`GPLv3 license <license>`.

This documentation is available on Read the Docs (`https://isogroup.readthedocs.io <https://isogroup.readthedocs.io/>`_) and can be downloaded as a `PDF file <https://readthedocs.org/projects/isocor/downloads/pdf/latest/>`_.


.. rubric:: Key features

* Group MS features into **isotopic clusters** without a priori knowledge (untargeted mode),
* Use in-house databases for **targeted annotation** of isotopic clusters,
* Generate a database with isotopologues for a given tracer element,
* Calculate **exact mass and retention time errors** for annotated clusters,
* Compatible with any tracer element,
* Handle multiple sample types (unlabeled, fully labeled, Pascal Triangle, etc)
* Open-source, free and easy to install everywhere where Python 3 and pip run,


.. seealso:: We strongly encourage you to read the :ref:`Tutorials` before using IsoGroup.


.. toctree::
   :maxdepth: 2
   :caption: Usage

   quickstart.rst
   tutorials_targeted.rst
   tutorials_untargeted.rst

.. toctree::
   :maxdepth: 1
   :caption: Miscellaneous

   definitions.rst
   license.rst
   library_doc.rst

.. todolist::
