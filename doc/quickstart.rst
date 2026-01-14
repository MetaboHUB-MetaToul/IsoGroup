..  _Quick start:

Quick start
********************************************************************************


Installation
------------------------------------------------

IsoGroup requires Python 3.7 or higher. If you do not have a Python environment
configured on your computer, we recommend that you follow the instructions
from `Anaconda <https://www.anaconda.com/download/>`_.

Then, open a terminal (e.g. run *Anaconda Prompt* if you have installed Anaconda) and type:

.. code-block:: bash

  pip install isogroup

You are now ready to start IsoGroup.

If this method does not work, you should ask your local system administrator or
the IT department "how to install a Python 3 package from PyPi" on your computer.

Alternatives & update
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you know that you do not have permission to install software systemwide,
you can install IsoGroup into your user directory using the :samp:`--user` flag:

.. code-block:: bash

  pip install --user isogroup


If you already have a previous version of IsoGroup installed, you can upgrade it to the latest version with:

.. code-block:: bash

  pip install --upgrade isogroup


Alternatively, you can also download all sources in a tarball from `GitHub <https://github.com/MetaboHUB-MetaToul/IsoGroup/tree/main/isogroup/base>`_,
but it will be more difficult to update IsoGroup later on.


Usage
------------------------------------------------

Command Line Interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IsoGroup can be executed from the command line interface (CLI) to perform either **untargeted isotopic clustering** or **targeted annotation**.

You can now use **two separate command-line executables** depending on your needs:

- **Untargeted isotopic clustering**:

.. code-block:: bash

  isogroup_untargeted [command line options]

Here after the available options with their full names are enumerated and detailed.

.. argparse::
   :module: isogroup.ui.cli
   :func: build_parser_untargeted
   :prog: isogroup_untargeted
   :nodescription:

- **Targeted annotation using a database**:

.. code-block:: bash

  isogroup_targeted [command line options ]

Here after the available options with their full names are enumerated and detailed.

.. argparse::
   :module: isogroup.ui.cli
   :func: build_parser_targeted
   :prog: isogroup_targeted
   :nodescription:


IsoGroup automatically carries out either untargeted isotopic clustering or targeted annotation of mass features

.. warning:: The annotation and clustering options must be carefully selected to ensure reliable interpretations of labeling data, as detailed in the Tutorials.

.. seealso:: Tutorial :ref:`First time using IsoGroup` has example data that you can use to test your installation.


Library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IsoGroup is also available as a library (a Python module) that you can import directly in your Python
scripts:

.. code-block:: python

  import isogroup

.. .. seealso::  Have a look at our :ref:`library showcase <Library documentation>` if you are interested into this experimental feature.
