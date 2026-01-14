# IsoGroup - **Iso**topic **Group**ing for mass spectrometry labeling experiments

[![PyPI version](https://badge.fury.io/py/IsoGroup.svg)](https://badge.fury.io/py/IsoGroup)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/isogroup.svg)](https://pypi.python.org/pypi/isogroup/)
[![Documentation Status](https://readthedocs.org/projects/isocor/badge/?version=latest)](http://isocor.readthedocs.io/?badge=latest) -->


## What is IsoGroup?
**IsoGroup is a scientific software dedicated to the processing of isotopic data from untargeted mass spectrometry (MS) labeling experiments**.
It is suitable for diverse applications including annotation, isotopic profiling, quantification, and fluxomics experiments.
IsoGroup supports two complementary modes of data processing: 
* In **targeted** mode, IsoGroup performs **annotation** of isotopic clusters based on a user-provided database of compounds.
  Mass features (m/z, retention time, intensity) are matched to known metabolites and isotopologues according to their exact mass and retention time.
  A database of isotopologues specific to the selected tracer is generated and used for the annotation.
* In **untargeted** mode, IsoGroup performs **clustering** by grouping MS features identified as :ref:`isotopologues <isotopologues>` into :ref:`isotopic clusters <isotopic clusters>` **without prior knowledge**, based solely on the characteristic mass differences induced by the tracer element.


The output of IsoGroup is a list of isotopic clusters with their m/z, retention time, intensity, and isotopologues. 

The code is open-source, and available under a GPLv3 license.

Detailed documentation can be found online at Read the Docs ([https://isogroup.readthedocs.io/](https://isogroup.readthedocs.io/)).
Check out the [Tutorials](https://isogroup.readthedocs.io/en/latest/tutorials.html) !

## Key features
* group MS features into **isotopic clusters**
* generate a database with isotopologues for a given tracer element,
* **annotate** isotopic clusters,
* calculate **exact mass and retention time errors** for annotated clusters,
* can be used with any tracer element,
* open-source, free and easy to install everywhere where Python 3 and pip run,

## Quick-start
IsoGroup requires Python 3.7 or higher and run on all platforms.
Please check [the documentation](https://isogroup.readthedocs.io/en/latest/quickstart.html) for complete
installation and usage instructions.

Use `pip` to **install IsoGroup from PyPi**:

```bash
$ pip install isogroup
```

<!-- Then, run IsoGroup in command line with:

```bash
$ isogroup -->

<!-- ``` -->

IsoGroup is also available directly from command-line and as a Python library.

## Bug and feature requests
If you have an idea on how we could improve IsoGroup please submit a new *issue*
to [our GitHub issue tracker](https://github.com/MetaboHUB-MetaToul/IsoGroup/issues).


## Developers guide
### Contributions
Contributions are very welcome! :heart:

Please work on your own fork,
follow [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide,
and make sure you pass all the tests before a pull request.

### Local install with pip
In development mode, do a `pip install -e /path/to/IsoGroup` to install
locally the development version.


### Build the documentation locally
Build the HTML documentation with:

```bash
$ cd doc
$ make html
```

The PDF documentation can be built locally by replacing `html` by `latexpdf`
in the command above. You will need a recent latex installation.


## Authors
Noémie Butin, Loïc Le Grégam, Pierre Millard, Rochelle Kouakou

## Contact
:email: noemie.butin@univ-perp.fr, loic.legregam@unige.ch, millard@insa-toulouse.fr, kouakou@insa-toulouse.fr
