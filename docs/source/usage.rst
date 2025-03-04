Usage
=====

.. _installation:

Installation
------------

To use CAD_to_OpenMC, first install it using either pip or a conda-environment:

If using pip, you will first need to install the python bindings to MOAB (short
for Mesh Oriented datABase), from a special location, since the MOAB-rpoject
historically uses conda-environments. There is an official PYPI-packge in the development
tree however so this sitaution should soon be remedied.

.. code-block:: console

(.venv) $ pip install https://github.com/shimwell/wheels/raw/refs/heads/main/moab/moab-wheels-ubuntu-latest/moab-5.5.1-cp<pyv>-cp<pyv>-manylinux_2_28_x86_64.whl
(.venv) $ pip install CAD_to_OpenMC

where you will have to replace <pyv> by your abbreviated python version. Eg. "312" for python 3.12. 3.9..3.12 are acceptable.

If you prefer a conda (or mamba) environments, you may simply use

.. code-block:: console
$ conda install cad-to-openmc


Creating surface meshes
----------------

To create your first mesh you may use the supplied example data file: ``7pin.step``
which may be found in the examples data directory.
to do so (in a python console)
.. code-block:: console
$ import CAD_to_OpenMC.assembly as ab
$ a=ab.Assembly(['examples/7pin.step'])
$ a.run(h5m_filename = '7pin.h5m')


