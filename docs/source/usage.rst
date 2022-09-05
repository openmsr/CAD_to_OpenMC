Usage
=====

.. _installation:

Installation
------------

To use CAD_to_OpenMC, first install it using pip:

.. code-block:: console

   (.venv) $ pip install CAD_to_OpenMC

Creating surface meshes
----------------

To create your first mesh you may use the supplied example data file: ``7pin.step``
which may be found in the examples data directory.
to do so (in a python console)

>>> import CAD_to_OpenMC.assebmly as ab
>>> a=ab.Assembly()
>>> a.stp_files=['examples/7pin.step']
>>> a.import_stp_files()
>>> a.stl2h5m()

