Welcome to The Documentation for CAD_to_OpenMC
=============================================
**CAD_to_OpenMC** is a python tool intened to ease the process of
creating OpenMC-ready surface meshes from a CAD-generated representation in the form
of a .step (or .stp) -file.

It borrows quite heavily from  from [Paramak](https://github.com/fusion-energy/paramak) and [cad_to_h5m](https://github.com/shimwell/step_to_h5m). Furthermore it relies on teh open source packages gmsh and cadquery2 for its functionality.

The mesh is created either by means of gmsh or as a simple .stl-export from cadquery. In the latter case the generated mesh may be refined using the mesh-refinement-tool [Mmg](https://www.mmgtools.org).

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   api
