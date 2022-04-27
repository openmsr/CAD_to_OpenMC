# step_to_h5m
This is intended to be a python package heavily inspired by [Paramak](https://github.com/fusion-energy/paramak).
It's raison d'etre is to enable a link between the online CAD-tool onShape and the neutron transport code OpenMC.

We will use cadQuery and it's links to OCCT to enable import and imprinting/merging algorithms.

The code structure we intend is a main class Assembly which contains the geometry in terms of a set of volumes.
We then use gmsh to generate a surface mesh which may be used by OpenMC through DAGMC/MOAB (i.e. .h5m-files)
