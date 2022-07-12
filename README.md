[![Python application](https://github.com/openmsr/step_to_h5m/actions/workflows/python-app.yml/badge.svg?branch=factory_meshclasses)](https://github.com/openmsr/step_to_h5m/actions/workflows/python-app.yml)

# CAD_to_OpenMC

This is intended to be a python package heavily inspired by [Paramak](https://github.com/fusion-energy/paramak), and borrows a lot from [step_to_h5m]( https://github.com/fusion-energy/step_to_h5m).
It's "raison d'etre" is to enable a link between the online CAD-tool onShape and the neutron transport code OpenMC.

We will use cadQuery and its links to OCCT to enable import and imprinting/merging algorithms.

The code structure we intend is a main class Assembly which contains the geometry in terms of a set of volumes.
We then use gmsh to generate a surface mesh which may be used by OpenMC through DAGMC/MOAB (i.e. .h5m-files)

# To install/set up in a virtual python environment
_replace \<name\> with an arbitrary name for your virtual environment_
1. In the directory where you want your environment to reside do: ```python -m venv <name>```
2. Activate the environment: ```source <name>/bin/activate```
3. Build and install moab (if not already installed). The moab team relies on conda for standard installation but are working on a pip-based solution. Once that is done moab would simply be added to the requirements-file instead.
  1. Clone the moab code-repository: e.g. ```git clone git@bitbucket.org:fathomteam/moab.git```
  2. Configure and build the code: ```mkdir build; cd build; cmake .. -DENABLE_PYMOAB=1 -DCMAKE_INSTALL_PREFIX=<path/to/venv/>; make; make install```
4. Install the package: ```pip install CAD_to_OpenMC```

# Run a test case:
The follwing code-snippet can now be run to explore the capabilities of Assembly/step_to_h5m
```
  import CAD_to_OpenMC.assembly as ab
  a=ab.Assembly()
  a.stp_files=["file.step"]
  a.import_stp_files()
  a.solids_to_h5m(backend='gmsh')
```
Unless anything else is specified this proceudre simply uses the default CAD_to_OpenMC parameters to create meshes.
E.g. for the "gmsh"-backend this means sampling curves 20 times, a minimum mesh-element size of 0.1, and a maximum mesh element size of 10.
This procedure will in turn call OCP and gmsh to produce a mesh with merged surfaces in the output file "dagmc.h5m"

The other available meshing backend is the stl-export from CadQuery (accessible by setting ```backend='stl'```) which uses the parameters ```stl_tol``` and ```stl_ang_tol``` to set meshing quality.

Parameters are changed by means of altering entries in the ```mesher_config```-dictionary. Like:
<code>
 ab.mesher_config['min_mesh_size']=0.2
</code>

Below is a list of the available parameters and their
meanings:

<dl>
<dt>min_mesh_size</dt>
<dd>Minimum mesh element size (gmsh backend)</dd>
<dt>max_mesh_size</dt>
<dd>Maximum mesh element size (gmsh backend)</dd>
<dt>curve_samples</dt>
<dd>The number of point samples used to approximate a curve, 1D-mesh. (gmsh backend)</dd>
<dt>mesh_algorithm</dt>
<dd>Integer specifying which mesh algorithm to use (gmsh backend) 1: Adaptive algorithm - generally the most robust choice.</dd>
<dt>vetoed</dt>
<dd>A python list of surfaces that should be excluded from meshin. Useful when some surface fails for whatever reason</dd>
<dt>threads</dt>
<dd>The number of threads to be used to speed up the meshing algorithms. Useful for multicore-computers.</dd>
<dt>tolerance</dt>
<ddRelative mesh tolerance (stl backend). Lower this to get a finer mesh.</dd>
<dt>angular_tolerance</dt>
<dd>Relative angular mesh tolerance (stl backend) Lower this to get a better angular resolution.</dd>
</dl>
