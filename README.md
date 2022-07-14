[![Python application](https://github.com/openmsr/step_to_h5m/actions/workflows/python-app.yml/badge.svg?branch=factory_meshclasses)](https://github.com/openmsr/step_to_h5m/actions/workflows/python-app.yml)

# CAD_to_OpenMC

This is intended to be a python package heavily inspired by [Paramak](https://github.com/fusion-energy/paramak), and borrows a lot from [step_to_h5m]( https://github.com/fusion-energy/step_to_h5m).
It's "raison d'etre" is to establish an open source-based link between CAD tools in general and the neutron transport code OpenMC.

Although most CAD-tools use some other internal and/or native representation for geometry, most, if not all, will be able to export to the STEP-file format. Therefore this is the format we use 

We will use cadQuery and its links to OCCT to enable import and imprinting/merging algorithms.

The code structure relies on a main class Assembly which maintains the geometry in terms of a list of instances of the subclass Entity.
A geometry is imported from a (set of) .step files into the entity-list. This list is passed on to a mesher-class which generates a meshed geometry.

# To install/set up in a virtual python environment
_replace \<name\> with an arbitrary name for your virtual environment_
1. In the directory where you want your environment to reside do: ```python -m venv <name>```
2. Activate the environment: ```source <name>/bin/activate```
3. Build and install moab (if not already installed). The moab team relies on conda for standard installation but are working on a pip-based solution. Once that is done moab would simply be added to the requirements-file instead.
    1. Clone the moab code-repository: e.g. ```git clone git@bitbucket.org:fathomteam/moab.git```
    2. Configure and build the code:
    ```bash
      mkdir build;
      cd build; cmake .. -DENABLE_PYMOAB=1 -DCMAKE_INSTALL_PREFIX=<path/to/venv/>;
      make;
      make install;
    ```    
    3. Additionally you will need to build the python interface layer.
    ```bash
      cd pymoab
      sudo python setup.py install
    ```
4. Install the main package: ```pip install CAD_to_OpenMC```. This will pip-install all the required python packages in the virtual environment. This ensures that no additional conflicts are introduced with the system python.

Should you wish to install the development version of this package you may do so by cloning this repository and replace the last command by: ```pip install <path/to/repo>```. This procedure will build and install the python package locally directly from source.

# Run a test case:
The follwing code-snippet can now be run to explore the capabilities of Assembly/step_to_h5m. We supply a couple of example .step-files in the examples directory. In this example we'll be using a geometry with a set of 7 pin-cells.

```python
  import CAD_to_OpenMC.assembly as ab
  a=ab.Assembly()
  a.stp_files=["examples/7pin.step"]
  a.import_stp_files()
  a.solids_to_h5m()
```

Unless anything else is specified this procedure simply uses the default CAD_to_OpenMC parameters to create meshes using the default choice of meshing backend - namely gmsh.
E.g. for the "gmsh"-backend this means sampling curves 20 times, a minimum mesh-element size of 0.1, and a maximum mesh element size of 10.
This procedure will in turn call OCP and gmsh to produce a mesh with merged surfaces in the output file "dagmc.h5m"

The other available meshing backend is the stl-export from CadQuery (accessible by setting ```backend='stl'```) which uses the parameters ```tolerance``` and ```angular_tolerance``` to set meshing quality.

Parameters are changed by means of altering entries in the ```mesher_config```-dictionary defined  in the assemlby modulde root namespace. Like:
```python
 ab.mesher_config['min_mesh_size']=0.2
 ```

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
