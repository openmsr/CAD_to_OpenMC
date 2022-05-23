# step_to_h5m
This is intended to be a python package heavily inspired by [Paramak](https://github.com/fusion-energy/paramak).
It's raison d'etre is to enable a link between the online CAD-tool onShape and the neutron transport code OpenMC.

We will use cadQuery and it's links to OCCT to enable import and imprinting/merging algorithms.

The code structure we intend is a main class Assembly which contains the geometry in terms of a set of volumes.
We then use gmsh to generate a surface mesh which may be used by OpenMC through DAGMC/MOAB (i.e. .h5m-files)

# To install/set up in a virtual python environment
_replace \<name\> with an arbitrary name for your virtual environment_
1. Clone the github repository as you'd normally do
2. In the directory where you want your environment to reside do: ```python -m venv <name>```
3. Activate the environment: ```source <name>/bin/activate```
4. install the requirements: ```pip install -r requirements.txt```
5. Build and install moab (if not already installed). The moab team relies on conda for standard installation but are working on a pip-based solution. Once that is done moab would simply be added to the requirements-file instead.
  1. Clone the moab code-repository: e.g. ```git clone git@bitbucket.org:fathomteam/moab.git```
  2. Configure and build the code: ```mkdir build; cd build; cmake .. -DENABLE_PYMOAB=1; make; make install```

# Run a test case:
The follwing code-snippet can now be run to explore the capabilities of Assembly/step_to_h5m
```
  import assembly as ab
  a=ab.Assembly()
  a.stp_files=["file.step"]
  a.import_stp_files()
  a.export_brep('file.brep')
  a.brep_to_h5m(brep_filename='file.brep',min_mesh_size=0.1, max_mesh_size=10, samples=20)
```

N.b. the last parameters to brep_to_h5m are simply echoing their default values.
This procedure will in turn call OCP and gmsh to prduce a mesh with merged surfaces in the output file "dagmc.h5m"

The ```export_brep```-step may be omitted, in which case a temporary file will be written.
