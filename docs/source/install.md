# To install/set up in a virtual python environment
_replace \<name\> with an arbitrary name for your virtual environment_
1. In the directory where you want your environment to reside do: ```python -m venv <name>```
2. Activate the environment: ```source <name>/bin/activate```
3. Build and install moab (if not already installed). The moab team relies on conda for standard installation but are working on a pip-based solution. Once that is done moab would simply be added to the requirements-file instead.
    1. Clone the moab code-repository: e.g. ```git clone git@bitbucket.org:fathomteam/moab.git```
    2. Configure and build the code:
    ```bash
      mkdir build;
      cd build; cmake .. -DENABLE_PYMOAB=1 -DENABLE_HDF5=1 -DCMAKE_INSTALL_PREFIX=<name>;
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

5. Optionally install the msh refinement tool mmg (https://www.mmgtools.org), which may be run in conjunction with the cq/stl-mesher backend to avoid the high aspect ratio triangles that this backend tends to produce. Arch-linux users may install this from the AUR, otherwise get the source (and build from that) or binary builds from the mmg-site.

# To install in a conda environment
_replace \<name\> with an arbitrary name for your virtual environment_

If instead you prefer to use a [conda-environment](https://docs.conda.io/projects/conda/en/stable/), this is now as simple as:
1. create an environment, e.g. ```conda create -n <name>```
2. activate it: ```conda activate <name>```
3. install CAD_to_OpenMC: ```conda install cad-to-openmc```

You may of course replace conda with [mamba/micromamba](https://mamba.readthedocs.io/en/latest/), should you prefer to do so.
