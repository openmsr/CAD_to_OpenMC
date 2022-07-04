import gmsh
import cadquery2 as cq
import os
import tempfile
import assemblymesher

class MesherGMSH:
  def __init__(self, min_mesh_size, max_mesh_size, curve_samples, default, mesh_algorithm, vetoed, threads, solids):
    self.IntermediateLayer='brep'
    self._min_mesh_size=min_mesh_size
    self._max_mesh_size=max_mesh_size
    self.curve_samples=curve_samples
    self.mesh_algorithm=mesh_algorithm
    self.default=default
    self.solids=solids
    self.vetoed=vetoed
    self.threads=threads
    self._gmsh_init()
    self._cq_solids_to_gmsh()

  def __del__(self):
    gmsh.finalize()

  def _gmsh_init(self):
      gmsh.initialize()
      if (self.verbose>1):
          gmsh.option.setNumber("General.Terminal",1)
      else:
          gmsh.option.setNumber("General.Terminal",0)

      if(not self.default):
        gmsh.option.setString("Geometry.OCCTargetUnit","CM")
        #do this by means of properties instead
        if(self.threads is not None):
          gmsh.option.setNumber("General.NumThreads",self.threads)

        gmsh.option.setNumber("Mesh.Algorithm", self.mesh_algorithm)
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.min_mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.max_mesh_size)

        gmsh.option.setNumber("Mesh.MaxRetries",3)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints",0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", self.curve_samples)

  def _cq_solids_to_gmsh(self):
      #should check of solids is in fact a compound?
      compound=cq.Compound.makeCompound(self.solids)

      with tempfile.TemporaryDirectory() as td:
        outpath=os.path.join(td,'export.',self.IntermediateLayer)
        with open(outpath,'w') as fp:
          compound.exportBrep(outpath)
        gmsh.model.occ.importShapes(outpath)
      gmsh.model.occ.synchronize()
      gmsh.model.add(f"model from assembly.py: {outpath}")

  def _generate_mesh(self):
      if(self.verbose>0):
          print("INFO: GMSH generate surface mesh")
      gmsh.model.mesh.generate(2)

  def generate_stls(self):
      """export all the optionally merged volumes as stl-files
      and returns the list of files. This list may subsequently be iterated upon
      to merge into a h5m-file. Expects that the geometry has been surface-mesh by gmsh
      so we have a list of volumes to operate on.
      We do this be greating gmsh physical groups and export 1 group at a time."""
      self._generate_mesh()
      stls=[]
      for dim,vid in self.volumes:
         if (dim!=3):
             #appears not to be a volume - skip
             continue
         ents = gmsh.model.getAdjacencies(dim,vid)
         #may want to skip some problematic surfaces
         if(self.vetoed):
            picked_ents=[e for e in ents[1] if e not in self.vetoed]
         else:
            picked_ents=ents[1]
         pg = gmsh.model.addPhysicalGroup(2,picked_ents)
         ps = gmsh.model.setPhysicalName(2,pg,f'surfaces_on_volume_{vid}')
         filename=f'volume_{vid}.stl'
         try:
            gmsh.write(filename)
            stls.append((vid,filename))
         except:
            print(f'WARNING: Could not write volume {vid}. Skipping')
         gmsh.model.removePhysicalGroups([]) # remove group again
      return stls


class MesherGMSHBuilder:
  def __init__(self):
    self._instance = None

  def __call__(self, min_mesh_size, max_mesh_size, curve_samples, default, fname, **_ignored):
    if not self._instance:
      self._instance = MesherGMSH(min_mesh_size,max_mesh_size, curve_samples,default)
    return self._instance
