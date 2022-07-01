import gmsh


class MesherGMSH:
  def __init__(self, min_mesh_size, max_mesh_size, curve_samples, default, mesh_algorithm, solids):
    self._min_mesh_size=min_mesh_size
    self._max_mesh_size=max_mesh_size
    self.curveSamples=curveSamples
    self.mesh_algorithm=mesh_algorithm
    self.solids=solids
    self._gmsh_init()

  def _gmsh_init(self):
      self._export_brep()

      gmsh.initialize()
      if (self.verbose>1):
          gmsh.option.setNumber("General.Terminal",1)
      else:
          gmsh.option.setNumber("General.Terminal",0)

      gmsh.model.add(f"model from Assembly.py {brep_fn}")
      if(not default):
        gmsh.option.setString("Geometry.OCCTargetUnit","CM")
        #do this by means of properties instead
        if(threads is not None):
          gmsh.option.setNumber("General.NumThreads",threads)

        gmsh.option.setNumber("Mesh.Algorithm", self.mesh_algorithm)
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.min_mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.max_mesh_size)

        gmsh.option.setNumber("Mesh.MaxRetries",3)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints",0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", self.curve_samples)
      self.volumes = gmsh.model.occ.importShapes(self.fname)
      gmsh.model.occ.synchronize()

  def _gmsh_deinit(self):
      gmsh.finalize()

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
         pg = gmsh.model.addPhysicalGroup(2,ents[1])
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
