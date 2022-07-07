import gmsh
import cadquery2 as cq
import os
import tempfile

class MesherGMSH:
  def __init__(self, min_mesh_size, max_mesh_size, curve_samples, default, mesh_algorithm, vetoed, threads, entities):
    self.IntermediateLayer='brep'
    self.min_mesh_size=min_mesh_size
    self.max_mesh_size=max_mesh_size
    self.curve_samples=curve_samples
    self.mesh_algorithm=mesh_algorithm
    self.default=default
    self.entities=entities
    self.vetoed=vetoed
    self.threads=threads
    self.verbose=2
    self._gmsh_init()
    self._cq_solids_to_gmsh()

  #def __del__(self):
  #  gmsh.finalize()

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
      import glob
      #create a compund cq solid from entities
      solids=[e.solid for e in self.entities]
      compound=cq.Compound.makeCompound(solids)
      with tempfile.TemporaryDirectory() as td:
        outpath=os.path.join(td,'export.'+self.IntermediateLayer)
        compound.exportBrep(outpath)
        vols=gmsh.model.occ.importShapes(outpath)
      gmsh.model.occ.synchronize()
      self._reorder()

  def _reorder(self):
      #in the process of exporting/importing, the ordering may have been jumbled
      brep_volume_dimtags=gmsh.model.occ.getEntities(3)
      tag_idx=[]
      for dimtag in brep_volume_dimtags:
        cms=gmsh.model.occ.get_center_of_mass(dimtag[0],dimtag[1])
        bb0=gmsh.model.occ.get_bounding_box(dimtag[0],dimtag[1])
        bb=[abs(bb0[0]-bb0[3]),abs(bb0[1]-bb0[4]),abs(bb0[2]-bb0[5])]
        vol=gmsh.model.occ.get_mass(dimtag[0],dimtag[1])
        j=self._find_similar(cms,bb,vol)
        tag_idx.append(j)
      #tag_idx is now the reordered list of indices
      #use that to reorder the entities to match the newly imported
      #list of volumes in the gmsh geometry representation
      print("tag_idx",tag_idx)
      self.entities=[self.entities[i] for i in tag_idx if i!=-1]

  def _find_similar(self,cms,bb,vol):
      closest=1e12
      iclose=-1
      for i,e in enumerate(self.entities):
          simi=e.similarity(cms,bb,vol,tolerance=1e5)
          print(cms,bb,vol)
          print(i,e.center,e.bb,e.volume,simi)
          if (simi<closest):
              iclose=i
              closest=simi
      return iclose

  def _generate_mesh(self):
      if(self.verbose>0):
          print("INFO: GMSH generate surface mesh")
      print(gmsh.model.getEntities(2))
      gmsh.model.mesh.generate(2)

  def generate_stls(self):
      """export all the optionally merged volumes as stl-files
      and returns the list of files. This list may subsequently be iterated upon
      to merge into a h5m-file. Expects that the geometry has been surface-mesh by gmsh
      so we have a list of volumes to operate on.
      We do this be greating gmsh physical groups and export 1 group at a time."""
      self._generate_mesh()
      stls=[]
      vols=gmsh.model.getEntities(3)
      for i,e in enumerate(self.entities):
        #gmsh volume ids run from 1
        dim,vid=vols[i]

        if (dim!=3):
            #appears not to be a volume - skip
            continue
        ents = gmsh.model.getAdjacencies(dim,vid)
        #may want to skip some problematic surfaces
        if(self.vetoed):
            picked_ents=[f for f in ents[1] if f not in self.vetoed]
        else:
            picked_ents=ents[1]
        pg = gmsh.model.addPhysicalGroup(2,picked_ents)
        ps = gmsh.model.setPhysicalName(2,pg,f'surfaces_on_volume_{vid}')
        filename=f'volume_{vid}.stl'
        try:
            gmsh.write(filename)
            e.stl=filename
        except:
            e.stl=None
            print(f'WARNING: Could not write volume {vid}. Skipping')
        gmsh.model.removePhysicalGroups([]) # remove group again

class MesherGMSHBuilder:
  def __init__(self):
    self._instance = None

  def __call__(self, min_mesh_size, max_mesh_size, curve_samples, default, mesh_algorithm, vetoed, threads, entities,**_ignored):
    if not self._instance:
      self._instance = MesherGMSH(min_mesh_size, max_mesh_size, curve_samples, default, mesh_algorithm, vetoed, threads, entities)
    return self._instance
