import h5py
import numpy as np
#The layout of the HDF5 file is roughly:
#
#  tstt/
#   |
#   +----elemtypes
#   |
#   +----nodes/
#   |       |
#   |       +---coordinates
#   |       |
#   |       +---tags/
#   |
#   +----elements/
#   |       |
#   |       +---"group"/
#   |       |      |
#   |       |      +---connectivity
#   |       |      |
#   |       |      +---tags/
#   |       |
#   |       +---"group"/
#   |       |      |
#   |       |      +---connectivity
#   |       |      |
#   |       |      +---tags/
#  ...     ...
#   |
#   +----sets/
#   |       |
#   |       +----list
#   |       |
#   |       +----contents
#   |       |
#   |       +----parents
#   |       |
#   |       +----children
#   |
#   +----tags/
#           |
#           +---"name"/
#           |      |
#           |      +----type
#           |      |
#           |      +----id_list
#           |      |
#           |      +----values
#           |
#           +---"name"/
#           |      |
#           |      +----type
#           |      |
#           |      +----id_list
#           |      |
#           |      +----values
#          ...
#
#
# The groups contain:
# nodes: all vertices as an n,d dataset of floats. n is the number of vertices and d is the number of coordinates per vertex == 3
# elements: application defined number of subgroups. Each subgroup defines one or more mesh elements that have same topology and length of connectivity
#           subgroup names are combos of type and connectivity: Edge2 Tri3 etc.

elemtype={"Edge": 0, "Tri": 1, "Quad": 2, "Polygon": 3, "Tet": 4, "Pyramid":5, "Prism":6, "Knife":7, "Hex":9, "Polyhedron":10}


#This is an interface to allow to output stl-based meshes in the h5m format.
class File:
  def __init__(self, filename="dagmc.h5m"):
    #open an interface to the file
    self.h5pyf=h5py.File(filename,"w")

    #set up the basic h5m structure
    #The root grou is always tstt
    self.root_group=self.h5pyf.create_group('tstt')
    for g in ['elements', 'nodes', 'sets', 'tags']:
      self.root_group.create_group(g)
    self.root_group['elemtypes']=h5py.enum_dtype(elemtype, basetype='u1')
    dt=h5py.string_dtype(encoding='ascii', length=None)
    self.root_group['history']=self.h5pyf.create_dataset('history', (4,), dtype=dt)
    self.root_group['history'][0]='h5m'
    self.root_group['history'][1]='1.0'
    dd=datetime.datetime.now()
    self.root_group['history'][1]=dd.strftime("%x")
    self.root_group['history'][1]=dd.strftime("%X.%f")

  def __del__(self):
    pass

  def add_vertices(vertices,n=None,d=None):
    if not n:
      n=vertices.shape[0]
    if not d:
      n=vertices.shape[1]
    ds=self.root_group['nodes'].create_dataset('coordinates', (n,d), dtype='f8')
    #assume we start at id 1
    ds.attr['start_id']=np.uint64(1)
    ds=vertices[:n,:d]

  def add_elements(connectivity,elemname):
    if elemname not in elemtype.keys():
      print(f'{elemname} is not a h5m-supported element type.')
    ds=self.root_group['elements'].create_dataset(f'{elemname}{connectivity.shape[1]}', elements.shape, dtype='f8')
    #assume we start at id 1
    ds.attr
    ds=vertices[:n,:d]
