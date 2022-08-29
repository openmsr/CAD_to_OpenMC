"""
  Set of utility functions that handle IO with the .mesh format
"""
import struct
import numpy
from .stl_utils import *

def write_dummy_dotsol(dest,nvertices:int=0):
  header="""MeshVersionFormatted
1

Dimension 3
SolAtVertices
"""
  with open(dest,"w") as fp:
    fp.write(header)
    fp.write(f'{nvertices}\n')
    fp.write('1 1\n')
    fp.write('0\n'*nvertices)
    fp.write('End\n')


def _print_vertex(v,label=None):
  if not label:
    print(f'{v[0]:24} {v[1]:24} {v[2]:24}  1')
  else:
    print(f'{v[0]:24} {v[1]:24} {v[2]:24}  {label}')

def _print_triangle(tri,label=None):
  if not label:
    print(f' {tri[0]} {tri[1]} {tri[2]} 1')
  else:
    print(f' {tri[0]} {tri[1]} {tri[2]} {label}')

def _write_triangle(fp,tri,label=None):
  if not label:
    fp.write(f' {tri[0]} {tri[1]} {tri[2]} 1\n')
  else:
    fp.write(f' {tri[0]} {tri[1]} {tri[2]} {label}\n')

def _write_vertex(fp,v,label=None):
  if not label:
    fp.write(f'{v[0]:24} {v[1]:24} {v[2]:24}  1\n')
  else:
    fp.write(f'{v[0]:24} {v[1]:24} {v[2]:24}  {label}\n')

def write_dotmesh(dest, vertices, triangles,vertex_labels=None, triangle_labels=None):
  """ vertices and triangles can be lists of lists or 2d arrays """
  header=""" MeshVersionFormatted
 2
 Dimension
 3
"""
  nv=vertices.shape[0]
  nt=triangles.shape[0]

  with open(dest,"w") as fp:
    fp.write(header)
    fp.write(' Vertices\n')
    fp.write(f' {nv}\n')
    for i in range(vertices.shape[0]):
      _write_vertex(fp,vertices[i])
    fp.write(' Triangles\n')
    fp.write(f' {nt}\n')
    for i in range(nt):
      _write_triangle(fp,triangles[i]+1)
    fp.write(' End\n')
