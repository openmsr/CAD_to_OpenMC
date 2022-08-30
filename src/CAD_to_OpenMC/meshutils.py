"""
  Set of utility functions that handle IO with the .mesh format
"""
import struct
import numpy
import pathlib as pl
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

def _write_edge(fp,e,label=None):
  if not label:
    fp.write(f'{e[0]} {e[1]} 1\n')
  else:
    fp.write(f'{e[0]} {e[1]} {label}\n')


def write_dotmesh(dest, vertices, triangles,vertex_labels=None, triangle_labels=None, edges=None, required_edges='all'):
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
      try:
        _write_vertex(fp,vertices[i],label=vertex_labels[i])
      except:
        _write_vertex(fp,vertices[i],label=None)
    fp.write(' Triangles\n')
    fp.write(f' {nt}\n')
    for i in range(nt):
      try:
        _write_triangle(fp,triangles[i]+1,label=triangle_labels[i])
      except:
        _write_triangle(fp,triangles[i]+1,label=None)
    if (edges is not None):
      ne=edges.shape[0]
      fp.write(f'Edges\n{ne}\n')
      for i in range(ne):
        try:
          _write_edge(fp,edges[i]+1,label=edge_labels[i])
        except:
          _write_edge(fp,edges[i]+1,label=None)
      if (required_edges=='all'):
        fp.write(f'RequiredEdges\n{ne}\n')
        for i in range(ne):
          fp.write(f'{i+1}\n')
    fp.write(' End\n')

def find_edges(triangles):
  """algorithm to find the edge lines in a set of triangles. I.e. the lines
  which are only present in a single triangle"""
  #loop over all triangles i = 0,...,N_t
  #loop over the set of 3 possible edges e\in{ [v_0,v_1],[v_1,v_2],[v_3,v_4] } in T_i, where j\in{0,1,2}
  #given an edge e_j, run through all triangles T_k, where i<k<=N_t
  #if that edge is found move to next
  #else append to list of edges
  edges=[]
  for i in range(triangles.shape[0]):
    Ti=triangles[i]
    ee_i=(sorted([Ti[0],Ti[1]]),sorted([Ti[1],Ti[2]]),sorted([Ti[2],Ti[0]]))
    for j in range(3):
      e_j=ee_i[j]
      for k in range(triangles.shape[0]):
        if (k==i):
          continue
        Tk=triangles[k]
        ee_k=(sorted([Tk[0],Tk[1]]),sorted([Tk[1],Tk[2]]),sorted([Tk[2],Tk[0]]))
        if(e_j in ee_k):
          break
      else:
        edges.append(e_j)
  return edges
