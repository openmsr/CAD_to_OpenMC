import gmsh
import sys
import re
import os
import logging

"""
Script that checks the entitities of a given file by importing them into gmsh and
meshes either volumes (--vols) or surfaces (--sfc) my meshing them 1 by 1, Reporting any
that do not mesh cleanly with the default gmsh-parameters. The meshing may run in parallel
by adding a 2nd argument which is taken as the number of threads to run.
"""
argv=sys.argv
stepfile=argv[1]
if '--vols' in argv:
  mode='vol'
  argv.remove('--vols')
if '--sfc' in argv:
  mode='sfc'
  argv.remove('--sfc')

sids=[]
for a in argv:
  m=re.match('--sid=(\d+)',a)
  if m:
    sids.append(int(m.group(1)))
if len(sids)>0:
  mode='sfc'

if not mode:
  print("Please specify a running mode. Either \'--vol\', \'--sfc\', or \'--sid=N\'.")
  exit(0)

default_opts=False
if '--default' in argv:
  default_opts=True
  argv.remove('--default')

size=1
try:
    size=int(argv[2])
except:
    pass

logging.basicConfig(filename=f'output_{os.getpid()}.log', level=logging.DEBUG)


gmsh.initialize()
gmsh.option.setNumber("General.Terminal",0)
gmsh.model.occ.importShapes(stepfile)
gmsh.model.occ.synchronize()
vols=gmsh.model.getEntities(3)
sfcs=gmsh.model.getEntities(2)
max_id=50000
max_vid=(vols[-1])[1]
max_sfcid=(sfcs[-1])[1]

print(f'INFO: Found {len(vols)} volumes in {stepfile}')
print(f'INFO: Found {len(sfcs)} surfaces in {stepfile}')
print(f'INFO: Output wil be in: output_{os.getpid()}.log')

gmsh.option.set_number("General.NumThreads",size)
if not default_opts:
  gmsh.option.set_number("Mesh.Algorithm",1)
  gmsh.option.set_number("Mesh.MeshSizeMin",0.1)
  gmsh.option.set_number("Mesh.MeshSizeMax",10)
  gmsh.option.set_number("Mesh.MeshSizeFromCurvature",0)
  gmsh.option.set_number("Mesh.MaxRetries",3)
  gmsh.option.setNumber("Mesh.MeshSizeFromPoints",0)
  gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
else:
  logging.info('Using default gmsh meshing options')

if mode=='vol':
  for vol,vid in vols:
    if int(vid) % 10 ==0:
      print(f'At volume {vid} of {len(vols)}')
    logging.info(f'Volume {vid}')
    exclude=[(v,i) for (v,i) in vols if i!=vid]
    gmsh.model.removeEntities(exclude,recursive=True)
    try:
      gmsh.model.mesh.generate(2)
    except:
      logging.error(f'Vol {vid} did not mesh, listing connected surfaces')
      (up,down)=gmsh.model.getAdjacencies(3,vid)
      for j in down:
        logging.error(f'#SID {j}')
    gmsh.model.occ.synchronize()

if mode=='sfc':
  for sfc,sid in sfcs:
    if len(sids)>0:
      if sid not in sids:
        continue
    logging.info(f'Check surface {sid}')
    exclude=[(s,i) for (s,i) in sfcs if i!=sid]
    (up,down)=gmsh.model.getAdjacencies(2,sid)
    s=gmsh.model.getEntityName(3,up[0])
    gmsh.model.removeEntities(vols)
    gmsh.model.removeEntities(exclude,recursive=True)
    try:
      gmsh.model.mesh.generate(2)
    except:
      logging.warning(f'WARNING: surface {sid} belonging to vol {up[0]},\"{s}\" did not mesh, listing connected curves')
      logging.warning(f'#SID {sid}')
      for j in down:
        logging.warning(f'#CID {j}')
    gmsh.model.occ.synchronize()

gmsh.finalize()
